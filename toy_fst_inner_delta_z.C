// Toy delta-z alignment solve for inner FST sensors only.
//
// Usage:
//   root -l -b -q 'toy_fst_inner_delta_z.C("align_test.root","fGeom.root")'
//
// This uses only branches already written to fwdAlign, plus fGeom.root to recover
// each sensor plane's U/V axes.  It is a first-order toy:
//
//   residual_u ~= slope_u_vs_z * deltaZ
//   residual_v ~= slope_v_vs_z * deltaZ
//
// with residual = measurement - prediction from GenFit.

#include "TCanvas.h"
#include "TFile.h"
#include "TGeoBBox.h"
#include "TGeoManager.h"
#include "TGeoMatrix.h"
#include "TGeoNode.h"
#include "TGeoShape.h"
#include "TGeoTube.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLine.h"
#include "TMath.h"
#include "TString.h"
#include "TStyle.h"
#include "TTree.h"
#include "TVector3.h"

#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <string>

namespace {
    constexpr int kFstNumDisks = 3;
    constexpr int kFstNumWedgePerDisk = 12;
    constexpr int kFstNumSensorsPerWedge = 3;
    constexpr int kFstNumSensors = kFstNumDisks * kFstNumWedgePerDisk * kFstNumSensorsPerWedge;
    constexpr double kInvalidAlignValue = -99999.0;

    struct SensorGeom {
        bool ok = false;
        TVector3 origin;
        TVector3 u;
        TVector3 v;
    };

    struct DeltaZAccumulator {
        int nRows = 0;
        int nTerms = 0;
        int nTerms0 = 0;
        int nTerms1 = 0;
        double numerator = 0.0;
        double denominator = 0.0;

        void add(double residual, double sigma, double slope, int dim) {
            if (!std::isfinite(residual) || !std::isfinite(sigma) || !std::isfinite(slope))
                return;
            if (std::abs(residual) > 0.9 * std::abs(kInvalidAlignValue))
                return;
            if (sigma <= 0.0 || std::abs(slope) < 1e-8)
                return;

            const double weight = 1.0 / (sigma * sigma);
            numerator += slope * residual * weight;
            denominator += slope * slope * weight;
            ++nTerms;
            if (dim == 0) ++nTerms0;
            if (dim == 1) ++nTerms1;
        }

        bool valid() const {
            return denominator > 0.0;
        }

        double value() const {
            return valid() ? numerator / denominator : std::numeric_limits<double>::quiet_NaN();
        }

        double error() const {
            return valid() ? 1.0 / std::sqrt(denominator) : std::numeric_limits<double>::quiet_NaN();
        }
    };

    bool hasBranch(TTree *tree, const char *name) {
        return tree && tree->GetBranch(name);
    }

    int globalSensorIndex(int disk, int wedge, int sensor) {
        return disk * kFstNumWedgePerDisk * kFstNumSensorsPerWedge
             + wedge * kFstNumSensorsPerWedge
             + sensor;
    }

    bool buildFstSensorGeometry(const char *geomFilename, std::array<SensorGeom, kFstNumSensors> &sensors) {
        static const int kEventSensorToFtusCopy[kFstNumSensorsPerWedge] = {3, 1, 2};
        static const int kElecToGeantWedge[kFstNumWedgePerDisk] =
            {2, 7, 1, 12, 6, 11, 5, 10, 4, 9, 3, 8};
        static const int kElecToGeantWedgeDisk2[kFstNumWedgePerDisk] =
            {7, 1, 12, 6, 11, 5, 10, 4, 9, 3, 8, 2};

        TFile geomFile(geomFilename, "READ");
        if (geomFile.IsZombie()) {
            std::cerr << "Cannot open geometry file: " << geomFilename << std::endl;
            return false;
        }

        TGeoManager *geom = dynamic_cast<TGeoManager*>(geomFile.Get("dyson"));
        if (!geom) {
            std::cerr << "Cannot find TGeoManager named dyson in " << geomFilename << std::endl;
            return false;
        }

        for (int disk = 0; disk < kFstNumDisks; ++disk) {
            const int fstdCopy = disk + 4;
            const int *wedgeMap = (fstdCopy == 5) ? kElecToGeantWedgeDisk2 : kElecToGeantWedge;

            for (int wedge = 0; wedge < kFstNumWedgePerDisk; ++wedge) {
                const int geantFstwCopy = wedgeMap[wedge];

                for (int sensor = 0; sensor < kFstNumSensorsPerWedge; ++sensor) {
                    const int ftusCopy = kEventSensorToFtusCopy[sensor];
                    const int globalSensor = globalSensorIndex(disk, wedge, sensor);
                    TString path = TString::Format(
                        "/HALL_1/CAVE_1/FSTM_1/FSTD_%d/FSTW_%d/FTUS_%d",
                        fstdCopy,
                        geantFstwCopy,
                        ftusCopy
                    );

                    if (!geom->cd(path)) {
                        std::cerr << "Cannot cd to " << path.Data() << std::endl;
                        continue;
                    }

                    TGeoMatrix *matrix = geom->GetCurrentMatrix();
                    TGeoNode *node = geom->GetCurrentNode();
                    if (!matrix || !node || !node->GetVolume())
                        continue;

                    const Double_t *translation = matrix->GetTranslation();
                    TVector3 origin(translation[0], translation[1], translation[2]);

                    TGeoShape *shape = node->GetVolume()->GetShape();
                    Double_t activeLocal[3] = {0.0, 0.0, 0.0};
                    bool haveActiveCenter = false;
                    if (auto tube = dynamic_cast<TGeoTubeSeg*>(shape)) {
                        const double rCenter = 0.5 * (tube->GetRmin() + tube->GetRmax());
                        const double phiCenter = 0.5 * (tube->GetPhi1() + tube->GetPhi2()) * TMath::DegToRad();
                        activeLocal[0] = rCenter * std::cos(phiCenter);
                        activeLocal[1] = rCenter * std::sin(phiCenter);
                        activeLocal[2] = 0.0;
                        haveActiveCenter = true;
                    } else if (auto box = dynamic_cast<TGeoBBox*>(shape)) {
                        const Double_t *boxOrigin = box->GetOrigin();
                        activeLocal[0] = boxOrigin[0];
                        activeLocal[1] = boxOrigin[1];
                        activeLocal[2] = boxOrigin[2];
                        haveActiveCenter = true;
                    }
                    if (haveActiveCenter) {
                        Double_t activeMaster[3] = {0.0, 0.0, 0.0};
                        matrix->LocalToMaster(activeLocal, activeMaster);
                        origin.SetXYZ(activeMaster[0], activeMaster[1], activeMaster[2]);
                    }

                    const Double_t *rot = matrix->GetRotationMatrix();
                    TVector3 u(rot[0], rot[3], rot[6]);
                    TVector3 v(rot[1], rot[4], rot[7]);
                    if (u.Cross(v).Z() < 0)
                        v = -v;

                    sensors[globalSensor].ok = true;
                    sensors[globalSensor].origin = origin;
                    sensors[globalSensor].u = u.Unit();
                    sensors[globalSensor].v = v.Unit();
                }
            }
        }

        return true;
    }

    void printResult(const char *label, const DeltaZAccumulator &acc) {
        if (!acc.valid()) {
            std::cout << label << ": no valid solve terms" << std::endl;
            return;
        }
        std::cout << Form(
            "%-18s  rows=%6d terms=%6d (u=%5d v=%5d)  deltaZ=% .6f cm  err=% .6f cm",
            label,
            acc.nRows,
            acc.nTerms,
            acc.nTerms0,
            acc.nTerms1,
            acc.value(),
            acc.error()
        ) << std::endl;
    }
}

void toy_fst_inner_delta_z(
    const char *alignFilename = "align_test.root",
    const char *geomFilename = "fGeom.root",
    const char *outputRoot = "toy_fst_inner_delta_z.root",
    const char *outputPdf = "toy_fst_inner_delta_z.pdf"
) {
    gStyle->SetOptStat(1110);

    std::array<SensorGeom, kFstNumSensors> sensorGeom;
    if (!buildFstSensorGeometry(geomFilename, sensorGeom))
        return;

    TFile alignFile(alignFilename, "READ");
    if (alignFile.IsZombie()) {
        std::cerr << "Cannot open alignment file: " << alignFilename << std::endl;
        return;
    }

    TTree *tree = dynamic_cast<TTree*>(alignFile.Get("fwdAlign"));
    if (!tree) {
        std::cerr << "Cannot find TTree fwdAlign in " << alignFilename << std::endl;
        return;
    }

    const char *requiredBranches[] = {
        "fstGlobalSensor", "fstDisk", "fstWedge", "fstSensor",
        "hasResidual", "residualDim", "fitConverged",
        "trackPx", "trackPy", "trackPz", "trackP", "trackEta",
        "resUnbiased0", "resUnbiased1",
        "resUnbiasedSigma0", "resUnbiasedSigma1"
    };
    for (const char *branch : requiredBranches) {
        if (!hasBranch(tree, branch)) {
            std::cerr << "Missing required branch: " << branch << std::endl;
            return;
        }
    }

    int fstGlobalSensor = -1;
    int fstDisk = -1;
    int fstWedge = -1;
    int fstSensor = -1;
    int hasResidual = 0;
    int residualDim = 0;
    int fitConverged = 0;
    int trackNHitsFit = 0;
    float trackPx = 0.0f;
    float trackPy = 0.0f;
    float trackPz = 0.0f;
    float trackP = 0.0f;
    float trackEta = 0.0f;
    float resUnbiased0 = kInvalidAlignValue;
    float resUnbiased1 = kInvalidAlignValue;
    float resUnbiasedSigma0 = kInvalidAlignValue;
    float resUnbiasedSigma1 = kInvalidAlignValue;

    tree->SetBranchAddress("fstGlobalSensor", &fstGlobalSensor);
    tree->SetBranchAddress("fstDisk", &fstDisk);
    tree->SetBranchAddress("fstWedge", &fstWedge);
    tree->SetBranchAddress("fstSensor", &fstSensor);
    tree->SetBranchAddress("hasResidual", &hasResidual);
    tree->SetBranchAddress("residualDim", &residualDim);
    tree->SetBranchAddress("fitConverged", &fitConverged);
    if (hasBranch(tree, "trackNHitsFit"))
        tree->SetBranchAddress("trackNHitsFit", &trackNHitsFit);
    tree->SetBranchAddress("trackPx", &trackPx);
    tree->SetBranchAddress("trackPy", &trackPy);
    tree->SetBranchAddress("trackPz", &trackPz);
    tree->SetBranchAddress("trackP", &trackP);
    tree->SetBranchAddress("trackEta", &trackEta);
    tree->SetBranchAddress("resUnbiased0", &resUnbiased0);
    tree->SetBranchAddress("resUnbiased1", &resUnbiased1);
    tree->SetBranchAddress("resUnbiasedSigma0", &resUnbiasedSigma0);
    tree->SetBranchAddress("resUnbiasedSigma1", &resUnbiasedSigma1);

    DeltaZAccumulator allInner;
    std::array<DeltaZAccumulator, kFstNumDisks> byDisk;
    std::array<DeltaZAccumulator, kFstNumSensors> bySensor;

    TH1D *hDzTerm0 = new TH1D("hDzTerm0", "Inner FST toy #Delta z candidates from residual 0;resUnbiased0/slope0 [cm];terms", 160, -20, 20);
    TH1D *hDzTerm1 = new TH1D("hDzTerm1", "Inner FST toy #Delta z candidates from residual 1;resUnbiased1/slope1 [cm];terms", 160, -20, 20);
    TH1D *hSlope0 = new TH1D("hSlope0", "Inner FST local slope 0;dmeas0/dz;hits", 120, -0.8, 0.8);
    TH1D *hSlope1 = new TH1D("hSlope1", "Inner FST local slope 1;dmeas1/dz;hits", 120, -0.8, 0.8);
    TH2D *hResVsSlope0 = new TH2D("hResVsSlope0", "Inner FST resUnbiased0 vs slope0;dmeas0/dz;resUnbiased0 [cm]", 120, -0.8, 0.8, 120, -1.0, 1.0);
    TH2D *hResVsSlope1 = new TH2D("hResVsSlope1", "Inner FST resUnbiased1 vs slope1;dmeas1/dz;resUnbiased1 [cm]", 120, -0.8, 0.8, 120, -1.0, 1.0);

    Long64_t nRead = 0;
    Long64_t nSelected = 0;
    Long64_t nNoGeom = 0;
    Long64_t nBadSlope = 0;

    const Long64_t nEntries = tree->GetEntries();
    for (Long64_t i = 0; i < nEntries; ++i) {
        tree->GetEntry(i);
        ++nRead;

        if (fstSensor != 0)
            continue;
        if (fstGlobalSensor < 0 || fstGlobalSensor >= kFstNumSensors)
            continue;
        if (fstDisk < 0 || fstDisk >= kFstNumDisks)
            continue;
        if (hasResidual != 1 || residualDim < 2 || fitConverged != 1)
            continue;
        if (!(trackEta >= 2.5 && trackEta <= 4.0))
            continue;
        if (!(trackP > 0.5))
            continue;
        if (hasBranch(tree, "trackNHitsFit") && trackNHitsFit > 0 && trackNHitsFit < 5)
            continue;
        if (!sensorGeom[fstGlobalSensor].ok) {
            ++nNoGeom;
            continue;
        }
        if (!std::isfinite(trackPz) || std::abs(trackPz) < 1e-8)
            continue;

        const TVector3 mom(trackPx, trackPy, trackPz);
        const double slope0 = mom.Dot(sensorGeom[fstGlobalSensor].u) / trackPz;
        const double slope1 = mom.Dot(sensorGeom[fstGlobalSensor].v) / trackPz;
        if (!std::isfinite(slope0) || !std::isfinite(slope1)) {
            ++nBadSlope;
            continue;
        }

        ++nSelected;
        ++allInner.nRows;
        ++byDisk[fstDisk].nRows;
        ++bySensor[fstGlobalSensor].nRows;

        allInner.add(resUnbiased0, resUnbiasedSigma0, slope0, 0);
        byDisk[fstDisk].add(resUnbiased0, resUnbiasedSigma0, slope0, 0);
        bySensor[fstGlobalSensor].add(resUnbiased0, resUnbiasedSigma0, slope0, 0);
        if (std::abs(slope0) > 1e-8)
            hDzTerm0->Fill(resUnbiased0 / slope0);
        hSlope0->Fill(slope0);
        hResVsSlope0->Fill(slope0, resUnbiased0);

        allInner.add(resUnbiased1, resUnbiasedSigma1, slope1, 1);
        byDisk[fstDisk].add(resUnbiased1, resUnbiasedSigma1, slope1, 1);
        bySensor[fstGlobalSensor].add(resUnbiased1, resUnbiasedSigma1, slope1, 1);
        if (std::abs(slope1) > 1e-8)
            hDzTerm1->Fill(resUnbiased1 / slope1);
        hSlope1->Fill(slope1);
        hResVsSlope1->Fill(slope1, resUnbiased1);
    }

    std::cout << "\nToy inner-FST delta-z solve" << std::endl;
    std::cout << "  input tree entries: " << nRead << std::endl;
    std::cout << "  selected inner FST rows: " << nSelected << std::endl;
    std::cout << "  skipped no geometry: " << nNoGeom << std::endl;
    std::cout << "  skipped bad slope: " << nBadSlope << std::endl;
    std::cout << "  selection: fstSensor==0, hasResidual, residualDim>=2, fitConverged, 2.5<=eta<=4.0, P>0.5, NHitsFit>=5 if available" << std::endl;
    std::cout << "  sign convention: residual = measurement - prediction, so deltaZ = residual / slope" << std::endl;
    printResult("all inner", allInner);
    for (int disk = 0; disk < kFstNumDisks; ++disk) {
        printResult(Form("disk %d", disk), byDisk[disk]);
    }

    TGraphErrors *gDisk = new TGraphErrors();
    gDisk->SetName("gDeltaZByDisk");
    gDisk->SetTitle("Toy inner-FST #Delta z by disk;fstDisk;#Delta z [cm]");
    for (int disk = 0; disk < kFstNumDisks; ++disk) {
        if (!byDisk[disk].valid())
            continue;
        const int ip = gDisk->GetN();
        gDisk->SetPoint(ip, disk, byDisk[disk].value());
        gDisk->SetPointError(ip, 0.0, byDisk[disk].error());
    }

    TGraphErrors *gSensor = new TGraphErrors();
    gSensor->SetName("gDeltaZByInnerGlobalSensor");
    gSensor->SetTitle("Toy inner-FST #Delta z by inner global sensor;fstGlobalSensor;#Delta z [cm]");
    for (int globalSensor = 0; globalSensor < kFstNumSensors; globalSensor += kFstNumSensorsPerWedge) {
        if (!bySensor[globalSensor].valid() || bySensor[globalSensor].nRows < 5)
            continue;
        const int ip = gSensor->GetN();
        gSensor->SetPoint(ip, globalSensor, bySensor[globalSensor].value());
        gSensor->SetPointError(ip, 0.0, bySensor[globalSensor].error());
    }

    TFile out(outputRoot, "RECREATE");
    hDzTerm0->Write();
    hDzTerm1->Write();
    hSlope0->Write();
    hSlope1->Write();
    hResVsSlope0->Write();
    hResVsSlope1->Write();
    gDisk->Write();
    gSensor->Write();
    out.Close();

    TCanvas c("cToyFstInnerDeltaZ", "Toy inner FST delta-z alignment", 1400, 1000);
    c.Divide(2, 2);
    c.cd(1);
    hDzTerm0->Draw();
    c.cd(2);
    hDzTerm1->Draw();
    c.cd(3);
    gDisk->SetMarkerStyle(20);
    gDisk->Draw("AP");
    {
        TLine zero(-0.5, 0.0, 2.5, 0.0);
        zero.SetLineStyle(2);
        zero.SetLineColor(kRed + 1);
        zero.Draw();
    }
    c.cd(4);
    gSensor->SetMarkerStyle(20);
    gSensor->Draw("AP");
    {
        TLine zero(-1.0, 0.0, 108.0, 0.0);
        zero.SetLineStyle(2);
        zero.SetLineColor(kRed + 1);
        zero.Draw();
    }
    c.SaveAs(outputPdf);

    std::cout << "Wrote " << outputRoot << std::endl;
    std::cout << "Wrote " << outputPdf << std::endl;
}
