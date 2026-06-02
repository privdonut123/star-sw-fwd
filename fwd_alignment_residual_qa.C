// ROOT macro for first-pass Forward alignment residual QA.
//
// Usage:
//   root4star -l -b -q 'fwd_alignment_residual_qa.C("align_test.root","fwd_align_qa")'
//
// Outputs:
//   fwd_align_qa.root
//   fwd_align_qa.pdf

#include "TCanvas.h"
#include "TCut.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TKey.h"
#include "TLegend.h"
#include "TLine.h"
#include "TProfile.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TTree.h"
#include "TPaveText.h"
#include "TString.h"

#include <iostream>
#include <vector>

bool fwdAlignHasBranch(TTree *tree, const char *name) {
    return tree && tree->GetBranch(name);
}

void fwdAlignDrawLine(double x1, double y1, double x2, double y2, int color = kRed + 1, int style = 2) {
    TLine *line = new TLine(x1, y1, x2, y2);
    line->SetLineColor(color);
    line->SetLineStyle(style);
    line->Draw();
}

void fwdAlignSavePage(TCanvas *canvas, const TString &pdf) {
    canvas->Print(pdf);
}

void fwd_alignment_residual_qa(
    const char *inputFilename = "align_test.root",
    const char *outputPrefix = "fwd_align_qa",
    bool requireFullyConverged = false,
    double residualRangeMicron = 2000.0
) {
    gROOT->SetBatch(kTRUE);
    gStyle->SetOptStat(1110);
    gStyle->SetOptFit(0);
    gStyle->SetPalette(1);

    TFile *inputFile = TFile::Open(inputFilename, "READ");
    if (!inputFile || inputFile->IsZombie()) {
        std::cerr << "Cannot open input file: " << inputFilename << std::endl;
        return;
    }

    TTree *tree = dynamic_cast<TTree*>(inputFile->Get("fwdAlign"));
    if (!tree) {
        std::cerr << "Cannot find TTree fwdAlign in " << inputFilename << std::endl;
        inputFile->ls();
        return;
    }

    const char *requiredBranches[] = {
        "run", "event", "trackIndex", "pointIndex", "detId", "hitId",
        "fstGlobalSensor", "fstDisk", "fstWedge", "fstSensor",
        "measurementDim", "residualDim", "hasResidual", "fitConverged",
        "fitConvergedFully", "sorting", "meas0", "meas1",
        "resBiased0", "resBiased1", "resUnbiased0", "resUnbiased1"
    };
    const int nRequiredBranches = sizeof(requiredBranches) / sizeof(requiredBranches[0]);
    bool missingBranch = false;
    for (int i = 0; i < nRequiredBranches; ++i) {
        if (!fwdAlignHasBranch(tree, requiredBranches[i])) {
            std::cerr << "Missing required fwdAlign branch: " << requiredBranches[i] << std::endl;
            missingBranch = true;
        }
    }
    if (missingBranch) {
        std::cerr << "Stopping because the alignment tree schema is not the expected one." << std::endl;
        return;
    }

    bool hasPullBranches =
        fwdAlignHasBranch(tree, "pullBiased0") &&
        fwdAlignHasBranch(tree, "pullBiased1") &&
        fwdAlignHasBranch(tree, "pullUnbiased0") &&
        fwdAlignHasBranch(tree, "pullUnbiased1") &&
        fwdAlignHasBranch(tree, "resBiasedSigma0") &&
        fwdAlignHasBranch(tree, "resBiasedSigma1") &&
        fwdAlignHasBranch(tree, "resUnbiasedSigma0") &&
        fwdAlignHasBranch(tree, "resUnbiasedSigma1");

    TString rootOutput = TString::Format("%s.root", outputPrefix);
    TString pdfOutput = TString::Format("%s.pdf", outputPrefix);
    TFile *outputFile = TFile::Open(rootOutput, "RECREATE");
    if (!outputFile || outputFile->IsZombie()) {
        std::cerr << "Cannot create output file: " << rootOutput.Data() << std::endl;
        return;
    }

    TString fitCut = requireFullyConverged ? "fitConvergedFully>0" : "fitConverged>0";
    TCut fstResidualCut = TString::Format(
        "detId==45&&hasResidual>0&&residualDim==2&&fstGlobalSensor>=0&&%s",
        fitCut.Data()
    ).Data();
    TCut fstPullCut = TString::Format(
        "detId==45&&hasResidual>0&&residualDim==2&&fstGlobalSensor>=0&&%s"
        "&&pullBiased0>-90000&&pullBiased1>-90000&&pullUnbiased0>-90000&&pullUnbiased1>-90000",
        fitCut.Data()
    ).Data();
    TCut fstAllCut = "detId==45&&fstGlobalSensor>=0";

    Long64_t nEntries = tree->GetEntries();
    Long64_t nFstRows = tree->GetEntries(fstAllCut);
    Long64_t nFstResidualRows = tree->GetEntries(fstResidualCut);

    std::cout << "Input file: " << inputFilename << std::endl;
    std::cout << "fwdAlign entries: " << nEntries << std::endl;
    std::cout << "FST rows: " << nFstRows << std::endl;
    std::cout << "FST residual rows after cuts: " << nFstResidualRows << std::endl;
    std::cout << "Base cut: " << TString(fstResidualCut).Data() << std::endl;
    std::cout << "Residual units in plots: microns (tree stores cm)" << std::endl;
    std::cout << "Pull branches available: " << (hasPullBranches ? "yes" : "no") << std::endl;

    TCanvas *canvas = new TCanvas("cFwdAlignQa", "Forward alignment QA", 1200, 900);
    canvas->Print(pdfOutput + "[");

    // Summary page.
    canvas->Clear();
    TPaveText *summary = new TPaveText(0.08, 0.14, 0.92, 0.88, "NDC");
    summary->SetFillColor(0);
    summary->SetTextAlign(12);
    summary->AddText("Forward alignment residual QA");
    summary->AddText(TString::Format("Input: %s", inputFilename));
    summary->AddText(TString::Format("Output ROOT: %s", rootOutput.Data()));
    summary->AddText(TString::Format("Output PDF: %s", pdfOutput.Data()));
    summary->AddText(TString::Format("Total fwdAlign rows: %lld", nEntries));
    summary->AddText(TString::Format("FST rows: %lld", nFstRows));
    summary->AddText(TString::Format("FST residual rows after cuts: %lld", nFstResidualRows));
    summary->AddText(TString::Format("Fit cut: %s", fitCut.Data()));
    summary->AddText("Residual plots use resUnbiased/resBiased * 10000, so units are microns.");
    summary->AddText(hasPullBranches ? "Pull plots are true GenFit residual/sigma pulls." : "Pull branches are not present in this input; rerun afterburner with updated maker.");
    summary->AddText("FST local coordinates are the GenFit planar-measurement coordinates.");
    summary->AddText("FTT is intentionally not analyzed in this macro.");
    summary->Draw();
    fwdAlignSavePage(canvas, pdfOutput);

    // Basic schema / hit-test page.
    TH1D *hDetId = new TH1D("hDetId", "Alignment rows by detector;detId;rows", 80, 0, 80);
    TH2D *hDim = new TH2D("hDim", "FST measurement and residual dimensions;measurementDim;residualDim", 6, -0.5, 5.5, 6, -0.5, 5.5);
    TH1D *hHasResidual = new TH1D("hHasResidual", "FST hasResidual flag;hasResidual;rows", 3, -0.5, 2.5);
    TH1D *hSortMinusSensor = new TH1D("hSortMinusSensor", "FST sorting - fstGlobalSensor, expect 1;sorting - fstGlobalSensor;rows", 21, -9.5, 11.5);
    tree->Draw("detId>>hDetId", "", "goff");
    tree->Draw("residualDim:measurementDim>>hDim", fstAllCut, "goff");
    tree->Draw("hasResidual>>hHasResidual", fstAllCut, "goff");
    tree->Draw("sorting-fstGlobalSensor>>hSortMinusSensor", fstAllCut, "goff");

    canvas->Clear();
    canvas->Divide(2, 2);
    canvas->cd(1); hDetId->Draw();
    canvas->cd(2); hDim->Draw("COLZ TEXT");
    canvas->cd(3); hHasResidual->Draw();
    canvas->cd(4); hSortMinusSensor->Draw();
    fwdAlignSavePage(canvas, pdfOutput);

    // Residual distributions.
    TH1D *hResU = new TH1D("hResU", "FST unbiased residual 0;resUnbiased0 [um];rows", 160, -residualRangeMicron, residualRangeMicron);
    TH1D *hResV = new TH1D("hResV", "FST unbiased residual 1;resUnbiased1 [um];rows", 160, -residualRangeMicron, residualRangeMicron);
    TH1D *hBiasedU = new TH1D("hBiasedU", "FST biased residual 0;resBiased0 [um];rows", 160, -residualRangeMicron, residualRangeMicron);
    TH1D *hBiasedV = new TH1D("hBiasedV", "FST biased residual 1;resBiased1 [um];rows", 160, -residualRangeMicron, residualRangeMicron);
    tree->Draw("resUnbiased0*10000>>hResU", fstResidualCut, "goff");
    tree->Draw("resUnbiased1*10000>>hResV", fstResidualCut, "goff");
    tree->Draw("resBiased0*10000>>hBiasedU", fstResidualCut, "goff");
    tree->Draw("resBiased1*10000>>hBiasedV", fstResidualCut, "goff");

    canvas->Clear();
    canvas->Divide(2, 2);
    canvas->cd(1); hResU->Draw(); fwdAlignDrawLine(0, 0, 0, hResU->GetMaximum());
    canvas->cd(2); hResV->Draw(); fwdAlignDrawLine(0, 0, 0, hResV->GetMaximum());
    canvas->cd(3); hBiasedU->Draw(); fwdAlignDrawLine(0, 0, 0, hBiasedU->GetMaximum());
    canvas->cd(4); hBiasedV->Draw(); fwdAlignDrawLine(0, 0, 0, hBiasedV->GetMaximum());
    fwdAlignSavePage(canvas, pdfOutput);

    if (hasPullBranches) {
        TH1D *hPullU = new TH1D("hPullU", "FST unbiased pull 0;pullUnbiased0;rows", 160, -10, 10);
        TH1D *hPullV = new TH1D("hPullV", "FST unbiased pull 1;pullUnbiased1;rows", 160, -10, 10);
        TH1D *hPullBiasedU = new TH1D("hPullBiasedU", "FST biased pull 0;pullBiased0;rows", 160, -10, 10);
        TH1D *hPullBiasedV = new TH1D("hPullBiasedV", "FST biased pull 1;pullBiased1;rows", 160, -10, 10);
        tree->Draw("pullUnbiased0>>hPullU", fstPullCut, "goff");
        tree->Draw("pullUnbiased1>>hPullV", fstPullCut, "goff");
        tree->Draw("pullBiased0>>hPullBiasedU", fstPullCut, "goff");
        tree->Draw("pullBiased1>>hPullBiasedV", fstPullCut, "goff");

        canvas->Clear();
        canvas->Divide(2, 2);
        canvas->cd(1); hPullU->Draw(); fwdAlignDrawLine(0, 0, 0, hPullU->GetMaximum());
        canvas->cd(2); hPullV->Draw(); fwdAlignDrawLine(0, 0, 0, hPullV->GetMaximum());
        canvas->cd(3); hPullBiasedU->Draw(); fwdAlignDrawLine(0, 0, 0, hPullBiasedU->GetMaximum());
        canvas->cd(4); hPullBiasedV->Draw(); fwdAlignDrawLine(0, 0, 0, hPullBiasedV->GetMaximum());
        fwdAlignSavePage(canvas, pdfOutput);

        TH1D *hSigmaU = new TH1D("hSigmaU", "FST unbiased residual sigma 0;#sigma_{res,0} [um];rows", 160, 0, residualRangeMicron);
        TH1D *hSigmaV = new TH1D("hSigmaV", "FST unbiased residual sigma 1;#sigma_{res,1} [um];rows", 160, 0, residualRangeMicron);
        TH1D *hSigmaBiasedU = new TH1D("hSigmaBiasedU", "FST biased residual sigma 0;#sigma_{res,0} [um];rows", 160, 0, residualRangeMicron);
        TH1D *hSigmaBiasedV = new TH1D("hSigmaBiasedV", "FST biased residual sigma 1;#sigma_{res,1} [um];rows", 160, 0, residualRangeMicron);
        tree->Draw("resUnbiasedSigma0*10000>>hSigmaU", fstPullCut, "goff");
        tree->Draw("resUnbiasedSigma1*10000>>hSigmaV", fstPullCut, "goff");
        tree->Draw("resBiasedSigma0*10000>>hSigmaBiasedU", fstPullCut, "goff");
        tree->Draw("resBiasedSigma1*10000>>hSigmaBiasedV", fstPullCut, "goff");

        canvas->Clear();
        canvas->Divide(2, 2);
        canvas->cd(1); hSigmaU->Draw();
        canvas->cd(2); hSigmaV->Draw();
        canvas->cd(3); hSigmaBiasedU->Draw();
        canvas->cd(4); hSigmaBiasedV->Draw();
        fwdAlignSavePage(canvas, pdfOutput);
    }

    // Occupancy and FST geometry grouping.
    TH1D *hFstDisk = new TH1D("hFstDisk", "FST residual rows by disk;fstDisk;rows", 3, -0.5, 2.5);
    TH1D *hFstSensor = new TH1D("hFstSensor", "FST residual rows by global sensor;fstGlobalSensor;rows", 108, -0.5, 107.5);
    TH2D *hWedgeDisk = new TH2D("hWedgeDisk", "FST residual occupancy;fstDisk;fstWedge", 3, -0.5, 2.5, 12, -0.5, 11.5);
    TH2D *hSensorDisk = new TH2D("hSensorDisk", "FST sensor-in-wedge residual occupancy;fstDisk;fstSensor", 3, -0.5, 2.5, 3, -0.5, 2.5);
    tree->Draw("fstDisk>>hFstDisk", fstResidualCut, "goff");
    tree->Draw("fstGlobalSensor>>hFstSensor", fstResidualCut, "goff");
    tree->Draw("fstWedge:fstDisk>>hWedgeDisk", fstResidualCut, "goff");
    tree->Draw("fstSensor:fstDisk>>hSensorDisk", fstResidualCut, "goff");

    canvas->Clear();
    canvas->Divide(2, 2);
    canvas->cd(1); hFstDisk->Draw();
    canvas->cd(2); hFstSensor->Draw();
    canvas->cd(3); hWedgeDisk->Draw("COLZ TEXT");
    canvas->cd(4); hSensorDisk->Draw("COLZ TEXT");
    fwdAlignSavePage(canvas, pdfOutput);

    // Mean residual profiles by detector grouping.
    TProfile *pResUByDisk = new TProfile("pResUByDisk", "Mean FST unbiased residual 0 by disk;fstDisk;<resUnbiased0> [um]", 3, -0.5, 2.5);
    TProfile *pResVByDisk = new TProfile("pResVByDisk", "Mean FST unbiased residual 1 by disk;fstDisk;<resUnbiased1> [um]", 3, -0.5, 2.5);
    TProfile *pResUByWedge = new TProfile("pResUByWedge", "Mean FST unbiased residual 0 by wedge;fstWedge;<resUnbiased0> [um]", 12, -0.5, 11.5);
    TProfile *pResVByWedge = new TProfile("pResVByWedge", "Mean FST unbiased residual 1 by wedge;fstWedge;<resUnbiased1> [um]", 12, -0.5, 11.5);
    tree->Draw("resUnbiased0*10000:fstDisk>>pResUByDisk", fstResidualCut, "goff");
    tree->Draw("resUnbiased1*10000:fstDisk>>pResVByDisk", fstResidualCut, "goff");
    tree->Draw("resUnbiased0*10000:fstWedge>>pResUByWedge", fstResidualCut, "goff");
    tree->Draw("resUnbiased1*10000:fstWedge>>pResVByWedge", fstResidualCut, "goff");

    canvas->Clear();
    canvas->Divide(2, 2);
    canvas->cd(1); pResUByDisk->Draw("E1"); fwdAlignDrawLine(-0.5, 0, 2.5, 0);
    canvas->cd(2); pResVByDisk->Draw("E1"); fwdAlignDrawLine(-0.5, 0, 2.5, 0);
    canvas->cd(3); pResUByWedge->Draw("E1"); fwdAlignDrawLine(-0.5, 0, 11.5, 0);
    canvas->cd(4); pResVByWedge->Draw("E1"); fwdAlignDrawLine(-0.5, 0, 11.5, 0);
    fwdAlignSavePage(canvas, pdfOutput);

    if (hasPullBranches) {
        TProfile *pPullUByDisk = new TProfile("pPullUByDisk", "Mean FST unbiased pull 0 by disk;fstDisk;<pullUnbiased0>", 3, -0.5, 2.5);
        TProfile *pPullVByDisk = new TProfile("pPullVByDisk", "Mean FST unbiased pull 1 by disk;fstDisk;<pullUnbiased1>", 3, -0.5, 2.5);
        TProfile *pPullUByWedge = new TProfile("pPullUByWedge", "Mean FST unbiased pull 0 by wedge;fstWedge;<pullUnbiased0>", 12, -0.5, 11.5);
        TProfile *pPullVByWedge = new TProfile("pPullVByWedge", "Mean FST unbiased pull 1 by wedge;fstWedge;<pullUnbiased1>", 12, -0.5, 11.5);
        tree->Draw("pullUnbiased0:fstDisk>>pPullUByDisk", fstPullCut, "goff");
        tree->Draw("pullUnbiased1:fstDisk>>pPullVByDisk", fstPullCut, "goff");
        tree->Draw("pullUnbiased0:fstWedge>>pPullUByWedge", fstPullCut, "goff");
        tree->Draw("pullUnbiased1:fstWedge>>pPullVByWedge", fstPullCut, "goff");

        canvas->Clear();
        canvas->Divide(2, 2);
        canvas->cd(1); pPullUByDisk->Draw("E1"); fwdAlignDrawLine(-0.5, 0, 2.5, 0);
        canvas->cd(2); pPullVByDisk->Draw("E1"); fwdAlignDrawLine(-0.5, 0, 2.5, 0);
        canvas->cd(3); pPullUByWedge->Draw("E1"); fwdAlignDrawLine(-0.5, 0, 11.5, 0);
        canvas->cd(4); pPullVByWedge->Draw("E1"); fwdAlignDrawLine(-0.5, 0, 11.5, 0);
        fwdAlignSavePage(canvas, pdfOutput);
    }

    TProfile *pResUByGlobalSensor = new TProfile("pResUByGlobalSensor", "Mean FST unbiased residual 0 by global sensor;fstGlobalSensor;<resUnbiased0> [um]", 108, -0.5, 107.5);
    TProfile *pResVByGlobalSensor = new TProfile("pResVByGlobalSensor", "Mean FST unbiased residual 1 by global sensor;fstGlobalSensor;<resUnbiased1> [um]", 108, -0.5, 107.5);
    TH2D *hResUByGlobalSensor = new TH2D("hResUByGlobalSensor", "FST unbiased residual 0 by global sensor;fstGlobalSensor;resUnbiased0 [um]", 108, -0.5, 107.5, 160, -residualRangeMicron, residualRangeMicron);
    TH2D *hResVByGlobalSensor = new TH2D("hResVByGlobalSensor", "FST unbiased residual 1 by global sensor;fstGlobalSensor;resUnbiased1 [um]", 108, -0.5, 107.5, 160, -residualRangeMicron, residualRangeMicron);
    tree->Draw("resUnbiased0*10000:fstGlobalSensor>>pResUByGlobalSensor", fstResidualCut, "goff");
    tree->Draw("resUnbiased1*10000:fstGlobalSensor>>pResVByGlobalSensor", fstResidualCut, "goff");
    tree->Draw("resUnbiased0*10000:fstGlobalSensor>>hResUByGlobalSensor", fstResidualCut, "goff");
    tree->Draw("resUnbiased1*10000:fstGlobalSensor>>hResVByGlobalSensor", fstResidualCut, "goff");

    canvas->Clear();
    canvas->Divide(1, 2);
    canvas->cd(1); pResUByGlobalSensor->Draw("E1"); fwdAlignDrawLine(-0.5, 0, 107.5, 0);
    canvas->cd(2); pResVByGlobalSensor->Draw("E1"); fwdAlignDrawLine(-0.5, 0, 107.5, 0);
    fwdAlignSavePage(canvas, pdfOutput);

    canvas->Clear();
    canvas->Divide(1, 2);
    canvas->cd(1); hResUByGlobalSensor->Draw("COLZ");
    canvas->cd(2); hResVByGlobalSensor->Draw("COLZ");
    fwdAlignSavePage(canvas, pdfOutput);

    if (hasPullBranches) {
        TProfile *pPullUByGlobalSensor = new TProfile("pPullUByGlobalSensor", "Mean FST unbiased pull 0 by global sensor;fstGlobalSensor;<pullUnbiased0>", 108, -0.5, 107.5);
        TProfile *pPullVByGlobalSensor = new TProfile("pPullVByGlobalSensor", "Mean FST unbiased pull 1 by global sensor;fstGlobalSensor;<pullUnbiased1>", 108, -0.5, 107.5);
        TH2D *hPullUByGlobalSensor = new TH2D("hPullUByGlobalSensor", "FST unbiased pull 0 by global sensor;fstGlobalSensor;pullUnbiased0", 108, -0.5, 107.5, 160, -10, 10);
        TH2D *hPullVByGlobalSensor = new TH2D("hPullVByGlobalSensor", "FST unbiased pull 1 by global sensor;fstGlobalSensor;pullUnbiased1", 108, -0.5, 107.5, 160, -10, 10);
        tree->Draw("pullUnbiased0:fstGlobalSensor>>pPullUByGlobalSensor", fstPullCut, "goff");
        tree->Draw("pullUnbiased1:fstGlobalSensor>>pPullVByGlobalSensor", fstPullCut, "goff");
        tree->Draw("pullUnbiased0:fstGlobalSensor>>hPullUByGlobalSensor", fstPullCut, "goff");
        tree->Draw("pullUnbiased1:fstGlobalSensor>>hPullVByGlobalSensor", fstPullCut, "goff");

        canvas->Clear();
        canvas->Divide(1, 2);
        canvas->cd(1); pPullUByGlobalSensor->Draw("E1"); fwdAlignDrawLine(-0.5, 0, 107.5, 0);
        canvas->cd(2); pPullVByGlobalSensor->Draw("E1"); fwdAlignDrawLine(-0.5, 0, 107.5, 0);
        fwdAlignSavePage(canvas, pdfOutput);

        canvas->Clear();
        canvas->Divide(1, 2);
        canvas->cd(1); hPullUByGlobalSensor->Draw("COLZ");
        canvas->cd(2); hPullVByGlobalSensor->Draw("COLZ");
        fwdAlignSavePage(canvas, pdfOutput);
    }

    // Local-coordinate dependence tests. These are useful for spotting rotations or scale effects.
    TH2D *hResUVsMeas0 = new TH2D("hResUVsMeas0", "FST residual 0 vs local meas0;meas0 [cm];resUnbiased0 [um]", 120, -40, 40, 160, -residualRangeMicron, residualRangeMicron);
    TH2D *hResUVsMeas1 = new TH2D("hResUVsMeas1", "FST residual 0 vs local meas1;meas1 [cm];resUnbiased0 [um]", 120, -40, 40, 160, -residualRangeMicron, residualRangeMicron);
    TH2D *hResVVsMeas0 = new TH2D("hResVVsMeas0", "FST residual 1 vs local meas0;meas0 [cm];resUnbiased1 [um]", 120, -40, 40, 160, -residualRangeMicron, residualRangeMicron);
    TH2D *hResVVsMeas1 = new TH2D("hResVVsMeas1", "FST residual 1 vs local meas1;meas1 [cm];resUnbiased1 [um]", 120, -40, 40, 160, -residualRangeMicron, residualRangeMicron);
    tree->Draw("resUnbiased0*10000:meas0>>hResUVsMeas0", fstResidualCut, "goff");
    tree->Draw("resUnbiased0*10000:meas1>>hResUVsMeas1", fstResidualCut, "goff");
    tree->Draw("resUnbiased1*10000:meas0>>hResVVsMeas0", fstResidualCut, "goff");
    tree->Draw("resUnbiased1*10000:meas1>>hResVVsMeas1", fstResidualCut, "goff");

    canvas->Clear();
    canvas->Divide(2, 2);
    canvas->cd(1); hResUVsMeas0->Draw("COLZ");
    canvas->cd(2); hResUVsMeas1->Draw("COLZ");
    canvas->cd(3); hResVVsMeas0->Draw("COLZ");
    canvas->cd(4); hResVVsMeas1->Draw("COLZ");
    fwdAlignSavePage(canvas, pdfOutput);

    if (hasPullBranches) {
        TH2D *hPullUVsMeas0 = new TH2D("hPullUVsMeas0", "FST pull 0 vs local meas0;meas0 [cm];pullUnbiased0", 120, -40, 40, 160, -10, 10);
        TH2D *hPullUVsMeas1 = new TH2D("hPullUVsMeas1", "FST pull 0 vs local meas1;meas1 [cm];pullUnbiased0", 120, -40, 40, 160, -10, 10);
        TH2D *hPullVVsMeas0 = new TH2D("hPullVVsMeas0", "FST pull 1 vs local meas0;meas0 [cm];pullUnbiased1", 120, -40, 40, 160, -10, 10);
        TH2D *hPullVVsMeas1 = new TH2D("hPullVVsMeas1", "FST pull 1 vs local meas1;meas1 [cm];pullUnbiased1", 120, -40, 40, 160, -10, 10);
        tree->Draw("pullUnbiased0:meas0>>hPullUVsMeas0", fstPullCut, "goff");
        tree->Draw("pullUnbiased0:meas1>>hPullUVsMeas1", fstPullCut, "goff");
        tree->Draw("pullUnbiased1:meas0>>hPullVVsMeas0", fstPullCut, "goff");
        tree->Draw("pullUnbiased1:meas1>>hPullVVsMeas1", fstPullCut, "goff");

        canvas->Clear();
        canvas->Divide(2, 2);
        canvas->cd(1); hPullUVsMeas0->Draw("COLZ");
        canvas->cd(2); hPullUVsMeas1->Draw("COLZ");
        canvas->cd(3); hPullVVsMeas0->Draw("COLZ");
        canvas->cd(4); hPullVVsMeas1->Draw("COLZ");
        fwdAlignSavePage(canvas, pdfOutput);
    }

    // Disk-resolved profiles; useful as first crude alignment constants.
    TProfile *pResUBySensorDisk0 = new TProfile("pResUBySensorDisk0", "Disk 0 mean residual 0 by sensor;fstGlobalSensor;<resUnbiased0> [um]", 36, -0.5, 35.5);
    TProfile *pResUBySensorDisk1 = new TProfile("pResUBySensorDisk1", "Disk 1 mean residual 0 by sensor;fstGlobalSensor;<resUnbiased0> [um]", 36, 35.5, 71.5);
    TProfile *pResUBySensorDisk2 = new TProfile("pResUBySensorDisk2", "Disk 2 mean residual 0 by sensor;fstGlobalSensor;<resUnbiased0> [um]", 36, 71.5, 107.5);
    TProfile *pResVBySensorDisk0 = new TProfile("pResVBySensorDisk0", "Disk 0 mean residual 1 by sensor;fstGlobalSensor;<resUnbiased1> [um]", 36, -0.5, 35.5);
    TProfile *pResVBySensorDisk1 = new TProfile("pResVBySensorDisk1", "Disk 1 mean residual 1 by sensor;fstGlobalSensor;<resUnbiased1> [um]", 36, 35.5, 71.5);
    TProfile *pResVBySensorDisk2 = new TProfile("pResVBySensorDisk2", "Disk 2 mean residual 1 by sensor;fstGlobalSensor;<resUnbiased1> [um]", 36, 71.5, 107.5);
    tree->Draw("resUnbiased0*10000:fstGlobalSensor>>pResUBySensorDisk0", fstResidualCut && "fstDisk==0", "goff");
    tree->Draw("resUnbiased0*10000:fstGlobalSensor>>pResUBySensorDisk1", fstResidualCut && "fstDisk==1", "goff");
    tree->Draw("resUnbiased0*10000:fstGlobalSensor>>pResUBySensorDisk2", fstResidualCut && "fstDisk==2", "goff");
    tree->Draw("resUnbiased1*10000:fstGlobalSensor>>pResVBySensorDisk0", fstResidualCut && "fstDisk==0", "goff");
    tree->Draw("resUnbiased1*10000:fstGlobalSensor>>pResVBySensorDisk1", fstResidualCut && "fstDisk==1", "goff");
    tree->Draw("resUnbiased1*10000:fstGlobalSensor>>pResVBySensorDisk2", fstResidualCut && "fstDisk==2", "goff");

    canvas->Clear();
    canvas->Divide(3, 2);
    canvas->cd(1); pResUBySensorDisk0->Draw("E1"); fwdAlignDrawLine(-0.5, 0, 35.5, 0);
    canvas->cd(2); pResUBySensorDisk1->Draw("E1"); fwdAlignDrawLine(35.5, 0, 71.5, 0);
    canvas->cd(3); pResUBySensorDisk2->Draw("E1"); fwdAlignDrawLine(71.5, 0, 107.5, 0);
    canvas->cd(4); pResVBySensorDisk0->Draw("E1"); fwdAlignDrawLine(-0.5, 0, 35.5, 0);
    canvas->cd(5); pResVBySensorDisk1->Draw("E1"); fwdAlignDrawLine(35.5, 0, 71.5, 0);
    canvas->cd(6); pResVBySensorDisk2->Draw("E1"); fwdAlignDrawLine(71.5, 0, 107.5, 0);
    fwdAlignSavePage(canvas, pdfOutput);

    if (hasPullBranches) {
        TProfile *pPullUBySensorDisk0 = new TProfile("pPullUBySensorDisk0", "Disk 0 mean pull 0 by sensor;fstGlobalSensor;<pullUnbiased0>", 36, -0.5, 35.5);
        TProfile *pPullUBySensorDisk1 = new TProfile("pPullUBySensorDisk1", "Disk 1 mean pull 0 by sensor;fstGlobalSensor;<pullUnbiased0>", 36, 35.5, 71.5);
        TProfile *pPullUBySensorDisk2 = new TProfile("pPullUBySensorDisk2", "Disk 2 mean pull 0 by sensor;fstGlobalSensor;<pullUnbiased0>", 36, 71.5, 107.5);
        TProfile *pPullVBySensorDisk0 = new TProfile("pPullVBySensorDisk0", "Disk 0 mean pull 1 by sensor;fstGlobalSensor;<pullUnbiased1>", 36, -0.5, 35.5);
        TProfile *pPullVBySensorDisk1 = new TProfile("pPullVBySensorDisk1", "Disk 1 mean pull 1 by sensor;fstGlobalSensor;<pullUnbiased1>", 36, 35.5, 71.5);
        TProfile *pPullVBySensorDisk2 = new TProfile("pPullVBySensorDisk2", "Disk 2 mean pull 1 by sensor;fstGlobalSensor;<pullUnbiased1>", 36, 71.5, 107.5);
        tree->Draw("pullUnbiased0:fstGlobalSensor>>pPullUBySensorDisk0", fstPullCut && "fstDisk==0", "goff");
        tree->Draw("pullUnbiased0:fstGlobalSensor>>pPullUBySensorDisk1", fstPullCut && "fstDisk==1", "goff");
        tree->Draw("pullUnbiased0:fstGlobalSensor>>pPullUBySensorDisk2", fstPullCut && "fstDisk==2", "goff");
        tree->Draw("pullUnbiased1:fstGlobalSensor>>pPullVBySensorDisk0", fstPullCut && "fstDisk==0", "goff");
        tree->Draw("pullUnbiased1:fstGlobalSensor>>pPullVBySensorDisk1", fstPullCut && "fstDisk==1", "goff");
        tree->Draw("pullUnbiased1:fstGlobalSensor>>pPullVBySensorDisk2", fstPullCut && "fstDisk==2", "goff");

        canvas->Clear();
        canvas->Divide(3, 2);
        canvas->cd(1); pPullUBySensorDisk0->Draw("E1"); fwdAlignDrawLine(-0.5, 0, 35.5, 0);
        canvas->cd(2); pPullUBySensorDisk1->Draw("E1"); fwdAlignDrawLine(35.5, 0, 71.5, 0);
        canvas->cd(3); pPullUBySensorDisk2->Draw("E1"); fwdAlignDrawLine(71.5, 0, 107.5, 0);
        canvas->cd(4); pPullVBySensorDisk0->Draw("E1"); fwdAlignDrawLine(-0.5, 0, 35.5, 0);
        canvas->cd(5); pPullVBySensorDisk1->Draw("E1"); fwdAlignDrawLine(35.5, 0, 71.5, 0);
        canvas->cd(6); pPullVBySensorDisk2->Draw("E1"); fwdAlignDrawLine(71.5, 0, 107.5, 0);
        fwdAlignSavePage(canvas, pdfOutput);
    }

    outputFile->Write();
    canvas->Print(pdfOutput + "]");
    outputFile->Close();
    inputFile->Close();

    std::cout << "Wrote " << rootOutput.Data() << std::endl;
    std::cout << "Wrote " << pdfOutput.Data() << std::endl;
}
