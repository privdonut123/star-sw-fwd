// Inspect FST FTUS geometry nodes used by Forward tracking.
//
// Usage:
//   root -l -b -q 'inspect_fst_ftus_geometry.C("fGeom.root")'
//   root -l -b -q 'inspect_fst_ftus_geometry.C("fGeom.root",0,0)'
//   root -l -b -q 'inspect_fst_ftus_geometry.C("fGeom.root",-1,-1,"fst_ftus_geometry.csv")'

#include "TFile.h"
#include "TGeoBBox.h"
#include "TGeoManager.h"
#include "TGeoMatrix.h"
#include "TGeoNode.h"
#include "TGeoShape.h"
#include "TGeoTube.h"
#include "TGeoVolume.h"
#include "TMath.h"
#include "TString.h"

#include <cmath>
#include <cstring>
#include <fstream>
#include <iostream>
#include <limits>

void inspect_fst_ftus_geometry(
    const char *geomFilename = "fGeom.root",
    int diskFilter = -1,
    int wedgeFilter = -1,
    const char *outputCsv = ""
) {
    static const int kNumDisks = 3;
    static const int kNumWedgesPerDisk = 12;
    static const int kNumSensorsPerWedge = 3;
    static const int kEventSensorToFtusCopy[kNumSensorsPerWedge] = {3, 1, 2};
    static const int kElecToGeantWedge[kNumWedgesPerDisk] = {2, 7, 1, 12, 6, 11, 5, 10, 4, 9, 3, 8};
    static const int kElecToGeantWedgeDisk2[kNumWedgesPerDisk] = {7, 1, 12, 6, 11, 5, 10, 4, 9, 3, 8, 2};

    TFile inputFile(geomFilename, "READ");
    if (inputFile.IsZombie()) {
        std::cerr << "Cannot open geometry file: " << geomFilename << std::endl;
        return;
    }

    TGeoManager *geom = dynamic_cast<TGeoManager*>(inputFile.Get("dyson"));
    if (!geom) {
        std::cerr << "Cannot find TGeoManager named dyson in " << geomFilename << std::endl;
        inputFile.ls();
        return;
    }

    std::ofstream csvFile;
    std::ostream *out = &std::cout;
    if (outputCsv && std::strlen(outputCsv) > 0) {
        csvFile.open(outputCsv);
        if (!csvFile.is_open()) {
            std::cerr << "Cannot open output CSV: " << outputCsv << std::endl;
            return;
        }
        out = &csvFile;
    }

    *out << "disk,wedge,eventSensor,geantFSTW,ftusCopy,"
         << "nodeX,nodeY,nodeZ,nodeR,nodePhiDeg,"
         << "activeX,activeY,activeZ,activeR,activePhiDeg,"
         << "uPhiDeg,vPhiDeg,crossUvZ,"
         << "shapeRMin,shapeRMax,shapePhi1Deg,shapePhi2Deg,path\n";

    for (int disk = 0; disk < kNumDisks; ++disk) {
        if (diskFilter >= 0 && disk != diskFilter)
            continue;

        const int fstdCopy = disk + 4;
        const int *wedgeMap = (fstdCopy == 5) ? kElecToGeantWedgeDisk2 : kElecToGeantWedge;

        for (int wedge = 0; wedge < kNumWedgesPerDisk; ++wedge) {
            if (wedgeFilter >= 0 && wedge != wedgeFilter)
                continue;

            const int geantFstwCopy = wedgeMap[wedge];

            for (int eventSensor = 0; eventSensor < kNumSensorsPerWedge; ++eventSensor) {
                const int ftusCopy = kEventSensorToFtusCopy[eventSensor];
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
                if (!matrix || !node || !node->GetVolume()) {
                    std::cerr << "Incomplete geometry node at " << path.Data() << std::endl;
                    continue;
                }

                const Double_t *nodeTranslation = matrix->GetTranslation();
                const double nodeX = nodeTranslation[0];
                const double nodeY = nodeTranslation[1];
                const double nodeZ = nodeTranslation[2];
                const double nodeR = std::hypot(nodeX, nodeY);
                const double nodePhiDeg = std::atan2(nodeY, nodeX) * TMath::RadToDeg();
                const Double_t *rotation = matrix->GetRotationMatrix();
                const double ux = rotation[0];
                const double uy = rotation[3];
                double vx = rotation[1];
                double vy = rotation[4];
                const double crossUvZ = ux * vy - uy * vx;
                if (crossUvZ < 0.0) {
                    vx = -vx;
                    vy = -vy;
                }
                const double uPhiDeg = std::atan2(uy, ux) * TMath::RadToDeg();
                const double vPhiDeg = std::atan2(vy, vx) * TMath::RadToDeg();

                double shapeRMin = std::numeric_limits<double>::quiet_NaN();
                double shapeRMax = std::numeric_limits<double>::quiet_NaN();
                double shapePhi1Deg = std::numeric_limits<double>::quiet_NaN();
                double shapePhi2Deg = std::numeric_limits<double>::quiet_NaN();
                Double_t activeLocal[3] = {0.0, 0.0, 0.0};
                bool haveActiveCenter = false;

                TGeoShape *shape = node->GetVolume()->GetShape();
                if (auto tube = dynamic_cast<TGeoTubeSeg*>(shape)) {
                    shapeRMin = tube->GetRmin();
                    shapeRMax = tube->GetRmax();
                    shapePhi1Deg = tube->GetPhi1();
                    shapePhi2Deg = tube->GetPhi2();
                    const double rCenter = 0.5 * (shapeRMin + shapeRMax);
                    const double phiCenter = 0.5 * (shapePhi1Deg + shapePhi2Deg) * TMath::DegToRad();
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

                Double_t activeMaster[3] = {
                    std::numeric_limits<double>::quiet_NaN(),
                    std::numeric_limits<double>::quiet_NaN(),
                    std::numeric_limits<double>::quiet_NaN()
                };
                if (haveActiveCenter)
                    matrix->LocalToMaster(activeLocal, activeMaster);

                const double activeR = std::hypot(activeMaster[0], activeMaster[1]);
                const double activePhiDeg = std::atan2(activeMaster[1], activeMaster[0]) * TMath::RadToDeg();

                *out << disk << ","
                     << wedge << ","
                     << eventSensor << ","
                     << geantFstwCopy << ","
                     << ftusCopy << ","
                     << nodeX << ","
                     << nodeY << ","
                     << nodeZ << ","
                     << nodeR << ","
                     << nodePhiDeg << ","
                     << activeMaster[0] << ","
                     << activeMaster[1] << ","
                         << activeMaster[2] << ","
                         << activeR << ","
                         << activePhiDeg << ","
                         << uPhiDeg << ","
                         << vPhiDeg << ","
                         << crossUvZ << ","
                         << shapeRMin << ","
                         << shapeRMax << ","
                     << shapePhi1Deg << ","
                     << shapePhi2Deg << ","
                     << path.Data()
                     << "\n";
            }
        }
    }

    if (csvFile.is_open())
        std::cout << "Wrote " << outputCsv << std::endl;
}
