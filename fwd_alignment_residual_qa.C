// ROOT macro for first-pass Forward alignment residual QA.
//
// Usage:
//   root4star -l -b -q
//   'fwd_alignment_residual_qa.C("align_test.root","fwd_align_qa")'
//
// Outputs:
//   fwd_align_qa.root
//   fwd_align_qa.pdf

#include "TCanvas.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLine.h"
#include "TPaveText.h"
#include "TProfile.h"
#include "TROOT.h"
#include "TString.h"
#include "TStyle.h"
#include "TTree.h"

#include <cmath>
#include <iostream>

///////////////////////////////////////////////////////////////////////////////
// USER EDIT SECTION
//
// If you are new to this macro, start here. These are the most common knobs.
// The rest of the file just books histograms, loops over fwdAlign, fills plots,
// and draws the PDF pages.

// Residual/pull track quality cuts.
const double kFwdQaTrackEtaMin = 2.5;
const double kFwdQaTrackEtaMax = 4.0;
const double kFwdQaTrackPMin = 1.0;  // GeV/c
const double kFwdQaTrackPtMin = 0.2; // GeV/c
const int kFwdQaDefaultMinTrackNFstHits = 3;

// Local measurement map binning.
const int kFwdAlignLocalMeas0Bins = 12;
const double kFwdAlignLocalMeas0Min = -6.0;
const double kFwdAlignLocalMeas0Max = 6.0;
const int kFwdAlignLocalMeas1Bins = 24;
const double kFwdAlignLocalMeas1Min = -12.0;
const double kFwdAlignLocalMeas1Max = 12.0;

// Raw FST and track-slope diagnostic binning.
const int kFwdQaRawRBins = 24;
const double kFwdQaRawRMin = 4.0;
const double kFwdQaRawRMax = 28.0;
const int kFwdQaEtaBins = 10;
const int kFwdQaPolarSlopeBins = 10;
const double kFwdQaPolarSlopeMin = 0.0;
const double kFwdQaPolarSlopeMax = 0.2;

// FST measurement closure diagnostic binning.
// Closure means:
//   fstClosure = (fstPlaneOrigin + meas0*U + meas1*V) - fstHitGlobal
const int kFwdQaClosureBins = 160;
const double kFwdQaCmToMicron = 10000.0;
const double kFwdQaClosureSignedMinMicron = -10000.0;
const double kFwdQaClosureSignedMaxMicron = 10000.0;
const double kFwdQaClosureMagMaxMicron = 20000.0;

///////////////////////////////////////////////////////////////////////////////
// Fixed FST constants. Usually do not edit these.

const int kFwdAlignNumFstDisks = 3;
const int kFwdAlignNumFstWedges = 12;
const int kFwdAlignNumFstSensorsPerWedge = 3;
const int kFwdAlignNumFstSensors = 108;
const int kFwdAlignFstDetId = 45;
const float kFwdAlignInvalid = -90000.f;

///////////////////////////////////////////////////////////////////////////////
// Small helper functions.

bool fwdAlignHasBranch(TTree *tree, const char *name) {
  return tree && tree->GetBranch(name);
}

bool fwdAlignBindBranch(TTree *tree, const char *name, void *address,
                        bool required = true) {
  if (!fwdAlignHasBranch(tree, name)) {
    if (required)
      std::cerr << "Missing required fwdAlign branch: " << name << std::endl;
    return false;
  }
  tree->SetBranchAddress(name, address);
  return true;
}

bool fwdAlignValid(float value) { return value > kFwdAlignInvalid; }

bool fwdAlignValidFstDisk(int disk) {
  return disk >= 0 && disk < kFwdAlignNumFstDisks;
}

bool fwdAlignValidFstWedge(int wedge) {
  return wedge >= 0 && wedge < kFwdAlignNumFstWedges;
}

bool fwdAlignValidFstSensor(int sensor) {
  return sensor >= 0 && sensor < kFwdAlignNumFstSensorsPerWedge;
}

bool fwdAlignValidFstGlobalSensor(int globalSensor) {
  return globalSensor >= 0 && globalSensor < kFwdAlignNumFstSensors;
}

void fwdAlignDrawLine(double x1, double y1, double x2, double y2,
                      int color = kRed + 1, int style = 2) {
  TLine *line = new TLine(x1, y1, x2, y2);
  line->SetLineColor(color);
  line->SetLineStyle(style);
  line->Draw();
}

void fwdAlignSavePage(TCanvas *canvas, const TString &pdf) {
  canvas->Print(pdf);
}

TH2D *fwdAlignMakeMeasurementMap(const TString &name, const TString &title) {
  TString fullTitle = title + ";meas0 [cm];meas1 [cm];rows";
  TH2D *hist = new TH2D(name.Data(), fullTitle.Data(), kFwdAlignLocalMeas0Bins,
                        kFwdAlignLocalMeas0Min, kFwdAlignLocalMeas0Max,
                        kFwdAlignLocalMeas1Bins, kFwdAlignLocalMeas1Min,
                        kFwdAlignLocalMeas1Max);
  hist->SetStats(false);
  return hist;
}

///////////////////////////////////////////////////////////////////////////////
// Branch bundle.
//
// This struct is just a container for one row of fwdAlign. The bind() method
// connects each member to the matching TTree branch once at the start.

struct FwdAlignBranches {
  int run;
  int event;
  int trackIndex;
  int pointIndex;
  int measurementIndex;
  int detId;
  int hitId;
  int fstGlobalSensor;
  int fstDisk;
  int fstWedge;
  int fstSensor;
  int measurementDim;
  int residualDim;
  int hasResidual;
  int nSeeds;
  int nFitTracks;
  int ndf;
  int fitConverged;
  int fitConvergedFully;
  int fitConvergedPartially;
  int trackNHitsFit;
  int trackNFstHits;

  float chi2;
  float pval;
  float trackPx;
  float trackPy;
  float trackPz;
  float trackP;
  float trackPt;
  float trackEta;
  float sorting;
  float meas0;
  float meas1;
  float meas2;
  float trackPred0;
  float trackPred1;
  float trackPred2;
  float fstRawR;
  float fstRawStripPhi;
  float fstMeanPhiStrip;
  float fstHitGlobalX;
  float fstHitGlobalY;
  float fstHitGlobalZ;
  float fstMeasGlobalX;
  float fstMeasGlobalY;
  float fstMeasGlobalZ;
  float fstClosureX;
  float fstClosureY;
  float fstClosureZ;
  float fstClosureU;
  float fstClosureV;
  float fstClosureMag;
  float resBiased0;
  float resBiased1;
  float resBiased2;
  float resBiasedSigma0;
  float resBiasedSigma1;
  float resBiasedSigma2;
  float pullBiased0;
  float pullBiased1;
  float pullBiased2;
  float resUnbiased0;
  float resUnbiased1;
  float resUnbiased2;
  float resUnbiasedSigma0;
  float resUnbiasedSigma1;
  float resUnbiasedSigma2;
  float pullUnbiased0;
  float pullUnbiased1;
  float pullUnbiased2;

  bool hasPullBranches;
  bool hasRawPredictionBranches;
  bool hasFstClosureBranches;
  bool hasTrackSlopeBranches;
  bool hasTrackNFstHitsBranch;

  FwdAlignBranches() {
    run = 0;
    event = 0;
    trackIndex = 0;
    pointIndex = 0;
    measurementIndex = 0;
    detId = 0;
    hitId = 0;
    fstGlobalSensor = -1;
    fstDisk = -1;
    fstWedge = -1;
    fstSensor = -1;
    measurementDim = 0;
    residualDim = 0;
    hasResidual = 0;
    nSeeds = 0;
    nFitTracks = 0;
    ndf = 0;
    fitConverged = 0;
    fitConvergedFully = 0;
    fitConvergedPartially = 0;
    trackNHitsFit = 0;
    trackNFstHits = 0;
    chi2 = 0;
    pval = 0;
    trackPx = 0;
    trackPy = 0;
    trackPz = 0;
    trackP = 0;
    trackPt = 0;
    trackEta = 0;
    sorting = 0;
    meas0 = 0;
    meas1 = 0;
    meas2 = 0;
    trackPred0 = kFwdAlignInvalid;
    trackPred1 = kFwdAlignInvalid;
    trackPred2 = kFwdAlignInvalid;
    fstRawR = kFwdAlignInvalid;
    fstRawStripPhi = kFwdAlignInvalid;
    fstMeanPhiStrip = kFwdAlignInvalid;
    fstHitGlobalX = kFwdAlignInvalid;
    fstHitGlobalY = kFwdAlignInvalid;
    fstHitGlobalZ = kFwdAlignInvalid;
    fstMeasGlobalX = kFwdAlignInvalid;
    fstMeasGlobalY = kFwdAlignInvalid;
    fstMeasGlobalZ = kFwdAlignInvalid;
    fstClosureX = kFwdAlignInvalid;
    fstClosureY = kFwdAlignInvalid;
    fstClosureZ = kFwdAlignInvalid;
    fstClosureU = kFwdAlignInvalid;
    fstClosureV = kFwdAlignInvalid;
    fstClosureMag = kFwdAlignInvalid;
    resBiased0 = kFwdAlignInvalid;
    resBiased1 = kFwdAlignInvalid;
    resBiased2 = kFwdAlignInvalid;
    resBiasedSigma0 = kFwdAlignInvalid;
    resBiasedSigma1 = kFwdAlignInvalid;
    resBiasedSigma2 = kFwdAlignInvalid;
    pullBiased0 = kFwdAlignInvalid;
    pullBiased1 = kFwdAlignInvalid;
    pullBiased2 = kFwdAlignInvalid;
    resUnbiased0 = kFwdAlignInvalid;
    resUnbiased1 = kFwdAlignInvalid;
    resUnbiased2 = kFwdAlignInvalid;
    resUnbiasedSigma0 = kFwdAlignInvalid;
    resUnbiasedSigma1 = kFwdAlignInvalid;
    resUnbiasedSigma2 = kFwdAlignInvalid;
    pullUnbiased0 = kFwdAlignInvalid;
    pullUnbiased1 = kFwdAlignInvalid;
    pullUnbiased2 = kFwdAlignInvalid;
    hasPullBranches = false;
    hasRawPredictionBranches = false;
    hasFstClosureBranches = false;
    hasTrackSlopeBranches = false;
    hasTrackNFstHitsBranch = false;
  }

  bool bind(TTree *tree) {
    bool ok = true;
    ok &= fwdAlignBindBranch(tree, "run", &run);
    ok &= fwdAlignBindBranch(tree, "event", &event);
    ok &= fwdAlignBindBranch(tree, "trackIndex", &trackIndex);
    ok &= fwdAlignBindBranch(tree, "pointIndex", &pointIndex);
    fwdAlignBindBranch(tree, "measurementIndex", &measurementIndex, false);
    ok &= fwdAlignBindBranch(tree, "detId", &detId);
    ok &= fwdAlignBindBranch(tree, "hitId", &hitId);
    ok &= fwdAlignBindBranch(tree, "fstGlobalSensor", &fstGlobalSensor);
    ok &= fwdAlignBindBranch(tree, "fstDisk", &fstDisk);
    ok &= fwdAlignBindBranch(tree, "fstWedge", &fstWedge);
    ok &= fwdAlignBindBranch(tree, "fstSensor", &fstSensor);
    ok &= fwdAlignBindBranch(tree, "measurementDim", &measurementDim);
    ok &= fwdAlignBindBranch(tree, "residualDim", &residualDim);
    ok &= fwdAlignBindBranch(tree, "hasResidual", &hasResidual);
    fwdAlignBindBranch(tree, "nSeeds", &nSeeds, false);
    fwdAlignBindBranch(tree, "nFitTracks", &nFitTracks, false);
    fwdAlignBindBranch(tree, "chi2", &chi2, false);
    fwdAlignBindBranch(tree, "ndf", &ndf, false);
    fwdAlignBindBranch(tree, "pval", &pval, false);
    ok &= fwdAlignBindBranch(tree, "fitConverged", &fitConverged);
    ok &= fwdAlignBindBranch(tree, "fitConvergedFully", &fitConvergedFully);
    fwdAlignBindBranch(tree, "fitConvergedPartially", &fitConvergedPartially,
                       false);
    fwdAlignBindBranch(tree, "trackNHitsFit", &trackNHitsFit, false);

    hasTrackNFstHitsBranch = fwdAlignHasBranch(tree, "trackNFstHits");
    if (hasTrackNFstHitsBranch)
      fwdAlignBindBranch(tree, "trackNFstHits", &trackNFstHits, false);

    fwdAlignBindBranch(tree, "trackPx", &trackPx, false);
    fwdAlignBindBranch(tree, "trackPy", &trackPy, false);
    fwdAlignBindBranch(tree, "trackPz", &trackPz, false);
    ok &= fwdAlignBindBranch(tree, "trackP", &trackP);
    fwdAlignBindBranch(tree, "trackPt", &trackPt, false);
    ok &= fwdAlignBindBranch(tree, "trackEta", &trackEta);
    ok &= fwdAlignBindBranch(tree, "sorting", &sorting);
    ok &= fwdAlignBindBranch(tree, "meas0", &meas0);
    ok &= fwdAlignBindBranch(tree, "meas1", &meas1);
    fwdAlignBindBranch(tree, "meas2", &meas2, false);
    ok &= fwdAlignBindBranch(tree, "resBiased0", &resBiased0);
    ok &= fwdAlignBindBranch(tree, "resBiased1", &resBiased1);
    fwdAlignBindBranch(tree, "resBiased2", &resBiased2, false);
    ok &= fwdAlignBindBranch(tree, "resUnbiased0", &resUnbiased0);
    ok &= fwdAlignBindBranch(tree, "resUnbiased1", &resUnbiased1);
    fwdAlignBindBranch(tree, "resUnbiased2", &resUnbiased2, false);

    hasPullBranches = fwdAlignHasBranch(tree, "pullBiased0") &&
                      fwdAlignHasBranch(tree, "pullBiased1") &&
                      fwdAlignHasBranch(tree, "pullUnbiased0") &&
                      fwdAlignHasBranch(tree, "pullUnbiased1") &&
                      fwdAlignHasBranch(tree, "resBiasedSigma0") &&
                      fwdAlignHasBranch(tree, "resBiasedSigma1") &&
                      fwdAlignHasBranch(tree, "resUnbiasedSigma0") &&
                      fwdAlignHasBranch(tree, "resUnbiasedSigma1");
    if (hasPullBranches) {
      fwdAlignBindBranch(tree, "resBiasedSigma0", &resBiasedSigma0, false);
      fwdAlignBindBranch(tree, "resBiasedSigma1", &resBiasedSigma1, false);
      fwdAlignBindBranch(tree, "resBiasedSigma2", &resBiasedSigma2, false);
      fwdAlignBindBranch(tree, "pullBiased0", &pullBiased0, false);
      fwdAlignBindBranch(tree, "pullBiased1", &pullBiased1, false);
      fwdAlignBindBranch(tree, "pullBiased2", &pullBiased2, false);
      fwdAlignBindBranch(tree, "resUnbiasedSigma0", &resUnbiasedSigma0, false);
      fwdAlignBindBranch(tree, "resUnbiasedSigma1", &resUnbiasedSigma1, false);
      fwdAlignBindBranch(tree, "resUnbiasedSigma2", &resUnbiasedSigma2, false);
      fwdAlignBindBranch(tree, "pullUnbiased0", &pullUnbiased0, false);
      fwdAlignBindBranch(tree, "pullUnbiased1", &pullUnbiased1, false);
      fwdAlignBindBranch(tree, "pullUnbiased2", &pullUnbiased2, false);
    }

    hasRawPredictionBranches = fwdAlignHasBranch(tree, "fstRawR") &&
                               fwdAlignHasBranch(tree, "fstMeanPhiStrip") &&
                               fwdAlignHasBranch(tree, "trackPred0") &&
                               fwdAlignHasBranch(tree, "trackPred1");
    if (hasRawPredictionBranches) {
      fwdAlignBindBranch(tree, "trackPred0", &trackPred0, false);
      fwdAlignBindBranch(tree, "trackPred1", &trackPred1, false);
      fwdAlignBindBranch(tree, "trackPred2", &trackPred2, false);
      fwdAlignBindBranch(tree, "fstRawR", &fstRawR, false);
      fwdAlignBindBranch(tree, "fstRawStripPhi", &fstRawStripPhi, false);
      fwdAlignBindBranch(tree, "fstMeanPhiStrip", &fstMeanPhiStrip, false);
    }

    hasFstClosureBranches = fwdAlignHasBranch(tree, "fstHitGlobalX") &&
                            fwdAlignHasBranch(tree, "fstHitGlobalY") &&
                            fwdAlignHasBranch(tree, "fstHitGlobalZ") &&
                            fwdAlignHasBranch(tree, "fstMeasGlobalX") &&
                            fwdAlignHasBranch(tree, "fstMeasGlobalY") &&
                            fwdAlignHasBranch(tree, "fstMeasGlobalZ") &&
                            fwdAlignHasBranch(tree, "fstClosureX") &&
                            fwdAlignHasBranch(tree, "fstClosureY") &&
                            fwdAlignHasBranch(tree, "fstClosureZ") &&
                            fwdAlignHasBranch(tree, "fstClosureU") &&
                            fwdAlignHasBranch(tree, "fstClosureV") &&
                            fwdAlignHasBranch(tree, "fstClosureMag");
    if (hasFstClosureBranches) {
      fwdAlignBindBranch(tree, "fstHitGlobalX", &fstHitGlobalX, false);
      fwdAlignBindBranch(tree, "fstHitGlobalY", &fstHitGlobalY, false);
      fwdAlignBindBranch(tree, "fstHitGlobalZ", &fstHitGlobalZ, false);
      fwdAlignBindBranch(tree, "fstMeasGlobalX", &fstMeasGlobalX, false);
      fwdAlignBindBranch(tree, "fstMeasGlobalY", &fstMeasGlobalY, false);
      fwdAlignBindBranch(tree, "fstMeasGlobalZ", &fstMeasGlobalZ, false);
      fwdAlignBindBranch(tree, "fstClosureX", &fstClosureX, false);
      fwdAlignBindBranch(tree, "fstClosureY", &fstClosureY, false);
      fwdAlignBindBranch(tree, "fstClosureZ", &fstClosureZ, false);
      fwdAlignBindBranch(tree, "fstClosureU", &fstClosureU, false);
      fwdAlignBindBranch(tree, "fstClosureV", &fstClosureV, false);
      fwdAlignBindBranch(tree, "fstClosureMag", &fstClosureMag, false);
    }

    hasTrackSlopeBranches = fwdAlignHasBranch(tree, "trackPt") &&
                            fwdAlignHasBranch(tree, "trackPz");
    return ok;
  }
};

void fwd_alignment_residual_qa(
    const char *inputFilename = "align_test.root",
    const char *outputPrefix = "fwd_align_qa",
    bool requireFullyConverged = false, double residualRangeMicron = 5000.0,
    int minTrackNFstHits = kFwdQaDefaultMinTrackNFstHits) {
  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(1110);
  gStyle->SetOptFit(0);
  gStyle->SetPalette(1);

  TFile *inputFile = TFile::Open(inputFilename, "READ");
  if (!inputFile || inputFile->IsZombie()) {
    std::cerr << "Cannot open input file: " << inputFilename << std::endl;
    return;
  }

  TTree *tree = dynamic_cast<TTree *>(inputFile->Get("fwdAlign"));
  if (!tree) {
    std::cerr << "Cannot find TTree fwdAlign in " << inputFilename << std::endl;
    inputFile->ls();
    return;
  }

  FwdAlignBranches row;
  if (!row.bind(tree)) {
    std::cerr
        << "Stopping because the alignment tree schema is not the expected one."
        << std::endl;
    return;
  }

  TString rootOutput = TString::Format("%s.root", outputPrefix);
  TString pdfOutput = TString::Format("%s.pdf", outputPrefix);
  TFile *outputFile = TFile::Open(rootOutput, "RECREATE");
  if (!outputFile || outputFile->IsZombie()) {
    std::cerr << "Cannot create output file: " << rootOutput.Data()
              << std::endl;
    return;
  }
  outputFile->cd();

  // 1. Build the cut strings printed in the terminal and summary page.
  TString fitCut =
      requireFullyConverged ? "fitConvergedFully>0" : "fitConverged>0";
  TString trackCut =
      TString::Format("trackEta>=%.1f&&trackEta<=%.1f&&trackP>=%.1f",
                      kFwdQaTrackEtaMin, kFwdQaTrackEtaMax, kFwdQaTrackPMin);
  if (minTrackNFstHits > 0 && row.hasTrackNFstHitsBranch)
    trackCut += TString::Format("&&trackNFstHits>=%d", minTrackNFstHits);

  // 2. Book histograms. The event loop below is the only place they are filled.
  TH1D *hDetId =
      new TH1D("hDetId", "Alignment rows by detector;detId;rows", 80, 0, 80);
  TH2D *hDim = new TH2D(
      "hDim",
      "FST measurement and residual dimensions;measurementDim;residualDim", 6,
      -0.5, 5.5, 6, -0.5, 5.5);
  TH1D *hHasResidual = new TH1D(
      "hHasResidual", "FST hasResidual flag;hasResidual;rows", 3, -0.5, 2.5);
  TH1D *hSortMinusSensor = new TH1D(
      "hSortMinusSensor",
      "FST sorting - fstGlobalSensor, expect 1;sorting - fstGlobalSensor;rows",
      21, -9.5, 11.5);

  TH2D *hMeas1VsMeas0 = fwdAlignMakeMeasurementMap(
      "hMeas1VsMeas0", "FST local measurement map, all sensors");
  TH2D *hMeas1VsMeas0Disk[kFwdAlignNumFstDisks] = {0};
  TH2D *hMeas1VsMeas0Wedge[kFwdAlignNumFstWedges] = {0};
  TH2D *hMeas1VsMeas0Sensor[kFwdAlignNumFstSensorsPerWedge] = {0};
  TH2D *hMeas1VsMeas0GlobalSensor[kFwdAlignNumFstSensors] = {0};
  for (int disk = 0; disk < kFwdAlignNumFstDisks; ++disk) {
    hMeas1VsMeas0Disk[disk] = fwdAlignMakeMeasurementMap(
        TString::Format("hMeas1VsMeas0Disk%d", disk),
        TString::Format("FST local measurement map, disk %d", disk));
  }
  for (int wedge = 0; wedge < kFwdAlignNumFstWedges; ++wedge) {
    hMeas1VsMeas0Wedge[wedge] = fwdAlignMakeMeasurementMap(
        TString::Format("hMeas1VsMeas0Wedge%d", wedge),
        TString::Format("FST local measurement map, wedge %d", wedge));
  }
  for (int sensor = 0; sensor < kFwdAlignNumFstSensorsPerWedge; ++sensor) {
    hMeas1VsMeas0Sensor[sensor] = fwdAlignMakeMeasurementMap(
        TString::Format("hMeas1VsMeas0Sensor%d", sensor),
        TString::Format("FST local measurement map, sensor-in-wedge %d",
                        sensor));
  }
  for (int globalSensor = 0; globalSensor < kFwdAlignNumFstSensors;
       ++globalSensor) {
    int disk = globalSensor / 36;
    int wedge = (globalSensor / 3) % 12;
    int sensor = globalSensor % 3;
    hMeas1VsMeas0GlobalSensor[globalSensor] = fwdAlignMakeMeasurementMap(
        TString::Format("hMeas1VsMeas0GlobalSensor%d", globalSensor),
        TString::Format(
            "FST local measurement map, global sensor %d (d%d w%d s%d)",
            globalSensor, disk, wedge, sensor));
  }

  TH1D *hResU =
      new TH1D("hResU", "FST unbiased residual 0;resUnbiased0 [um];rows", 160,
               -residualRangeMicron, residualRangeMicron);
  TH1D *hResV =
      new TH1D("hResV", "FST unbiased residual 1;resUnbiased1 [um];rows", 160,
               -residualRangeMicron, residualRangeMicron);
  TH1D *hBiasedU =
      new TH1D("hBiasedU", "FST biased residual 0;resBiased0 [um];rows", 160,
               -residualRangeMicron, residualRangeMicron);
  TH1D *hBiasedV =
      new TH1D("hBiasedV", "FST biased residual 1;resBiased1 [um];rows", 160,
               -residualRangeMicron, residualRangeMicron);

  TH1D *hPullU = 0;
  TH1D *hPullV = 0;
  TH1D *hPullBiasedU = 0;
  TH1D *hPullBiasedV = 0;
  TH1D *hSigmaU = 0;
  TH1D *hSigmaV = 0;
  TH1D *hSigmaBiasedU = 0;
  TH1D *hSigmaBiasedV = 0;
  if (row.hasPullBranches) {
    hPullU = new TH1D("hPullU", "FST unbiased pull 0;pullUnbiased0;rows", 160,
                      -10, 10);
    hPullV = new TH1D("hPullV", "FST unbiased pull 1;pullUnbiased1;rows", 160,
                      -10, 10);
    hPullBiasedU = new TH1D("hPullBiasedU",
                            "FST biased pull 0;pullBiased0;rows", 160, -10, 10);
    hPullBiasedV = new TH1D("hPullBiasedV",
                            "FST biased pull 1;pullBiased1;rows", 160, -10, 10);
    hSigmaU = new TH1D("hSigmaU",
                       "FST unbiased residual sigma 0;#sigma_{res,0} [um];rows",
                       160, 0, residualRangeMicron);
    hSigmaV = new TH1D("hSigmaV",
                       "FST unbiased residual sigma 1;#sigma_{res,1} [um];rows",
                       160, 0, residualRangeMicron);
    hSigmaBiasedU = new TH1D(
        "hSigmaBiasedU", "FST biased residual sigma 0;#sigma_{res,0} [um];rows",
        160, 0, residualRangeMicron);
    hSigmaBiasedV = new TH1D(
        "hSigmaBiasedV", "FST biased residual sigma 1;#sigma_{res,1} [um];rows",
        160, 0, residualRangeMicron);
  }

  TH1D *hFstDisk = new TH1D(
      "hFstDisk", "FST residual rows by disk;fstDisk;rows", 3, -0.5, 2.5);
  TH1D *hFstSensor = new TH1D(
      "hFstSensor", "FST residual rows by global sensor;fstGlobalSensor;rows",
      108, -0.5, 107.5);
  TH2D *hWedgeDisk =
      new TH2D("hWedgeDisk", "FST residual occupancy;fstDisk;fstWedge", 3, -0.5,
               2.5, 12, -0.5, 11.5);
  TH2D *hSensorDisk = new TH2D(
      "hSensorDisk", "FST sensor-in-wedge residual occupancy;fstDisk;fstSensor",
      3, -0.5, 2.5, 3, -0.5, 2.5);

  TProfile *pResUByDisk = new TProfile(
      "pResUByDisk",
      "Mean FST unbiased residual 0 by disk;fstDisk;<resUnbiased0> [um]", 3,
      -0.5, 2.5);
  TProfile *pResVByDisk = new TProfile(
      "pResVByDisk",
      "Mean FST unbiased residual 1 by disk;fstDisk;<resUnbiased1> [um]", 3,
      -0.5, 2.5);
  TProfile *pResUByWedge = new TProfile(
      "pResUByWedge",
      "Mean FST unbiased residual 0 by wedge;fstWedge;<resUnbiased0> [um]", 12,
      -0.5, 11.5);
  TProfile *pResVByWedge = new TProfile(
      "pResVByWedge",
      "Mean FST unbiased residual 1 by wedge;fstWedge;<resUnbiased1> [um]", 12,
      -0.5, 11.5);
  TProfile *pResUByGlobalSensor =
      new TProfile("pResUByGlobalSensor",
                   "Mean FST unbiased residual 0 by global "
                   "sensor;fstGlobalSensor;<resUnbiased0> [um]",
                   108, -0.5, 107.5);
  TProfile *pResVByGlobalSensor =
      new TProfile("pResVByGlobalSensor",
                   "Mean FST unbiased residual 1 by global "
                   "sensor;fstGlobalSensor;<resUnbiased1> [um]",
                   108, -0.5, 107.5);
  TH2D *hResUByGlobalSensor = new TH2D(
      "hResUByGlobalSensor",
      "FST unbiased residual 0 by global sensor;fstGlobalSensor;resUnbiased0 "
      "[um]",
      108, -0.5, 107.5, 160, -residualRangeMicron, residualRangeMicron);
  TH2D *hResVByGlobalSensor = new TH2D(
      "hResVByGlobalSensor",
      "FST unbiased residual 1 by global sensor;fstGlobalSensor;resUnbiased1 "
      "[um]",
      108, -0.5, 107.5, 160, -residualRangeMicron, residualRangeMicron);

  TProfile *pPullUByDisk = 0;
  TProfile *pPullVByDisk = 0;
  TProfile *pPullUByWedge = 0;
  TProfile *pPullVByWedge = 0;
  TProfile *pPullUByGlobalSensor = 0;
  TProfile *pPullVByGlobalSensor = 0;
  TH2D *hPullUByGlobalSensor = 0;
  TH2D *hPullVByGlobalSensor = 0;
  if (row.hasPullBranches) {
    pPullUByDisk =
        new TProfile("pPullUByDisk",
                     "Mean FST unbiased pull 0 by disk;fstDisk;<pullUnbiased0>",
                     3, -0.5, 2.5);
    pPullVByDisk =
        new TProfile("pPullVByDisk",
                     "Mean FST unbiased pull 1 by disk;fstDisk;<pullUnbiased1>",
                     3, -0.5, 2.5);
    pPullUByWedge = new TProfile(
        "pPullUByWedge",
        "Mean FST unbiased pull 0 by wedge;fstWedge;<pullUnbiased0>", 12, -0.5,
        11.5);
    pPullVByWedge = new TProfile(
        "pPullVByWedge",
        "Mean FST unbiased pull 1 by wedge;fstWedge;<pullUnbiased1>", 12, -0.5,
        11.5);
    pPullUByGlobalSensor =
        new TProfile("pPullUByGlobalSensor",
                     "Mean FST unbiased pull 0 by global "
                     "sensor;fstGlobalSensor;<pullUnbiased0>",
                     108, -0.5, 107.5);
    pPullVByGlobalSensor =
        new TProfile("pPullVByGlobalSensor",
                     "Mean FST unbiased pull 1 by global "
                     "sensor;fstGlobalSensor;<pullUnbiased1>",
                     108, -0.5, 107.5);
    hPullUByGlobalSensor = new TH2D(
        "hPullUByGlobalSensor",
        "FST unbiased pull 0 by global sensor;fstGlobalSensor;pullUnbiased0",
        108, -0.5, 107.5, 160, -10, 10);
    hPullVByGlobalSensor = new TH2D(
        "hPullVByGlobalSensor",
        "FST unbiased pull 1 by global sensor;fstGlobalSensor;pullUnbiased1",
        108, -0.5, 107.5, 160, -10, 10);
  }

  TH2D *hResUVsMeas0 = new TH2D(
      "hResUVsMeas0",
      "FST residual 0 vs local meas0;meas0 [cm];resUnbiased0 [um]",
      kFwdAlignLocalMeas0Bins, kFwdAlignLocalMeas0Min, kFwdAlignLocalMeas0Max,
      160, -residualRangeMicron, residualRangeMicron);
  TH2D *hResUVsMeas1 = new TH2D(
      "hResUVsMeas1",
      "FST residual 0 vs local meas1;meas1 [cm];resUnbiased0 [um]",
      kFwdAlignLocalMeas1Bins, kFwdAlignLocalMeas1Min, kFwdAlignLocalMeas1Max,
      160, -residualRangeMicron, residualRangeMicron);
  TH2D *hResVVsMeas0 = new TH2D(
      "hResVVsMeas0",
      "FST residual 1 vs local meas0;meas0 [cm];resUnbiased1 [um]",
      kFwdAlignLocalMeas0Bins, kFwdAlignLocalMeas0Min, kFwdAlignLocalMeas0Max,
      160, -residualRangeMicron, residualRangeMicron);
  TH2D *hResVVsMeas1 = new TH2D(
      "hResVVsMeas1",
      "FST residual 1 vs local meas1;meas1 [cm];resUnbiased1 [um]",
      kFwdAlignLocalMeas1Bins, kFwdAlignLocalMeas1Min, kFwdAlignLocalMeas1Max,
      160, -residualRangeMicron, residualRangeMicron);

  TH2D *hPullUVsMeas0 = 0;
  TH2D *hPullUVsMeas1 = 0;
  TH2D *hPullVVsMeas0 = 0;
  TH2D *hPullVVsMeas1 = 0;
  if (row.hasPullBranches) {
    hPullUVsMeas0 = new TH2D(
        "hPullUVsMeas0", "FST pull 0 vs local meas0;meas0 [cm];pullUnbiased0",
        kFwdAlignLocalMeas0Bins, kFwdAlignLocalMeas0Min, kFwdAlignLocalMeas0Max,
        160, -10, 10);
    hPullUVsMeas1 = new TH2D(
        "hPullUVsMeas1", "FST pull 0 vs local meas1;meas1 [cm];pullUnbiased0",
        kFwdAlignLocalMeas1Bins, kFwdAlignLocalMeas1Min, kFwdAlignLocalMeas1Max,
        160, -10, 10);
    hPullVVsMeas0 = new TH2D(
        "hPullVVsMeas0", "FST pull 1 vs local meas0;meas0 [cm];pullUnbiased1",
        kFwdAlignLocalMeas0Bins, kFwdAlignLocalMeas0Min, kFwdAlignLocalMeas0Max,
        160, -10, 10);
    hPullVVsMeas1 = new TH2D(
        "hPullVVsMeas1", "FST pull 1 vs local meas1;meas1 [cm];pullUnbiased1",
        kFwdAlignLocalMeas1Bins, kFwdAlignLocalMeas1Min, kFwdAlignLocalMeas1Max,
        160, -10, 10);
  }

  TH2D *hRes0VsRawR[kFwdAlignNumFstSensorsPerWedge] = {0};
  TH2D *hRes0VsMeanPhi[kFwdAlignNumFstSensorsPerWedge] = {0};
  TH2D *hTrackPred0VsMeas0[kFwdAlignNumFstSensorsPerWedge] = {0};
  TH2D *hTrackPred1VsMeas1[kFwdAlignNumFstSensorsPerWedge] = {0};
  TProfile *pRes0VsRawR[kFwdAlignNumFstSensorsPerWedge] = {0};
  TProfile *pRes0VsMeanPhi[kFwdAlignNumFstSensorsPerWedge] = {0};
  TProfile *pRes0VsRawRDiskSensor[kFwdAlignNumFstDisks]
                                 [kFwdAlignNumFstSensorsPerWedge] = {{0}};
  if (row.hasRawPredictionBranches) {
    for (int sensor = 0; sensor < kFwdAlignNumFstSensorsPerWedge; ++sensor) {
      hRes0VsRawR[sensor] = new TH2D(
          TString::Format("hRes0VsRawR_sensor%d", sensor),
          TString::Format("FST sensor %d residual 0 vs raw R;fstRawR "
                          "[cm];resUnbiased0 [cm]",
                          sensor),
          kFwdQaRawRBins, kFwdQaRawRMin, kFwdQaRawRMax, 160, -20.0, 20.0);
      hRes0VsMeanPhi[sensor] =
          new TH2D(TString::Format("hRes0VsMeanPhi_sensor%d", sensor),
                   TString::Format("FST sensor %d residual 0 vs mean phi "
                                   "strip;fstMeanPhiStrip;resUnbiased0 [cm]",
                                   sensor),
                   128, -0.5, 127.5, 160, -20.0, 20.0);
      hTrackPred0VsMeas0[sensor] =
          new TH2D(TString::Format("hTrackPred0VsMeas0_sensor%d", sensor),
                   TString::Format("FST sensor %d track prediction 0 vs "
                                   "meas0;meas0 [cm];trackPred0 [cm]",
                                   sensor),
                   60, -10.0, 10.0, 120, -20.0, 20.0);
      hTrackPred1VsMeas1[sensor] =
          new TH2D(TString::Format("hTrackPred1VsMeas1_sensor%d", sensor),
                   TString::Format("FST sensor %d track prediction 1 vs "
                                   "meas1;meas1 [cm];trackPred1 [cm]",
                                   sensor),
                   60, -12.0, 12.0, 120, -12.0, 12.0);
      pRes0VsRawR[sensor] = new TProfile(
          TString::Format("pRes0VsRawR_sensor%d", sensor),
          TString::Format("Mean FST sensor %d residual 0 vs raw "
                          "R;fstRawR [cm];<resUnbiased0> [cm]",
                          sensor),
          kFwdQaRawRBins, kFwdQaRawRMin, kFwdQaRawRMax, -20.0, 20.0);
      pRes0VsMeanPhi[sensor] = new TProfile(
          TString::Format("pRes0VsMeanPhi_sensor%d", sensor),
          TString::Format("Mean FST sensor %d residual 0 vs mean phi "
                          "strip;fstMeanPhiStrip;<resUnbiased0> [cm]",
                          sensor),
          64, -0.5, 127.5, -20.0, 20.0);
      for (int disk = 0; disk < kFwdAlignNumFstDisks; ++disk) {
        pRes0VsRawRDiskSensor[disk][sensor] = new TProfile(
            TString::Format("pRes0VsRawR_disk%d_sensor%d", disk, sensor),
            TString::Format("Mean FST disk %d sensor %d residual 0 vs raw "
                            "R;fstRawR [cm];<resUnbiased0> [cm]",
                            disk, sensor),
            kFwdQaRawRBins, kFwdQaRawRMin, kFwdQaRawRMax, -20.0, 20.0);
      }
    }
  }

  TH2D *hRes0VsEta[kFwdAlignNumFstSensorsPerWedge] = {0};
  TH2D *hRes0VsPolarSlope[kFwdAlignNumFstSensorsPerWedge] = {0};
  TProfile *pRes0VsEta[kFwdAlignNumFstSensorsPerWedge] = {0};
  TProfile *pRes0VsPolarSlope[kFwdAlignNumFstSensorsPerWedge] = {0};
  TProfile *pRes0VsEtaDiskSensor[kFwdAlignNumFstDisks]
                                [kFwdAlignNumFstSensorsPerWedge] = {{0}};
  TProfile *pRes0VsPolarSlopeDiskSensor[kFwdAlignNumFstDisks]
                                       [kFwdAlignNumFstSensorsPerWedge] = {{0}};
  if (row.hasTrackSlopeBranches) {
    for (int sensor = 0; sensor < kFwdAlignNumFstSensorsPerWedge; ++sensor) {
      hRes0VsEta[sensor] =
          new TH2D(TString::Format("hRes0VsEta_sensor%d", sensor),
                   TString::Format("FST sensor %d residual 0 vs track "
                                   "eta;trackEta;resUnbiased0 [cm]",
                                   sensor),
                   kFwdQaEtaBins, kFwdQaTrackEtaMin, kFwdQaTrackEtaMax, 160,
                   -20.0, 20.0);
      hRes0VsPolarSlope[sensor] = new TH2D(
          TString::Format("hRes0VsPolarSlope_sensor%d", sensor),
          TString::Format("FST sensor %d residual 0 vs "
                          "trackPt/trackPz;trackPt/trackPz;resUnbiased0 [cm]",
                          sensor),
          kFwdQaPolarSlopeBins, kFwdQaPolarSlopeMin, kFwdQaPolarSlopeMax, 160,
          -20.0, 20.0);
      pRes0VsEta[sensor] = new TProfile(
          TString::Format("pRes0VsEta_sensor%d", sensor),
          TString::Format("Mean FST sensor %d residual 0 vs track "
                          "eta;trackEta;<resUnbiased0> [cm]",
                          sensor),
          kFwdQaEtaBins, kFwdQaTrackEtaMin, kFwdQaTrackEtaMax, -20.0, 20.0);
      pRes0VsPolarSlope[sensor] = new TProfile(
          TString::Format("pRes0VsPolarSlope_sensor%d", sensor),
          TString::Format("Mean FST sensor %d residual 0 vs "
                          "trackPt/trackPz;trackPt/trackPz;<resUnbiased0> [cm]",
                          sensor),
          kFwdQaPolarSlopeBins, kFwdQaPolarSlopeMin, kFwdQaPolarSlopeMax, -20.0,
          20.0);
      for (int disk = 0; disk < kFwdAlignNumFstDisks; ++disk) {
        pRes0VsEtaDiskSensor[disk][sensor] = new TProfile(
            TString::Format("pRes0VsEta_disk%d_sensor%d", disk, sensor),
            TString::Format("Mean FST disk %d sensor %d residual 0 vs track "
                            "eta;trackEta;<resUnbiased0> [cm]",
                            disk, sensor),
            kFwdQaEtaBins, kFwdQaTrackEtaMin, kFwdQaTrackEtaMax, -20.0, 20.0);
        pRes0VsPolarSlopeDiskSensor[disk][sensor] = new TProfile(
            TString::Format("pRes0VsPolarSlope_disk%d_sensor%d", disk, sensor),
            TString::Format("Mean FST disk %d sensor %d residual 0 vs "
                            "trackPt/trackPz;trackPt/trackPz;<resUnbiased0> "
                            "[cm]",
                            disk, sensor),
            kFwdQaPolarSlopeBins, kFwdQaPolarSlopeMin, kFwdQaPolarSlopeMax,
            -20.0, 20.0);
      }
    }
  }

  TH1D *hFstClosureX = 0;
  TH1D *hFstClosureY = 0;
  TH1D *hFstClosureZ = 0;
  TH1D *hFstClosureU = 0;
  TH1D *hFstClosureV = 0;
  TH1D *hFstClosureMag = 0;
  TH2D *hFstClosureUByGlobalSensor = 0;
  TH2D *hFstClosureVByGlobalSensor = 0;
  TProfile *pFstClosureUByGlobalSensor = 0;
  TProfile *pFstClosureVByGlobalSensor = 0;
  TProfile *pFstClosureMagByGlobalSensor = 0;
  TProfile *pFstClosureUVsRawR[kFwdAlignNumFstSensorsPerWedge] = {0};
  TProfile *pFstClosureVVsRawR[kFwdAlignNumFstSensorsPerWedge] = {0};
  TProfile *pFstClosureUVsMeanPhi[kFwdAlignNumFstSensorsPerWedge] = {0};
  TProfile *pFstClosureVVsMeanPhi[kFwdAlignNumFstSensorsPerWedge] = {0};
  if (row.hasFstClosureBranches) {
    hFstClosureX = new TH1D(
        "hFstClosureX",
        "FST closure global X;fstMeasGlobalX - fstHitGlobalX [um];rows",
        kFwdQaClosureBins, kFwdQaClosureSignedMinMicron,
        kFwdQaClosureSignedMaxMicron);
    hFstClosureY = new TH1D(
        "hFstClosureY",
        "FST closure global Y;fstMeasGlobalY - fstHitGlobalY [um];rows",
        kFwdQaClosureBins, kFwdQaClosureSignedMinMicron,
        kFwdQaClosureSignedMaxMicron);
    hFstClosureZ = new TH1D(
        "hFstClosureZ",
        "FST closure global Z;fstMeasGlobalZ - fstHitGlobalZ [um];rows",
        kFwdQaClosureBins, kFwdQaClosureSignedMinMicron,
        kFwdQaClosureSignedMaxMicron);
    hFstClosureU = new TH1D(
        "hFstClosureU",
        "FST closure local U;(meas global - hit global) dot U [um];rows",
        kFwdQaClosureBins, kFwdQaClosureSignedMinMicron,
        kFwdQaClosureSignedMaxMicron);
    hFstClosureV = new TH1D(
        "hFstClosureV",
        "FST closure local V;(meas global - hit global) dot V [um];rows",
        kFwdQaClosureBins, kFwdQaClosureSignedMinMicron,
        kFwdQaClosureSignedMaxMicron);
    hFstClosureMag =
        new TH1D("hFstClosureMag",
                 "FST closure magnitude;|meas global - hit global| [um];rows",
                 kFwdQaClosureBins, 0.0, kFwdQaClosureMagMaxMicron);
    hFstClosureUByGlobalSensor = new TH2D(
        "hFstClosureUByGlobalSensor",
        "FST closure U by global sensor;fstGlobalSensor;closureU [um]", 108,
        -0.5, 107.5, kFwdQaClosureBins, kFwdQaClosureSignedMinMicron,
        kFwdQaClosureSignedMaxMicron);
    hFstClosureVByGlobalSensor = new TH2D(
        "hFstClosureVByGlobalSensor",
        "FST closure V by global sensor;fstGlobalSensor;closureV [um]", 108,
        -0.5, 107.5, kFwdQaClosureBins, kFwdQaClosureSignedMinMicron,
        kFwdQaClosureSignedMaxMicron);
    pFstClosureUByGlobalSensor = new TProfile(
        "pFstClosureUByGlobalSensor",
        "Mean FST closure U by global sensor;fstGlobalSensor;<closureU> [um]",
        108, -0.5, 107.5, kFwdQaClosureSignedMinMicron,
        kFwdQaClosureSignedMaxMicron);
    pFstClosureVByGlobalSensor = new TProfile(
        "pFstClosureVByGlobalSensor",
        "Mean FST closure V by global sensor;fstGlobalSensor;<closureV> [um]",
        108, -0.5, 107.5, kFwdQaClosureSignedMinMicron,
        kFwdQaClosureSignedMaxMicron);
    pFstClosureMagByGlobalSensor = new TProfile(
        "pFstClosureMagByGlobalSensor",
        "Mean FST closure magnitude by global sensor;fstGlobalSensor;"
        "<closure magnitude> [um]",
        108, -0.5, 107.5, 0.0, kFwdQaClosureMagMaxMicron);

    for (int sensor = 0; sensor < kFwdAlignNumFstSensorsPerWedge; ++sensor) {
      pFstClosureUVsRawR[sensor] = new TProfile(
          TString::Format("pFstClosureUVsRawR_sensor%d", sensor),
          TString::Format("Mean FST sensor %d closure U vs raw R;fstRawR "
                          "[cm];<closureU> [um]",
                          sensor),
          kFwdQaRawRBins, kFwdQaRawRMin, kFwdQaRawRMax,
          kFwdQaClosureSignedMinMicron, kFwdQaClosureSignedMaxMicron);
      pFstClosureVVsRawR[sensor] = new TProfile(
          TString::Format("pFstClosureVVsRawR_sensor%d", sensor),
          TString::Format("Mean FST sensor %d closure V vs raw R;fstRawR "
                          "[cm];<closureV> [um]",
                          sensor),
          kFwdQaRawRBins, kFwdQaRawRMin, kFwdQaRawRMax,
          kFwdQaClosureSignedMinMicron, kFwdQaClosureSignedMaxMicron);
      pFstClosureUVsMeanPhi[sensor] = new TProfile(
          TString::Format("pFstClosureUVsMeanPhi_sensor%d", sensor),
          TString::Format("Mean FST sensor %d closure U vs mean phi "
                          "strip;fstMeanPhiStrip;<closureU> [um]",
                          sensor),
          64, -0.5, 127.5, kFwdQaClosureSignedMinMicron,
          kFwdQaClosureSignedMaxMicron);
      pFstClosureVVsMeanPhi[sensor] = new TProfile(
          TString::Format("pFstClosureVVsMeanPhi_sensor%d", sensor),
          TString::Format("Mean FST sensor %d closure V vs mean phi "
                          "strip;fstMeanPhiStrip;<closureV> [um]",
                          sensor),
          64, -0.5, 127.5, kFwdQaClosureSignedMinMicron,
          kFwdQaClosureSignedMaxMicron);
    }
  }

  TProfile *pResUBySensorDisk[kFwdAlignNumFstDisks] = {0};
  TProfile *pResVBySensorDisk[kFwdAlignNumFstDisks] = {0};
  TProfile *pPullUBySensorDisk[kFwdAlignNumFstDisks] = {0};
  TProfile *pPullVBySensorDisk[kFwdAlignNumFstDisks] = {0};
  for (int disk = 0; disk < kFwdAlignNumFstDisks; ++disk) {
    const double xMin = disk * 36 - 0.5;
    const double xMax = (disk + 1) * 36 - 0.5;
    pResUBySensorDisk[disk] = new TProfile(
        TString::Format("pResUBySensorDisk%d", disk),
        TString::Format("Disk %d mean residual 0 by sensor;fstGlobalSensor;"
                        "<resUnbiased0> [um]",
                        disk),
        36, xMin, xMax);
    pResVBySensorDisk[disk] = new TProfile(
        TString::Format("pResVBySensorDisk%d", disk),
        TString::Format("Disk %d mean residual 1 by sensor;fstGlobalSensor;"
                        "<resUnbiased1> [um]",
                        disk),
        36, xMin, xMax);
    if (row.hasPullBranches) {
      pPullUBySensorDisk[disk] = new TProfile(
          TString::Format("pPullUBySensorDisk%d", disk),
          TString::Format("Disk %d mean pull 0 by sensor;fstGlobalSensor;"
                          "<pullUnbiased0>",
                          disk),
          36, xMin, xMax);
      pPullVBySensorDisk[disk] = new TProfile(
          TString::Format("pPullVBySensorDisk%d", disk),
          TString::Format("Disk %d mean pull 1 by sensor;fstGlobalSensor;"
                          "<pullUnbiased1>",
                          disk),
          36, xMin, xMax);
    }
  }

  // 3. Main event loop. Every selection is a named boolean so it is easy to
  // inspect or temporarily loosen during debugging.
  Long64_t nEntries = tree->GetEntries();
  Long64_t nFstRows = 0;
  Long64_t nFstMeasurementRows = 0;
  Long64_t nFstClosureRows = 0;
  Long64_t nFstResidualRows = 0;

  for (Long64_t entry = 0; entry < nEntries; ++entry) {
    tree->GetEntry(entry);

    const bool isFst = row.detId == kFwdAlignFstDetId;
    const bool hasFstSensor = row.fstGlobalSensor >= 0;
    const bool validGlobalSensor =
        fwdAlignValidFstGlobalSensor(row.fstGlobalSensor);
    const bool validDisk = fwdAlignValidFstDisk(row.fstDisk);
    const bool validWedge = fwdAlignValidFstWedge(row.fstWedge);
    const bool validSensor = fwdAlignValidFstSensor(row.fstSensor);
    const bool fstAll = isFst && hasFstSensor;
    const bool fstMeasurement = fstAll && row.measurementDim == 2 &&
                                fwdAlignValid(row.meas0) &&
                                fwdAlignValid(row.meas1);
    const bool fitOk = requireFullyConverged ? row.fitConvergedFully > 0
                                             : row.fitConverged > 0;
    const bool trackOk =
        row.trackEta >= kFwdQaTrackEtaMin &&
        row.trackEta <= kFwdQaTrackEtaMax && row.trackP >= kFwdQaTrackPMin &&
        row.trackPt >= kFwdQaTrackPtMin &&
        (minTrackNFstHits <= 0 || !row.hasTrackNFstHitsBranch ||
         row.trackNFstHits >= minTrackNFstHits);
    const bool residualBase = fstAll && row.hasResidual > 0 &&
                              row.residualDim == 2 && fitOk && trackOk;
    const bool residualOk = residualBase && fwdAlignValid(row.resBiased0) &&
                            fwdAlignValid(row.resBiased1) &&
                            fwdAlignValid(row.resUnbiased0) &&
                            fwdAlignValid(row.resUnbiased1);
    const bool pullOk =
        residualOk && row.hasPullBranches && fwdAlignValid(row.pullBiased0) &&
        fwdAlignValid(row.pullBiased1) && fwdAlignValid(row.pullUnbiased0) &&
        fwdAlignValid(row.pullUnbiased1) && row.resBiasedSigma0 > 0 &&
        row.resBiasedSigma1 > 0 && row.resUnbiasedSigma0 > 0 &&
        row.resUnbiasedSigma1 > 0;
    const bool rawPredictionOk =
        residualOk && row.hasRawPredictionBranches &&
        fwdAlignValid(row.fstRawR) && fwdAlignValid(row.fstMeanPhiStrip) &&
        fwdAlignValid(row.trackPred0) && fwdAlignValid(row.trackPred1);
    const bool closureOk =
        fstMeasurement && row.hasFstClosureBranches &&
        fwdAlignValid(row.fstClosureX) && fwdAlignValid(row.fstClosureY) &&
        fwdAlignValid(row.fstClosureZ) && fwdAlignValid(row.fstClosureU) &&
        fwdAlignValid(row.fstClosureV) && fwdAlignValid(row.fstClosureMag);
    const bool closureRawOk =
        closureOk && row.hasRawPredictionBranches && validSensor &&
        fwdAlignValid(row.fstRawR) && fwdAlignValid(row.fstMeanPhiStrip);
    const bool slopeOk =
        residualOk && row.hasTrackSlopeBranches && row.trackPz != 0;

    hDetId->Fill(row.detId);
    if (fstAll) {
      ++nFstRows;
      hDim->Fill(row.measurementDim, row.residualDim);
      hHasResidual->Fill(row.hasResidual);
      hSortMinusSensor->Fill(row.sorting - row.fstGlobalSensor);
    }

    if (fstMeasurement) {
      ++nFstMeasurementRows;
      hMeas1VsMeas0->Fill(row.meas0, row.meas1);
      if (validDisk)
        hMeas1VsMeas0Disk[row.fstDisk]->Fill(row.meas0, row.meas1);
      if (validWedge)
        hMeas1VsMeas0Wedge[row.fstWedge]->Fill(row.meas0, row.meas1);
      if (validSensor)
        hMeas1VsMeas0Sensor[row.fstSensor]->Fill(row.meas0, row.meas1);
      if (validGlobalSensor)
        hMeas1VsMeas0GlobalSensor[row.fstGlobalSensor]->Fill(row.meas0,
                                                             row.meas1);
    }

    if (closureOk) {
      ++nFstClosureRows;
      const double closureXUm = row.fstClosureX * kFwdQaCmToMicron;
      const double closureYUm = row.fstClosureY * kFwdQaCmToMicron;
      const double closureZUm = row.fstClosureZ * kFwdQaCmToMicron;
      const double closureUUm = row.fstClosureU * kFwdQaCmToMicron;
      const double closureVUm = row.fstClosureV * kFwdQaCmToMicron;
      const double closureMagUm = row.fstClosureMag * kFwdQaCmToMicron;
      hFstClosureX->Fill(closureXUm);
      hFstClosureY->Fill(closureYUm);
      hFstClosureZ->Fill(closureZUm);
      hFstClosureU->Fill(closureUUm);
      hFstClosureV->Fill(closureVUm);
      hFstClosureMag->Fill(closureMagUm);
      if (validGlobalSensor) {
        hFstClosureUByGlobalSensor->Fill(row.fstGlobalSensor,
                                         closureUUm);
        hFstClosureVByGlobalSensor->Fill(row.fstGlobalSensor,
                                         closureVUm);
        pFstClosureUByGlobalSensor->Fill(row.fstGlobalSensor,
                                         closureUUm);
        pFstClosureVByGlobalSensor->Fill(row.fstGlobalSensor,
                                         closureVUm);
        pFstClosureMagByGlobalSensor->Fill(row.fstGlobalSensor,
                                           closureMagUm);
      }
      if (closureRawOk) {
        pFstClosureUVsRawR[row.fstSensor]->Fill(row.fstRawR,
                                                closureUUm);
        pFstClosureVVsRawR[row.fstSensor]->Fill(row.fstRawR,
                                                closureVUm);
        pFstClosureUVsMeanPhi[row.fstSensor]->Fill(row.fstMeanPhiStrip,
                                                   closureUUm);
        pFstClosureVVsMeanPhi[row.fstSensor]->Fill(row.fstMeanPhiStrip,
                                                   closureVUm);
      }
    }

    if (!residualOk)
      continue;

    ++nFstResidualRows;
    const double resUUm = row.resUnbiased0 * 10000.0;
    const double resVUm = row.resUnbiased1 * 10000.0;
    const double biasedUUm = row.resBiased0 * 10000.0;
    const double biasedVUm = row.resBiased1 * 10000.0;

    hResU->Fill(resUUm);
    hResV->Fill(resVUm);
    hBiasedU->Fill(biasedUUm);
    hBiasedV->Fill(biasedVUm);
    hFstDisk->Fill(row.fstDisk);
    hFstSensor->Fill(row.fstGlobalSensor);
    hWedgeDisk->Fill(row.fstDisk, row.fstWedge);
    hSensorDisk->Fill(row.fstDisk, row.fstSensor);

    pResUByDisk->Fill(row.fstDisk, resUUm);
    pResVByDisk->Fill(row.fstDisk, resVUm);
    pResUByWedge->Fill(row.fstWedge, resUUm);
    pResVByWedge->Fill(row.fstWedge, resVUm);
    pResUByGlobalSensor->Fill(row.fstGlobalSensor, resUUm);
    pResVByGlobalSensor->Fill(row.fstGlobalSensor, resVUm);
    hResUByGlobalSensor->Fill(row.fstGlobalSensor, resUUm);
    hResVByGlobalSensor->Fill(row.fstGlobalSensor, resVUm);
    if (validDisk) {
      pResUBySensorDisk[row.fstDisk]->Fill(row.fstGlobalSensor, resUUm);
      pResVBySensorDisk[row.fstDisk]->Fill(row.fstGlobalSensor, resVUm);
    }

    hResUVsMeas0->Fill(row.meas0, resUUm);
    hResUVsMeas1->Fill(row.meas1, resUUm);
    hResVVsMeas0->Fill(row.meas0, resVUm);
    hResVVsMeas1->Fill(row.meas1, resVUm);

    if (pullOk) {
      hPullU->Fill(row.pullUnbiased0);
      hPullV->Fill(row.pullUnbiased1);
      hPullBiasedU->Fill(row.pullBiased0);
      hPullBiasedV->Fill(row.pullBiased1);
      hSigmaU->Fill(row.resUnbiasedSigma0 * 10000.0);
      hSigmaV->Fill(row.resUnbiasedSigma1 * 10000.0);
      hSigmaBiasedU->Fill(row.resBiasedSigma0 * 10000.0);
      hSigmaBiasedV->Fill(row.resBiasedSigma1 * 10000.0);
      pPullUByDisk->Fill(row.fstDisk, row.pullUnbiased0);
      pPullVByDisk->Fill(row.fstDisk, row.pullUnbiased1);
      pPullUByWedge->Fill(row.fstWedge, row.pullUnbiased0);
      pPullVByWedge->Fill(row.fstWedge, row.pullUnbiased1);
      pPullUByGlobalSensor->Fill(row.fstGlobalSensor, row.pullUnbiased0);
      pPullVByGlobalSensor->Fill(row.fstGlobalSensor, row.pullUnbiased1);
      hPullUByGlobalSensor->Fill(row.fstGlobalSensor, row.pullUnbiased0);
      hPullVByGlobalSensor->Fill(row.fstGlobalSensor, row.pullUnbiased1);
      hPullUVsMeas0->Fill(row.meas0, row.pullUnbiased0);
      hPullUVsMeas1->Fill(row.meas1, row.pullUnbiased0);
      hPullVVsMeas0->Fill(row.meas0, row.pullUnbiased1);
      hPullVVsMeas1->Fill(row.meas1, row.pullUnbiased1);
      if (validDisk) {
        pPullUBySensorDisk[row.fstDisk]->Fill(row.fstGlobalSensor,
                                              row.pullUnbiased0);
        pPullVBySensorDisk[row.fstDisk]->Fill(row.fstGlobalSensor,
                                              row.pullUnbiased1);
      }
    }

    if (rawPredictionOk && validSensor) {
      hRes0VsRawR[row.fstSensor]->Fill(row.fstRawR, row.resUnbiased0);
      hRes0VsMeanPhi[row.fstSensor]->Fill(row.fstMeanPhiStrip,
                                          row.resUnbiased0);
      hTrackPred0VsMeas0[row.fstSensor]->Fill(row.meas0, row.trackPred0);
      hTrackPred1VsMeas1[row.fstSensor]->Fill(row.meas1, row.trackPred1);
      pRes0VsRawR[row.fstSensor]->Fill(row.fstRawR, row.resUnbiased0);
      pRes0VsMeanPhi[row.fstSensor]->Fill(row.fstMeanPhiStrip,
                                          row.resUnbiased0);
      if (validDisk)
        pRes0VsRawRDiskSensor[row.fstDisk][row.fstSensor]->Fill(
            row.fstRawR, row.resUnbiased0);
    }

    if (slopeOk && validSensor) {
      const double polarSlope = row.trackPt / row.trackPz;
      hRes0VsEta[row.fstSensor]->Fill(row.trackEta, row.resUnbiased0);
      hRes0VsPolarSlope[row.fstSensor]->Fill(polarSlope, row.resUnbiased0);
      pRes0VsEta[row.fstSensor]->Fill(row.trackEta, row.resUnbiased0);
      pRes0VsPolarSlope[row.fstSensor]->Fill(polarSlope, row.resUnbiased0);
      if (validDisk) {
        pRes0VsEtaDiskSensor[row.fstDisk][row.fstSensor]->Fill(
            row.trackEta, row.resUnbiased0);
        pRes0VsPolarSlopeDiskSensor[row.fstDisk][row.fstSensor]->Fill(
            polarSlope, row.resUnbiased0);
      }
    }
  }

  // 4. Terminal summary.
  std::cout << "Input file: " << inputFilename << std::endl;
  std::cout << "fwdAlign entries: " << nEntries << std::endl;
  std::cout << "FST rows: " << nFstRows << std::endl;
  std::cout << "FST measurement rows after cuts: " << nFstMeasurementRows
            << std::endl;
  std::cout << "FST closure rows: " << nFstClosureRows << std::endl;
  std::cout << "FST residual rows after cuts: " << nFstResidualRows
            << std::endl;
  std::cout << "Residual/pull track cut: " << trackCut.Data() << std::endl;
  std::cout << "Residual/pull base cut: detId==45&&hasResidual>0"
            << "&&residualDim==2&&fstGlobalSensor>=0&&" << fitCut.Data() << "&&"
            << trackCut.Data() << "&&resBiased0>-90000&&resBiased1>-90000"
            << "&&resUnbiased0>-90000&&resUnbiased1>-90000" << std::endl;
  std::cout << TString::Format("Residual/pull vs local-meas binning: meas0 "
                               "%.1f..%.1f cm, meas1 %.1f..%.1f cm, 1 cm bins",
                               kFwdAlignLocalMeas0Min, kFwdAlignLocalMeas0Max,
                               kFwdAlignLocalMeas1Min, kFwdAlignLocalMeas1Max)
                   .Data()
            << std::endl;
  std::cout << "Residual units in plots: microns (tree stores cm)" << std::endl;
  std::cout << "Closure units in plots: microns (tree stores cm)" << std::endl;
  std::cout << "Pull branches available: "
            << (row.hasPullBranches ? "yes" : "no") << std::endl;
  std::cout << "Raw-strip/prediction diagnostic branches available: "
            << (row.hasRawPredictionBranches ? "yes" : "no") << std::endl;
  std::cout << "FST measurement closure branches available: "
            << (row.hasFstClosureBranches ? "yes" : "no") << std::endl;
  std::cout << "Track-slope diagnostic branches available: "
            << (row.hasTrackSlopeBranches ? "yes" : "no") << std::endl;
  std::cout << "trackNFstHits branch available: "
            << (row.hasTrackNFstHitsBranch ? "yes" : "no") << std::endl;
  if (minTrackNFstHits > 0 && !row.hasTrackNFstHitsBranch) {
    std::cout << "Requested trackNFstHits>=" << minTrackNFstHits
              << ", but this input tree does not have trackNFstHits; skipping "
                 "that cut."
              << std::endl;
  }

  // 5. Draw PDF pages. Histograms were filled above; this section only controls
  // page layout.
  TCanvas *canvas =
      new TCanvas("cFwdAlignQa", "Forward alignment QA", 1200, 900);
  canvas->Print(pdfOutput + "[");

  canvas->Clear();
  TPaveText *summary = new TPaveText(0.08, 0.14, 0.92, 0.88, "NDC");
  summary->SetFillColor(0);
  summary->SetTextAlign(12);
  summary->AddText("Forward alignment residual QA");
  summary->AddText("Implementation: explicit branch loop, no TTree::Draw");
  summary->AddText(TString::Format("Input: %s", inputFilename));
  summary->AddText(TString::Format("Output ROOT: %s", rootOutput.Data()));
  summary->AddText(TString::Format("Output PDF: %s", pdfOutput.Data()));
  summary->AddText(TString::Format("Total fwdAlign rows: %lld", nEntries));
  summary->AddText(TString::Format("FST rows: %lld", nFstRows));
  summary->AddText(TString::Format("FST measurement rows after cuts: %lld",
                                   nFstMeasurementRows));
  summary->AddText(TString::Format("FST closure rows: %lld",
                                   nFstClosureRows));
  summary->AddText(
      TString::Format("FST residual rows after cuts: %lld", nFstResidualRows));
  summary->AddText(TString::Format("Fit cut: %s", fitCut.Data()));
  summary->AddText(
      TString::Format("Residual/pull track cut: %s", trackCut.Data()));
  if (minTrackNFstHits > 0 && !row.hasTrackNFstHitsBranch) {
    summary->AddText(TString::Format(
        "trackNFstHits>=%d requested but branch is missing", minTrackNFstHits));
  }
  summary->AddText(
      TString::Format("Residual/pull vs local-meas bins: meas0 %.1f..%.1f cm, "
                      "meas1 %.1f..%.1f cm, 1 cm bins.",
                      kFwdAlignLocalMeas0Min, kFwdAlignLocalMeas0Max,
                      kFwdAlignLocalMeas1Min, kFwdAlignLocalMeas1Max));
  summary->AddText("Measurement maps do not apply the trackP/trackEta cut.");
  summary->AddText("Residual plots use resUnbiased/resBiased * 10000, so units "
                   "are microns.");
  summary->AddText("Closure plots use fstClosure* * 10000, so units are "
                   "microns.");
  summary->AddText(row.hasPullBranches
                       ? "Pull plots are true GenFit residual/sigma pulls."
                       : "Pull branches are not present in this input; rerun "
                         "afterburner with updated maker.");
  summary->AddText(row.hasRawPredictionBranches
                       ? "Raw-strip and trackPred diagnostic plots are enabled."
                       : "Raw-strip or trackPred branches are missing; "
                         "raw-strip diagnostics skipped.");
  summary->AddText(row.hasFstClosureBranches
                       ? "FST measurement closure plots are enabled."
                       : "FST closure branches are missing; rerun afterburner "
                         "with updated maker.");
  summary->AddText(row.hasTrackSlopeBranches
                       ? "Track eta/slope diagnostic plots are enabled."
                       : "Track slope branches are missing; slope diagnostics "
                         "skipped.");
  summary->AddText(
      "FST local coordinates are the GenFit planar-measurement coordinates.");
  summary->AddText("FTT is intentionally not analyzed in this macro.");
  summary->Draw();
  fwdAlignSavePage(canvas, pdfOutput);

  canvas->Clear();
  canvas->Divide(2, 2);
  canvas->cd(1);
  hDetId->Draw();
  canvas->cd(2);
  hDim->Draw("COLZ TEXT");
  canvas->cd(3);
  hHasResidual->Draw();
  canvas->cd(4);
  hSortMinusSensor->Draw();
  fwdAlignSavePage(canvas, pdfOutput);

  canvas->Clear();
  hMeas1VsMeas0->Draw("COLZ");
  fwdAlignSavePage(canvas, pdfOutput);

  canvas->Clear();
  canvas->Divide(3, 1);
  for (int disk = 0; disk < kFwdAlignNumFstDisks; ++disk) {
    canvas->cd(disk + 1);
    hMeas1VsMeas0Disk[disk]->Draw("COLZ");
  }
  fwdAlignSavePage(canvas, pdfOutput);

  canvas->Clear();
  canvas->Divide(4, 3);
  for (int wedge = 0; wedge < kFwdAlignNumFstWedges; ++wedge) {
    canvas->cd(wedge + 1);
    hMeas1VsMeas0Wedge[wedge]->Draw("COLZ");
  }
  fwdAlignSavePage(canvas, pdfOutput);

  canvas->Clear();
  canvas->Divide(3, 1);
  for (int sensor = 0; sensor < kFwdAlignNumFstSensorsPerWedge; ++sensor) {
    canvas->cd(sensor + 1);
    hMeas1VsMeas0Sensor[sensor]->Draw("COLZ");
  }
  fwdAlignSavePage(canvas, pdfOutput);

  for (int firstGlobalSensor = 0; firstGlobalSensor < kFwdAlignNumFstSensors;
       firstGlobalSensor += 12) {
    canvas->Clear();
    canvas->Divide(4, 3);
    for (int offset = 0; offset < 12; ++offset) {
      canvas->cd(offset + 1);
      hMeas1VsMeas0GlobalSensor[firstGlobalSensor + offset]->Draw("COLZ");
    }
    fwdAlignSavePage(canvas, pdfOutput);
  }

  if (row.hasFstClosureBranches) {
    canvas->Clear();
    canvas->Divide(3, 2);
    canvas->cd(1);
    hFstClosureX->Draw();
    fwdAlignDrawLine(0, 0, 0, hFstClosureX->GetMaximum());
    canvas->cd(2);
    hFstClosureY->Draw();
    fwdAlignDrawLine(0, 0, 0, hFstClosureY->GetMaximum());
    canvas->cd(3);
    hFstClosureZ->Draw();
    fwdAlignDrawLine(0, 0, 0, hFstClosureZ->GetMaximum());
    canvas->cd(4);
    hFstClosureU->Draw();
    fwdAlignDrawLine(0, 0, 0, hFstClosureU->GetMaximum());
    canvas->cd(5);
    hFstClosureV->Draw();
    fwdAlignDrawLine(0, 0, 0, hFstClosureV->GetMaximum());
    canvas->cd(6);
    hFstClosureMag->Draw();
    fwdAlignDrawLine(0, 0, 0, hFstClosureMag->GetMaximum());
    fwdAlignSavePage(canvas, pdfOutput);

    canvas->Clear();
    canvas->Divide(1, 3);
    canvas->cd(1);
    pFstClosureUByGlobalSensor->Draw("E1");
    fwdAlignDrawLine(-0.5, 0.0, 107.5, 0.0);
    canvas->cd(2);
    pFstClosureVByGlobalSensor->Draw("E1");
    fwdAlignDrawLine(-0.5, 0.0, 107.5, 0.0);
    canvas->cd(3);
    pFstClosureMagByGlobalSensor->Draw("E1");
    fwdAlignSavePage(canvas, pdfOutput);

    canvas->Clear();
    canvas->Divide(1, 2);
    canvas->cd(1);
    hFstClosureUByGlobalSensor->Draw("COLZ");
    fwdAlignDrawLine(-0.5, 0.0, 107.5, 0.0);
    canvas->cd(2);
    hFstClosureVByGlobalSensor->Draw("COLZ");
    fwdAlignDrawLine(-0.5, 0.0, 107.5, 0.0);
    fwdAlignSavePage(canvas, pdfOutput);

    if (row.hasRawPredictionBranches) {
      canvas->Clear();
      canvas->Divide(3, 2);
      for (int sensor = 0; sensor < kFwdAlignNumFstSensorsPerWedge; ++sensor) {
        canvas->cd(sensor + 1);
        pFstClosureUVsRawR[sensor]->Draw("E1");
        fwdAlignDrawLine(kFwdQaRawRMin, 0.0, kFwdQaRawRMax, 0.0);
        canvas->cd(sensor + 4);
        pFstClosureVVsRawR[sensor]->Draw("E1");
        fwdAlignDrawLine(kFwdQaRawRMin, 0.0, kFwdQaRawRMax, 0.0);
      }
      fwdAlignSavePage(canvas, pdfOutput);

      canvas->Clear();
      canvas->Divide(3, 2);
      for (int sensor = 0; sensor < kFwdAlignNumFstSensorsPerWedge; ++sensor) {
        canvas->cd(sensor + 1);
        pFstClosureUVsMeanPhi[sensor]->Draw("E1");
        fwdAlignDrawLine(-0.5, 0.0, 127.5, 0.0);
        canvas->cd(sensor + 4);
        pFstClosureVVsMeanPhi[sensor]->Draw("E1");
        fwdAlignDrawLine(-0.5, 0.0, 127.5, 0.0);
      }
      fwdAlignSavePage(canvas, pdfOutput);
    }
  }

  canvas->Clear();
  canvas->Divide(2, 2);
  canvas->cd(1);
  hResU->Draw();
  fwdAlignDrawLine(0, 0, 0, hResU->GetMaximum());
  canvas->cd(2);
  hResV->Draw();
  fwdAlignDrawLine(0, 0, 0, hResV->GetMaximum());
  canvas->cd(3);
  hBiasedU->Draw();
  fwdAlignDrawLine(0, 0, 0, hBiasedU->GetMaximum());
  canvas->cd(4);
  hBiasedV->Draw();
  fwdAlignDrawLine(0, 0, 0, hBiasedV->GetMaximum());
  fwdAlignSavePage(canvas, pdfOutput);

  if (row.hasPullBranches) {
    canvas->Clear();
    canvas->Divide(2, 2);
    canvas->cd(1);
    hPullU->Draw();
    fwdAlignDrawLine(0, 0, 0, hPullU->GetMaximum());
    canvas->cd(2);
    hPullV->Draw();
    fwdAlignDrawLine(0, 0, 0, hPullV->GetMaximum());
    canvas->cd(3);
    hPullBiasedU->Draw();
    fwdAlignDrawLine(0, 0, 0, hPullBiasedU->GetMaximum());
    canvas->cd(4);
    hPullBiasedV->Draw();
    fwdAlignDrawLine(0, 0, 0, hPullBiasedV->GetMaximum());
    fwdAlignSavePage(canvas, pdfOutput);

    canvas->Clear();
    canvas->Divide(2, 2);
    canvas->cd(1);
    hSigmaU->Draw();
    canvas->cd(2);
    hSigmaV->Draw();
    canvas->cd(3);
    hSigmaBiasedU->Draw();
    canvas->cd(4);
    hSigmaBiasedV->Draw();
    fwdAlignSavePage(canvas, pdfOutput);
  }

  canvas->Clear();
  canvas->Divide(2, 2);
  canvas->cd(1);
  hFstDisk->Draw();
  canvas->cd(2);
  hFstSensor->Draw();
  canvas->cd(3);
  hWedgeDisk->Draw("COLZ TEXT");
  canvas->cd(4);
  hSensorDisk->Draw("COLZ TEXT");
  fwdAlignSavePage(canvas, pdfOutput);

  canvas->Clear();
  canvas->Divide(2, 2);
  canvas->cd(1);
  pResUByDisk->Draw("E1");
  fwdAlignDrawLine(-0.5, 0, 2.5, 0);
  canvas->cd(2);
  pResVByDisk->Draw("E1");
  fwdAlignDrawLine(-0.5, 0, 2.5, 0);
  canvas->cd(3);
  pResUByWedge->Draw("E1");
  fwdAlignDrawLine(-0.5, 0, 11.5, 0);
  canvas->cd(4);
  pResVByWedge->Draw("E1");
  fwdAlignDrawLine(-0.5, 0, 11.5, 0);
  fwdAlignSavePage(canvas, pdfOutput);

  if (row.hasPullBranches) {
    canvas->Clear();
    canvas->Divide(2, 2);
    canvas->cd(1);
    pPullUByDisk->Draw("E1");
    fwdAlignDrawLine(-0.5, 0, 2.5, 0);
    canvas->cd(2);
    pPullVByDisk->Draw("E1");
    fwdAlignDrawLine(-0.5, 0, 2.5, 0);
    canvas->cd(3);
    pPullUByWedge->Draw("E1");
    fwdAlignDrawLine(-0.5, 0, 11.5, 0);
    canvas->cd(4);
    pPullVByWedge->Draw("E1");
    fwdAlignDrawLine(-0.5, 0, 11.5, 0);
    fwdAlignSavePage(canvas, pdfOutput);
  }

  canvas->Clear();
  canvas->Divide(1, 2);
  canvas->cd(1);
  pResUByGlobalSensor->Draw("E1");
  fwdAlignDrawLine(-0.5, 0, 107.5, 0);
  canvas->cd(2);
  pResVByGlobalSensor->Draw("E1");
  fwdAlignDrawLine(-0.5, 0, 107.5, 0);
  fwdAlignSavePage(canvas, pdfOutput);

  canvas->Clear();
  canvas->Divide(1, 2);
  canvas->cd(1);
  hResUByGlobalSensor->Draw("COLZ");
  canvas->cd(2);
  hResVByGlobalSensor->Draw("COLZ");
  fwdAlignSavePage(canvas, pdfOutput);

  if (row.hasPullBranches) {
    canvas->Clear();
    canvas->Divide(1, 2);
    canvas->cd(1);
    pPullUByGlobalSensor->Draw("E1");
    fwdAlignDrawLine(-0.5, 0, 107.5, 0);
    canvas->cd(2);
    pPullVByGlobalSensor->Draw("E1");
    fwdAlignDrawLine(-0.5, 0, 107.5, 0);
    fwdAlignSavePage(canvas, pdfOutput);

    canvas->Clear();
    canvas->Divide(1, 2);
    canvas->cd(1);
    hPullUByGlobalSensor->Draw("COLZ");
    canvas->cd(2);
    hPullVByGlobalSensor->Draw("COLZ");
    fwdAlignSavePage(canvas, pdfOutput);
  }

  canvas->Clear();
  canvas->Divide(2, 2);
  canvas->cd(1);
  hResUVsMeas0->Draw("COLZ");
  canvas->cd(2);
  hResUVsMeas1->Draw("COLZ");
  canvas->cd(3);
  hResVVsMeas0->Draw("COLZ");
  canvas->cd(4);
  hResVVsMeas1->Draw("COLZ");
  fwdAlignSavePage(canvas, pdfOutput);

  if (row.hasPullBranches) {
    canvas->Clear();
    canvas->Divide(2, 2);
    canvas->cd(1);
    hPullUVsMeas0->Draw("COLZ");
    canvas->cd(2);
    hPullUVsMeas1->Draw("COLZ");
    canvas->cd(3);
    hPullVVsMeas0->Draw("COLZ");
    canvas->cd(4);
    hPullVVsMeas1->Draw("COLZ");
    fwdAlignSavePage(canvas, pdfOutput);
  }

  if (row.hasRawPredictionBranches) {
    canvas->Clear();
    canvas->Divide(3, 2);
    for (int sensor = 0; sensor < kFwdAlignNumFstSensorsPerWedge; ++sensor) {
      canvas->cd(sensor + 1);
      hRes0VsRawR[sensor]->Draw("COLZ");
      fwdAlignDrawLine(kFwdQaRawRMin, 0.0, kFwdQaRawRMax, 0.0);
      canvas->cd(sensor + 4);
      hRes0VsMeanPhi[sensor]->Draw("COLZ");
      fwdAlignDrawLine(-0.5, 0.0, 127.5, 0.0);
    }
    fwdAlignSavePage(canvas, pdfOutput);

    canvas->Clear();
    canvas->Divide(3, 2);
    for (int sensor = 0; sensor < kFwdAlignNumFstSensorsPerWedge; ++sensor) {
      canvas->cd(sensor + 1);
      hTrackPred0VsMeas0[sensor]->Draw("COLZ");
      fwdAlignDrawLine(-10.0, -10.0, 10.0, 10.0, kRed + 1, 2);
      canvas->cd(sensor + 4);
      hTrackPred1VsMeas1[sensor]->Draw("COLZ");
      fwdAlignDrawLine(-12.0, -12.0, 12.0, 12.0, kRed + 1, 2);
    }
    fwdAlignSavePage(canvas, pdfOutput);

    canvas->Clear();
    canvas->Divide(3, 2);
    for (int sensor = 0; sensor < kFwdAlignNumFstSensorsPerWedge; ++sensor) {
      canvas->cd(sensor + 1);
      pRes0VsRawR[sensor]->Draw("E1");
      fwdAlignDrawLine(kFwdQaRawRMin, 0.0, kFwdQaRawRMax, 0.0);
      canvas->cd(sensor + 4);
      pRes0VsMeanPhi[sensor]->Draw("E1");
      fwdAlignDrawLine(-0.5, 0.0, 127.5, 0.0);
    }
    fwdAlignSavePage(canvas, pdfOutput);

    for (int disk = 0; disk < kFwdAlignNumFstDisks; ++disk) {
      canvas->Clear();
      canvas->Divide(3, 1);
      for (int sensor = 0; sensor < kFwdAlignNumFstSensorsPerWedge; ++sensor) {
        canvas->cd(sensor + 1);
        pRes0VsRawRDiskSensor[disk][sensor]->Draw("E1");
        fwdAlignDrawLine(kFwdQaRawRMin, 0.0, kFwdQaRawRMax, 0.0);
      }
      fwdAlignSavePage(canvas, pdfOutput);
    }
  }

  if (row.hasTrackSlopeBranches) {
    canvas->Clear();
    canvas->Divide(3, 2);
    for (int sensor = 0; sensor < kFwdAlignNumFstSensorsPerWedge; ++sensor) {
      canvas->cd(sensor + 1);
      hRes0VsEta[sensor]->Draw("COLZ");
      fwdAlignDrawLine(kFwdQaTrackEtaMin, 0.0, kFwdQaTrackEtaMax, 0.0);
      canvas->cd(sensor + 4);
      hRes0VsPolarSlope[sensor]->Draw("COLZ");
      fwdAlignDrawLine(kFwdQaPolarSlopeMin, 0.0, kFwdQaPolarSlopeMax, 0.0);
    }
    fwdAlignSavePage(canvas, pdfOutput);

    canvas->Clear();
    canvas->Divide(3, 2);
    for (int sensor = 0; sensor < kFwdAlignNumFstSensorsPerWedge; ++sensor) {
      canvas->cd(sensor + 1);
      pRes0VsEta[sensor]->Draw("E1");
      fwdAlignDrawLine(kFwdQaTrackEtaMin, 0.0, kFwdQaTrackEtaMax, 0.0);
      canvas->cd(sensor + 4);
      pRes0VsPolarSlope[sensor]->Draw("E1");
      fwdAlignDrawLine(kFwdQaPolarSlopeMin, 0.0, kFwdQaPolarSlopeMax, 0.0);
    }
    fwdAlignSavePage(canvas, pdfOutput);

    for (int disk = 0; disk < kFwdAlignNumFstDisks; ++disk) {
      canvas->Clear();
      canvas->Divide(3, 2);
      for (int sensor = 0; sensor < kFwdAlignNumFstSensorsPerWedge; ++sensor) {
        canvas->cd(sensor + 1);
        pRes0VsEtaDiskSensor[disk][sensor]->Draw("E1");
        fwdAlignDrawLine(kFwdQaTrackEtaMin, 0.0, kFwdQaTrackEtaMax, 0.0);
        canvas->cd(sensor + 4);
        pRes0VsPolarSlopeDiskSensor[disk][sensor]->Draw("E1");
        fwdAlignDrawLine(kFwdQaPolarSlopeMin, 0.0, kFwdQaPolarSlopeMax, 0.0);
      }
      fwdAlignSavePage(canvas, pdfOutput);
    }
  }

  canvas->Clear();
  canvas->Divide(3, 2);
  for (int disk = 0; disk < kFwdAlignNumFstDisks; ++disk) {
    const double xMin = disk * 36 - 0.5;
    const double xMax = (disk + 1) * 36 - 0.5;
    canvas->cd(disk + 1);
    pResUBySensorDisk[disk]->Draw("E1");
    fwdAlignDrawLine(xMin, 0, xMax, 0);
    canvas->cd(disk + 4);
    pResVBySensorDisk[disk]->Draw("E1");
    fwdAlignDrawLine(xMin, 0, xMax, 0);
  }
  fwdAlignSavePage(canvas, pdfOutput);

  if (row.hasPullBranches) {
    canvas->Clear();
    canvas->Divide(3, 2);
    for (int disk = 0; disk < kFwdAlignNumFstDisks; ++disk) {
      const double xMin = disk * 36 - 0.5;
      const double xMax = (disk + 1) * 36 - 0.5;
      canvas->cd(disk + 1);
      pPullUBySensorDisk[disk]->Draw("E1");
      fwdAlignDrawLine(xMin, 0, xMax, 0);
      canvas->cd(disk + 4);
      pPullVBySensorDisk[disk]->Draw("E1");
      fwdAlignDrawLine(xMin, 0, xMax, 0);
    }
    fwdAlignSavePage(canvas, pdfOutput);
  }

  outputFile->Write();
  canvas->Print(pdfOutput + "]");
  outputFile->Close();
  inputFile->Close();

  std::cout << "Wrote " << rootOutput.Data() << std::endl;
  std::cout << "Wrote " << pdfOutput.Data() << std::endl;
}
