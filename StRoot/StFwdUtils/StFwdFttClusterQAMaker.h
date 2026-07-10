#ifndef ST_FWD_FTT_CLUSTER_QA_MAKER_H
#define ST_FWD_FTT_CLUSTER_QA_MAKER_H

// ----------------------------------------------------------------------------
// StFwdFttClusterQAMaker
//
// Diagnostic maker that quantifies how often the three independent sTGC strip
// directions in a quadrant overlap so that a 2D "point" can be formed:
//     H (horizontal strips)            measure y
//     V (vertical strips)              measure x
//     D (diagonal 45deg strips,        measure x+y
//        kFttDiagonalH + kFttDiagonalV)
// Each direction is an INDEPENDENT avalanche, so a wrong H+V combination is a
// "ghost". Any 2 of {H,V,D} determine a point; the 3rd confirms / rejects it.
// Headline question: in real data do >=2 directions come together in valid
// strip-group regions, or do strips mostly fire independently (only 1 of 3)?
//
// This maker reads the StEvent StFttCollection clusters (which carry row(),
// orientation(), local x() in mm, plane/quadrant), so the rigorous row-aware
// strip-group geometry can be applied. The overlap logic mirrors the rigorous
// point maker:
//   - H x V validity  : is_Group1..8 (copied from
//                       StRoot/StFttClusterPointMaker/StFttClusterPointMaker.h)
//   - diagonal match  : |(x+y) - clu_d->x()*sqrt(2)| < 1.60*3 mm, with the
//                       DiagV/DiagH precedence of StFttPointMakerGroups.cxx.
//
// Add to a chain after the FTT cluster maker / fwdTrack, e.g.:
//     StFwdFttClusterQAMaker* q = new StFwdFttClusterQAMaker();
//     q->setOutputFilename("StFwdFttClusterQA.root");
//     chain->AddAfter("fwdTrack", q);
// ----------------------------------------------------------------------------

#include <map>
#include <vector>

#include "StChain/StMaker.h"
#include "TString.h"
#include "TH1.h"
#include "TH2.h"

class StEvent;
class StFttCollection;
class StFttDb;

class StFwdFttClusterQAMaker : public StMaker {
  public:
    StFwdFttClusterQAMaker();
    ~StFwdFttClusterQAMaker() {/* nada */};

    int Init();
    int Make();
    int Finish();

    void setOutputFilename(TString f) { mOutputFilename = f; }
    void setDebug(bool d = true) { mDebug = d; }

  protected:
    // ---- per-cluster info (plain fields) ----
    struct CluInfo {
        int    row;
        double x;
    };

    // number of "plane" history slots: planes 0..3 plus index 4 = all planes
    static const int kNP = 5;

    // ---- histogram bookkeeping (mirrors StFwdFitQAMaker) ----
    std::map<TString, TH1*> mHists;

    TH1* addHist(TH1* h) {
        mHists[h->GetName()] = h;
        h->SetDirectory(0);
        return h;
    }
    TH1* getHist(TString n) {
        if (mHists.count(n)) return mHists[n];
        LOG_ERROR << "Missing histogram: " << n.Data() << endm;
        return new TH1F("NULL", "NULL", 1, 0, 1);
    }
    // fill the per-plane histogram and its "all-planes" companion
    void fillPA(TString base, int plane, double v) {
        getHist(Form("%s_p%d", base.Data(), plane))->Fill(v);
        getHist(Form("%s_all", base.Data()))->Fill(v);
    }
    TH2* h2(TString n) { return (TH2*)getHist(n); }

    void bookHistos();
    void processEvent();
    void printSummary();

    // returns index of closest diagonal within +/- kDiagWin, or -1; sets residual
    int matchDiagonal(double intercept, const std::vector<CluInfo>& diag, double& residual);

    // ---- rigorous strip-group geometry, copied verbatim from
    // StRoot/StFttClusterPointMaker/StFttClusterPointMaker.h:46-62
    // (row_x = V cluster row, row_y = H cluster row, x = V x, y = H x; mm) ----
    bool is_Group1(int row_x, int row_y, double x, double y) const { return ( (14.60 <= x && x <= 172.29) && (14.60 <= y && y <= 172.29) && (row_x == 0) && (row_y == 0) ); }
    bool is_Group2(int row_x, int row_y, double x, double y) const { return ( (172.29 <= x && x <= 360.09) && (14.60 <= y && y <= 172.29) && (row_x == 0) && (row_y == 1)); }
    bool is_Group3(int row_x, int row_y, double x, double y) const {
        return ( ( ( (360.09 <= x && x <= 504.2) && (14.60 <= y && y <= 172.29) ) || ( (504.2<= x && x <= 548.3) && (14.60 <= y && y <= 216.89) ) ) && (row_x == 0) && (row_y == 2) );
    }
    bool is_Group4(int row_x, int row_y, double x, double y) const { return ( (14.60 <= x && x <= 172.29) && (172.29 <= y && y <= 360.09) && (row_x == 1) && (row_y == 0) ); }
    bool is_Group5(int row_x, int row_y, double x, double y) const {
        return ( ( ((172.29 <= x && x <= 315.4) && (172.29 <= y && y <= 360.09)) || ((315.4 <= x && x <= 360.09) && (172.29 <= y && y <= 410.9)) || ((360.09 <= x && x <= 410.9) && (315.4 <= y && y <= 410.9)) ) && (row_x == 1) && (row_y == 1) );
    }
    bool is_Group6(int row_x, int row_y, double x, double y) const { return ((360.09 <= x && x <= 504.2) && (172.29 <= y && y <= 315.4)) && (row_x == 1) && (row_y == 2); }
    bool is_Group7(int row_x, int row_y, double x, double y) const {
        return ( ( ( (360.09 <= y && y <= 504.2) && (14.60 <= x && x <= 172.29) ) || ( (504.2<= y && y <= 548.3) && (14.60 <= x && x <= 216.89) ) ) && (row_x == 2) && (row_y == 0));
    }
    bool is_Group8(int row_x, int row_y, double x, double y) const { return ( ((360.09 <= y && y <= 504.2) && (172.29 <= x && x <= 315.4)) && (row_x == 2) && (row_y == 1)); }

    // reproduces GhostHitRejection_StripGroup: valid H x V crossing region
    bool inAnyStripGroup(int row_x, int row_y, double x, double y) const {
        return is_Group1(row_x,row_y,x,y) || is_Group2(row_x,row_y,x,y) ||
               is_Group3(row_x,row_y,x,y) || is_Group4(row_x,row_y,x,y) ||
               is_Group5(row_x,row_y,x,y) || is_Group6(row_x,row_y,x,y) ||
               is_Group7(row_x,row_y,x,y) || is_Group8(row_x,row_y,x,y);
    }
    // For H x D / V x D pairings only one strip row is known; scan the partner
    // rows 0..2 honoring the one known row. knownIsV==true -> knownRow is V row.
    bool inValidStripGroup(double x, double y, int knownRow, bool knownIsV) const {
        for (int r = 0; r <= 2; ++r) {
            if (knownIsV) { if (inAnyStripGroup(knownRow, r, x, y)) return true; }
            else          { if (inAnyStripGroup(r, knownRow, x, y)) return true; }
        }
        return false;
    }

    // ---- members ----
    StFttDb* mFttDb;
    bool     mDebug;
    TString  mOutputFilename;

    // diagonal match window, matches active path in StFttPointMakerGroups.cxx
    double kDiagWin;

    // per-event cluster store [rob 0..15][orientation 0..4]
    std::vector<CluInfo> mClu[16][5];

    // summary counters ([p]=plane 0..3, [4]=all)
    Long64_t cnt_nDir[kNP][4];
    double   sum_nH[kNP], sum_nV[kNP], sum_nD[kNP];
    Long64_t sum_cHV[kNP], sum_cHD[kNP], sum_cVD[kNP], sum_cHVconf[kNP], sum_cHVghost[kNP];
    Long64_t cnt_diagV, cnt_diagH;
    Long64_t nRobInstances, nEventsProcessed;

    ClassDef(StFwdFttClusterQAMaker, 0);
};

#endif
