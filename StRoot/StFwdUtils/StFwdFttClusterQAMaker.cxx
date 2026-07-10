#include "StFwdUtils/StFwdFttClusterQAMaker.h"

#include "TFile.h"
#include "TMath.h"

#include "StEvent/StEvent.h"
#include "StEvent/StEnumerations.h"
#include "StEvent/StFttCollection.h"
#include "StEvent/StFttCluster.h"

#include "StFttDbMaker/StFttDb.h"

//______________________________________________________________________________
StFwdFttClusterQAMaker::StFwdFttClusterQAMaker()
    : StMaker("fwdFttClusterQA"),
      mFttDb(nullptr),
      mDebug(false),
      mOutputFilename("StFwdFttClusterQA.root"),
      kDiagWin(1.60 * 3),   // 4.8 mm, matches StFttPointMakerGroups.cxx
      cnt_diagV(0), cnt_diagH(0),
      nRobInstances(0), nEventsProcessed(0) {
    for (int p = 0; p < kNP; ++p) {
        for (int d = 0; d < 4; ++d) cnt_nDir[p][d] = 0;
        sum_nH[p] = sum_nV[p] = sum_nD[p] = 0;
        sum_cHV[p] = sum_cHD[p] = sum_cVD[p] = sum_cHVconf[p] = sum_cHVghost[p] = 0;
    }
}

//______________________________________________________________________________
int StFwdFttClusterQAMaker::Init() {
    bookHistos();
    return kStOK;
}

//______________________________________________________________________________
void StFwdFttClusterQAMaker::bookHistos() {
    const char* lab[8] = {"none","H","V","HV","D","HD","VD","HVD"};
    for (int p = 0; p < kNP; ++p) {
        TString s = (p < 4) ? Form("p%d", p) : "all";

        addHist(new TH1I(Form("nDir_%s", s.Data()),   Form("# directions present per ROB (%s);n directions;ROBs", s.Data()), 4, -0.5, 3.5));
        TH1I* hs = new TH1I(Form("subset_%s", s.Data()), Form("direction subset per ROB (%s);;ROBs", s.Data()), 8, -0.5, 7.5);
        for (int b = 0; b < 8; ++b) hs->GetXaxis()->SetBinLabel(b + 1, lab[b]);
        addHist(hs);

        addHist(new TH1I(Form("nH_%s", s.Data()),  Form("# H clusters per ROB (%s);nH;ROBs", s.Data()),  21, -0.5, 20.5));
        addHist(new TH1I(Form("nV_%s", s.Data()),  Form("# V clusters per ROB (%s);nV;ROBs", s.Data()),  21, -0.5, 20.5));
        addHist(new TH1I(Form("nD_%s", s.Data()),  Form("# D clusters per ROB (%s);nD;ROBs", s.Data()),  21, -0.5, 20.5));
        addHist(new TH1I(Form("nDH_%s", s.Data()), Form("# DiagH clusters per ROB (%s);nDH;ROBs", s.Data()), 21, -0.5, 20.5));
        addHist(new TH1I(Form("nDV_%s", s.Data()), Form("# DiagV clusters per ROB (%s);nDV;ROBs", s.Data()), 21, -0.5, 20.5));

        addHist(new TH1I(Form("cHV_%s", s.Data()),      Form("HxV candidate points per ROB (%s);cHV;ROBs", s.Data()), 31, -0.5, 30.5));
        addHist(new TH1I(Form("cHD_%s", s.Data()),      Form("HxD candidate points per ROB (%s);cHD;ROBs", s.Data()), 31, -0.5, 30.5));
        addHist(new TH1I(Form("cVD_%s", s.Data()),      Form("VxD candidate points per ROB (%s);cVD;ROBs", s.Data()), 31, -0.5, 30.5));
        addHist(new TH1I(Form("cHVconf_%s", s.Data()),  Form("HxVxD confirmed points per ROB (%s);cHVconf;ROBs", s.Data()), 31, -0.5, 30.5));
        addHist(new TH1I(Form("cHVghost_%s", s.Data()), Form("HxV unconfirmed (ghost) per ROB (%s);cHVghost;ROBs", s.Data()), 31, -0.5, 30.5));

        addHist(new TH1F(Form("residual_%s", s.Data()),      Form("(x+y)-d*sqrt2, matched diagonal (%s);residual [mm];pairs", s.Data()), 200, -20, 20));
        addHist(new TH1F(Form("residualWide_%s", s.Data()),  Form("(x+y)-d*sqrt2, all diagonals (%s);residual [mm];pairs", s.Data()), 400, -100, 100));
    }

    addHist(new TH1I("cHV_evt",      "HxV candidates per event;cHV;events", 201, -0.5, 200.5));
    addHist(new TH1I("cHD_evt",      "HxD candidates per event;cHD;events", 201, -0.5, 200.5));
    addHist(new TH1I("cVD_evt",      "VxD candidates per event;cVD;events", 201, -0.5, 200.5));
    addHist(new TH1I("cHVconf_evt",  "HxVxD confirmed per event;cHVconf;events", 201, -0.5, 200.5));
    addHist(new TH1I("cHVghost_evt", "HxV ghost per event;cHVghost;events", 201, -0.5, 200.5));

    TH1I* hd = new TH1I("diagType", "confirming diagonal type;;count", 2, -0.5, 1.5);
    hd->GetXaxis()->SetBinLabel(1, "DiagV");
    hd->GetXaxis()->SetBinLabel(2, "DiagH");
    addHist(hd);

    for (int p = 0; p < 4; ++p) {
        addHist(new TH2F(Form("pointmap_p%d", p),   Form("HxVxD confirmed points, plane %d;x [mm];y [mm]", p), 140, 0, 560, 140, 0, 560));
        addHist(new TH2F(Form("pointmapHV_p%d", p), Form("HxV valid crossings, plane %d;x [mm];y [mm]", p),     140, 0, 560, 140, 0, 560));
        addHist(new TH2F(Form("pointmapHD_p%d", p), Form("HxD candidate points, plane %d;x [mm];y [mm]", p),    140, 0, 560, 140, 0, 560));
        addHist(new TH2F(Form("pointmapVD_p%d", p), Form("VxD candidate points, plane %d;x [mm];y [mm]", p),    140, 0, 560, 140, 0, 560));
    }
}

//______________________________________________________________________________
int StFwdFttClusterQAMaker::Make() {
    StEvent* event = (StEvent*)GetInputDS("StEvent");
    if (!event) {
        LOG_WARN << "StFwdFttClusterQAMaker::Make - no StEvent" << endm;
        return kStOK;
    }
    StFttCollection* col = event->fttCollection();
    if (!col) {
        if (mDebug) LOG_INFO << "StFwdFttClusterQAMaker::Make - no StFttCollection" << endm;
        return kStOK;
    }

    // StFttDb gives the canonical rob mapping; fall back to quad + plane*4.
    mFttDb = static_cast<StFttDb*>(GetDataSet("fttDb"));
    if (!mFttDb) mFttDb = static_cast<StFttDb*>(GetDataSet("fttDbMkr"));

    LOG_INFO << "FTTQADBG numberOfClusters=" << col->numberOfClusters()
             << " numberOfPoints=" << col->numberOfPoints()
             << " numberOfRawHits=" << col->numberOfRawHits() << endm;

    // clear and (re)fill the per-event cluster store
    for (int r = 0; r < 16; ++r)
        for (int o = 0; o < 5; ++o) mClu[r][o].clear();

    for (StFttCluster* clu : col->clusters()) {
        if (!clu) continue;
        int o = (int)clu->orientation();
        if (o < 0 || o > 3) continue;                 // skip unknown orientation
        int rob = mFttDb ? (int)mFttDb->rob(clu)
                         : (int)clu->quadrant() + (int)clu->plane() * 4;
        if (rob < 0 || rob > 15) continue;
        CluInfo info;
        info.row = (int)clu->row();
        info.x   = clu->x();
        mClu[rob][o].push_back(info);
    }

    processEvent();
    nEventsProcessed++;
    return kStOK;
}

//______________________________________________________________________________
int StFwdFttClusterQAMaker::matchDiagonal(double intercept, const std::vector<CluInfo>& diag, double& residual) {
    int    best = -1;
    double bestDist = 1e9;
    for (size_t j = 0; j < diag.size(); ++j) {
        double dval = diag[j].x * TMath::Sqrt2();
        double diff = intercept - dval;
        if (TMath::Abs(diff) < kDiagWin && TMath::Abs(diff) < bestDist) {
            bestDist = TMath::Abs(diff);
            best = (int)j;
            residual = diff;
        }
    }
    return best;
}

//______________________________________________________________________________
void StFwdFttClusterQAMaker::processEvent() {
    Long64_t evtHV = 0, evtHD = 0, evtVD = 0, evtHVconf = 0, evtHVghost = 0;

    for (int rob = 0; rob < 16; ++rob) {
        int plane = rob / 4;
        nRobInstances++;

        std::vector<CluInfo>& Hs  = mClu[rob][kFttHorizontal];
        std::vector<CluInfo>& Vs  = mClu[rob][kFttVertical];
        std::vector<CluInfo>& DHs = mClu[rob][kFttDiagonalH];
        std::vector<CluInfo>& DVs = mClu[rob][kFttDiagonalV];

        int nH = (int)Hs.size(), nV = (int)Vs.size();
        int nDH = (int)DHs.size(), nDV = (int)DVs.size();
        int nD = nDH + nDV;

        // merged diagonal list for the H/V x D pairings
        std::vector<CluInfo> Ds;
        Ds.reserve(nD);
        Ds.insert(Ds.end(), DHs.begin(), DHs.end());
        Ds.insert(Ds.end(), DVs.begin(), DVs.end());

        bool hasH = nH > 0, hasV = nV > 0, hasD = nD > 0;
        int nDir = (int)hasH + (int)hasV + (int)hasD;
        int subset = (hasH ? 1 : 0) | (hasV ? 2 : 0) | (hasD ? 4 : 0);

        fillPA("nDir", plane, nDir);
        fillPA("subset", plane, subset);
        fillPA("nH", plane, nH);
        fillPA("nV", plane, nV);
        fillPA("nD", plane, nD);
        fillPA("nDH", plane, nDH);
        fillPA("nDV", plane, nDV);

        cnt_nDir[plane][nDir]++; cnt_nDir[4][nDir]++;
        sum_nH[plane] += nH; sum_nH[4] += nH;
        sum_nV[plane] += nV; sum_nV[4] += nV;
        sum_nD[plane] += nD; sum_nD[4] += nD;

        int cHV = 0, cHVconf = 0, cHVghost = 0, cHD = 0, cVD = 0;

        // ---- H x V (+ diagonal confirmation = triple) ----
        for (int iv = 0; iv < nV; ++iv) {
            double x = Vs[iv].x; int rowx = Vs[iv].row;
            if (x < 1e-5 || x > 1e5) continue;
            for (int ih = 0; ih < nH; ++ih) {
                double y = Hs[ih].x; int rowy = Hs[ih].row;
                if (y < 1e-5 || y > 1e5) continue;
                if (!inAnyStripGroup(rowx, rowy, x, y)) continue;
                cHV++;
                h2(Form("pointmapHV_p%d", plane))->Fill(x, y);

                double intercept = x + y;
                double residual = 0;
                int j = -1, diagKind = -1; // 0=DV, 1=DH
                if (x > y) {
                    j = matchDiagonal(intercept, DVs, residual); diagKind = 0;
                    if (j < 0) { j = matchDiagonal(intercept, DHs, residual); diagKind = 1; }
                } else {
                    j = matchDiagonal(intercept, DHs, residual); diagKind = 1;
                    if (j < 0) { j = matchDiagonal(intercept, DVs, residual); diagKind = 0; }
                }
                if (j >= 0) {
                    cHVconf++;
                    fillPA("residual", plane, residual);
                    h2(Form("pointmap_p%d", plane))->Fill(x, y);
                    getHist("diagType")->Fill(diagKind == 0 ? 0 : 1);
                    if (diagKind == 0) cnt_diagV++; else cnt_diagH++;
                } else {
                    cHVghost++;
                }

                // residual to ALL diagonals (window visualization)
                for (int d = 0; d < nD; ++d)
                    fillPA("residualWide", plane, intercept - Ds[d].x * TMath::Sqrt2());
            }
        }

        // ---- H x D ----
        for (int ih = 0; ih < nH; ++ih) {
            double y = Hs[ih].x; int rowy = Hs[ih].row;
            if (y < 1e-5 || y > 1e5) continue;
            for (int d = 0; d < nD; ++d) {
                double I = Ds[d].x * TMath::Sqrt2();
                double x = I - y;
                if (x < 1e-5 || x > 1e5) continue;
                if (inValidStripGroup(x, y, rowy, /*knownIsV=*/false)) {
                    cHD++;
                    h2(Form("pointmapHD_p%d", plane))->Fill(x, y);
                }
            }
        }

        // ---- V x D ----
        for (int iv = 0; iv < nV; ++iv) {
            double x = Vs[iv].x; int rowx = Vs[iv].row;
            if (x < 1e-5 || x > 1e5) continue;
            for (int d = 0; d < nD; ++d) {
                double I = Ds[d].x * TMath::Sqrt2();
                double y = I - x;
                if (y < 1e-5 || y > 1e5) continue;
                if (inValidStripGroup(x, y, rowx, /*knownIsV=*/true)) {
                    cVD++;
                    h2(Form("pointmapVD_p%d", plane))->Fill(x, y);
                }
            }
        }

        fillPA("cHV", plane, cHV);
        fillPA("cHD", plane, cHD);
        fillPA("cVD", plane, cVD);
        fillPA("cHVconf", plane, cHVconf);
        fillPA("cHVghost", plane, cHVghost);

        sum_cHV[plane] += cHV; sum_cHV[4] += cHV;
        sum_cHD[plane] += cHD; sum_cHD[4] += cHD;
        sum_cVD[plane] += cVD; sum_cVD[4] += cVD;
        sum_cHVconf[plane] += cHVconf; sum_cHVconf[4] += cHVconf;
        sum_cHVghost[plane] += cHVghost; sum_cHVghost[4] += cHVghost;

        evtHV += cHV; evtHD += cHD; evtVD += cVD;
        evtHVconf += cHVconf; evtHVghost += cHVghost;
    }

    getHist("cHV_evt")->Fill(evtHV);
    getHist("cHD_evt")->Fill(evtHD);
    getHist("cVD_evt")->Fill(evtVD);
    getHist("cHVconf_evt")->Fill(evtHVconf);
    getHist("cHVghost_evt")->Fill(evtHVghost);
}

//______________________________________________________________________________
void StFwdFttClusterQAMaker::printSummary() {
    LOG_INFO << "=================== StFwdFttClusterQAMaker summary ===================" << endm;
    LOG_INFO << Form("Events processed      : %lld", nEventsProcessed) << endm;
    LOG_INFO << Form("ROB instances (16/evt): %lld", nRobInstances) << endm;
    LOG_INFO << " plane |   <nH>   <nV>   <nD> | f(>=1) f(>=2) f(=3) | <cHV> <cHD> <cVD>/evt | HVconf HVghost ghostFrac" << endm;
    double evt = (nEventsProcessed > 0) ? (double)nEventsProcessed : 1.0;
    for (int p = 0; p < kNP; ++p) {
        Long64_t robs = 0; for (int d = 0; d < 4; ++d) robs += cnt_nDir[p][d];
        double denom = (robs > 0) ? (double)robs : 1.0;
        double f1 = (cnt_nDir[p][1] + cnt_nDir[p][2] + cnt_nDir[p][3]) / denom;
        double f2 = (cnt_nDir[p][2] + cnt_nDir[p][3]) / denom;
        double f3 = (cnt_nDir[p][3]) / denom;
        double meanNH = sum_nH[p] / denom, meanNV = sum_nV[p] / denom, meanND = sum_nD[p] / denom;
        Long64_t hvtot = sum_cHVconf[p] + sum_cHVghost[p];
        double ghostFrac = (hvtot > 0) ? (double)sum_cHVghost[p] / (double)hvtot : 0.0;
        TString s = (p < 4) ? Form("  %d  ", p) : " ALL ";
        LOG_INFO << Form(" %s | %6.3f %6.3f %6.3f | %5.1f%% %5.1f%% %4.1f%% | %5.2f %5.2f %5.2f | %6lld %6lld   %5.1f%%",
                         s.Data(), meanNH, meanNV, meanND, 100*f1, 100*f2, 100*f3,
                         sum_cHV[p]/evt, sum_cHD[p]/evt, sum_cVD[p]/evt,
                         sum_cHVconf[p], sum_cHVghost[p], 100*ghostFrac) << endm;
    }
    Long64_t robs = 0; for (int d = 0; d < 4; ++d) robs += cnt_nDir[4][d];
    double denom = (robs > 0) ? (double)robs : 1.0;
    LOG_INFO << "Key metric (independence test, ALL planes):" << endm;
    LOG_INFO << Form("  ROBs with >=2 of {H,V,D} : %6.2f%%   (overlap points possible)", 100*(cnt_nDir[4][2]+cnt_nDir[4][3])/denom) << endm;
    LOG_INFO << Form("  ROBs with exactly 1 dir  : %6.2f%%   (points impossible)", 100*cnt_nDir[4][1]/denom) << endm;
    LOG_INFO << Form("  ROBs with >=1 dir  : %6.2f%%   (points possible)", 100*(cnt_nDir[4][1]+cnt_nDir[4][2]+cnt_nDir[4][3])/denom) << endm;
    LOG_INFO << Form("  ROBs with 0 directions   : %6.2f%%", 100*cnt_nDir[4][0]/denom) << endm;
    Long64_t hvtot = sum_cHVconf[4] + sum_cHVghost[4];
    double dHV = (hvtot > 0) ? (double)hvtot : 1.0;
    LOG_INFO << "Ghost estimate (is_GroupN-valid HxV pairs, ALL planes):" << endm;
    LOG_INFO << Form("  total valid HxV pairs : %lld", hvtot) << endm;
    LOG_INFO << Form("  confirmed by diagonal : %lld (%.1f%%)", sum_cHVconf[4], 100*sum_cHVconf[4]/dHV) << endm;
    LOG_INFO << Form("  unconfirmed (ghosts)  : %lld (%.1f%%)", sum_cHVghost[4], 100*sum_cHVghost[4]/dHV) << endm;
    Long64_t dtot = cnt_diagV + cnt_diagH;
    double dd = (dtot > 0) ? (double)dtot : 1.0;
    LOG_INFO << Form("  confirming diag type  : DiagV %.1f%%, DiagH %.1f%%", 100*cnt_diagV/dd, 100*cnt_diagH/dd) << endm;
    LOG_INFO << "=====================================================================" << endm;
}

//______________________________________________________________________________
int StFwdFttClusterQAMaker::Finish() {
    printSummary();

    TDirectory* prevDir = gDirectory;
    TFile* fOut = new TFile(mOutputFilename.Data(), "RECREATE");
    fOut->cd();
    for (auto& nh : mHists) {
        nh.second->SetDirectory(gDirectory);
        nh.second->Write();
    }
    fOut->Close();
    gDirectory = prevDir;

    LOG_INFO << "StFwdFttClusterQAMaker wrote " << mOutputFilename.Data() << endm;
    return kStOK;
}
