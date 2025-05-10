#include "StFwdUtils/StFwdFitQAMaker.h"
#include "StFwdTrackMaker/Common.h"

#include "TMath.h"
#include "TVector3.h"

#include <limits>
#include <map>
#include <string>
#include <string>
#include <vector>

#include "StBFChain/StBFChain.h"

#include "StEvent/StEvent.h"
#include "StEvent/StGlobalTrack.h"
#include "StEvent/StHelixModel.h"
#include "StEvent/StPrimaryTrack.h"
#include "StEvent/StRnDHit.h"
#include "StEvent/StRnDHitCollection.h"
#include "StEvent/StTrack.h"
#include "StEvent/StTrackGeometry.h"
#include "StEvent/StTrackNode.h"
#include "StEvent/StPrimaryVertex.h"
#include "StEvent/StEnumerations.h"
#include "StEvent/StTrackDetectorInfo.h"
#include "StEvent/StFttPoint.h"
#include "StEvent/StFcsHit.h"
#include "StEvent/StFcsCluster.h"
#include "StEvent/StFttCollection.h"
#include "StEvent/StFcsCollection.h"
#include "StEvent/StTriggerData.h"
#include "StEvent/StFstHitCollection.h"
#include "StEvent/StFstHit.h"
#include "StEvent/StFwdTrackCollection.h"
#include "StChain/StChainOpt.h"

#include "StEventUtilities/StEventHelper.h"

#include "StRoot/StEpdUtil/StEpdGeom.h"

#include "StMuDSTMaker/COMMON/StMuDstMaker.h"
#include "StMuDSTMaker/COMMON/StMuDst.h"
#include "StMuDSTMaker/COMMON/StMuFwdTrack.h"
#include "StMuDSTMaker/COMMON/StMuFwdTrackCollection.h"


#include "tables/St_g2t_fts_hit_Table.h"
#include "tables/St_g2t_track_Table.h"
#include "tables/St_g2t_vertex_Table.h"
#include "tables/St_g2t_event_Table.h"

#include "StarMagField/StarMagField.h"

#include "St_base/StMessMgr.h"
#include "StarClassLibrary/StPhysicalHelix.hh"
#include "StarClassLibrary/SystemOfUnits.h"


#include "TROOT.h"
#include "TLorentzVector.h"
#include "StEvent/StFwdTrack.h"
#include "StFcsDbMaker/StFcsDb.h"

//________________________________________________________________________
StFwdFitQAMaker::StFwdFitQAMaker() : StMaker("fwdFitQA"), mOutputFilename("StFwdFitQAMaker.root"){};
int StFwdFitQAMaker::Finish() {

    auto prevDir = gDirectory;

    // output file name
    string name = "StFwdFitQAMaker.root";
    TFile *fOutput = new TFile(mOutputFilename.Data(), "RECREATE");
    fOutput->cd();


    // accumulated analysis
    addHist( (TH1*)(getHist( "RcMatchedMcEta" )->Clone( "EffEta" )) );
    getHist("EffEta")->Divide( (TH1*)getHist( "McEta" ) );

    addHist( (TH1*)(getHist( "RcMatchedMcPt" )->Clone( "EffPt" )) );
    getHist("EffPt")->Divide( (TH1*)getHist( "McPt" ) );

    addHist( (TH1*)(getHist( "RcMatchedMcPhi" )->Clone( "EffPhi" )) );
    getHist("EffPhi")->Divide( (TH1*)getHist( "McPhi" ) );


        // 3 FST hits on MC Track
    addHist( (TH1*)(getHist( "RcMatched3FSTMcEtaPrimary" )->Clone( "Eff3FSTEtaPrimary" )) );
    getHist("Eff3FSTEtaPrimary")->Divide( (TH1*)getHist( "McEta3FST" ) );

    addHist( (TH1*)(getHist( "RcMatched3FSTMcPtPrimary" )->Clone( "Eff3FSTPtPrimary" )) );
    getHist("Eff3FSTPtPrimary")->Divide( (TH1*)getHist( "McPt3FST" ) );

    addHist( (TH1*)(getHist( "RcMatched3FSTMcPhiPrimary" )->Clone( "Eff3FSTPhiPrimary" )) );
    getHist("Eff3FSTPhiPrimary")->Divide( (TH1*)getHist( "McPhi3FST" ) );

    addHist( (TH1*)(getHist( "RcMatched3FSTMcEtaGlobal" )->Clone( "Eff3FSTEtaGlobal" )) );
    getHist("Eff3FSTEtaGlobal")->Divide( (TH1*)getHist( "McEta3FST" ) );

    addHist( (TH1*)(getHist( "RcMatched3FSTMcPtGlobal" )->Clone( "Eff3FSTPtGlobal" )) );
    getHist("Eff3FSTPtGlobal")->Divide( (TH1*)getHist( "McPt3FST" ) );

    addHist( (TH1*)(getHist( "RcMatched3FSTMcPhiGlobal" )->Clone( "Eff3FSTPhiGlobal" )) );
    getHist("Eff3FSTPhiGlobal")->Divide( (TH1*)getHist( "McPhi3FST" ) );

    getHist( "TrackStats")->GetXaxis()->SetBinLabel( 1, "MC" );
    getHist( "TrackStats")->GetXaxis()->SetBinLabel( 2, "MCFWD" );
    getHist( "TrackStats")->GetXaxis()->SetBinLabel( 3, "MCAcc" );
    getHist( "TrackStats")->GetXaxis()->SetBinLabel( 4, "RC" );
    getHist( "TrackStats")->GetXaxis()->SetBinLabel( 5, "GL" );
    getHist( "TrackStats")->GetXaxis()->SetBinLabel( 6, "GLMC" );
    getHist( "TrackStats")->GetXaxis()->SetBinLabel( 7, "PR" );
    getHist( "TrackStats")->GetXaxis()->SetBinLabel( 8, "PRMC" );
    getHist( "TrackStats")->GetXaxis()->SetBinLabel( 9, "BEST" );


    for (auto nh : mHists) {
        nh.second->SetDirectory(gDirectory);
        nh.second->Write();
    }

    // restore previous directory
    gDirectory = prevDir;

    LOG_INFO << "Writing StFwdFitQAMaker output" << endm;

    return kStOk;
}

//________________________________________________________________________
int StFwdFitQAMaker::Init() {
    LOG_DEBUG << "StFwdFitQAMaker::Init" << endm;

    addHist( new TH1F("McQ", ";McQ; counts", 3, -2, 1) );
    addHist( new TH1F("McNumFst", ";NumFST; counts", 10, 0, 10) );
    addHist( new TH1F("McNumFtt", ";NumFTT; counts", 10, 0, 10) );
    addHist( new TH1F("McNumHits", ";NumHits; counts", 10, 0, 10) );
    addHist( new TH1F("McPt", ";McPt; counts", 100, 0, 5) );
    addHist( new TH1F("McEta", ";McEta; counts", 50, 0, 5) );
    addHist( new TH1F("McPhi", ";McPhi; counts", 32, -2*3.1415926, 2*3.1415926) );
    addHist( new TH1F("McStartVtx", ";StartVtx; counts", 10, 0, 10) );
    addHist( new TH1F("McPID", ";GEANT PID; counts", 20, 0, 20) );

    addHist( new TH1F("McPt3FST", ";McPt; counts", 100, 0, 5) );
    addHist( new TH1F("McEta3FST", ";McEta; counts", 50, 0, 5) );
    addHist( new TH1F("McPhi3FST", ";McPhi; counts", 32, -2*3.1415926, 2*3.1415926) );

    addHist( new TH1F("RcQ", ";RcQ; counts", 11, -5.5, 5.5) );
    addHist( new TH1F("RcPt", ";RcPt; counts", 100, 0, 5) );
    addHist( new TH1F("RcEta", ";RcEta; counts", 50, 0, 5) );
    addHist( new TH1F("RcPhi", ";RcPhi; counts", 32, -2*3.1415926, 2*3.1415926) );
    addHist( new TH1F("RcPID", "RC;GEANT PID; counts", 20, 0, 20) );

    addHist( new TH1F("RcMatchedMcPt", ";RcMatchedMcPt; counts", 100, 0, 5) );
    addHist( new TH1F("RcMatchedMcEta", ";RcMatchedMcEta; counts", 50, 0, 5) );
    addHist( new TH1F("RcMatchedMcPhi", ";RcMatchedMcPhi; counts", 32, -2*3.1415926, 2*3.1415926) );

    addHist( new TH1F("RcMatched3FSTMcPtPrimary", ";RcMatched3FSTMcPt; counts", 100, 0, 5) );
    addHist( new TH1F("RcMatched3FSTMcEtaPrimary", ";RcMatched3FSTMcEta; counts", 50, 0, 5) );
    addHist( new TH1F("RcMatched3FSTMcPhiPrimary", ";RcMatched3FSTMcPhi; counts", 32, -2*3.1415926, 2*3.1415926) );

    addHist( new TH1F("RcMatched3FSTMcPtGlobal", ";RcMatched3FSTMcPt; counts", 100, 0, 5) );
    addHist( new TH1F("RcMatched3FSTMcEtaGlobal", ";RcMatched3FSTMcEta; counts", 50, 0, 5) );
    addHist( new TH1F("RcMatched3FSTMcPhiGlobal", ";RcMatched3FSTMcPhi; counts", 32, -2*3.1415926, 2*3.1415926) );

    addHist( new TH1F("RcQOverPt", ";RcQ/RcPt; counts", 2000, -1, 1) );

    addHist( new TH1F("RcIdTruth", ";RcIdTruth; counts", 102, -2, 100) );
    addHist( new TH1F( "deltaRelPt", ";(pT_{RC} - pT_{MC})/pT_{MC}; count", 400, -2, 2 ) );
    addHist( new TH1F( "deltaRelInvPt", ";(pT_{RC}^{-1} - pT_{MC}^{-1})/pT_{MC}^{-1}; count", 400, -2, 2 ) );
    addHist( new TH1F( "curveResolutionPrim", ";(C_{RC} - C_{MC})/C_{MC}; count", 400, -2, 2 ) );
    addHist( new TH1F( "curveResolutionGlobal", ";(C_{RC} - C_{MC})/C_{MC}; count", 400, -2, 2 ) );
    addHist( new TH2F( "deltaRelInvPtVsPt", ";pT_{MC};(pT_{RC}^{-1} - pT_{MC}^{-1})/pT_{MC}^{-1}; count", 500, 0, 5, 400, -1, 4 ) );

    addHist( new TH2F("RcQMcQ", ";McQ; RcQ", 5, -2.5, 2.5, 5, -2.5, 2.5) );
    addHist( new TH2F("RcPtMcPt", ";McPt; RcPt", 100, 0.0, 1.0, 100, 0.0, 1.0) );
    addHist( new TH2F("RcPtMcPtPrim", "Primary;McPt; RcPt", 100, 0.0, 1.0, 100, 0.0, 1.0) );
    addHist( new TH2F("RcPtMcPtGlobal", "Global;McPt; RcPt", 100, 0.0, 1.0, 100, 0.0, 1.0) );
    addHist( new TH2F( "curveResolutionVsPtPrim", ";Pt_{MC};(C_{RC} - C_{MC})/C_{MC}; ", 100, 0, 5.0, 400, -2, 2 ) );
    addHist( new TH2F( "curveResolutionVsPtGlobal", ";Pt_{MC};(C_{RC} - C_{MC})/C_{MC}; ", 100, 0, 5.0, 400, -2, 2 ) );

    addHist( new TH1F( "TrackStats", "Track Stats; Cat;Counts", 10, 0, 10 )  );
    addHist( new TH1F( "GlobalMult", "Global Mult; N;Counts", 50, 0, 50 )  );
    addHist( new TH1F( "BeamlineMult", "Beamline Mult; N;Counts", 50, 0, 50 )  );
    addHist( new TH1F( "PrimaryMult", "Primary Mult; N;Counts", 50, 0, 50 )  );
    addHist( new TH1F( "SecondaryMult", "Secondary Mult; N;Counts", 50, 0, 50 )  );
    
    addHist( new TH1F( "GoodGlobalMult", "Global Mult; N;Counts", 50, 0, 50 )  );
    addHist( new TH1F( "GoodBeamlineMult", "Beamline Mult; N;Counts", 50, 0, 50 )  );
    addHist( new TH1F( "GoodPrimaryMult", "Primary Mult; N;Counts", 50, 0, 50 )  );
    addHist( new TH1F( "GoodSecondaryMult", "Secondary Mult; N;Counts", 50, 0, 50 )  );

    addHist( new TH2F("QidVsPt", "QMatrix; Pt; Qid;", 100, 0, 2.0, 2, -0.5, 1.5 ) );
    addHist( new TH2F("QidVsPtPrim", "QMatrix; Pt; Qid;", 100, 0, 5.0, 2, -0.5, 1.5 ) );
    addHist( new TH2F("QidVsPtGlobal", "QMatrix; Pt; Qid;", 100, 0, 5.0, 2, -0.5, 1.5 ) );
    addHist( new TH2F("QMisidVsPtPrim", "QMatrix; Pt; Qid;", 100, 0, 5.0, 2, -0.5, 1.5 ) );
    addHist( new TH2F("QMisidVsPtGlobal", "QMatrix; Pt; Qid;", 100, 0, 5.0, 2, -0.5, 1.5 ) );

    // EPD
    addHist( new TH2F("epdId", ";EPD ID (pp); EPD ID (tt)", 200, -0.5, 199.5, 200, -0.5, 199.5) );
    addHist( new TH1F("epdEnergy", ";EPD nMIPs; counts", 200, 0, 10) );
    addHist( new TH2F("epdEnergyPP", ";EPD ID (pp); EPD nMIPs;",14, -0.5, 13.5, 200, 0, 10) );
    addHist( new TH2F("epdEnergyTT", ";EPD ID (tt); EPD nMIPs;",32, -0.5, 31.5, 200, 0, 10) );

    addHist( new TH2F("epdTrackdXdY", ";dX [cm]; dY [cm]", 400, -100, 100, 400, -100, 100) );
    addHist( new TH2F("epdTrackdRdPhi", ";dX [cm]; dY [cm]", 500, 0, 250, 120, -TMath::Pi(), TMath::Pi()) );

    addHist( new TH2F("epdTrackdXdYMixed", ";dX [cm]; dY [cm]", 400, -100, 100, 400, -100, 100) );
    addHist( new TH2F("epdTrackdRdPhiMixed", ";dX [cm]; dY [cm]", 500, 0, 250, 120, -TMath::Pi(), TMath::Pi()) );

    // EPD vs. Global
    addHist( new TH2F("epdGlobaldXdY", ";dX [cm]; dY [cm]", 400, -100, 100, 400, -100, 100) );
    addHist( new TH2F("epdGlobaldRdPhi", ";dX [cm]; dY [cm]", 500, 0, 250, 120, -TMath::Pi(), TMath::Pi()) );

    addHist( new TH2F("epdGlobaldXdYMixed", ";dX [cm]; dY [cm]", 400, -100, 100, 400, -100, 100) );
    addHist( new TH2F("epdGlobaldRdPhiMixed", ";dX [cm]; dY [cm]", 500, 0, 250, 120, -TMath::Pi(), TMath::Pi()) );

    // EPD vs. Beamline
    addHist( new TH2F("epdBeamlinedXdY", ";dX [cm]; dY [cm]", 400, -100, 100, 400, -100, 100) );
    addHist( new TH2F("epdBeamlinedRdPhi", ";dX [cm]; dY [cm]", 500, 0, 250, 120, -TMath::Pi(), TMath::Pi()) );

    addHist( new TH2F("epdBeamlinedXdYMixed", ";dX [cm]; dY [cm]", 400, -100, 100, 400, -100, 100) );
    addHist( new TH2F("epdBeamlinedRdPhiMixed", ";dX [cm]; dY [cm]", 500, 0, 250, 120, -TMath::Pi(), TMath::Pi()) );

    // EPD vs. Primary
    addHist( new TH2F("epdPrimarydXdY", ";dX [cm]; dY [cm]", 400, -100, 100, 400, -100, 100) );
    addHist( new TH2F("epdPrimarydRdPhi", ";dX [cm]; dY [cm]", 500, 0, 250, 120, -TMath::Pi(), TMath::Pi()) );

    addHist( new TH2F("epdPrimarydXdYMixed", ";dX [cm]; dY [cm]", 400, -100, 100, 400, -100, 100) );
    addHist( new TH2F("epdPrimarydRdPhiMixed", ";dX [cm]; dY [cm]", 500, 0, 250, 120, -TMath::Pi(), TMath::Pi()) );
    return kStOK;
}
//________________________________________________________________________
int StFwdFitQAMaker::Make() {
    LOG_DEBUG << "StFwdFitQAMaker::Make" << endm;
    long long itStart = FwdTrackerUtils::nowNanoSecond();
    if (!mAnalyzeMuDst)
        ProcessFwdTracks();
    else
        ProcessFwdMuTracks();
    
    ProcessData();
    LOG_DEBUG << "Processing Fwd Track Fit QA took: " << (FwdTrackerUtils::nowNanoSecond() - itStart) * 1e6 << " ms" << endm;
    return kStOK;
} // Make
//________________________________________________________________________
void StFwdFitQAMaker::Clear(const Option_t *opts) { LOG_DEBUG << "StFwdFitQAMaker::CLEAR" << endm; }
//________________________________________________________________________
void StFwdFitQAMaker::ProcessData(  ){
    // This is an example of how to process fwd track collection
    LOG_DEBUG << "StFwdFitQAMaker::ProcessData" << endm;
    StEvent *stEvent = static_cast<StEvent *>(GetInputDS("StEvent"));
    if (!stEvent)
        return;
    

    auto mFcsDb = static_cast<StFcsDb*>(GetDataSet("fcsDb"));
    if (!mFcsDb)
        return;
    StEpdGeom epdgeo;

    std::vector<TVector3> mFcsPreHits;

    // LOAD PRESHOWER HITS (EPD)
    for ( int det = 4; det < 6; det ++ ) {

        StSPtrVecFcsHit& hits = stEvent->fcsCollection()->hits(det);
        int nh=stEvent->fcsCollection()->numberOfHits(det);
        LOG_INFO << "Det " << det << " has " << nh << " hits" << endm;
        for ( int i = 0; i < nh; i++ ){
            StFcsHit* hit=hits[i];

            if(det==kFcsPresNorthDetId || det==kFcsPresSouthDetId){ //EPD
                double zepd=375.0;
                int pp,tt,n;
                double x[5],y[5];

                mFcsDb->getEPDfromId(det,hit->id(),pp,tt);
                getHist2( "epdId" )->Fill( pp,tt );
                getHist ( "epdEnergy" )->Fill( hit->energy() );
                getHist2( "epdEnergyPP" )->Fill( pp, hit->energy() );
                getHist2( "epdEnergyTT" )->Fill( tt, hit->energy() );
                if ( hit->energy() < 0.2 ) continue;

                // Get STAR position of EPD tile corners
                epdgeo.GetCorners(100*pp+tt,&n,x,y);
                double x0 = (x[0] + x[1] + x[2] + x[3]) / 4.0;
                double y0 = (y[0] + y[1] + y[2] + y[3]) / 4.0;
                mFcsPreHits.push_back( TVector3( x0, y0, zepd ) );
            } // if det
        } // for i
    } // for det


    StFwdTrackCollection * ftc = stEvent->fwdTrackCollection();
    if (!ftc)
        return;

    size_t nGlobal = 0;
    size_t nBeamline = 0;
    size_t nPrimary = 0;
    size_t nSecondary = 0;

    size_t nGoodGlobal = 0;
    size_t nGoodBeamline = 0;
    size_t nGoodPrimary = 0;
    size_t nGoodSecondary = 0;

    for ( auto fwdTrack : ftc->tracks() ){

        if (fwdTrack->isGlobalTrack() ){
            nGlobal++;
            if ( fwdTrack->chi2() > 0.001 && fwdTrack->chi2() < 6000 ) nGoodGlobal++;
        }
        if (fwdTrack->isBeamLineConstrainedTrack() ){
            nBeamline++;
            if ( fwdTrack->chi2() > 0.001 && fwdTrack->chi2() < 6000 ) nGoodBeamline++;
        }
        if (fwdTrack->isPrimaryTrack() ){
            nPrimary++;
            if ( fwdTrack->chi2() > 0.001 && fwdTrack->chi2() < 6000 ) nGoodPrimary++;
        }
        if (fwdTrack->isFwdVertexConstrainedTrack() ){
            nSecondary++;
            if ( fwdTrack->chi2() > 0.001 && fwdTrack->chi2() < 6000 ) nGoodSecondary++;            
        }

        // get EPD projection
        auto proj = fwdTrack->getProjectionFor(kFcsPresId);
        
        if ( fwdTrack->chi2() < 0.01 ) continue;
        // now loop on EPD hits
        for ( auto epdHit : mFcsPreHits){
            TVector3 track(proj.mXYZ.x(),proj.mXYZ.y(),proj.mXYZ.z());
            double dx = track.X() - epdHit.X();
            double dy = track.Y() - epdHit.Y();
            double dR = TMath::Sqrt(dx*dx + dy*dy);
            double dPhi = track.DeltaPhi(epdHit);

            getHist2( "epdTrackdXdY" )->Fill( dx, dy );
            getHist2( "epdTrackdRdPhi" )->Fill( dR, dPhi );

            if ( fwdTrack->isGlobalTrack() ){
                getHist2( "epdGlobaldXdY" )->Fill( dx, dy );
                getHist2( "epdGlobaldRdPhi" )->Fill( dR, dPhi );
            }
            if ( fwdTrack->isBeamLineConstrainedTrack() ){
                getHist2( "epdBeamlinedXdY" )->Fill( dx, dy );
                getHist2( "epdBeamlinedRdPhi" )->Fill( dR, dPhi );
            }
            if ( fwdTrack->isPrimaryTrack() ){
                getHist2( "epdPrimarydXdY" )->Fill( dx, dy );
                getHist2( "epdPrimarydRdPhi" )->Fill( dR, dPhi );
            }
        }

        // Mixed event
        for ( auto epdHit : mFcsPreHitsLastEvent){
            TVector3 track(proj.mXYZ.x(),proj.mXYZ.y(),proj.mXYZ.z());
            double dx = track.X() - epdHit.X();
            double dy = track.Y() - epdHit.Y();
            double dR = TMath::Sqrt(dx*dx + dy*dy);
            double dPhi = track.DeltaPhi(epdHit);

            getHist2( "epdTrackdXdYMixed" )->Fill( dx, dy );
            getHist2( "epdTrackdRdPhiMixed" )->Fill( dR, dPhi );
        }
        
    }

    mFcsPreHitsLastEvent.clear();
    for ( auto epdHit : mFcsPreHits){
        mFcsPreHitsLastEvent.push_back(epdHit);
    }

    getHist( "GlobalMult" )->Fill( nGlobal );
    getHist( "BeamlineMult" )->Fill( nBeamline );
    getHist( "PrimaryMult" )->Fill( nPrimary );
    getHist( "SecondaryMult" )->Fill( nSecondary );

    getHist( "GoodGlobalMult" )->Fill( nGoodGlobal );
    getHist( "GoodBeamlineMult" )->Fill( nGoodBeamline );
    getHist( "GoodPrimaryMult" )->Fill( nGoodPrimary );
    getHist( "GoodSecondaryMult" )->Fill( nGoodSecondary );

}
//________________________________________________________________________
void StFwdFitQAMaker::ProcessFwdTracks(  ){
    // This is an example of how to process fwd track collection
    LOG_DEBUG << "StFwdFitQAMaker::ProcessFwdTracks" << endm;
    StEvent *stEvent = static_cast<StEvent *>(GetInputDS("StEvent"));
    if (!stEvent)
        return;
    StFwdTrackCollection * ftc = stEvent->fwdTrackCollection();
    if (!ftc)
        return;

    // Get geant tracks
    St_g2t_track *g2t_track = (St_g2t_track *)GetDataSet("geant/g2t_track");

    if (!g2t_track)
        return;

    mcTracks.clear();
    LOG_DEBUG << g2t_track->GetNRows() << " mc tracks in geant/g2t_track " << endm;
    for (int irow = 0; irow < g2t_track->GetNRows(); irow++) {
        g2t_track_st *track = (g2t_track_st *)g2t_track->At(irow);

        if (0 == track)
            continue;

        McFwdTrack mct;

        mct.id = track->id;
        mct.p.SetXYZ( track->p[0], track->p[1], track->p[2] );
        mct.q = track->charge;
        mct.g2track = track;
        mct.pid = (int)track->ge_pid;

        mcTracks.push_back( mct );

        getHist( "McPID" )->Fill( mct.pid );
        getHist( "McEta" )->Fill( mct.p.Eta() );
        getHist( "McPt" )->Fill( mct.p.Pt() );
        getHist( "McPhi" )->Fill( mct.p.Phi() );
        getHist( "McQ" )->Fill( mct.p.Pt() );
        getHist( "TrackStats" )->Fill( 0 );
    } // get MC tracks

    /***********************/
    St_g2t_fts_hit *g2t_fsi_hits = (St_g2t_fts_hit *)GetDataSet("geant/g2t_fsi_hit");
    if ( !g2t_fsi_hits ){
        LOG_DEBUG << "No g2t_fts_hits, cannot load FST hits from GEANT" << endm;
        return;
    }
    int nfsi = g2t_fsi_hits->GetNRows();

    for (int i = 0; i < nfsi; i++) {
        g2t_fts_hit_st *git = (g2t_fts_hit_st *)g2t_fsi_hits->At(i);
        if (0 == git) continue; // geant hit
        int track_id = git->track_p - 1;
        if ( track_id >= 0 && track_id < mcTracks.size() ){
            mcTracks[track_id].numFST++;
        }
    } // Get FST hits on MC tracks

    St_g2t_fts_hit *g2t_stg_hits = (St_g2t_fts_hit *)GetDataSet("geant/g2t_stg_hit");
    if (!g2t_stg_hits){
        LOG_WARN << "geant/g2t_stg_hit is empty" << endm;
        return;
    }
    int nstg = g2t_stg_hits->GetNRows();
    for (int i = 0; i < nstg; i++) {

        g2t_fts_hit_st *git = (g2t_fts_hit_st *)g2t_stg_hits->At(i);
        if (0 == git) continue; // geant hit
        int track_id = git->track_p - 1;
        int volume_id = git->volume_id;

        // only use the hits on the front modules
        if ( volume_id % 2 ==0 )
            continue;
        if ( track_id >= 0 && track_id < mcTracks.size() ){
            mcTracks[track_id].numFTT++;

        }
    }
    for ( auto mct : mcTracks ){
        getHist( "McNumFst" )->Fill( mct.numFST );
        getHist( "McNumFtt" )->Fill( mct.numFTT );
        getHist( "McNumHits" )->Fill( mct.numFTT + mct.numFST );

        if ( mct.numFST >= 1 || mct.numFTT >= 1 ){
            getHist( "TrackStats" )->Fill( 1 );
        }
        if ( mct.numFST >= 3 ){
            getHist( "McEta3FST" )->Fill( mct.p.Eta() );
            getHist( "McPt3FST" )->Fill( mct.p.Pt() );
            getHist( "McPhi3FST" )->Fill( mct.p.Phi() );

        }
        if ( mct.numFST >=3 || mct.numFTT >= 4){
            getHist( "TrackStats" )->Fill( 2 );
        }
    }

    /**************************************************************************
    */
    for ( auto fwdTrack : ftc->tracks() ){
        getHist( "TrackStats" )->Fill( 3 );
        int idx = fwdTrack->idTruth() - 1;
        if ( fwdTrack->vertexIndex() > 100 ){
            getHist( "TrackStats" )->Fill( 4 );
            if ( idx >= 0 && idx < mcTracks.size() )
                getHist( "TrackStats" )->Fill( 5 );
        } else {
            getHist( "TrackStats" )->Fill( 6 );
                if ( idx >= 0 && idx < mcTracks.size() )
                    getHist( "TrackStats" )->Fill( 7 );
        }

        getHist( "RcIdTruth" )->Fill( fwdTrack->idTruth() );
        getHist( "RcEta" )->Fill( fwdTrack->momentum().pseudoRapidity() );
        getHist( "RcPt" )->Fill( fwdTrack->momentum().perp() );
        getHist( "RcPhi" )->Fill( fwdTrack->momentum().phi() );



        LOG_INFO << "idx = " << idx << endm;
        if ( idx < 0 || idx >= mcTracks.size() ) continue;
        McFwdTrack mct = mcTracks[idx];
        // getHist( "TrackStats" )->Fill( 8 );

        getHist( "RcMatchedMcEta" )->Fill( mct.p.Eta() );
        getHist( "RcMatchedMcPt" )->Fill( mct.p.Pt() );
        getHist( "RcMatchedMcPhi" )->Fill( mct.p.Phi() );

        getHist( "RcPtMcPt" )->Fill( mct.p.Pt(), fwdTrack->momentum().perp() );

        if ( fwdTrack->vertexIndex() < 100 ){
            getHist( "RcPtMcPtPrim" )->Fill( mct.p.Pt(), fwdTrack->momentum().perp() );
        } else {
            getHist( "RcPtMcPtGlobal" )->Fill( mct.p.Pt(), fwdTrack->momentum().perp() );
        }

        if ( mct.numFST < 1 || mct.numFTT < 1 ) continue;
        // if ( mct.numFTT < 4 ) continue;
        getHist( "TrackStats" )->Fill( 8 );


        if ( fwdTrack->trackType() == StFwdTrack::kPrimaryVertexConstrained ){
            getHist( "RcMatched3FSTMcEtaPrimary" ) ->Fill( mct.p.Eta() );
            getHist( "RcMatched3FSTMcPtPrimary" )  ->Fill( mct.p.Pt() );
            getHist( "RcMatched3FSTMcPhiPrimary" ) ->Fill( mct.p.Phi() );
        } else {
            getHist( "RcMatched3FSTMcEtaGlobal" ) ->Fill( mct.p.Eta() );
            getHist( "RcMatched3FSTMcPtGlobal" )  ->Fill( mct.p.Pt() );
            getHist( "RcMatched3FSTMcPhiGlobal" ) ->Fill( mct.p.Phi() );
        }
        getHist( "RcPID" )              ->Fill( mct.pid );
        ((TH2*)getHist( "RcQMcQ" ))     ->Fill( mct.q, fwdTrack->charge() );
        double qop = fwdTrack->charge() / fwdTrack->momentum().perp();
        getHist( "RcQOverPt" )           ->Fill( qop );
        getHist( "RcQ" )                ->Fill( fwdTrack->charge() );

        double dEta  = fwdTrack->momentum().pseudoRapidity() - mct.p.Eta();
        // double dPhi  = fwdTrack->momentum().DeltaPhi( mct.p );
        // double dR    = fwdTrack->momentum().DeltaR( mct.p );
        double dPt   = fwdTrack->momentum().perp() - mct.p.Pt();
        double dRPt  = (fwdTrack->momentum().perp() - mct.p.Pt()) / mct.p.Pt();
        double iPt   = 1.0 / fwdTrack->momentum().perp();
        double curve   = 1.0 / fwdTrack->momentum().perp();
        double curveMc = 1.0 / mct.p.Pt();
        double iMcPt = 1.0 / mct.p.Pt();
        double dRIPt = (iPt - iMcPt) / iMcPt;
        double dCurve = (curve - curveMc) / curveMc;

        getHist( "deltaRelPt" ) ->Fill( dRPt );
        getHist( "deltaRelInvPt" )->Fill( dRIPt );
        ((TH2*)getHist( "deltaRelInvPtVsPt" ))->Fill( mct.p.Pt(), dRIPt );
        if ( mct.g2track )
            getHist( "McStartVtx" )->Fill( mct.g2track->start_vertex_p );
        
        if ( fwdTrack->trackType() == StFwdTrack::kPrimaryVertexConstrained ){
            getHist( "curveResolutionPrim" )->Fill( dCurve );
            getHist( "curveResolutionVsPtPrim" )->Fill( mct.p.Pt(), dCurve );
            getHist("QidVsPtPrim")->Fill( mct.p.Pt(), mct.q == fwdTrack->charge() ? 1 : 0 );
            getHist("QMisidVsPtPrim")->Fill( mct.p.Pt(), mct.q == fwdTrack->charge() ? 0 : 1 );
        } else {
            getHist( "curveResolutionGlobal" )->Fill( dCurve );
            getHist( "curveResolutionVsPtGlobal" )->Fill( mct.p.Pt(), dCurve );
            getHist("QidVsPtGlobal")->Fill( mct.p.Pt(), mct.q == fwdTrack->charge() ? 1 : 0 );
            getHist("QMisidVsPtGlobal")->Fill( mct.p.Pt(), mct.q == fwdTrack->charge() ? 0 : 1 );
        }
    }


}
