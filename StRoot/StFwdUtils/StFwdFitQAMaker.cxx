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
StFwdFitQAMaker::StFwdFitQAMaker() : StMaker("fwdFitQA"){};
int StFwdFitQAMaker::Finish() {

    auto prevDir = gDirectory;

    // output file name
    string name = "StFwdFitQAMaker.root";
    TFile *fOutput = new TFile(name.c_str(), "RECREATE");
    fOutput->cd();


    // accumulated analysis
    addHist( (TH1*)(getHist( "RcMatchedMcEta" )->Clone( "EffEta" )) );
    getHist("EffEta")->Divide( (TH1*)getHist( "McEta" ) );

    addHist( (TH1*)(getHist( "RcMatchedMcPt" )->Clone( "EffPt" )) );
    getHist("EffPt")->Divide( (TH1*)getHist( "McPt" ) );

    addHist( (TH1*)(getHist( "RcMatchedMcPhi" )->Clone( "EffPhi" )) );
    getHist("EffPhi")->Divide( (TH1*)getHist( "McPhi" ) );


        // 3 FST hits on MC Track
    addHist( (TH1*)(getHist( "RcMatched3FSTMcEta" )->Clone( "Eff3FSTEta" )) );
    getHist("Eff3FSTEta")->Divide( (TH1*)getHist( "McEta3FST" ) );

    addHist( (TH1*)(getHist( "RcMatched3FSTMcPt" )->Clone( "Eff3FSTPt" )) );
    getHist("Eff3FSTPt")->Divide( (TH1*)getHist( "McPt3FST" ) );

    addHist( (TH1*)(getHist( "RcMatched3FSTMcPhi" )->Clone( "Eff3FSTPhi" )) );
    getHist("Eff3FSTPhi")->Divide( (TH1*)getHist( "McPhi3FST" ) );

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

    addHist( new TH1F("RcMatched3FSTMcPt", ";RcMatched3FSTMcPt; counts", 100, 0, 5) );
    addHist( new TH1F("RcMatched3FSTMcEta", ";RcMatched3FSTMcEta; counts", 50, 0, 5) );
    addHist( new TH1F("RcMatched3FSTMcPhi", ";RcMatched3FSTMcPhi; counts", 32, -2*3.1415926, 2*3.1415926) );

    addHist( new TH1F("RcQOverP", ";RcQ/RcP; counts", 2000, -1, 1) );

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

    addHist( new TH1F( "TrackStats", "Track Stats; Cat;Counts", 10, 0, 10 )  );
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
    LOG_DEBUG << "Processing Fwd Track Fit QA took: " << (FwdTrackerUtils::nowNanoSecond() - itStart) * 1e6 << " ms" << endm;
    return kStOK;
} // Make
//________________________________________________________________________
void StFwdFitQAMaker::Clear(const Option_t *opts) { LOG_DEBUG << "StFwdFitQAMaker::CLEAR" << endm; }
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

        if ( mct.numFST < 3 && mct.numFTT < 4 ) continue;
        // if ( mct.numFTT < 4 ) continue;
        getHist( "TrackStats" )->Fill( 8 );


        getHist( "RcMatched3FSTMcEta" ) ->Fill( mct.p.Eta() );
        getHist( "RcMatched3FSTMcPt" )  ->Fill( mct.p.Pt() );
        getHist( "RcMatched3FSTMcPhi" ) ->Fill( mct.p.Phi() );
        getHist( "RcPID" )              ->Fill( mct.pid );
        ((TH2*)getHist( "RcQMcQ" ))     ->Fill( mct.q, fwdTrack->charge() );
        double qop = fwdTrack->charge() / fwdTrack->momentum().mag();
        getHist( "RcQOverP" )           ->Fill( qop );
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
        
        if ( fwdTrack->isPrimary() ){
            getHist( "curveResolutionPrim" )->Fill( dCurve );
        } else {
            getHist( "curveResolutionGlobal" )->Fill( dCurve );
        }
    }


}