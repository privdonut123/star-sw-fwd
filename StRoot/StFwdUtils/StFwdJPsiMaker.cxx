#include "StFwdUtils/StFwdJPsiMaker.h"
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
StFwdJPsiMaker::StFwdJPsiMaker() : StMaker("fwdJPsi"){};
int StFwdJPsiMaker::Finish() { 
    
    auto prevDir = gDirectory;
        
    // output file name
    string name = "StFwdJPsiMaker.root";
    TFile *fOutput = new TFile(name.c_str(), "RECREATE");
    fOutput->cd();
    for (auto nh : mHists) {
        nh.second->SetDirectory(gDirectory);
        nh.second->Write();
    }

    // restore previous directory
    gDirectory = prevDir;

    LOG_INFO << "Writing StFwdJPsiMaker output" << endm;

    return kStOk; 
}

//________________________________________________________________________
int StFwdJPsiMaker::Init() { 
    LOG_DEBUG << "StFwdJPsiMaker::Init" << endm; 

    addHist( new TH1F("fwdMultFailed", ";N_{ch}^{FWD}; counts", 100, 0, 100) );

    addHist( new TH1F("Mll", ";Mll; counts/[10 MeV]", 200, 2.0, 4.0) );
    addHist( new TH1F("rc0Mll", ";rc Mll; counts/[10 MeV]", 200, 2.0, 4.0) );
    addHist( new TH1F("rc1Mll", ";rc Mll; counts/[10 MeV]", 200, 2.0, 4.0) );
    addHist( new TH1F("rc2Mll", ";rc Mll; counts/[10 MeV]", 200, 2.0, 4.0) );

    return kStOK;
}
//________________________________________________________________________
int StFwdJPsiMaker::Make() {
    LOG_DEBUG << "StFwdJPsiMaker::Make" << endm;
    StEvent *event = (StEvent *)GetDataSet("StEvent");
    if (event){
        StFttCollection *fttCol = event->fttCollection();
        if (fttCol){
            LOG_INFO << "The Ftt Collection has " << fttCol->numberOfPoints() << " points" << endm;
        }
    }
    long long itStart = FwdTrackerUtils::nowNanoSecond();
    
    ProcessFwdMcTracks();

    // if (!mAnalyzeMuDst)
        ProcessFwdTracks();
    // else 
    //     ProcessFwdMuTracks();
    LOG_DEBUG << "Processing Fwd Tracks took: " << (FwdTrackerUtils::nowNanoSecond() - itStart) * 1e6 << " ms" << endm;
    return kStOK;
} // Make
//________________________________________________________________________
void StFwdJPsiMaker::Clear(const Option_t *opts) { LOG_DEBUG << "StFwdJPsiMaker::CLEAR" << endm; }

void StFwdJPsiMaker::ProcessFwdMcTracks( ){
    // Get geant tracks
    St_g2t_track *g2t_track = (St_g2t_track *)GetDataSet("geant/g2t_track");

    if (!g2t_track)
        return;

    double electronMass = 0.00051099895000;
    double muonMass = 0.1056583755;

    LOG_DEBUG << g2t_track->GetNRows() << " mc tracks in geant/g2t_track " << endm;

    TLorentzVector lv1, lv2, lv;
    for (int i = 0; i < g2t_track->GetNRows(); i++) {
        g2t_track_st *t1 = (g2t_track_st *)g2t_track->At(i);
        if (0 == t1)
            continue;

        int q1 = t1->charge;
        int pid1 = t1->ge_pid;
        float t1pt = t1->pt;
        float t1eta = t1->eta;
        double m1 = electronMass;
        if ( pid1 == 5 || pid1 == 6)
            m1 = muonMass;
        float t1phi = std::atan2(t1->p[1], t1->p[0]); //track->phi;
        lv1.SetPtEtaPhiM( t1pt, t1eta, t1phi, m1 );
        LOG_INFO << TString::Format( "JPSI track pt=%f, eta=%f, phi=%f, id=%d", t1->pt, t1eta, t1phi, t1->id ) << endm;
        

        for (int j = 0; j < g2t_track->GetNRows(); j++) {
            g2t_track_st *t2 = (g2t_track_st *)g2t_track->At(j);
            if (0 == t2 || i == j)
                continue;
            
            
            int q2 = t2->charge; 
            // LOG_INFO << "q1="<<q1 << ", q2="<<q2<<endm;
            
            
            if (q1 == q2) continue;   
            int pid2 = t2->ge_pid;

            // LOG_INFO << "pid1="<<pid1 << ", pid2="<<pid2<<endm;
            // if ( pid1 == 2 && pid2 == 3 || pid1 == 3 && pid2 == 2 ){}
            // else continue;
            if (t1->start_vertex_p != t2->start_vertex_p) continue;

            double m2 = electronMass;
            if ( pid2 == 5 || pid2 == 6)
                m2 = muonMass;

            float t2pt = t2->pt;
            float t2eta = t2->eta;

            if (t1->eta < 0 || t2->eta < 0) continue;
            float t2phi = std::atan2(t2->p[1], t2->p[0]); //track->phi;

            // LOG_INFO << TString::Format( "t1 pt=%f, eta=%f, phi=%f, id=%d", t1->pt, t1eta, t1phi, t1->id ) << endm;
            // LOG_INFO << TString::Format( "t2 pt=%f, eta=%f, phi=%f, id=%d", t2->pt, t2eta, t2phi, t2->id ) << endm;
            lv2.SetPtEtaPhiM( t2pt, t2eta, t2phi, m2 );

            lv = lv1 + lv2;
            // LOG_INFO << "M=" << lv.M() << endm;
            getHist( "Mll" )->Fill( lv.M() ) ;

        }

        LOG_INFO << "JPSI event complete" << endm;

        // int track_id = track->id;
        // float pt2 = track->p[0] * track->p[0] + track->p[1] * track->p[1];
        // float pt = std::sqrt(pt2);
        // float eta = track->eta;
        // float phi = std::atan2(track->p[1], track->p[0]); //track->phi;
        // int q = track->charge;
    }
} // ProcessFwdMcTracks

//________________________________________________________________________
void StFwdJPsiMaker::ProcessFwdTracks(  ){
    // This is an example of how to process fwd track collection
    LOG_DEBUG << "StFwdJPsiMaker::ProcessFwdTracks" << endm;
    StEvent *stEvent = static_cast<StEvent *>(GetInputDS("StEvent"));
    if (!stEvent)
        return;
    StFwdTrackCollection * ftc = stEvent->fwdTrackCollection();
    if (!ftc)
        return;

    LOG_DEBUG << "Checking FcsCollection" << endm;
    StFcsCollection *fcs = stEvent->fcsCollection();
    if (!fcs) return;

    StFcsDb *mFcsDb = static_cast<StFcsDb *>(GetDataSet("fcsDb"));


    LOG_INFO << "FwdTrackCollection has: " << ftc->tracks().size() << " tracks" << endm;
    TLorentzVector lv, lv1, lv2;
    double muonMass = 0.1056583755;
    for ( auto fwdTrack1 : ftc->tracks() ){
        if ( !fwdTrack1->didFitConvergeFully() ) {        
            continue;
        }

        StThreeVectorD p1 = fwdTrack1->momentum();
        lv1.SetXYZM( p1.x(), p1.y(), p1.z(), muonMass );
        for ( auto fwdTrack2 : ftc->tracks() ){
            if ( !fwdTrack2->didFitConvergeFully() ) {        
            continue;
            }
            StThreeVectorD p2 = fwdTrack2->momentum();
            lv2.SetXYZM( p2.x(), p2.y(), p2.z(), muonMass );
            lv = lv1+lv2;

            getHist( "rc0Mll" )->Fill( lv.M() );

            if ( fwdTrack1->charge() == fwdTrack2->charge() ) continue;
            getHist( "rc1Mll" )->Fill( lv.M() );

            if ( fwdTrack1->numberOfFitPoints() < 8 || fwdTrack2->numberOfFitPoints() < 8 ) continue;
            getHist( "rc2Mll" )->Fill( lv.M() );
        }
    } // Loop ftc->tracks()
} // ProcessFwdTracks

//________________________________________________________________________
void StFwdJPsiMaker::ProcessFwdMuTracks(  ){
    // This is an example of how to process fwd track collection
    LOG_DEBUG << "StFwdJPsiMaker::ProcessFwdMuTracks" << endm;
    StMuDstMaker *mMuDstMaker = (StMuDstMaker *)GetMaker("MuDst");
    if(!mMuDstMaker) {
        LOG_WARN << " No MuDstMaker ... bye-bye" << endm;
        return;
    }
    StMuDst *mMuDst = mMuDstMaker->muDst();
    if(!mMuDst) {
        LOG_WARN << " No MuDst ... bye-bye" << endm;
        return;
    }
    StMuFwdTrackCollection * ftc = mMuDst->muFwdTrackCollection();
    if (!ftc) return;

    StMuFcsCollection *fcs = mMuDst->muFcsCollection();
    if (!fcs) return;

    LOG_INFO << "Number of StMuFwdTracks: " << ftc->numberOfFwdTracks() << endl;

    StFcsDb *mFcsDb = static_cast<StFcsDb *>(GetDataSet("fcsDb"));

    size_t fwdMultFST = 0;
    size_t fwdMultEcalMatch = 0;
    size_t fwdMultHcalMatch = 0;

    for ( size_t iTrack = 0; iTrack < ftc->numberOfFwdTracks(); iTrack++ ){
        StMuFwdTrack * muFwdTrack = ftc->getFwdTrack( iTrack );
        // LOG_DEBUG << TString::Format("StMuFwdTrack[ nProjections=%lu, nFTTSeeds=%lu, nFSTSeeds=%lu, mPt=%f ]", muFwdTrack->mProjections.size(), muFwdTrack->mFTTPoints.size(), muFwdTrack->mFSTPoints.size(), muFwdTrack->momentum().Pt()) << endm;

        LOG_DEBUG << "StMuFwdTrack has " << muFwdTrack->mEcalClusters.GetEntries() << " Ecal matched" << endm;
        LOG_DEBUG << "StMuFwdTrack has " << muFwdTrack->mHcalClusters.GetEntries() << " Hcal matched" << endm;
    } // iTrack
}
