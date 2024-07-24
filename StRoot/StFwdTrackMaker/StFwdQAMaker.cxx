#include "StFwdTrackMaker/StFwdQAMaker.h"
#include "StFwdQAMaker.h"
#include "St_base/StMessMgr.h"
#include "StBFChain/StBFChain.h"
#include "StFwdTrackMaker/StFwdTrackMaker.h"

#include "StFwdTrackMaker/include/Tracker/FwdTracker.h"
#include "StFwdTrackMaker/include/Tracker/ObjExporter.h"
// StEvent includes
#include "StEvent/StBTofCollection.h"
#include "StEvent/StBTofHeader.h"
#include "StEvent/StEvent.h"
#include "StEvent/StFttCluster.h"
#include "StEvent/StFttCollection.h"
#include "StEvent/StFcsCluster.h"
#include "StEvent/StFcsCollection.h"
#include "StFcsDbMaker/StFcsDb.h"
#include "StRoot/StEpdUtil/StEpdGeom.h"
#include "StEvent/StFwdTrackCollection.h"
#include "StEvent/StFwdTrack.h"


#include "StMuDSTMaker/COMMON/StMuDstMaker.h"
#include "StMuDSTMaker/COMMON/StMuDst.h"
#include "StMuDSTMaker/COMMON/StMuFstCollection.h"
#include "StMuDSTMaker/COMMON/StMuFstHit.h"
#include "StMuDSTMaker/COMMON/StMuPrimaryVertex.h"
#include "StMuDSTMaker/COMMON/StMuFwdTrack.h"
#include "StMuDSTMaker/COMMON/StMuFwdTrackCollection.h"
#include "StMuDSTMaker/COMMON/StMuFcsCollection.h"
#include "StMuDSTMaker/COMMON/StMuFcsCluster.h"
#include "StMuDSTMaker/COMMON/StMuFcsHit.h"
#include "StMuDSTMaker/COMMON/StMuFttCluster.h"
#include "StMuDSTMaker/COMMON/StMuFttPoint.h"

/** Clear the FwdTreeData from one event to next */
void FwdTreeData::clear(){
    header.clear();
    fttSeeds.reset();
    fttPoints.reset();
    fttClusters.reset();
    fstSeeds.reset();
    reco.reset();
    seeds.reset();
    wcal.reset();
    hcal.reset();
}



StFwdQAMaker::StFwdQAMaker() : StMaker("fwdQAMaker"), mTreeFile(nullptr), mTree(nullptr) {

}

int StFwdQAMaker::Init() {

    mTreeFile = new TFile("fwdtree.root", "RECREATE");
    mTree = new TTree("fwd", "fwd tracking tree");

    mTree->Branch("header",           &mTreeData. header, 3200, 99 );
    mTree->Branch("nSeedTracks",      &mTreeData.nSeedTracks, "nSeedTracks/I");
    mTreeData.fstSeeds.createBranch(mTree, "fst");
    mTreeData.fttSeeds.createBranch(mTree, "ftt");
    mTreeData.fttPoints.createBranch(mTree, "fttPoints");
    mTreeData.fttClusters.createBranch(mTree, "fttClusters");
    mTreeData.wcal.createBranch(mTree, "wcalClusters");
    mTreeData.hcal.createBranch(mTree, "hcalClusters");

    mTreeData.wcalHits.createBranch(mTree, "wcalHits");
    mTreeData.hcalHits.createBranch(mTree, "hcalHits");

    mTreeData.reco.createBranch(mTree, "reco");
    mTreeData.seeds.createBranch(mTree, "seeds");
    return kStOk;
}
int StFwdQAMaker::Finish() {

    if ( mTreeFile && mTree ){
        mTreeFile->cd();
        mTree->Write();
        mTreeFile->Write();
        LOG_DEBUG << "StFwdQA File written" << endm;
    }
    return kStOk;
}
int StFwdQAMaker::Make() {
    LOG_INFO << "SETUP START" << endm;
    // setup the datasets / makers
    mMuDstMaker = (StMuDstMaker *)GetMaker("MuDst");
    if(mMuDstMaker) {
        mMuDst = mMuDstMaker->muDst();
        mMuForwardTrackCollection = mMuDst->muFwdTrackCollection();
        mMuFcsCollection = mMuDst->muFcsCollection();
        if (mMuForwardTrackCollection){
            LOG_DEBUG << "Number of StMuFwdTracks: " << mMuForwardTrackCollection->numberOfFwdTracks() << endm;
        }
    } else {
        LOG_DEBUG << "No StMuDstMaker found: " << mMuDstMaker << endm;
    }
    mFcsDb = static_cast<StFcsDb *>(GetDataSet("fcsDb"));

    mFwdTrackMaker = (StFwdTrackMaker*) GetMaker( "fwdTrack" );
    if (!mFwdTrackMaker) {
        LOG_WARN << "No StFwdTrackMaker found, skipping StFwdQAMaker" << endm;
        // return kStOk;
    }
    // Event header info from stevent
    StEvent *mStEvent = static_cast<StEvent *>(GetInputDS("StEvent"));
    LOG_DEBUG << "SETUP COMPLETE" << endm;
    // Get the primary vertex used by the FWD Tracker
    auto eventPV = mFwdTrackMaker->GetEventPrimaryVertex();
    LOG_DEBUG << "HEADER COMPLETE" << endm;

    FillTracks();
    FillFttClusters();
    FillFcsStMuDst();
    mTree->Fill();
    return kStOk;
}
void StFwdQAMaker::Clear(const Option_t *opts) {
    mTreeData.clear();
    return;
}

void StFwdQAMaker::FillTracks() {
    if ( mMuForwardTrackCollection ){
        LOG_DEBUG << "Adding " << mMuForwardTrackCollection->numberOfFwdTracks() << " FwdTracks (MuDst)" << endm;
        for ( size_t iTrack = 0; iTrack < mMuForwardTrackCollection->numberOfFwdTracks(); iTrack++ ){
            auto muTrack = mMuForwardTrackCollection->getFwdTrack(iTrack);
            mTreeData.reco.add( muTrack );
            if ( iTrack > 5000 ) {
                LOG_WARN << "Truncating to 5000 tracks" << endm;
                break;
            }
        }
    }
    LOG_DEBUG << "TRACKS COMPLETE" << endm;
}

void StFwdQAMaker::FillFcsStMuDst( ) {

    if ( !mMuDst ){
        LOG_DEBUG << "No mMuDst found, skipping StFwdQAMaker::FillFcsStEvent" << endm;
        return;
    }
    StMuFcsCollection* fcs = mMuDst->muFcsCollection();
    if ( !fcs ){
        LOG_DEBUG << "No muFcsCollection found, skipping StFwdQAMaker::FillFcsStEvent" << endm;
        return;
    }

    StEpdGeom epdgeo;
    // LOAD ECAL / HCAL CLUSTERS
    LOG_INFO << "MuDst has #fcs clusters: " << fcs->numberOfClusters() << endm;
    for( size_t i = 0; i < fcs->numberOfClusters(); i++){
        StMuFcsCluster * clu = fcs->getCluster(i);
        if ( clu->detectorId() == kFcsEcalNorthDetId || clu->detectorId() == kFcsEcalSouthDetId ){
            LOG_DEBUG << "Adding WCAL Cluster to FwdTree" << endm;
            mTreeData.wcal.add( clu );
        } else if ( clu->detectorId() == kFcsHcalNorthDetId || clu->detectorId() == kFcsHcalSouthDetId ){
            LOG_DEBUG << "Adding HCAL Cluster to FwdTree" << endm;
            mTreeData.hcal.add( clu );
        }
    }

    // LOAD ECAL / HCAL CLUSTERS
    LOG_INFO << "MuDst has #fcs hits: " << fcs->numberOfHits() << endm;
    for( size_t i = 0; i < fcs->numberOfHits(); i++){
        StMuFcsHit * hit = fcs->getHit(i);
        if ( hit->detectorId() == kFcsEcalNorthDetId || hit->detectorId() == kFcsEcalSouthDetId ){
            LOG_DEBUG << "Adding WCAL Cluster to FwdTree" << endm;
            mTreeData.wcalHits.add( hit );
        } else if ( hit->detectorId() == kFcsHcalNorthDetId || hit->detectorId() == kFcsHcalSouthDetId ){
            LOG_DEBUG << "Adding HCAL Cluster to FwdTree" << endm;
            mTreeData.hcalHits.add( hit );
        }
    }
}


void StFwdQAMaker::FillFttClusters(){

    auto muFttCollection = mMuDst->muFttCollection();
    if ( muFttCollection ){
        LOG_DEBUG << "MuDst has #ftt clusters: " << muFttCollection->numberOfClusters() << endm;
        for ( size_t i = 0; i < muFttCollection->numberOfClusters(); i++ ){
            StMuFttCluster * c = muFttCollection->getCluster(i);
            mTreeData.fttClusters.add( c );
        }

        for ( size_t i = 0; i < muFttCollection->numberOfPoints(); i++ ){
            StMuFttPoint * c = muFttCollection->getPoint(i);
            mTreeData.fttPoints.add( c );
        }


    }

    

    // float pz[] = {280.90499, 303.70498, 326.60501, 349.40499};
    // float SCALE = 1.0;
    // TVector3 cp;
    // FwdTreeHit fh;
    // FwdTreeFttCluster ftc;
    // StEvent *event = static_cast<StEvent *>(GetInputDS("StEvent"));
    // if ( !event || !event->fttCollection() ) return;

    // LOG_DEBUG << "FTT RawHits: " << event->fttCollection()->numberOfRawHits() << endm;
    // LOG_DEBUG << "FTT Clusters: " << event->fttCollection()->numberOfClusters() << endm;

    // for ( size_t i = 0; i < event->fttCollection()->numberOfClusters(); i++ ){
    //     StFttCluster* c = event->fttCollection()->clusters()[i];
    //     if ( c->nStrips() < 1 ) continue;
    //     float dw = 0.05, dlh = 60.0, dlv = 60.0;
    //     float mx = 0.0, my = 0.0;
    //     float sx = 1.0, sy = 1.0;


    //     if ( c->quadrant() == kFttQuadrantA ){
    //         mx = 0; my = 0;
    //         sx = 1.0; sy = 1.0;
    //     } else if ( c->quadrant() == kFttQuadrantB ){
    //         mx = 10.16*SCALE; my = 0.0*SCALE;
    //         sy = -1;
    //         dlv = -dlv;

    //     } else if ( c->quadrant() == kFttQuadrantC ){
    //         mx = -10.16*SCALE ; my = -00.0*SCALE;
    //         sx = -1.0; sy = -1.0;
    //         dlh = -dlh; dlv = -dlv;

    //     } else if ( c->quadrant() == kFttQuadrantD ){
    //         sx = -1;
    //         dlh = -dlh;
    //     }

    //     cp.SetZ( -pz[ c->plane() ] * SCALE );
    //     if ( c->orientation() == kFttHorizontal ){
    //         cp.SetY( my + sy * c->x()/10.0 * SCALE );
    //         cp.SetX( mx );
    //     } else if ( c->orientation() == kFttVertical ){
    //         cp.SetX( mx + sx * c->x()/10.0 * SCALE );
    //         cp.SetY( my );
    //     }
    //     // fh.set( cp.X(), cp.Y(), cp.Z(), c->nStrips(), c->quadrant() );
    //     ftc.pos.SetXYZ( cp.X(), cp.Y(), cp.Z() );
    //     ftc.mQuadrant = c->quadrant();
    //     ftc.mRow = c->row();
    //     ftc.mOrientation = c->orientation();
    //     ftc.mPlane = c->plane();
    //     ftc.mId = c->id();
    //     ftc.mNStrips = c->nStrips();
    //     ftc.mX = c->x();
    //     ftc.mSumAdc = c->sumAdc();

    //     mTreeData.fttClusters.add( ftc );
    // }
}
