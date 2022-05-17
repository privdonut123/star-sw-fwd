// \class StFcsTrackMatchMaker
// \author Akio Ogawa
//
//  $Id: StFcsTrackMatchMaker.cxx,v 1.1 2021/03/30 13:34:15 akio Exp $
//  $Log: StFcsTrackMatchMaker.cxx,v $
//  Revision 1.1  2021/03/30 13:34:15  akio
//  Moved from $CVSROOT/offline/upgrade/akio/ to $CVSROOT/StRoot/StSpinPool/

#include "StFcsTrackMatchMaker.h"
#include "StEvent/StFwdTrack.h"
#include "StEvent/StEnumerations.h"
#include "StMessMgr.h"
#include "Stypes.h"
#include "StEventTypes.h"

#include "StThreeVectorF.hh"
#include "StFcsDbMaker/StFcsDb.h"
#include "StRoot/StEpdUtil/StEpdGeom.h"

//#include "TApplication.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"

ClassImp(StFcsTrackMatchMaker)

StFcsTrackMatchMaker::StFcsTrackMatchMaker(const char* name): StMaker(name) {}

StFcsTrackMatchMaker::~StFcsTrackMatchMaker(){}

int StFcsTrackMatchMaker::Init(){  
    mFcsDb=static_cast<StFcsDb*>(GetDataSet("fcsDb"));  
    if(!mFcsDb){
	LOG_ERROR  << "StFcsTrackMatchMaker::InitRun Failed to get StFcsDb" << endm;
	return kStFatal;
    }
    mEpdgeo=new StEpdGeom;
    if(mFilename){
	mFile = new TFile(mFilename,"RECREATE");
	mHdx[0] = new TH1F("dx_EcalTrk","dx Ecal-Track",100,-100,100);
	mHdy[0] = new TH1F("dy_EcalTrk","dy Ecal-Track",100,-100,100);
	mHdr[0] = new TH1F("dr_EcalTrk","dr Ecal-Track",100,-100,100);
	mNtrk[0]= new TH1F("NTrk_Ecal","NTrk_Ecal",10,0.0,10.0);
	mNclu[0]= new TH1F("NEcalClu_Trk","NEcalClu_Trk",10,0.0,10.0);
	mHdx[1] = new TH1F("dx_HcalTrk","dx Hcal-Track",100,-100,100);
	mHdy[1] = new TH1F("dy_HcalTrk","dy Hcal-Track",100,-100,100);	
	mHdr[1] = new TH1F("dr_EcalTrk","dr Ecal-Track",100,-100,100);
	mNtrk[1]= new TH1F("NTrk_Hcal","NTrk_Hcal",10,0.0,10.0);
	mNclu[1]= new TH1F("NHcalClu_Trk","NHcalClu_Trk",10,0.0,10.0);
    }
    return kStOK;
}

int StFcsTrackMatchMaker::Finish(){
    if(mFile){
	LOG_INFO << "Closing "<<mFilename<<endm;
	mFile->Close();
    }
    return kStOK;
}

int StFcsTrackMatchMaker::Make(){
    StEvent* event = (StEvent*)GetInputDS("StEvent");
    if(!event) {
	LOG_ERROR << "StFcsTrackMatchMaker::Make did not find StEvent"<<endm; 
	return kStErr;
    }
    mFcsColl = event->fcsCollection();
    if(!mFcsColl) {
	LOG_ERROR << "StFcsTrackMatchMaker::Make did not find StEvent->StFcsCollection"<<endm; 
	return kStErr;
    }
    mFwdTrkColl = event->fwdTrackCollection();
    if(mFwdTrkColl) {
	LOG_ERROR << "StFcsTrackMatchMaker::Make did not find StEvent->fwdTrackCollection"<<endm; 
	return kStErr;
    }
    
    int ntrk=mFwdTrkColl->numberOfTracks();
    int nMatch[2]={0,0};
    for(int itrk=0; itrk<ntrk; itrk++){
	StFwdTrack* trk=mFwdTrkColl->tracks()[itrk];
	//north or south from track
	int ns=0;
	if(trk->mProjections[1].XYZ.x()>0.0) ns=1;
	//Look for a Ecal & Hcal match for a track
	for(int ehp=0; ehp<2; ehp++){
	    const int iproj[3]={1,2,0}; //ecal, hcal, pres
	    StThreeVectorF proj = trk->mProjections[iproj[ehp]].XYZ;
	    int det = ns*2+ehp;
	    int nclu  = mFcsColl->numberOfClusters(det);
	    for(int iclu=0; iclu<nclu; iclu++){
		StFcsCluster* clu=mFcsColl->clusters(det)[iclu];
		float energy=clu->energy();
		if(energy>mMinEnergy){
		    StThreeVectorD xyz=mFcsDb->getStarXYZfromColumnRow(det,clu->x(),clu->y());
		    double dx = xyz.x() - proj.x();
		    double dy = xyz.y() - proj.y();
		    double dr = sqrt(dx*dx + dy*dy);
		    if(mFile){
			if(fabs(dy)<mMaxDistance) mHdx[ehp]->Fill(dx);
			if(fabs(dx)<mMaxDistance) mHdy[ehp]->Fill(dy);
			mHdr[ehp]->Fill(dr);
		    }
		    if(dr<mMaxDistance){
			// if(ehp==0){trk->addEcalCluster(clu);}
			// if(ehp==1){trk->addHcalCluster(clu);}
			// clu->addTrack(trk);
			nMatch[ehp]++;
		    }
		}
	    }
	}
    }

    // mNtrk[1]
    // mNclu[1]=

    LOG_INFO << Form("NTrack=%3d NEcalMatch=%3d NHcalMatch=%3d",ntrk,nMatch[0],nMatch[1])<<endm;
    return kStOK;
}
