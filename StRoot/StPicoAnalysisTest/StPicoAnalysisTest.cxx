#include "StPicoAnalysisTest.h"
#ifndef ST_PICO_ANALYSIS_TEST_h
#include "StRoot/StPicoEvent/StPicoEvent.h"
#include "StRoot/StPicoDstMaker/StPicoDstMaker.h"
#include "StRoot/StPicoEvent/StPicoDst.h"
#endif

#include "TDataSetIter.h"
#include "StDAQMaker/StDAQReader.h"

//#include "StRoot/StEvent/StFcsCollection.h" //have to replace this functionality

#include "StRoot/St_base/StMessMgr.h"
#include "StRoot/StPicoEvent/StPicoFcsHit.h"
#include "StRoot/StPicoEvent/StPicoFcsCluster.h"
#include "StRoot/StPicoEvent/StPicoFwdTrack.h"
#include "StRoot/StFcsDbMaker/StFcsDb.h"

#include "StLorentzVectorD.hh"

#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TString.h"
#include "TFile.h"

#include <string.h>
#include <time.h>

ClassImp(StPicoAnalysisTest)

StPicoAnalysisTest::StPicoAnalysisTest(const char *name, StPicoDstMaker *picoMaker):StMaker(name){
    mPicoDstMaker = picoMaker;
    mPicoDst = 0;

}

StPicoAnalysisTest::~StPicoAnalysisTest(){};
 
Int_t StPicoAnalysisTest::Init(){  
    /*
  mFcsDb = static_cast<StFcsDb*>(GetDataSet("fcsDb"));
  if(!mFcsDb){
    LOG_WARN << "StPicoAnalysisTest::Init cannot get fcsDb" << endm;
    return kStErr;
  }
*/
  if(mFilename){
    LOG_INFO << "StPicoAnalysisTest::Init - Opening "<<mFilename<<endm;
    mFile=new TFile(mFilename,"RECREATE");
  }

  const char* nameCut[mNCut] = {"All","ETOT","HTOT","Cone","SigmaMax","TrackMatch","ChargeSign"};
  
  for(int cut=0; cut<mNCut; cut++){
    /*mhetest[cut]  = new TH1F(Form("Rhetest_%s", nameCut[cut]),Form("hetest_%s",nameCut[cut]),100,0.0,100.0);
    mtrpytest[cut] = new TH1F(Form("Rtrpytest_%s",nameCut[cut]),Form("htrpytest_%s",nameCut[cut]),100,-10.0,10.0);
    mmasstest[cut] = new TH1F(Form("Rmasstest_%s",nameCut[cut]),Form("hmasstest_%s",nameCut[cut]),200,0,20.0);
*/

    mETot[cut]  = new TH1F(Form("RETot_%s", nameCut[cut]),Form("Epair/ETOT_%s",nameCut[cut]),50,0.0,1.1);
    mHTot[cut]  = new TH1F(Form("RHTot_%s", nameCut[cut]),Form("Epair/HTOT_%s",nameCut[cut]),50,0.0,3.0);
    mCone[cut]  = new TH1F(Form("RCone_%s", nameCut[cut]),Form("Epair/Cone_%s (R=%3.1f)",nameCut[cut],mConeR),50,0.0,1.1);
    mSigmax[cut]= new TH1F(Form("Sigmax_%s",nameCut[cut]),Form("SigmaMax_%s",  nameCut[cut]),50,0.0,2.0);
    mPToverET[cut] = new TH1F(Form("EToverPT_%s",nameCut[cut]),Form("EcalET/TrackPT_%s",nameCut[cut]),50,0.0,3.0);
    mChargeSum[cut]= new TH1F(Form("ChargeSum_%s",nameCut[cut]),Form("ChargeSum_%s",nameCut[cut]),5,-2.5,2.5); 

    mET[cut]   = new TH1F(Form("ET_%s"  ,nameCut[cut]),Form("ET_%s"   ,nameCut[cut]),50,0.0,5.0);
    mEZ[cut]   = new TH1F(Form("EZ_%s"  ,nameCut[cut]),Form("EZ_%s"   ,nameCut[cut]),50,0.0,120.0);
    mM[cut]    = new TH1F(Form("M_%s"   ,nameCut[cut]),Form("Mass_%s" ,nameCut[cut]),50,0.0,10.0);

    mZ[cut]    = new TH1F(Form("Z_%s"   ,nameCut[cut]),Form("Zll_%s"  ,nameCut[cut]),50,0.0,1.0);
    mCosT[cut] = new TH1F(Form("CosT_%s",nameCut[cut]),Form("CosT_%s" ,nameCut[cut]),50,-1.0,1.0);
    mPhi[cut]  = new TH1F(Form("Phi_%s" ,nameCut[cut]),Form("Phi_%s"  ,nameCut[cut]),50,-M_PI,M_PI);

    mXFPT[cut] = new TH2F(Form("XFPT_%s",nameCut[cut]),Form("XFPT12_%s;xF;ET",nameCut[cut]),50,0.0,0.5,50,0.0,8.0);
    mET12[cut] = new TH2F(Form("ET12_%s",nameCut[cut]),Form("ET12_%s;ET1;ET2",nameCut[cut]),50,0.0,8.0,50,0.0,8.0);
    mXY[cut]   = new TH2F(Form("XY_%s"  ,nameCut[cut]),Form("XY12_%s;X;Y"    ,nameCut[cut]),50,-130,130,50,-110,110);
    mPTET[cut] = new TH2F(Form("PTET_%s",nameCut[cut]),Form("ETvsPT_%s; ET(Ecal); TrkPT",nameCut[cut]),50,0,8,50,0,8);

  }
  return kStOk;
}

Double_t StPicoAnalysisTest::Pt(Double_t Px, Double_t Py) {
  //cout << "Pt " << sqrt(Px*Px + Py*Py) << endl;
  return sqrt(Px*Px + Py*Py);
}

Double_t StPicoAnalysisTest::Pt(StPicoFcsCluster *clu) {
  return StPicoAnalysisTest::Pt(clu->fourMomentum()[0],clu->fourMomentum()[1]);
}

Double_t StPicoAnalysisTest::Pt(StPicoFwdTrack *trk) {
  return StPicoAnalysisTest::Pt(trk->momentum()[0],trk->momentum()[1]);
}

Double_t StPicoAnalysisTest::Eta(StPicoFwdTrack *trk) {
  //cout << "Eta " << TMath::ATanH(trk->momentum()[2]/(Double_t) trk->momentum().Mag()) << endl;
  return TMath::ATanH(trk->momentum()[2]/(Double_t) trk->momentum().Mag());
}

Double_t StPicoAnalysisTest::Phi(StPicoFwdTrack *trk) {
  //cout << "Phi " << TMath::ACos(trk->momentum()[0]/(Double_t) StPicoAnalysisTest::Pt(trk)) << endl;
  return TMath::ACos(trk->momentum()[0]/(Double_t) StPicoAnalysisTest::Pt(trk));
}

Int_t StPicoAnalysisTest::Make(){
    /*
    TH1::SetDefaultSumw2(true);
    TH2::SetDefaultSumw2(true);
    TH3::SetDefaultSumw2(true);
*/
    if(!mPicoDstMaker) {
        LOG_WARN << " No PicoDstMaker! Skip! " << endm;
        return kStWarn;
    }
    mPicoDst = mPicoDstMaker->picoDst();
    if(!mPicoDst) {
        LOG_WARN << " No PicoDst! Skip! " << endm;
        return kStWarn;
    }
    mPicoEvent = static_cast<StPicoEvent*>(mPicoDst->event());
   //mPicoEvent = dynamic_cast<StPicoEvent *>(GetInputDS("StPicoEvent"));
    if (mPicoEvent) {
        //LOG_WARN << "HOORAY!" << endm;
        //return kStWarn;
    }

    if(!mPicoEvent) {
        LOG_WARN << " No PicoEvent! Skip! " << endm;
        return kStWarn;
    }
    
    cout << "PV: " << mPicoDst->event()->primaryVertex().z() << endl;

    //if the above doesn't work, look at this old line: StEvent *stEvent = static_cast<StEvent *>(GetInputDS("StEvent"));

//loop over tracks with : n = mPicoDst->numberOfTracks();
//in the loop, do mtrack=mPicoDst->track(i); and e.g. mPt=mtrack->pMom().Perp();

 // mFcsCollection=0;
  
  //Isaac's test:                                                                                                                    
  //StEvent *event = (StEvent *)GetDataSet("StEvent");//redundant/unnecessary since unused, but leaving for now.
  /*
  //Isaac's test:
  StEvent *stEvent = static_cast<StEvent *>(GetInputDS("StEvent"));
  if (!stEvent) {
    LOG_INFO << "No StEvent found" << endm;
    return kStErr;
  }
  */
  /*
  StFwdTrackCollection * ftc = stEvent->fwdTrackCollection();
  if (!ftc)
    return;
  */
 /* 
  LOG_INFO << "Checking FcsCollection" << endm;
  mFcsCollection = stEvent->fcsCollection();
  if (!mFcsCollection) {
    LOG_INFO << "No StFcsCollection found" << endm;
    return kStErr;
  }
*/
  //Isaac test:
  if (mPicoDst->numberOfFcsHits()!=0) {
    //LOG_INFO << "have " <<  mPicoDst->numberOfFcsHits() << " clusters" /*<< " " << mPicoDst->numberOfFcsClusters(1)*/ << endm;
  }
  else{
    //LOG_INFO << "no hits" << endm;
  }

  //No clusters found
  if(mPicoDst->numberOfFcsHits()==0) return kStOK;

/*
  //Find highest ET clusters for north and south
  StPicoFcsCluster* candidates[2]={0,0};
  //for(int ns=0; ns<2; ns++){//IM: NEED TO ADD THIS BACK IN
    //ISAAC TEST! //uncomment this if you ever comment out the "no clusters found" check above
    //if (mFcsCollection->numberOfClusters(ns)==0) {continue;} 
   // StPicoFcsCluster& ecal= mPicoDst->FcsCluster();//ns);      //IM: NEED TO ADD THIS BACK IN
    //sort by ET
    std::sort(ecal.begin(), ecal.end(), [](StPicoFcsCluster* a, StPicoFcsCluster* b) {
        return b->fourMomentum().perp() < a->fourMomentum().perp();
      });    
    
    //keep highest
    //comment out "if...mETCut)", if you want no Et cut for tests
    if(ecal[0]->fourMomentum().perp() > mETCut) candidates[ns] = ecal[0];
          
  //}  
  */

//cout << "HITS!" << endl;
 //new
 /*
for (int i = 0; i < mPicoDst->numberOfFcsHits(); ++ i) {
    StPicoFcsHit* clu = mPicoDst->fcsHit(i);
    //cout << clu->fourMomentum()[2] << endl;//pz
   //cout << clu->energy() << endl;
    mhetest[0]->Fill(clu->energy());//clu->fourMomentum()[2]);
}
*/
//cout << "CLUSTERS!" << endl;
 //new

////////TClonesArray* clusters = (TClonesArray*) StPicoDst::picoArray (StPicoArrays::FcsCluster);

// std::sort(clusters->begin(), clusters->end(), [](StPicoFcsCluster* a, StPicoFcsCluster* b) {
//    return b->fourMomentum()[1] < a->fourMomentum()[1];
//  });    
//TMath::Sort(clusters->GetEntriesFast(), clusters->GetArray(), compareClustersByMomentum);
//clusters->Sort(compareClustersByMomentum);

if (mPicoDst->numberOfFcsClusters() < 2) {return kStOK;} // later will probably want to count what fraction of the time this happens

StPicoFcsCluster* candidates[2] = {nullptr, nullptr};// = {mPicoDst->fcsCluster(0),mPicoDst->fcsCluster(1)};
/*
for (int i = 0; i < mPicoDst->numberOfFcsClusters(); ++ i) {
    StPicoFcsCluster* clu = mPicoDst->fcsCluster(i);
    cout << StPicoAnalysisTest::Pt(clu->fourMomentum()[0],clu->fourMomentum()[1]) << endl;
    if (!candidates[0]) {
      candidates[0] = clu;
      continue;
    }
    if (StPicoAnalysisTest::Pt(clu->fourMomentum()[0],clu->fourMomentum()[1]) > StPicoAnalysisTest::Pt(candidates[0]->fourMomentum()[0],candidates[0]->fourMomentum()[1])) {
      candidates[1] = candidates[0];
      candidates[0] = clu;
    }
    else if (!candidates[1]) {
      candidates[1] = clu;
    }
    else if (StPicoAnalysisTest::Pt(clu->fourMomentum()[0],clu->fourMomentum()[1]) > StPicoAnalysisTest::Pt(candidates[1]->fourMomentum()[0],candidates[1]->fourMomentum()[1])) {
      candidates[1] = clu;
    }
    //mhetest[0]->Fill(clu->energy());//clu->fourMomentum()[2]);
}
*/
//cout << StPicoAnalysisTest::Pt(candidates[0]->fourMomentum()[0],candidates[0]->fourMomentum()[1]) << " " << StPicoAnalysisTest::Pt(candidates[1]->fourMomentum()[0],candidates[1]->fourMomentum()[1]) << endl;
/*
for (int i = 0; i < mPicoDst->numberOfFcsHits(); ++ i) {
  if (mPicoDst->fcsHit(i)->fourMomentum().Eta() != 0) {
    cout << "REAL NITTY GRITTY" << endl;
    cout << mPicoDst->fcsHit(i)->fourMomentum().Eta() << " " << mPicoDst->fcsHit(i)->detectorId() << endl;
    cout << mPicoDst->fcsHit(i)->fourMomentum()[0] << " " << mPicoDst->fcsHit(i)->fourMomentum()[1] << " " << mPicoDst->fcsHit(i)->fourMomentum()[2] << " " << mPicoDst->fcsHit(i)->fourMomentum()[3] << endl;
    cout << "~~~" << endl;
  }
}
*/
for (int i = 0; i < mPicoDst->numberOfFcsClusters(); ++ i) {
    StPicoFcsCluster* clu = mPicoDst->fcsCluster(i);
    Double_t Pt = StPicoAnalysisTest::Pt(clu->fourMomentum()[0],clu->fourMomentum()[1]);
    //cout << Pt << " " << clu->detectorId() << endl;
    if (clu->detectorId() == 0) {//north ECal
      if (!candidates[0]) {
        candidates[0] = clu;
        continue;
      }
      if (Pt > StPicoAnalysisTest::Pt(candidates[0]->fourMomentum()[0],candidates[0]->fourMomentum()[1])) {
        candidates[0] = clu;
      }
    }
    if (clu->detectorId() == 1) { //south ECal
      if (!candidates[1]) {
        candidates[1] = clu;
        continue;
      }
      if (Pt > StPicoAnalysisTest::Pt(candidates[1]->fourMomentum()[0],candidates[1]->fourMomentum()[1])) {
        candidates[1] = clu;
      }
    }
    //mhetest[0]->Fill(clu->energy());//clu->fourMomentum()[2]);
}

if (!candidates[0] || !candidates[1]) {return kStOK;} // later will probably want to count what fraction of the time this happens

//cout << StPicoAnalysisTest::Pt(candidates[0]->fourMomentum()[0],candidates[0]->fourMomentum()[1]) << " " << StPicoAnalysisTest::Pt(candidates[1]->fourMomentum()[0],candidates[1]->fourMomentum()[1]) << endl;


//cout << "TRACKS!" << endl;
/*
for (int i = 0; i < mPicoDst->numberOfFwdTracks(); ++ i) {
    StPicoFwdTrack* trk = mPicoDst->fwdTrack(i);
    //cout << trk->momentum()[1] << endl; //py
    //mtrpytest[0]->Fill(trk->momentum()[1]);
}
*/

  
  //  LOG_INFO << "candidates[0] = " << candidates[0] << " candidates[1] " << candidates[1] << endm;
  
  //No lepton pair candidates found above mETCut
  if(StPicoAnalysisTest::Pt(candidates[0]->fourMomentum()[0],candidates[0]->fourMomentum()[1]) < mETCut || StPicoAnalysisTest::Pt(candidates[1]->fourMomentum()[0],candidates[1]->fourMomentum()[1]) < mETCut) return kStOK; 
/*//don't want to use the database with the picos
  //3 vectors for lepton candidates
  StThreeVectorD VN = mFcsDb->getStarXYZ(0,candidates[0]->x(),candidates[0]->y());//,float FcsZ=-1.0, float zVertex=0.0)
  StThreeVectorD VS = mFcsDb->getStarXYZ(1,candidates[1]->x(),candidates[1]->y());
  */

  //Getting TOT & Cone for isolation cut
  float tot[2] = {0,0}; //eh
  float cone[2]= {0,0}; //ns
  double eta[2],phi[2];
  eta[0]=candidates[0]->fourMomentum().Eta(); //VN.pseudoRapidity();
  phi[0]=candidates[0]->fourMomentum().Phi(); //VN.phi();
  eta[1]=candidates[1]->fourMomentum().Eta(); //VS.pseudoRapidity();
  phi[1]=candidates[1]->fourMomentum().Phi(); //VS.phi();
  
  //for(int eh=0; eh<2; eh++){    
  //  for(int ns=0; ns<2; ns++){
      //int det=eh*2 + ns;
      //StSPtrVecFcsHit& hits = mFcsCollection->hits(det);    
      int n = mPicoDst->numberOfFcsHits(); //mFcsCollection->numberOfHits(det);
      for(int i=0; i<n; i++) {
        StPicoFcsHit* hit = mPicoDst->fcsHit(i);
        if (!hit) {continue;}
        int det = hit->detectorId();
        int ns = det % 2;
        int eh = det / 2;
        if (eh > 1) {continue;} //preshower north and south IDs are 4 and 5.
        //cout << "TEST: " << hit->fourMomentum().E() << " " << hit->fourMomentum().Eta() << " " << hit->fourMomentum().Phi() << " " << StPicoAnalysisTest::Pt(hit->fourMomentum()[0], hit->fourMomentum()[1]) << endl;
	      tot[eh] += hit->fourMomentum().E();
	      //StThreeVectorD v = mFcsDb->getStarXYZ(hits[i]);
	      double e= hit->fourMomentum().Eta(); //v.pseudoRapidity();
	      double p= hit->fourMomentum().Phi(); //v.phi();
	      double deta = e-eta[ns];
	      double dphi = p-phi[ns];
	      while(dphi> M_PI) {dphi -= 2*M_PI;}
	      while(dphi<-M_PI) {dphi += 2*M_PI;}
	      double dr = sqrt(deta*deta + dphi*dphi);
	      if(dr < mConeR) cone[ns] += hit->fourMomentum().E(); //hits[i]->energy();
      }
  //  }
  // }


  //2 body decay kinematics

  StLorentzVectorD ln(candidates[0]->fourMomentum()[0], candidates[0]->fourMomentum()[1], candidates[0]->fourMomentum()[2], candidates[0]->fourMomentum()[3]);
  StLorentzVectorD ls(candidates[1]->fourMomentum()[0], candidates[1]->fourMomentum()[1], candidates[1]->fourMomentum()[2], candidates[1]->fourMomentum()[3]);
  StLorentzVectorD di = ln + ls;
  StLorentzVectorD bln = ln.boost(-di);
  StLorentzVectorD bls = ls.boost(-di);
  double EN  = ln.e();
  double ES  = ls.e();
  double E   = di.e();
  double ETN = ln.perp();
  double ETS = ls.perp();
  double ET  = di.perp();
  double EZ  = di.pz();
  double M   = di.m();
  double Z   = abs(EN-ES)/(EN+ES);
  double CosTN = bln.cosTheta();
  double CosTS = bls.cosTheta();
  double CosT  = CosTN; //take north one for now... When we have tracking, take positive charged
  double Phi   = di.phi(); 
  LOG_DEBUG << Form("AAA CosTheta N=%7.4f S=%7.4f",CosTN,CosTS) << endm;

  //cout << "M: " << M << endl;
  //mmasstest[0]->Fill(M);


  //Ecal cluster SigmaMax
  double SigmaMaxN = candidates[0]->sigmaMax();
  double SigmaMaxS = candidates[1]->sigmaMax();
  
  //Ratio of DiLepton candidate to TOT  
  double ratioETOT = E/tot[0];
  double ratioHTOT = 9.99;
  if(tot[1]>0) ratioHTOT=E/tot[1];

  //Ratio of DiLepton candidates to cone
  double ratioConeN = EN/cone[0];
  double ratioConeS = ES/cone[1];
  
  //top pT associated track

  //TEMP! For now don't have a list of associated tracks. Will just have to track match by hand now:
  double max_delR = /*99999*/ 0.1; //0.1 is random guess right now
  int track_counter0 = 0, track_counter1 = 0;
  double best_pt0 = -9999, best_pt1 = -9999;
  for (int i = 0; i < mPicoDst->numberOfFwdTracks(); ++ i) {
    //if (i == 0) { cout << "have tracks!" << endl;}
    StPicoFwdTrack* trk = mPicoDst->fwdTrack(i);
    //cout << eta[0] << " " << StPicoAnalysisTest::Eta(trk) << " " << phi[0] << " " << StPicoAnalysisTest::Phi(trk) << endl;
    double delR0 = sqrt(pow(eta[0] - StPicoAnalysisTest::Eta(trk),2) + pow(phi[0] - StPicoAnalysisTest::Phi(trk),2)); //DeltaR to candidate 0
    double delR1 = sqrt(pow(eta[1] - StPicoAnalysisTest::Eta(trk),2) + pow(phi[1] - StPicoAnalysisTest::Phi(trk),2)); //DeltaR to candidate 1
    //cout << delR0 << " " << delR1 << endl;
    if (delR0 < max_delR) { //found "matched" track
      //cout << "have a track < 0.4!0" << endl;
      //max_delR = delR; //update minimum distance
      if (StPicoAnalysisTest::Pt(mPicoDst->fwdTrack(i)) > best_pt0) {
        //cout << "have a matched track!0" << endl;
        track_counter0 = i;//this is now the highest pT matched track
        best_pt0 = StPicoAnalysisTest::Pt(mPicoDst->fwdTrack(track_counter0));
      }
    }
    if (delR1 < max_delR) { //found "matched" track
      //cout << "have a track < 0.4!0" << endl;
      //max_delR = delR; //update minimum distance
      if (StPicoAnalysisTest::Pt(mPicoDst->fwdTrack(i)) > best_pt1) {
        //cout << "have a matched track!1" << endl;
        track_counter1 = i;//this is now the highest pT matched track
        best_pt1 = StPicoAnalysisTest::Pt(mPicoDst->fwdTrack(track_counter1));
      }
    }
  }

  float pt0=0,pt1=0,r0=0,r1=0;
  int cg0=0, cg1=0;
  if (track_counter0 > 0 || (mPicoDst->numberOfFwdTracks() > 0 && track_counter0 == 0 && sqrt(pow(eta[0] - StPicoAnalysisTest::Eta(mPicoDst->fwdTrack(0)),2) + pow(phi[0] - StPicoAnalysisTest::Phi(mPicoDst->fwdTrack(0)),2)) < max_delR)) { //either it's been advanced which means there was a match, or it hasn't because the match was on the first one.
    //cout << "have a final matched track!0" << endl;
    pt0 = StPicoAnalysisTest::Pt(mPicoDst->fwdTrack(track_counter0));
    r0 = ETN/(float) pt0;
    cg0 = mPicoDst->fwdTrack(track_counter0)->charge();
  }
  if (track_counter1 > 0 || (mPicoDst->numberOfFwdTracks() > 0 && track_counter1 == 0 && sqrt(pow(eta[1] - StPicoAnalysisTest::Eta(mPicoDst->fwdTrack(0)),2) + pow(phi[1] - StPicoAnalysisTest::Phi(mPicoDst->fwdTrack(0)),2)) < max_delR)) { //either it's been advanced which means there was a match, or it hasn't because the match was on the first one.
    //cout << "have a final matched track!1" << endl; 
    pt1 = StPicoAnalysisTest::Pt(mPicoDst->fwdTrack(track_counter1));
    r1 = ETS/(float) pt1;
    cg1 = mPicoDst->fwdTrack(track_counter1)->charge();
  }

  /*
  StPicoFwdTrack *trk1=0, *trk2=0;
  float pt1=0,pt2=0,r1=0,r2=0;
  int cg1=0, cg2=0;
  if(candidates[0]->tracks().size()>0) {
    trk1=candidates[0]->tracks()[0]; 
    pt1=trk1->momentum().perp(); 
    r1=ETN/pt1;
    cg1=trk1->charge();
  }
  if(candidates[1]->tracks().size()>0) {
    trk2=candidates[1]->tracks()[0]; 
    pt2=trk2->momentum().perp(); 
    r2=ETS/pt2;
    cg2=trk2->charge();
  } 
  LOG_INFO << Form("trk1 pt=%6.2f et=%6.2f R=%6.4f cg=%2d",pt1,ETN,r1,cg1)<<endm;
  LOG_INFO << Form("trk2 pt=%6.2f et=%6.2f R=%6.4f cg=%2d",pt2,ETS,r2,cg2)<<endm;
*/

  for(int cut=0; cut<mNCut; cut++){
    if(cut==1 && ratioETOT<mETotCut) break;
    if(cut==2 && ratioHTOT<mHTotCut) break;
    if(cut==3 && (ratioConeN<mConeCut || ratioConeS<mConeCut)) break;    
    if(cut==4 && (SigmaMaxN > mSigmaMaxCut || SigmaMaxS > mSigmaMaxCut) ) break;
    if(cut==5 && (r0 < mETPTCutLow || r1 < mETPTCutLow || r0 > mETPTCutHigh || r1 > mETPTCutHigh )) break;    
    if(cut==6 && cg0 + cg1 != 0) break;

    //ecal[0]->fourMomentum().perp());//ratioETOT);

    mETot[cut]->Fill(ratioETOT);
    mHTot[cut]->Fill(ratioHTOT);
    mCone[cut]->Fill(ratioConeN);
    mCone[cut]->Fill(ratioConeS);
    mSigmax[cut]->Fill(SigmaMaxN); // 2nd moment
    mSigmax[cut]->Fill(SigmaMaxS); 
    /*if(trk1)*/ mPToverET[cut]->Fill(r0);//think later about how to apply a similar selection on existence of trk1, trk2 as before.
    /*if(trk2)*/ mPToverET[cut]->Fill(r1);
    /*if(trk1 && trk2)*/ mChargeSum[cut]->Fill(cg0 + cg1);

    mET  [cut]->Fill(ET);
    mEZ  [cut]->Fill(EZ);
    mM   [cut]->Fill(M);
    mZ   [cut]->Fill(Z);
    mCosT[cut]->Fill(CosT);
    mPhi [cut]->Fill(Phi);

    mET12[cut]->Fill(ETN,ETS);
    mXFPT[cut]->Fill(EN/255.0,ETN);
    mXFPT[cut]->Fill(ES/255.0,ETS);
    mXY[cut]->Fill(candidates[0]->x(), candidates[0]->y());// Mean x ("center of gravity") in local grid coordinate (1st moment).
    mXY[cut]->Fill(candidates[1]->x(),candidates[1]->y());
    /*if(trk1)*/ mPTET[cut]->Fill(ETN,pt0);
    /*if(trk2)*/ mPTET[cut]->Fill(ETS,pt1);

  }

  return kStOK; 
}
  
Int_t StPicoAnalysisTest::Finish(){
  mFile->Write();
  mFile->Close();
  printf("StPicoAnalysisTest::Finish - Closing %s\n",mFilename);
  return kStOK;
};

ClassImp(StPicoAnalysisTest);
