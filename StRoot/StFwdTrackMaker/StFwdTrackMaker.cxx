#include "StFwdTrackMaker/StFwdTrackMaker.h"
#include "StFwdTrackMaker/include/Tracker/FwdHit.h"
#include "StFwdTrackMaker/include/Tracker/FwdTracker.h"
#include "StFwdTrackMaker/include/Tracker/TrackFitter.h"
#include "StFwdTrackMaker/include/Tracker/FwdGeomUtils.h"
#include "StFwdTrackMaker/include/Tracker/ObjExporter.h"

#include "KiTrack/IHit.h"
#include "GenFit/Track.h"
#include "GenFit/GFRaveVertexFactory.h"

#include "TMath.h"

#include <limits>
#include <map>
#include <string>
#include <string>
#include <vector>

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

#include "StEventUtilities/StEventHelper.h"

#include "tables/St_g2t_fts_hit_Table.h"
#include "tables/St_g2t_track_Table.h"
#include "tables/St_g2t_vertex_Table.h"
#include "tables/St_g2t_event_Table.h"

#include "StarMagField/StarMagField.h"

#include "St_base/StMessMgr.h"
#include "StarClassLibrary/StPhysicalHelix.hh"
#include "StarClassLibrary/SystemOfUnits.h"

#include <SystemOfUnits.h>
#include <exception>

#include "TROOT.h"
#include "TLorentzVector.h"

#include "StRoot/StEpdUtil/StEpdGeom.h"
#include "StFcsDbMaker/StFcsDb.h"
#include "StFstUtil/StFstCollection.h"

#include "StEvent/StFwdTrack.h"

#define DLOG(...) printf(__VA_ARGS__);

FwdSystem* FwdSystem::sInstance = nullptr;
TMVA::Reader * BDTCrit2::reader = nullptr;
float BDTCrit2::Crit2_RZRatio = -999;
float BDTCrit2::Crit2_DeltaRho = -999;
float BDTCrit2::Crit2_DeltaPhi = -999;
float BDTCrit2::Crit2_StraightTrackRatio = -999;

size_t StarFieldAdaptor::sNCalls;

//_______________________________________________________________________________________
class GenfitUtils{
    public:

    // For now, accept anything we are passed, no matter what it is or how bad it is
    template<typename T> static bool accept( T ) { return true; }

    
}; // GenfitUtils


// 


// Basic sanity cuts on genfit tracks
template<> bool GenfitUtils::accept( genfit::Track *track )
{
    // This also gets rid of failed fits (but may need to explicitly
    // for fit failure...)
    if (track->getNumPoints() <= 0 ) return false; // fit may have failed

    auto cardinal = track->getCardinalRep();

    // Check that the track fit converged
    auto status = track->getFitStatus( cardinal );
    if ( !status->isFitConverged() ) {
    return false;
    }


    // Next, check that all points on the track have fitter info
    // (may be another indication of a failed fit?)
    for ( auto point : track->getPoints() ) {
    if ( !point->hasFitterInfo(cardinal) ) {
    return false;
    }
    }

    // Following line fails with an exception, because some tracks lack 
    //   forward update, or prediction in fitter info at the first point
    //
    // genfit::KalmanFitterInfo::getFittedState(bool) const of 
    //                         GenFit/fitters/src/KalmanFitterInfo.cc:250

    // Fitted state at the first point
    // const auto &atFirstPoint = track->getFittedState();

    // Getting the fitted state from a track occasionally fails, because
    // the first point on the fit doesn't have forward/backward fit
    // information.  So we want the first point with fit info...

    genfit::TrackPoint* first = nullptr;
    unsigned int ipoint = 0;
    for ( ipoint = 0; ipoint < track->getNumPoints(); ipoint++ ) {
    first = track->getPointWithFitterInfo( ipoint );
    if ( first ) break;
    } 

    // No points on the track have fit information
    if ( !first ) {
    LOG_INFO << "No fit information on track" << endm;
    return false;
    }

    auto& fittedState= track->getFittedState(ipoint);

    TVector3 momentum = fittedState.getMom();
    double   pt       = momentum.Perp();

    if (pt < 0.10 ) return false; // below this

    return true;

};


//______________________________________________________________________________________

class SiRasterizer {
  public:
    SiRasterizer() {}
    SiRasterizer(FwdTrackerConfig &_cfg) { setup(_cfg); }
    ~SiRasterizer() {}
    void setup(FwdTrackerConfig &_cfg) {
        cfg = _cfg;
        mRasterR = cfg.get<double>("SiRasterizer:r", 3.0);
        mRasterPhi = cfg.get<double>("SiRasterizer:phi", 0.1);
    }

    bool active() {
        return cfg.get<bool>("SiRasterizer:active", false);
    }

    TVector3 raster(TVector3 p0) {
        TVector3 p = p0;
        double r = p.Perp();
        double phi = p.Phi();
        const double minR = 5.0;
        // 5.0 is the r minimum of the Si
        p.SetPerp(minR + (std::floor((r - minR) / mRasterR) * mRasterR + mRasterR / 2.0));
        p.SetPhi(-TMath::Pi() + (std::floor((phi + TMath::Pi()) / mRasterPhi) * mRasterPhi + mRasterPhi / 2.0));
        return p;
    }

    FwdTrackerConfig cfg;
    double mRasterR, mRasterPhi;
};

//  Wrapper class around the forward tracker
class ForwardTracker : public ForwardTrackMaker {
  public:
    // Replaces original initialization.  Config file and hitloader
    // will be provided by the maker.
    void initialize( bool genHistograms ) {
        LOG_INFO << "ForwardTracker::initialize()" << endm;
        nEvents = 1; // only process single event

        // Create the forward system...
        FwdSystem::sInstance = new FwdSystem();

        // make our quality plotter
        mQualityPlotter = new QualityPlotter(mConfig);
        mQualityPlotter->makeHistograms(mConfig.get<size_t>("TrackFinder:nIterations", 1));

        // initialize the track fitter
        mTrackFitter = new TrackFitter(mConfig);
        mTrackFitter->setGenerateHistograms(genHistograms);
        mTrackFitter->setup();

        ForwardTrackMaker::initialize( genHistograms );
    }

    void finish() {

        if ( mGenHistograms ){
            mQualityPlotter->finish();
            writeEventHistograms();
        }

        if (FwdSystem::sInstance){
            delete FwdSystem::sInstance;
            FwdSystem::sInstance = 0;
        }
        if (mQualityPlotter){
            delete mQualityPlotter;
            mQualityPlotter = 0;
        }
        if (mTrackFitter){
            delete mTrackFitter;
            mTrackFitter= 0;
        }
    }
};



//________________________________________________________________________
StFwdTrackMaker::StFwdTrackMaker() : StMaker("fwdTrack"), mGenHistograms(false), mGenTree(false), mForwardTracker(nullptr), mForwardData(nullptr){
    SetAttr("useFtt",1);                 // Default Ftt on 
    SetAttr("useFst",1);                 // Default Fst on
    SetAttr("useFcs",1);                 // Default Fcs on
    SetAttr("config", "config.xml");     // Default configuration file (user may override before Init())
    SetAttr("fillEvent",1); // fill StEvent
};

int StFwdTrackMaker::Finish() {
    
    auto prevDir = gDirectory;
    if ( mGenHistograms ) {
        
        // output file name
        string name = mFwdConfig.get<string>("Output:url", "fwdTrackerOutput.root");
        TFile *fOutput = new TFile(name.c_str(), "RECREATE");
        fOutput->cd();

        fOutput->mkdir("StFwdTrackMaker");
        fOutput->cd("StFwdTrackMaker");
        for (auto nh : mHistograms) {
            nh.second->SetDirectory(gDirectory);
            nh.second->Write();
        }
        fOutput->cd("");
    }

    mForwardTracker->finish();

    gDirectory = prevDir;

    if (mGenTree) {
        mTreeFile->cd();
        mTree->Write();
        mTreeFile->Write();
    }
    return kStOk;
}

//________________________________________________________________________
int StFwdTrackMaker::Init() {

    // Initialize configuration file
    std::string configFile = SAttr("config");
    if (mConfigFile.length() > 4) {
        configFile = mConfigFile;
        LOG_INFO << "Forward Tracker is using config file : " <<  mConfigFile << endm;
    }

    mFwdConfig.load( configFile );

    if (mGenTree) {
        mTreeFile = new TFile("mltree.root", "RECREATE");
        mTree = new TTree("Stg", "stg hits");
        mTree->Branch("hN",         &mTreeData. hN, "hN/I");
        mTree->Branch("hX",         &mTreeData. hX, "hX/F");
        mTree->Branch("hY",         &mTreeData. hY, "hY/F");
        mTree->Branch("hZ",         &mTreeData. hZ, "hZ/F");
        
        mTree->Branch("hTrackId",   &mTreeData. hTrackId, "hTrackId/I");
        mTree->Branch("hVolumeId",  &mTreeData. hVolumeId, "hVolumeId/I");
        mTree->Branch("hPt",        &mTreeData. hPt, "hPt/F");
        mTree->Branch("hVertexId",  &mTreeData. hVertexId, "hVertexId/I");

        // mc tracks
        mTree->Branch("mcN",        &mTreeData. mcN, "mcN/I");
        mTree->Branch("mcPt",       &mTreeData. mcPt, "mcPt/F");
        mTree->Branch("mcEta",      &mTreeData. mcEta, "mcEta/F");
        mTree->Branch("mcPhi",      &mTreeData. mcPhi, "mcPhi/F");
        mTree->Branch("mcCharge",   &mTreeData. mcCharge, "mcCharge/I");
        mTree->Branch("mcVertexId", &mTreeData. mcVertexId, "mcVertexId/I");

        // mcverts
        mTree->Branch("vmcN",       &mTreeData. vmcN, "vmcN/I");
        mTree->Branch("vmcX",       &mTreeData. vmcX, "vmcX/F");
        mTree->Branch("vmcY",       &mTreeData. vmcY, "vmcY/F");
        mTree->Branch("vmcZ",       &mTreeData. vmcZ, "vmcZ/F");

        // rcverts
        mTree->Branch("vrcN",       &mTreeData. vrcN, "vrcN/I");
        mTree->Branch("vrcX",       &mTreeData. vrcX, "vrcX/F");
        mTree->Branch("vrcY",       &mTreeData. vrcY, "vrcY/F");
        mTree->Branch("vrcZ",       &mTreeData. vrcZ, "vrcZ/F");

        // rc tracks
        mTree->Branch("rcN",        &mTreeData. rcN, "rcN/I");
        mTree->Branch("rcPt",       &mTreeData. rcPt, "rcPt/F");
        mTree->Branch("rcEta",      &mTreeData. rcEta, "rcEta/F");
        mTree->Branch("rcPhi",      &mTreeData. rcPhi, "rcPhi/F");
        mTree->Branch("rcCharge",   &mTreeData. rcCharge, "rcCharge/I");
        mTree->Branch("rcTrackId",  &mTreeData. rcTrackId, "rcTrackId/I");
        mTree->Branch("rcNumFst",   &mTreeData. rcNumFst, "rcNumFst/I");
        mTree->Branch("rcQuality",  &mTreeData. rcQuality, "rcQuality/F");


    //     std::string path = "TrackFinder.Iteration[0].SegmentBuilder";
    //     std::vector<string> paths = mFwdConfig.childrenOf(path);

    //     for (string p : paths) {
    //         string name = mFwdConfig.get<string>(p + ":name", "");
    //         mTreeCrits[name]; // create the entry
    //         mTree->Branch(name.c_str(), &mTreeCrits[name]);
    //         mTree->Branch((name + "_trackIds").c_str(), &mTreeCritTrackIds[name]);

    //         if ( name == "Crit2_RZRatio" ){
    //             string n = name + "_x1";
    //             mTreeCrits[(n)]; mTree->Branch(n.c_str(), &mTreeCrits[n]);

    //             n = name + "_y1";
    //             mTreeCrits[(n)]; mTree->Branch(n.c_str(), &mTreeCrits[n]);

    //             n = name + "_z1";
    //             mTreeCrits[(n)]; mTree->Branch(n.c_str(), &mTreeCrits[n]);

    //             n = name + "_x2";
    //             mTreeCrits[(n)]; mTree->Branch(n.c_str(), &mTreeCrits[n]);

    //             n = name + "_y2";
    //             mTreeCrits[(n)]; mTree->Branch(n.c_str(), &mTreeCrits[n]);

    //             n = name + "_z2";
    //             mTreeCrits[(n)]; mTree->Branch(n.c_str(), &mTreeCrits[n]);

    //             n = name + "_h1";
    //             mTreeCritTrackIds[(n)]; mTree->Branch(n.c_str(), &mTreeCritTrackIds[n]);
    //             n = name + "_h2";
    //             mTreeCritTrackIds[(n)]; mTree->Branch(n.c_str(), &mTreeCritTrackIds[n]);
    //             n = name + "_h3";
    //             mTreeCritTrackIds[(n)]; mTree->Branch(n.c_str(), &mTreeCritTrackIds[n]);
    //         }

    //         if ( name == "Crit2_BDT" ){
    //             string n = name + "_DeltaPhi";
    //             mTreeCrits[(n)]; mTree->Branch(n.c_str(), &mTreeCrits[n]);
    //             n = name + "_DeltaRho";
    //             mTreeCrits[(n)]; mTree->Branch(n.c_str(), &mTreeCrits[n]);
    //             n = name + "_RZRatio";
    //             mTreeCrits[(n)]; mTree->Branch(n.c_str(), &mTreeCrits[n]);
    //             n = name + "_StraightTrackRatio";
    //             mTreeCrits[(n)]; mTree->Branch(n.c_str(), &mTreeCrits[n]);
    //         }
    //     }

    //     // Three hit criteria
    //     path = "TrackFinder.Iteration[0].ThreeHitSegments";
    //     paths = mFwdConfig.childrenOf(path);

    //     for (string p : paths) {
    //         string name = mFwdConfig.get<string>(p + ":name", "");
    //         mTreeCrits[name]; // create the entry
    //         mTree->Branch(name.c_str(), &mTreeCrits[name]);
    //         mTree->Branch((name + "_trackIds").c_str(), &mTreeCritTrackIds[name]);
    //     }

    //     mTree->SetAutoFlush(0);
    } // gen tree


    // create an SiRasterizer in case we need it 
    mSiRasterizer = std::shared_ptr<SiRasterizer>( new SiRasterizer(mFwdConfig));
    mForwardTracker = std::shared_ptr<ForwardTracker>(new ForwardTracker());
    mForwardTracker->setConfig(mFwdConfig);

    // only save criteria values if we are generating a tree.
    mForwardTracker->setSaveCriteriaValues(mGenTree);

    mForwardData = std::shared_ptr<FwdDataSource>(new FwdDataSource());
    mForwardTracker->setData(mForwardData);
    mForwardTracker->initialize( mGenHistograms );

    if ( mGenHistograms ){
        mHistograms["fwdVertexZ"] = new TH1D("fwdVertexZ", "FWD Vertex (RAVE);z", 1000, -50, 50);
        mHistograms["fwdVertexXY"] = new TH2D("fwdVertexXY", "FWD Vertex (RAVE);x;y", 100, -1, 1, 100, -1, 1);
        mHistograms["fwdVertexDeltaZ"] = new TH2D("fwdVertexDeltaZ", "FWD Vertex - MC Vertex;#Delta z", 100, -1, 1, 100, -1, 1);
        
        mHistograms["McEventEta"] = new TH1D("McEventEta", ";MC Track Eta", 1000, -5, 5);
        mHistograms["McEventPt"] = new TH1D("McEventPt", ";MC Track Pt (GeV/c)", 1000, 0, 10);
        mHistograms["McEventPhi"] = new TH1D("McEventPhi", ";MC Track Phi", 1000, 0, 6.2831852);

        // these are tracks within 2.5 < eta < 4.0
        mHistograms["McEventFwdEta"] = new TH1D("McEventFwdEta", ";MC Track Eta", 1000, -5, 5);
        mHistograms["McEventFwdPt"] = new TH1D("McEventFwdPt", ";MC Track Pt (GeV/c)", 1000, 0, 10);
        mHistograms["McEventFwdPhi"] = new TH1D("McEventFwdPhi", ";MC Track Phi", 1000, 0, 6.2831852);

        // create mHistograms
        mHistograms["nMcTracks"] = new TH1I("nMcTracks", ";# MC Tracks/Event", 1000, 0, 1000);
        mHistograms["nMcTracksFwd"] = new TH1I("nMcTracksFwd", ";# MC Tracks/Event", 1000, 0, 1000);
        mHistograms["nMcTracksFwdNoThreshold"] = new TH1I("nMcTracksFwdNoThreshold", ";# MC Tracks/Event", 1000, 0, 1000);

        mHistograms["nHitsSTGC"] = new TH1I("nHitsSTGC", ";# STGC Hits/Event", 1000, 0, 1000);
        mHistograms["nHitsFSI"] = new TH1I("nHitsFSI", ";# FSIT Hits/Event", 1000, 0, 1000);

        mHistograms["stgc_volume_id"] = new TH1I("stgc_volume_id", ";stgc_volume_id", 50, 0, 50);
        mHistograms["fsi_volume_id"] = new TH1I("fsi_volume_id", ";fsi_volume_id", 50, 0, 50);

        mHistograms["fsiHitDeltaR"] = new TH1F("fsiHitDeltaR", "FSI; delta r (cm); ", 500, -5, 5);
        mHistograms["fsiHitDeltaPhi"] = new TH1F("fsiHitDeltaPhi", "FSI; delta phi; ", 500, -5, 5);

        // there are 4 stgc stations
        for (int i = 0; i < 4; i++) {
            mHistograms[TString::Format("stgc%dHitMap", i).Data()] = new TH2F(TString::Format("stgc%dHitMap", i), TString::Format("STGC Layer %d; x (cm); y(cm)", i), 200, -100, 100, 200, -100, 100);

            mHistograms[TString::Format("stgc%dHitMapPrim", i).Data()] = new TH2F(TString::Format("stgc%dHitMapPrim", i), TString::Format("STGC Layer %d; x (cm); y(cm)", i), 200, -100, 100, 200, -100, 100);
            mHistograms[TString::Format("stgc%dHitMapSec", i).Data()] = new TH2F(TString::Format("stgc%dHitMapSec", i), TString::Format("STGC Layer %d; x (cm); y(cm)", i), 200, -100, 100, 200, -100, 100);
        }

        // There are 3 silicon stations
        for (int i = 0; i < 3; i++) {
            mHistograms[TString::Format("fsi%dHitMap", i).Data()] = new TH2F(TString::Format("fsi%dHitMap", i), TString::Format("FSI Layer %d; x (cm); y(cm)", i), 200, -100, 100, 200, -100, 100);
            mHistograms[TString::Format("fsi%dHitMapZ", i).Data()] = new TH2F(TString::Format("fsi%dHitMapZ", i), TString::Format("FSI Layer %d; x (cm); y(cm)", i), 200, -100, 100, 200, -100, 100);

            mHistograms[TString::Format("fsi%dHitMapR", i).Data()] = new TH1F(TString::Format("fsi%dHitMapR", i), TString::Format("FSI Layer %d; r (cm); ", i), 500, 0, 50);
            mHistograms[TString::Format("fsi%dHitMapPhi", i).Data()] = new TH1F(TString::Format("fsi%dHitMapPhi", i), TString::Format("FSI Layer %d; phi; ", i), 320, 0, TMath::Pi() * 2 + 0.1);
        }

    } // mGenHistograms
    LOG_INFO << "StFwdTrackMaker::Init" << endm;
    return kStOK;
};

TMatrixDSym makeSiCovMat(TVector3 hit, FwdTrackerConfig &xfg) {
    // we can calculate the CovMat since we know the det info, but in future we should probably keep this info in the hit itself

    float rSize = xfg.get<float>("SiRasterizer:r", 3.0);
    float phiSize = xfg.get<float>("SiRasterizer:phi", 0.004);

    // measurements on a plane only need 2x2
    // for Si geom we need to convert from cylindrical to cartesian coords
    TMatrixDSym cm(2);
    TMatrixD T(2, 2);
    TMatrixD J(2, 2);
    const float x = hit.X();
    const float y = hit.Y();
    const float R = sqrt(x * x + y * y);
    const float cosphi = x / R;
    const float sinphi = y / R;
    const float sqrt12 = sqrt(12.);

    const float dr = rSize / sqrt12;
    const float dphi = (phiSize) / sqrt12;

    // Setup the Transposed and normal Jacobian transform matrix;
    // note, the si fast sim did this wrong
    // row col
    T(0, 0) = cosphi;
    T(0, 1) = -R * sinphi;
    T(1, 0) = sinphi;
    T(1, 1) = R * cosphi;

    J(0, 0) = cosphi;
    J(0, 1) = sinphi;
    J(1, 0) = -R * sinphi;
    J(1, 1) = R * cosphi;

    TMatrixD cmcyl(2, 2);
    cmcyl(0, 0) = dr * dr;
    cmcyl(1, 1) = dphi * dphi;

    TMatrixD r = T * cmcyl * J;

    // note: float sigmaX = sqrt(r(0, 0));
    // note: float sigmaY = sqrt(r(1, 1));

    cm(0, 0) = r(0, 0);
    cm(1, 1) = r(1, 1);
    cm(0, 1) = r(0, 1);
    cm(1, 0) = r(1, 0);

    TMatrixDSym tamvoc(3);
    tamvoc( 0, 0 ) = cm(0, 0); tamvoc( 0, 1 ) = cm(0, 1); tamvoc( 0, 2 ) = 0.0;
    tamvoc( 1, 0 ) = cm(1, 0); tamvoc( 1, 1 ) = cm(1, 1); tamvoc( 1, 2 ) = 0.0;
    tamvoc( 2, 0 ) = 0.0;      tamvoc( 2, 1 ) = 0.0; tamvoc( 2, 2 )      = 0.01*0.01;


    return tamvoc;
}

void StFwdTrackMaker::loadStgcHits( FwdDataSource::McTrackMap_t &mcTrackMap, FwdDataSource::HitMap_t &hitMap, int count ){

    // Get the StEvent handle to see if the rndCollection is available
    StEvent *event = (StEvent *)this->GetDataSet("StEvent");

    StFttCollection *col = event->fttCollection();
    string fttFromSource = mFwdConfig.get<string>( "Source:ftt", "" );
    
    mTreeData.hN = 0;

    if ( col && col->numberOfPoints() > 0 || "DATA" == fttFromSource ){
        TMatrixDSym hitCov3(3);
        const double sigXY = 10;
        hitCov3(0, 0) = sigXY * sigXY;
        hitCov3(1, 1) = sigXY * sigXY;
        hitCov3(2, 2) = 20; // unused since they are loaded as points on plane
        for ( auto point : col->points() ){
            // LOG_INFO << "Adding FTT Point" << endm;
            int plane_id = 0;
            FwdHit *hit = new FwdHit(count++, point->xyz().x()/10.0, point->xyz().y()/10.0, point->xyz().z(), -point->plane(), 0, hitCov3, nullptr);
            mFttHits.push_back( TVector3( point->xyz().x()/10.0, point->xyz().y()/10.0, point->xyz().z() )  );
            if ( mGenHistograms ) {
                this->mHistograms[TString::Format("stgc%dHitMapSec", point->plane()).Data()]->Fill(point->xyz().x()/10.0, point->xyz().y()/10.0);
            }
            // Add the hit to the hit map
            hitMap[hit->getSector()].push_back(hit);

            // LOG_INFO << TString::Format( "(%f, %f, %f)", point->xyz().x(), point->xyz().y(), point->xyz().z() );

            if (mGenTree && mTreeData.hN < MAX_TREE_ELEMENTS) {
                
                mTreeData.hX.push_back( point->xyz().x()/10.0 );
                mTreeData.hY.push_back( point->xyz().y()/10.0 );
                mTreeData.hZ.push_back( point->xyz().z() );
                mTreeData.hTrackId.push_back( 0 );
                mTreeData.hVolumeId.push_back( point->plane() );
                mTreeData.hPt.push_back( 0 );
                mTreeData.hVertexId.push_back( 0 );
                mTreeData.hN++;
            }
        }

        return;
    }


    StRnDHitCollection *rndCollection = nullptr;
    if ( event) {
        rndCollection = event->rndHitCollection();
    }

    

    if ( !rndCollection || "GEANT" == fttFromSource ){
        LOG_DEBUG << "Loading sTGC hits directly from GEANT hits" << endm;
        loadStgcHitsFromGEANT( mcTrackMap, hitMap, count );
    } else {
        LOG_DEBUG << "loading sTGC from StEvent" << endm;
        loadStgcHitsFromStEvent( mcTrackMap, hitMap, count );
    }
} // loadStgcHits

void StFwdTrackMaker::loadStgcHitsFromGEANT( FwdDataSource::McTrackMap_t &mcTrackMap, FwdDataSource::HitMap_t &hitMap, int count ){
    /************************************************************/
    // STGC Hits
    St_g2t_fts_hit *g2t_stg_hits = (St_g2t_fts_hit *)GetDataSet("geant/g2t_stg_hit");


    this->mTreeData.hN = 0;
    if (!g2t_stg_hits){
        LOG_WARN << "geant/g2t_stg_hit is empty" << endm; 
        return;
    }

    // make the Covariance Matrix once and then reuse
    TMatrixDSym hitCov3(3);
    const double sigXY = 0.01;
    hitCov3(0, 0) = sigXY * sigXY;
    hitCov3(1, 1) = sigXY * sigXY;
    hitCov3(2, 2) = 0.0; // unused since they are loaded as points on plane

    int nstg = g2t_stg_hits->GetNRows();

    LOG_DEBUG << "This event has " << nstg << " stg hits in geant/g2t_stg_hit " << endm;
    if ( mGenHistograms ) {
        this->mHistograms["nHitsSTGC"]->Fill(nstg);
    }
    

    bool filterGEANT = mFwdConfig.get<bool>( "Source:fttFilter", false );
    for (int i = 0; i < nstg; i++) {

        g2t_fts_hit_st *git = (g2t_fts_hit_st *)g2t_stg_hits->At(i);
        if (0 == git)
            continue; // geant hit
        int track_id = git->track_p;
        int volume_id = git->volume_id;
        int plane_id = (volume_id - 1) / 4;           // from 1 - 16. four chambers per station
        float x = git->x[0] + gRandom->Gaus(0, sigXY); // 100 micron blur according to approx sTGC reso
        float y = git->x[1] + gRandom->Gaus(0, sigXY); // 100 micron blur according to approx sTGC reso
        float z = git->x[2];


        
        if (mGenTree && mTreeData.hN < MAX_TREE_ELEMENTS) {
            mTreeData.hX.push_back( x );
            mTreeData.hY.push_back( y );
            mTreeData.hZ.push_back( z );
            mTreeData.hTrackId.push_back( track_id );
            mTreeData.hVolumeId.push_back( plane_id );
            mTreeData.hPt.push_back( mcTrackMap[track_id]->mPt );
            mTreeData.hVertexId.push_back( mcTrackMap[track_id]->mStartVertex );
            mTreeData.hN++;
        } else if ( mGenTree ){
            LOG_WARN << "Truncating hits in TTree output" << endm;
        }

        if ( mGenHistograms ){
            this->mHistograms["stgc_volume_id"]->Fill(volume_id);
        }

        if (plane_id < 4 && plane_id >= 0) {
            if ( mGenHistograms ){
                this->mHistograms[TString::Format("stgc%dHitMap", plane_id).Data()]->Fill(x, y);
            }
        } else {
            continue;
        }

        // this rejects GEANT hits with eta -999 - do we understand this effect?
        if ( filterGEANT ) {
            if ( mcTrackMap[track_id] && fabs(mcTrackMap[track_id]->mEta) > 5.0 ){
                
                if ( mGenHistograms ) 
                    this->mHistograms[TString::Format("stgc%dHitMapSec", plane_id).Data()]->Fill(x, y);
                continue;
            } else if ( mcTrackMap[track_id] && fabs(mcTrackMap[track_id]->mEta) < 5.0 ){
                if ( mGenHistograms ) this->mHistograms[TString::Format("stgc%dHitMapPrim", plane_id).Data()]->Fill(x, y);
            }
        }

        FwdHit *hit = new FwdHit(count++, x, y, z, -plane_id, track_id, hitCov3, mcTrackMap[track_id]);

        // Add the hit to the hit map
        hitMap[hit->getSector()].push_back(hit);
        mFttHits.push_back( TVector3( x, y, z )  );

        // Add hit pointer to the track
        if (mcTrackMap[track_id])
            mcTrackMap[track_id]->addHit(hit);
    } // loop on hits

    if (mGenTree){
        LOG_INFO << "Saving " << mTreeData.hN << " hits in Tree" << endm;
    }
} // loadStgcHits

void StFwdTrackMaker::loadStgcHitsFromStEvent( FwdDataSource::McTrackMap_t &mcTrackMap, FwdDataSource::HitMap_t &hitMap, int count ){

    // Get the StEvent handle
    StEvent *event = (StEvent *)this->GetDataSet("StEvent");
    if (!event)
        return;

    StRnDHitCollection *rndCollection = event->rndHitCollection();
    if (!rndCollection) {
        LOG_INFO << "No StRnDHitCollection found" << endm;
        return;
    }

    const StSPtrVecRnDHit &hits = rndCollection->hits();

    // we will reuse this to hold the cov mat
    TMatrixDSym hitCov3(3);
    
    for (unsigned int stgchit_index = 0; stgchit_index < hits.size(); stgchit_index++) {
        StRnDHit *hit = hits[stgchit_index];

        if ( hit->layer() <= 6 ){
            // skip FST hits here
            continue;
        }

        int layer = hit->layer() - 9;

        const StThreeVectorF pos = hit->position();

        StMatrixF covmat = hit->covariantMatrix();

        // copy covariance matrix from StMatrixF
        std::copy(&covmat(0,0), &covmat(0,0) + 9, hitCov3.GetMatrixArray());

        shared_ptr<McTrack> mct = nullptr;
        if ( hit->idTruth() > 0 && mcTrackMap.count( hit->idTruth() ) ){
            mct = mcTrackMap[hit->idTruth()];
        }
        FwdHit *fhit = new FwdHit(count++, hit->position().x(), hit->position().y(), hit->position().z(), -layer, hit->idTruth(), hitCov3, mct);

        // Add the hit to the hit map
        hitMap[fhit->getSector()].push_back(fhit);

        // Add hit pointer to the track
        if ( hit->idTruth() > 0 && mcTrackMap[hit->idTruth()])
            mcTrackMap[hit->idTruth()]->addHit(fhit);
    }
} //loadStgcHitsFromStEvent

void StFwdTrackMaker::loadFstHits( FwdDataSource::McTrackMap_t &mcTrackMap, FwdDataSource::HitMap_t &hitMap, int count ){

    // Get the StEvent handle to see if the rndCollection is available
    StEvent *event = (StEvent *)this->GetDataSet("StEvent");


    LOG_INFO << "LOAD FST HITS" << endm;
    // StFstHitCollection *fstHitCollection = event->fstHitCollection();


    // TObjectSet *fstDataSet = (TObjectSet *) GetDataSet("fstRawHitAndCluster");

    // if (!fstDataSet) {
    //     LOG_WARN << "Make() - fstRawHitAndCluster dataset not found. No FST hits will be available for tracking" << endm;
    //     return;
    // }

    // StFstCollection *fstCollectionPtr = (StFstCollection *) fstDataSet->GetObject();

    // if ( !fstCollectionPtr ) {
    //     LOG_WARN << "Make() - StFstCollection not found. No FST hits will be available for tracking" << endm;
    //     return;
    // }

    // for (unsigned char wedgeIdx = 0; wedgeIdx < kFstNumWedges; wedgeIdx++) {
    //     StFstClusterCollection *clusterCollectionPtr = fstCollectionPtr->getClusterCollection(wedgeIdx );


    //     if ( clusterCollectionPtr ) {
    //         unsigned int numClusters = clusterCollectionPtr->getNumClusters();
    //         LOG_INFO << "Make() - Number of clusters found in wedge " << (int)(wedgeIdx + 1) << ": " << numClusters << endm;

    //         unsigned short idTruth = 0;
    //         unsigned char  nRawHits = -1, nRawHitsR = -1, nRawHitsPhi = -1;
    //         unsigned char  disk = -1, wedge = -1, sensor = -1, apv = -1;
    //         int  meanRStrip = -1, meanPhiStrip = -1;
    //                     float  charge = 0., chargeErr = 0.;
    //         unsigned char  maxTb = -1;
    //         int  key = -1;
    //         float phiInner = -999.9, phiOuter = -999.9;

    //         for (std::vector< StFstCluster * >::iterator clusterIter = clusterCollectionPtr->getClusterVec().begin(); clusterIter != clusterCollectionPtr->getClusterVec().end(); ++clusterIter) {
    //             LOG_INFO << "clusterIter = " << (*clusterIter) << endm;
    //             idTruth         = (*clusterIter)->getIdTruth();
    //             key             = (*clusterIter)->getKey();
    //             disk            = (*clusterIter)->getDisk();
    //             wedge           = (*clusterIter)->getWedge();
    //             sensor          = (*clusterIter)->getSensor();
    //             apv             = (*clusterIter)->getApv();
    //             meanRStrip      = (*clusterIter)->getMeanRStrip();
    //             meanPhiStrip    = (*clusterIter)->getMeanPhiStrip();
    //             maxTb           = (*clusterIter)->getMaxTimeBin();
    //             charge          = (*clusterIter)->getTotCharge();
    //             chargeErr       = (*clusterIter)->getTotChargeErr();
    //             nRawHits        = (*clusterIter)->getNRawHits();
    //             nRawHitsR       = (*clusterIter)->getNRawHitsR();
    //             nRawHitsPhi     = (*clusterIter)->getNRawHitsPhi();
    //             // nClusteringType = (*clusterIter)->getClusteringType();

    //             LOG_INFO << "Making new Fst Hit on disk: " << disk << endm;
    //             StFstHit *newHit = new StFstHit(disk, wedge, sensor, apv, charge, chargeErr, maxTb, meanRStrip, meanPhiStrip, nRawHits, nRawHitsR, nRawHitsPhi);
    //             newHit->setId(key);
    //             newHit->setIdTruth(idTruth);
    //             LOG_INFO << "MADE " << newHit << endm;


    //             int moduleIdx;
    //             if(disk == 1) // Disk 1
    //                 moduleIdx = wedge;
    //             else if(disk == 2)// Disk 2
    //                 moduleIdx = wedge-12;
    //             else if(disk == 3)// Disk 3
    //                 moduleIdx = wedge-24;
    //             // The simple transformation will be updated with the geomtry table in database later
    //             if(disk == 1 || disk == 3)
    //             {// Disk 1 & 3
    //                 phiInner = kFstphiStart[moduleIdx-1]*TMath::Pi()/6.0 + 0.5*kFstzDirct[moduleIdx-1]*kFstStripPitchPhi;
    //                 phiOuter = kFstphiStop[moduleIdx-1]*TMath::Pi()/6.0  - 0.5*kFstzDirct[moduleIdx-1]*kFstStripPitchPhi;
    //             }
    //             else if(disk == 2)
    //             { // Disk 2
    //                 phiInner = kFstphiStop[moduleIdx-1]*TMath::Pi()/6.0  - 0.5*kFstzDirct[moduleIdx-1]*kFstStripPitchPhi;
    //                 phiOuter = kFstphiStart[moduleIdx-1]*TMath::Pi()/6.0 + 0.5*kFstzDirct[moduleIdx-1]*kFstStripPitchPhi;
    //             }
    //             double local[3];
    //             if(meanRStrip < kFstNumRStripsPerWedge/2)
    //             { // inner
    //                 local[0] = kFstrStart[meanRStrip] + 0.5*kFstStripPitchR;
    //                 local[1] = phiInner + kFstzFilp[disk-1]*kFstzDirct[moduleIdx-1]*meanPhiStrip*kFstStripPitchPhi;
    //             }
    //             else
    //             {// outer
    //                 if(sensor == 1){
    //                     local[0] = kFstrStart[meanRStrip] + 0.5*kFstStripPitchR;
    //                     local[1] = phiOuter - kFstzFilp[disk-1]*kFstzDirct[moduleIdx-1]*meanPhiStrip*kFstStripPitchPhi - kFstzFilp[disk-1]*kFstzDirct[moduleIdx-1]*0.5*kFstStripGapPhi;
    //                 }
    //                 if(sensor == 2){
    //                     local[0] = kFstrStart[meanRStrip] + 0.5*kFstStripPitchR;
    //                     local[1] = phiOuter - kFstzFilp[disk-1]*kFstzDirct[moduleIdx-1]*meanPhiStrip*kFstStripPitchPhi + kFstzFilp[disk-1]*kFstzDirct[moduleIdx-1]*0.5*kFstStripGapPhi;
    //                 }
    //             }

    //             LOG_INFO << "Setting Fst Hit Position (local)" << endm;
    //             if(disk == 1) local[2] = 151.750; //unit: cm
    //             else if(disk == 2) local[2] = 165.248; //unit: cm
    //             else if(disk == 3) local[2] = 178.781; //unit: cm
    //             newHit->setLocalPosition(local[0], local[1], local[2]); //set local position on sensor

    //             float x0 = local[0] * cos( local[1] );
    //             float y0 = local[0] * sin( local[1] );
    //             mFstHits.push_back( TVector3( x0, y0, local[2] ) );

    //             // LOG_INFO << "Adding to fstHitCollection = " << fstHitCollection << endm;
    //             // fstHitCollection->addHit(newHit);

    //             // LOG_INFO << "DONE Adding to fstHitCollection = " << fstHitCollection << endm;
    //         } //cluster loop over
    //     }//end clusterCollectionPtr



    // }



    // // if ( fstHitCollection ){

    // //     LOG_INFO << "LOAD FST HITS, NOT NULL" << endm;
    // //     for ( unsigned int iw = 0; iw < kFstNumWedges; iw++ ){
    // //         StFstWedgeHitCollection * wc = fstHitCollection->wedge( iw );

    // //         for ( unsigned int is = 0; is < kFstNumSensorsPerWedge; is++ ){
    // //             StFstSensorHitCollection * sc = wc->sensor( is );

    // //             std::vector<StFstHit*> fsthits = sc->hits();

    // //             LOG_INFO << "fsthits.size() == " << fsthits.size() << endm;
    // //             for ( int ih = 0; ih < fsthits.size(); ih++ ){
    // //                 float vR = fsthits[ih]->localPosition(0);
    // //                 float vPhi = fsthits[ih]->localPosition(1);
    // //                 float vZ = fsthits[ih]->localPosition(2);

    // //                 LOG_INFO << "FST HIT: " << vR << ", " << vPhi << ", " << vZ << endm;                
    // //                 mFstHits.push_back( TVector3( 0, 0, 0)  );
    // //             }

    // //         }

    // //     }
        
    // // }

    // LOG_INFO << " FOUND " << mFstHits.size() << " FST HITS" << endm;
    // return ;



    StRnDHitCollection *rndCollection = nullptr;
    if (nullptr != event) {
        rndCollection = event->rndHitCollection();
    }
    bool siRasterizer = mFwdConfig.get<bool>( "SiRasterizer:active", false );
    if ( siRasterizer || rndCollection == nullptr ){
        LOG_DEBUG << "Loading Fst hits from GEANT with SiRasterizer" << endm;
        loadFstHitsFromGEANT( mcTrackMap, hitMap, count );
    } else {
        LOG_DEBUG << "Loading Fst hits from StEvent" << endm;
        loadFstHitsFromStEvent( mcTrackMap, hitMap, count );
    }
} // loadFstHits

void StFwdTrackMaker::loadFstHitsFromStEvent( FwdDataSource::McTrackMap_t &mcTrackMap, FwdDataSource::HitMap_t &hitMap, int count ){

    // Get the StEvent handle
    StEvent *event = (StEvent *)this->GetDataSet("StEvent");
    if (!event) 
        return;


    StRnDHitCollection *rndCollection = event->rndHitCollection();

    const StSPtrVecRnDHit &hits = rndCollection->hits();

    // we will reuse this to hold the cov mat
    TMatrixDSym hitCov3(3);
    
    for (unsigned int fsthit_index = 0; fsthit_index < hits.size(); fsthit_index++) {
        StRnDHit *hit = hits[fsthit_index];
    
        if ( hit->layer() > 6 ){
            // skip sTGC hits here
            continue;
        }

        const StThreeVectorF pos = hit->position();

        StMatrixF covmat = hit->covariantMatrix();

        // copy covariance matrix element by element from StMatrixF
        hitCov3(0,0) = covmat[0][0]; hitCov3(0,1) = covmat[0][1]; hitCov3(0,2) = covmat[0][2];
        hitCov3(1,0) = covmat[1][0]; hitCov3(1,1) = covmat[1][1]; hitCov3(1,2) = covmat[1][2];
        hitCov3(2,0) = covmat[2][0]; hitCov3(2,1) = covmat[2][1]; hitCov3(2,2) = covmat[2][2];

        FwdHit *fhit = new FwdHit(count++, hit->position().x(), hit->position().y(), hit->position().z(), hit->layer(), hit->idTruth(), hitCov3, mcTrackMap[hit->idTruth()]);
        // LOG_INFO << "FST HIT( " << hit->position().x() << ", " << hit->position().y() << ", " << hit->position().z() << ", " << hit->layer() << endm;
        size_t index = hit->layer()-4;
        // LOG_INFO << "FST Layer = " << hit->layer() << ", " << TString::Format("fsi%luHitMapZ", index).Data() << endm;
        if (mGenHistograms && index < 3 ){

            ((TH2*)mHistograms[TString::Format("fsi%luHitMapZ", index).Data()]) -> Fill( hit->position().x(), hit->position().y(), hit->position().z() );
        }

        // Add the hit to the hit map
        hitMap[fhit->getSector()].push_back(fhit);
        mFstHits.push_back( TVector3( hit->position().x(), hit->position().y(), hit->position().z())  );

    }
} //loadFstHitsFromStEvent

void StFwdTrackMaker::loadFstHitsFromGEANT( FwdDataSource::McTrackMap_t &mcTrackMap, FwdDataSource::HitMap_t &hitMap, int count ){
    /************************************************************/
    // Load FSI Hits from GEANT
    St_g2t_fts_hit *g2t_fsi_hits = (St_g2t_fts_hit *)GetDataSet("geant/g2t_fsi_hit");

    if ( !g2t_fsi_hits )
        return;

    int nfsi = g2t_fsi_hits->GetNRows();
    

    // reuse this to store cov mat
    TMatrixDSym hitCov3(3);
    
    if ( mGenHistograms ) this->mHistograms["nHitsFSI"]->Fill(nfsi);
    // LOG_INFO << "# fsi hits = " << nfsi << endm;

    for (int i = 0; i < nfsi; i++) {

        g2t_fts_hit_st *git = (g2t_fts_hit_st *)g2t_fsi_hits->At(i);
        
        if (0 == git)
            continue; // geant hit
        
        int track_id = git->track_p;
        int volume_id = git->volume_id;  // 4, 5, 6
        int d = volume_id / 1000;        // disk id
        int w = (volume_id % 1000) / 10; // wedge id
        int s = volume_id % 10;          // sensor id
        
        int plane_id = d - 4;
        float x = git->x[0];
        float y = git->x[1];
        float z = git->x[2];

        if (mSiRasterizer->active()) {
            TVector3 rastered = mSiRasterizer->raster(TVector3(git->x[0], git->x[1], git->x[2]));
            
            if ( mGenHistograms ) {
                this->mHistograms["fsiHitDeltaR"]->Fill(std::sqrt(x * x + y * y) - rastered.Perp());
                this->mHistograms["fsiHitDeltaPhi"]->Fill(std::atan2(y, x) - rastered.Phi());
            }
            x = rastered.X();
            y = rastered.Y();
        }


        if ( mGenHistograms ) this->mHistograms["fsi_volume_id"]->Fill(d);

        if (plane_id < 3 && plane_id >= 0) {

            if ( mGenHistograms ) {
                this->mHistograms[TString::Format("fsi%dHitMap", plane_id).Data()]->Fill(x, y);
                this->mHistograms[TString::Format("fsi%dHitMapR", plane_id).Data()]->Fill(std::sqrt(x * x + y * y));
                this->mHistograms[TString::Format("fsi%dHitMapPhi", plane_id).Data()]->Fill(std::atan2(y, x) + TMath::Pi());
            }
        } else {
            continue;
        }

        hitCov3 = makeSiCovMat( TVector3( x, y, z ), mFwdConfig );
        FwdHit *hit = new FwdHit(count++, x, y, z, d, track_id, hitCov3, mcTrackMap[track_id]);
        mFstHits.push_back( TVector3( x, y, z )  );

        // Add the hit to the hit map
        hitMap[hit->getSector()].push_back(hit);
    }
} // loadFstHitsFromGEANT

size_t StFwdTrackMaker::loadMcTracks( FwdDataSource::McTrackMap_t &mcTrackMap ){
    

    LOG_INFO << "Looking for GEANT sim vertex info" << endm;
    St_g2t_vertex *g2t_vertex = (St_g2t_vertex *)GetDataSet("geant/g2t_vertex");

    if ( g2t_vertex != nullptr ) {
        // Set the MC Vertex for track fitting
        g2t_vertex_st *vert = (g2t_vertex_st*)g2t_vertex->At(0);
        mForwardTracker->setEventVertex( TVector3( vert->ge_x[0], vert->ge_x[1], vert->ge_x[2] ) );
    }

    // Get geant tracks
    St_g2t_track *g2t_track = (St_g2t_track *)GetDataSet("geant/g2t_track");

    if (!g2t_track)
        return 0;

    size_t nShowers = 0;

    mTreeData.mcN = 1;
    LOG_DEBUG << g2t_track->GetNRows() << " mc tracks in geant/g2t_track " << endm;
    if ( mGenHistograms ) this->mHistograms["nMcTracks"]->Fill(g2t_track->GetNRows());

    for (int irow = 0; irow < g2t_track->GetNRows(); irow++) {
        g2t_track_st *track = (g2t_track_st *)g2t_track->At(irow);

        if (0 == track)
            continue;

        int track_id = track->id;
        float pt2 = track->p[0] * track->p[0] + track->p[1] * track->p[1];
        float pt = std::sqrt(pt2);
        float eta = track->eta;
        float phi = std::atan2(track->p[1], track->p[0]); //track->phi;
        int q = track->charge;

        if (!mcTrackMap[track_id] ) 
            mcTrackMap[track_id] = shared_ptr<McTrack>(new McTrack(pt, eta, phi, q, track->start_vertex_p));
        
        if (mGenTree && mTreeData.mcN < MAX_TREE_ELEMENTS) {
            mTreeData.mcPt.push_back( pt );
            mTreeData.mcEta.push_back( eta );
            mTreeData.mcPhi.push_back( phi );
            mTreeData.mcCharge.push_back( q );
            mTreeData.mcVertexId.push_back( track->start_vertex_p );

            // LOG_INFO << TString::Format( "MC Track: (pt=%f, eta=%f, phi=%f, start=%d, e=%f )", pt, eta, phi, track->start_vertex_p, track->e ) << endm;

            if (track->is_shower)
                nShowers++;

            mTreeData.mcN++;
        } else if ( mGenTree ) {
            LOG_WARN << "Truncating Mc tracks in TTree output" << endm;
        }

    } // loop on track (irow)


    // now check the Mc tracks against the McEvent filter
    size_t nForwardTracks = 0;
    size_t nForwardTracksNoThreshold = 0;
    for (auto mctm : mcTrackMap ){
        if ( mctm.second == nullptr ) continue;

        if ( mGenHistograms ){
            mHistograms[ "McEventPt" ] ->Fill( mctm.second->mPt );
            mHistograms[ "McEventEta" ] ->Fill( mctm.second->mEta );
            mHistograms[ "McEventPhi" ] ->Fill( mctm.second->mPhi );
        }

        if ( mctm.second->mEta > 2.5 && mctm.second->mEta < 4.0 ){
            
            if ( mGenHistograms ){
                mHistograms[ "McEventFwdPt" ] ->Fill( mctm.second->mPt );
                mHistograms[ "McEventFwdEta" ] ->Fill( mctm.second->mEta );
                mHistograms[ "McEventFwdPhi" ] ->Fill( mctm.second->mPhi );
            }

            nForwardTracksNoThreshold++;
            if ( mctm.second->mPt > 0.05  )
                nForwardTracks++;
        }
    } // loop on mcTrackMap

    if ( mGenHistograms ) {
        mHistograms[ "nMcTracksFwd" ]->Fill( nForwardTracks );
        mHistograms[ "nMcTracksFwdNoThreshold" ]->Fill( nForwardTracksNoThreshold );
    }


    return nForwardTracks;
} // loadMcTracks

void StFwdTrackMaker::loadFcs( ) {
    DLOG( "Loading FCS Hits and EPD hits \n" );

    StEvent *stEvent = static_cast<StEvent *>(GetInputDS("StEvent"));
    StFcsDb* fcsDb=static_cast<StFcsDb*>(GetDataSet("fcsDb"));
    if ( !stEvent || !fcsDb ){
        DLOG( "NO FCS DB!\n" );
        return;
    }
    StFcsCollection* fcsCol = stEvent->fcsCollection();
    StEpdGeom* epdgeo=new StEpdGeom;


    // LOAD ECAL / HCAL CLUSTERS
    
    for ( int idet = 0; idet  < 4; idet++ ){
        LOG_INFO << "Load FCS Clusters from Det " << idet << endm;
        StSPtrVecFcsCluster& clusters = fcsCol->clusters(idet);
        
        int nc=fcsCol->numberOfClusters(idet);

        DLOG( "FOUND %d clusters \n", nc )
        for ( int i = 0; i < nc; i++ ){
            LOG_INFO << "\tcluster[" << i << "]" << endm;
            StFcsCluster* clu = clusters[i];
            StThreeVectorD xyz = fcsDb->getStarXYZfromColumnRow(clu->detectorId(),clu->x(),clu->y());
            mFcsClusters.push_back( TVector3( xyz.x(), xyz.y(), xyz.z() ) );
            DLOG( "LOADING ECAL/HCAL Cluster\n" );
        }
    }

    // LOAD PRESHOWER HITS (EPD)
    for ( int det = 4; det < 6; det ++ ) {
        // LOG_INFO << "Load FCS Hits from Det " << det << endm;
        StSPtrVecFcsHit& hits = stEvent->fcsCollection()->hits(det);
        int nh=fcsCol->numberOfHits(det);
        for ( int i = 0; i < nh; i++ ){
            StFcsHit* hit=hits[i];
            // LOG_INFO << "\thit[" << i << "]" << endm;
            if(det==kFcsPresNorthDetId || det==kFcsPresSouthDetId){ //EPD
                 double zepd=375.0;
                 int pp,tt,id,n;
                 double x[5],y[5];
                 // LOG_INFO << "EPD Hit: " << hit->adcSum()  << ", energy: " << hit->energy() << endm;
                 if ( hit->energy() < 0.2 ) continue;
                 fcsDb->getEPDfromId(det,hit->id(),pp,tt);
                 epdgeo->GetCorners(100*pp+tt,&n,x,y);
                 double x0 = (x[0] + x[1] + x[2] + x[3]) / 4.0;
                 double y0 = (y[0] + y[1] + y[2] + y[3]) / 4.0;
                 mFcsPreHits.push_back( TVector3( x0, y0, zepd ) );
            }
        }
    }
    

    delete epdgeo;

}


//________________________________________________________________________
int StFwdTrackMaker::Make() {
    // START time for measuring tracking
    long long itStart = FwdTrackerUtils::nowNanoSecond();

    // Access forward Tracker maps
    FwdDataSource::McTrackMap_t &mcTrackMap = mForwardData->getMcTracks();
    FwdDataSource::HitMap_t &hitMap = mForwardData->getFttHits();
    FwdDataSource::HitMap_t &fsiHitMap = mForwardData->getFstHits();
    
    // clear vectors for visualization OBJ hits
    mFttHits.clear();
    mFstHits.clear();
    mFcsPreHits.clear();
    mFcsClusters.clear();
    mFwdTracks.clear();

    // default event vertex
    mForwardTracker->setEventVertex( TVector3( 0, 0, 0 ) );

    /**********************************************************************/
    // Load MC tracks
    size_t nForwardTracks = loadMcTracks( mcTrackMap );
    size_t maxForwardTracks = mFwdConfig.get<size_t>( "McEvent.Mult:max", 10000 );
    if ( nForwardTracks > maxForwardTracks ){
        LOG_DEBUG << "Skipping event with more than " << maxForwardTracks << " forward tracks" << endm;
        return kStOk;
    }

    /**********************************************************************/
    // Load sTGC 
    if ( IAttr("useFtt") ) {
        loadStgcHits( mcTrackMap, hitMap );
    }
    

    /**********************************************************************/
    // Load FST
    if ( IAttr("useFst") ) {
        loadFstHits( mcTrackMap, fsiHitMap );
    }

    /**********************************************************************/
    // Load FCS
    if ( IAttr("useFcs") ) {
        loadFcs();
    }

    /**********************************************************************/
    // Run Track finding + fitting
    LOG_INFO << ">>START Event Forward Tracking" << endm;
    mForwardTracker->doEvent();
    LOG_INFO << "<<FINISH Event Forward Tracking" << endm;

    


    /**********************************************************************/
    // Run Track finding + fitting
    
    const auto &genfitTracks = mForwardTracker -> globalTracks();
    if ( mVisualize  && genfitTracks.size() > 0 && genfitTracks.size() < 20 ) {
        const auto &seed_tracks = mForwardTracker -> getRecoTracks();

        ObjExporter woe;
        woe.output( 
            TString::Format( "ev%lu", eventIndex ).Data(), 
            seed_tracks, genfitTracks, mRaveVertices, 
            mFttHits, mFstHits, mFcsPreHits, mFcsClusters );
        eventIndex++;
        LOG_INFO << "Done Writing OBJ " << endm;
    }

    // mForwardTracker->getTrackFitter()->projectToFst( 0, genfitTracks[0] );

    // fill the ttree if we have it turned on (mGenTree)
    FillTTree();

    LOG_INFO << "Forward tracking on this event took " << (FwdTrackerUtils::nowNanoSecond() - itStart) * 1e-6 << " ms" << endm;


    StEvent *stEvent = static_cast<StEvent *>(GetInputDS("StEvent"));

    if ( IAttr("fillEvent") ) {

        if (!stEvent) {
            LOG_WARN << "No StEvent found. Forward tracks will not be saved" << endm;
            return kStWarn;
        }

        // Now fill StEvent
        // FillEvent();
        StFwdTrackCollection * ftc = stEvent->fwdTrackCollection();
        if ( !ftc ){
            ftc = new StFwdTrackCollection();
            stEvent->setFwdTrackCollection( ftc );
        }

        for ( auto *genfitTrack : genfitTracks ) {

            LOG_INFO << "Adding FwdTrack to StEvent" << endm;
            StFwdTrack *fwdTrack = new StFwdTrack( genfitTrack );
            ftc->addTrack( fwdTrack );
            // provide the projections to EPD, ECAL, HCAL
            float cov[9];
            auto tv3 = ObjExporter::trackPosition( genfitTrack, 375.0, cov );
            fwdTrack->mProjections.push_back( StFwdTrackProjection( StThreeVectorF( tv3.X(), tv3.Y(), tv3.Z() ), cov) ); // EPD

            tv3 = ObjExporter::trackPosition( genfitTrack, 715.0, cov );
            fwdTrack->mProjections.push_back( StFwdTrackProjection( StThreeVectorF( tv3.X(), tv3.Y(), tv3.Z() ), cov) ); // ECAL

            tv3 = ObjExporter::trackPosition( genfitTrack, 807.0, cov );
            fwdTrack->mProjections.push_back( StFwdTrackProjection( StThreeVectorF( tv3.X(), tv3.Y(), tv3.Z() ), cov) ); // HCAL
        }

        LOG_INFO << "StFwdTrackCollection has " << ftc->numberOfTracks() << " tracks now" << endm;
        ProcessFwdTracks();
        LOG_INFO << "DONE ProcessFwdTracks" << endm;

    } // IAttr FillEvent

    return kStOK;
} // Make

void StFwdTrackMaker::FitVertex(){
    const auto &seed_tracks = mForwardTracker -> getRecoTracks();
    const auto &genfitTracks = mForwardTracker -> globalTracks();

    if ( genfitTracks.size() > 2 ){
        genfit::GFRaveVertexFactory gfrvf;

        TMatrixDSym bscm(3);
        bscm(0, 0) = 1.1*1.1;
        bscm(1, 1) = 1.1*1.1;
        bscm(2, 2) = 10.5 * 10.5;
        gfrvf.setBeamspot( TVector3( 0, 0, 0 ), bscm );
        // std::vector< genfit::GFRaveVertex * > vertices;
        const auto &genfitTracks = mForwardTracker -> globalTracks();
        mRaveVertices.clear();
        gfrvf.findVertices( &mRaveVertices, genfitTracks, false );
        
        LOG_INFO << "mRaveVertices.size() = " << mRaveVertices.size() << endm;
        // printf( "EVENT VERTEX @(%f, %f, %f)\n", vert->ge_x[0], vert->ge_x[1], vert->ge_x[2] );
        for ( auto vert : mRaveVertices ){
            printf( "RAVE vertex @(%f, %f, %f)\n\n", vert->getPos().X(), vert->getPos().Y(), vert->getPos().Z() );
        }
    }
}


void StFwdTrackMaker::FillTTree(){

    St_g2t_vertex *g2t_vertex = (St_g2t_vertex *)GetDataSet("geant/g2t_vertex");
    if (mGenTree) {
        
        // VERTICES
        if ( g2t_vertex ){
            mTreeData.vmcN = g2t_vertex->GetNRows();
            if ( mTreeData.vmcN >= MAX_TREE_ELEMENTS ) mTreeData.vmcN = MAX_TREE_ELEMENTS;
            LOG_INFO << "Saving " << mTreeData.vmcN << " vertices in TTree" << endm;
            for ( int i = 0; i < mTreeData.vmcN; i++ ){
                g2t_vertex_st *vert = (g2t_vertex_st*)g2t_vertex->At(i);
                mTreeData.vmcX.push_back( vert->ge_x[0] );
                mTreeData.vmcY.push_back( vert->ge_x[1] );
                mTreeData.vmcZ.push_back( vert->ge_x[2] );
            }
        }

        // RAVE RECO VERTICES
        mTreeData.vrcN = mRaveVertices.size();
        if ( mTreeData.vrcN >= MAX_TREE_ELEMENTS ) mTreeData.vrcN = MAX_TREE_ELEMENTS;
        LOG_INFO << "Saving " << mTreeData.vrcN << " RAVE vertices in TTree" << endm;
        for ( int i = 0; i < mTreeData.vrcN; i++ ) {
            auto vert = mRaveVertices[i];
            mTreeData.vrcX.push_back( vert->getPos().X() );
            mTreeData.vrcY.push_back( vert->getPos().Y() );
            mTreeData.vrcZ.push_back( vert->getPos().Z() );
        }

        // if (mForwardTracker->getSaveCriteriaValues()) {
        //     for (auto crit : mForwardTracker->getTwoHitCriteria()) {
        //         string name = crit->getName();

        //         // special, save all hit info for this one
                

        //         if ( name == "Crit2_BDT" ){
        //             mTreeCrits["Crit2_BDT_DeltaPhi"].clear(); 
        //             mTreeCrits["Crit2_BDT_DeltaRho"].clear(); 
        //             mTreeCrits["Crit2_BDT_RZRatio"].clear(); 
        //             mTreeCrits["Crit2_BDT_StraightTrackRatio"].clear(); 

        //             for (auto kv : mForwardTracker->getCriteriaAllValues(name)) {
        //                 mTreeCrits["Crit2_BDT_DeltaPhi"].push_back( kv["Crit2_BDT_DeltaPhi"] );
        //                 mTreeCrits["Crit2_BDT_DeltaRho"].push_back( kv["Crit2_BDT_DeltaRho"] );
        //                 mTreeCrits["Crit2_BDT_RZRatio"].push_back( kv["Crit2_BDT_RZRatio"] );
        //                 mTreeCrits["Crit2_BDT_StraightTrackRatio"].push_back( kv["Crit2_BDT_StraightTrackRatio"] );
        //             }

        //         }

        //         if ( name == "Crit2_RZRatio" ){
        //             LOG_INFO << "allValues.size() = " << mForwardTracker->getCriteriaAllValues(name).size() << " == " << mForwardTracker->getCriteriaTrackIds(name).size() << endm;
        //             assert( mForwardTracker->getCriteriaAllValues(name).size() == mForwardTracker->getCriteriaTrackIds(name).size() && " Crit lengths must be equal" );
        //             mTreeCrits["Crit2_RZRatio_x1"].clear();
        //             mTreeCrits["Crit2_RZRatio_y1"].clear();
        //             mTreeCrits["Crit2_RZRatio_z1"].clear();
        //             mTreeCrits["Crit2_RZRatio_x2"].clear();
        //             mTreeCrits["Crit2_RZRatio_y2"].clear();
        //             mTreeCrits["Crit2_RZRatio_z2"].clear();

        //             mTreeCritTrackIds["Crit2_RZRatio_h1"].clear();
        //             mTreeCritTrackIds["Crit2_RZRatio_h2"].clear();
        //             mTreeCritTrackIds["Crit2_RZRatio_h3"].clear();
                    

        //             for (auto kv : mForwardTracker->getCriteriaAllValues(name)) {
        //                 mTreeCrits["Crit2_RZRatio_x1"].push_back( kv["Crit2_RZRatio_x1"] );
        //                 mTreeCrits["Crit2_RZRatio_y1"].push_back( kv["Crit2_RZRatio_y1"] );
        //                 mTreeCrits["Crit2_RZRatio_z1"].push_back( kv["Crit2_RZRatio_z1"] );

        //                 mTreeCrits["Crit2_RZRatio_x2"].push_back( kv["Crit2_RZRatio_x2"] );
        //                 mTreeCrits["Crit2_RZRatio_y2"].push_back( kv["Crit2_RZRatio_y2"] );
        //                 mTreeCrits["Crit2_RZRatio_z2"].push_back( kv["Crit2_RZRatio_z2"] );

        //                 mTreeCritTrackIds["Crit2_RZRatio_h1"].push_back( kv["Crit2_RZRatio_h1"] );
        //                 mTreeCritTrackIds["Crit2_RZRatio_h2"].push_back( kv["Crit2_RZRatio_h2"] );
        //                 mTreeCritTrackIds["Crit2_RZRatio_h3"].push_back( -1 );
        //             }
        //         }


        //         LOG_DEBUG << "Saving Criteria values from " << name << " in TTree" << endm;
        //         mTreeCrits[name].clear();
        //         mTreeCritTrackIds[name].clear();
        //         // copy by value so ROOT doesnt get lost (uses pointer to vector)
        //         for (float v : mForwardTracker->getCriteriaValues(name)) {
        //             mTreeCrits[name].push_back(v);
        //         }
        //         for (int v : mForwardTracker->getCriteriaTrackIds(name)) {
        //             mTreeCritTrackIds[name].push_back(v);
        //         }
        //     }

        //     // three hit criteria
        //     for (auto crit : mForwardTracker->getThreeHitCriteria()) {
        //         string name = crit->getName();

        //         // special, save all hit info for this one
        //         if ( name == "Crit2_RZRatio" ){
        //             LOG_INFO << "allValues.size() = " << mForwardTracker->getCriteriaAllValues(name).size() << " == " << mForwardTracker->getCriteriaTrackIds(name).size() << endm;
        //             assert( mForwardTracker->getCriteriaAllValues(name).size() == mForwardTracker->getCriteriaTrackIds(name).size() && " Crit lengths must be equal" );

        //             mTreeCritTrackIds["Crit2_RZRatio_h1"].clear();
        //             mTreeCritTrackIds["Crit2_RZRatio_h2"].clear();
        //             mTreeCritTrackIds["Crit2_RZRatio_h3"].clear();

        //             for (auto kv : mForwardTracker->getCriteriaAllValues(name)) {
        //                 mTreeCritTrackIds["Crit2_RZRatio_h1"].push_back( kv["Crit2_RZRatio_h1"] );
        //                 mTreeCritTrackIds["Crit2_RZRatio_h2"].push_back( kv["Crit2_RZRatio_h2"] );
        //                 mTreeCritTrackIds["Crit2_RZRatio_h3"].push_back( kv["Crit2_RZRatio_h3"] );
        //             }
        //         }


        //         LOG_DEBUG << "Saving Criteria values from " << name << " in TTree" << endm;
        //         mTreeCrits[name].clear();
        //         mTreeCritTrackIds[name].clear();
        //         // copy by value so ROOT doesnt get lost (uses pointer to vector)
        //         for (float v : mForwardTracker->getCriteriaValues(name)) {
        //             mTreeCrits[name].push_back(v);
        //         }
        //         for (int v : mForwardTracker->getCriteriaTrackIds(name)) {
        //             mTreeCritTrackIds[name].push_back(v);
        //         }
        //     }

        //     // clear them 
        //     mForwardTracker->clearSavedCriteriaValues();
        // }

        // SAVE RECO tracks

        mTreeData.rcN = 0;
        // Track seeds
        const auto &seedTracks = mForwardTracker -> getRecoTracks();
        const auto &fitMomenta = mForwardTracker -> getFitMomenta();
        const auto &numFstHits = mForwardTracker -> getNumFstHits();
        const auto &fitStatus  = mForwardTracker -> getFitStatus();

        LOG_INFO << "Saving " << seedTracks.size() << " seed tracks in TTree" << endm;
        for ( size_t i = 0; i < seedTracks.size(); i++ ){
            if ( i >= MAX_TREE_ELEMENTS ){
                LOG_WARN << "Truncating Reco tracks in TTree output" << endm;
                break;
            }

            int idt = 0;
            double qual = 0;
            idt = MCTruthUtils::dominantContribution(seedTracks[i], qual);

            
            mTreeData.rcQuality.push_back( qual );
            mTreeData.rcTrackId.push_back( idt );

            if ( seedTracks.size() == numFstHits.size() ){
                mTreeData.rcNumFst.push_back( numFstHits[i] );
                mTreeData.rcCharge.push_back( fitStatus[i].getCharge() );
            }

            if ( seedTracks.size() == fitMomenta.size() ){
                mTreeData.rcPt.push_back( fitMomenta[i].Pt() );
                mTreeData.rcEta.push_back( fitMomenta[i].Eta() );
                mTreeData.rcPhi.push_back( fitMomenta[i].Phi() );
            }

            mTreeData.rcN ++;
        }

        // only warn once, instead of on every iteration
        if ( seedTracks.size() != fitMomenta.size() ){                
            LOG_WARN << "Size mismatch between track seeds and track fits" << endm;
        }

        mTree->Fill();
    } // if mGenTree
}

//________________________________________________________________________
void StFwdTrackMaker::Clear(const Option_t *opts) {
    LOG_DEBUG << "StFwdTrackMaker::CLEAR" << endm;
    mForwardData->clear();

    if (mGenTree){
        mTreeData.hN = mTreeData.rcN = mTreeData.mcN = mTreeData.vmcN = mTreeData.vrcN = 0;
        mTreeData.hX.clear();
        mTreeData.hY.clear();
        mTreeData.hZ.clear();
    }
}
//________________________________________________________________________
void StFwdTrackMaker::ProcessFwdTracks(  ){

    StEvent *stEvent = static_cast<StEvent *>(GetInputDS("StEvent"));
    StFwdTrackCollection * ftc = stEvent->fwdTrackCollection();
    for ( auto fwdTrack : ftc->tracks() ){
        auto atEPD  = fwdTrack->mProjections[0];
        auto atECAL = fwdTrack->mProjections[1];
        auto atHCAL = fwdTrack->mProjections[2];

        LOG_INFO << "Track crosses EPD @ " << TString::Format( "(%0.2f+/-%0.2f, %0.2f+/-%0.2f, %0.2f+/-%0.2f)", atEPD.XYZ.x(), atEPD.dx(), atEPD.XYZ.y(), atEPD.dy(), atEPD.XYZ.z(), atEPD.dz() ) << endm;
    }
}