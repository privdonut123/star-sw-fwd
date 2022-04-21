#ifndef ST_FWD_TRACK_MAKER_H
#define ST_FWD_TRACK_MAKER_H

#include "StChain/StMaker.h"

#ifndef __CINT__
#include "GenFit/Track.h"
#endif

#include "FwdTrackerConfig.h"
#include "TVector3.h"

namespace KiTrack {
class IHit;
};

namespace genfit {
  class Track;
  class GFRaveVertex;
}

class ForwardTracker;
class FwdDataSource;
class FwdHit;
class StarFieldAdaptor;

class StGlobalTrack;
class StRnDHitCollection;
class StTrack;
class StTrackDetectorInfo;
class SiRasterizer;
class McTrack;

// ROOT includes
#include "TNtuple.h"
#include "TTree.h"
// STL includes
#include <vector>
#include <memory>



class StFwdTrack;

const size_t MAX_TREE_ELEMENTS = 4000;
struct FwdTreeData {


    // hits;
    int hN;
    vector<float> hX, hY, hZ, hPt;
    vector<int> hTrackId, hVolumeId, hVertexId;

    // RC tracks
    int rcN;
    vector<float> rcPt, rcEta, rcPhi, rcQuality;
    vector<int> rcTrackId, rcNumFst, rcCharge;

    // MC Tracks
    int mcN;
    vector<float> mcPt, mcEta, mcPhi;
    vector<int> mcVertexId, mcCharge;

    int vmcN;
    vector<float> vmcX, vmcY, vmcZ;
    int vrcN;
    vector<float> vrcX, vrcY, vrcZ;

    std::map<string, std::vector<float>> Crits;
    std::map<string, std::vector<int>> CritTrackIds;

};


class StFwdTrackMaker : public StMaker {

    ClassDef(StFwdTrackMaker, 0);

  public:
    StFwdTrackMaker();
    ~StFwdTrackMaker(){/* nada */};

    int Init();
    int Finish();
    int Make();
    void Clear(const Option_t *opts = "");

    enum { kInnerGeometry,
           kOuterGeometry };

    void SetConfigFile(std::string n) {
        mConfigFile = n;
    }
    void SetGenerateHistograms( bool _genHisto ){ mGenHistograms = _genHisto; }
    void SetGenerateTree(bool _genTree) { mGenTree = _genTree; }
    void SetVisualize( bool _viz ) { mVisualize = _viz; }

    vector<StFwdTrack*> mFwdTracks;

  private:
  protected:

    // Track Seed typdef 
    typedef std::vector<KiTrack::IHit *> Seed_t;

    
    // for Wavefront OBJ export
    size_t eventIndex = 0;
    

    bool mGenHistograms = false;
    bool mGenTree = false;
    std::string mConfigFile;


    std::map<std::string, TH1 *> mHistograms;
    TFile *mTreeFile = nullptr;
    TTree *mTree     = nullptr;
    FwdTreeData mTreeData;

    bool mVisualize = true;
    vector<TVector3> mFttHits;
    vector<TVector3> mFstHits;
    vector<TVector3> mFcsClusters;
    vector<TVector3> mFcsPreHits;

    std::vector< genfit::GFRaveVertex * > mRaveVertices;

    // // elements used only if the mGenTree = true
    // float mTreeX[MAX_TREE_ELEMENTS], mTreeY[MAX_TREE_ELEMENTS], mTreeZ[MAX_TREE_ELEMENTS], mTreeHPt[MAX_TREE_ELEMENTS];
    // int mTreeN, mTreeTID[MAX_TREE_ELEMENTS], mTreeVID[MAX_TREE_ELEMENTS], mTreeHSV[MAX_TREE_ELEMENTS];

    // int mTreeNTracks, mTreeRNTracks, mTreeRTID[MAX_TREE_ELEMENTS], mTreeVertID[MAX_TREE_ELEMENTS];
    // unsigned short mTreeRNumFst[MAX_TREE_ELEMENTS];
    // short mTreeQ[MAX_TREE_ELEMENTS], mTreeRQ[MAX_TREE_ELEMENTS];
    // float mTreePt[MAX_TREE_ELEMENTS], mTreeEta[MAX_TREE_ELEMENTS], mTreePhi[MAX_TREE_ELEMENTS];
    // float mTreeRPt[MAX_TREE_ELEMENTS], mTreeREta[MAX_TREE_ELEMENTS], mTreeRPhi[MAX_TREE_ELEMENTS], mTreeRQual[MAX_TREE_ELEMENTS];

    // // MC EVENT Vertices (Seed)
    // int mTreeNVert;
    // float mTreeVertX[MAX_TREE_ELEMENTS], mTreeVertY[MAX_TREE_ELEMENTS], mTreeVertZ[MAX_TREE_ELEMENTS];

    // // RC RAVE Verts
    // int mTreeNRave;
    // float mTreeRaveX[MAX_TREE_ELEMENTS], mTreeRaveY[MAX_TREE_ELEMENTS], mTreeRaveZ[MAX_TREE_ELEMENTS];
    // std::vector< genfit::GFRaveVertex * > mRaveVertices;

    // std::map<string, std::vector<float>> mTreeCrits;
    // std::map<string, std::vector<int>> mTreeCritTrackIds;

    void ProcessFwdTracks();

    // I could not get the library generation to succeed with these.
    // so I have removed them
    #ifndef __CINT__
        std::shared_ptr<SiRasterizer> mSiRasterizer;
        FwdTrackerConfig mFwdConfig;
        std::shared_ptr<ForwardTracker> mForwardTracker;
        std::shared_ptr<FwdDataSource> mForwardData;
        size_t loadMcTracks( std::map<int, std::shared_ptr<McTrack>> &mcTrackMap );
        void loadFcs();
        void loadStgcHits( std::map<int, std::shared_ptr<McTrack>> &mcTrackMap, std::map<int, std::vector<KiTrack::IHit *>> &hitMap, int count = 0 );
        void loadStgcHitsFromGEANT( std::map<int, std::shared_ptr<McTrack>> &mcTrackMap, std::map<int, std::vector<KiTrack::IHit *>> &hitMap, int count = 0 );
        void loadStgcHitsFromStEvent( std::map<int, std::shared_ptr<McTrack>> &mcTrackMap, std::map<int, std::vector<KiTrack::IHit *>> &hitMap, int count = 0 );
        void loadFstHits( std::map<int, std::shared_ptr<McTrack>> &mcTrackMap, std::map<int, std::vector<KiTrack::IHit *>> &hitMap, int count = 0 );
        void loadFstHitsFromGEANT( std::map<int, std::shared_ptr<McTrack>> &mcTrackMap, std::map<int, std::vector<KiTrack::IHit *>> &hitMap, int count = 0 );
        void loadFstHitsFromStEvent( std::map<int, std::shared_ptr<McTrack>> &mcTrackMap, std::map<int, std::vector<KiTrack::IHit *>> &hitMap, int count = 0 );
    #endif

    void FillTTree(); // if debugging ttree is turned on (mGenTree)
    void FitVertex();
    // Fill StEvent
    void FillEvent();
    void FillDetectorInfo(StTrackDetectorInfo *info, const genfit::Track *track, bool increment);
    void FillTrack(StTrack *otrack, const genfit::Track *itrack, const Seed_t &iseed, StTrackDetectorInfo *info);
    void FillTrackFlags(StTrack *otrack, const genfit::Track *itrack);
    void FillTrackGeometry(StTrack *otrack, const genfit::Track *itrack, double zplane, int io);
    void FillTrackDcaGeometry ( StGlobalTrack    *otrack, const genfit::Track *itrack );
    void FillTrackFitTraits(StTrack *otrack, const genfit::Track *itrack);
    void FillTrackMatches(StTrack *otrack, const genfit::Track *itrack);
};

#endif
