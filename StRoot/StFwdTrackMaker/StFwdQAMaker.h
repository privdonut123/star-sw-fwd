#ifndef ST_FWD_TREE_MAKER_H
#define ST_FWD_TREE_MAKER_H

#include "TClonesArray.h"
#ifndef __CINT__
#include "GenFit/Track.h"
#include "StFwdTrackMaker/include/Tracker/FwdHit.h"
#include "StMuDSTMaker/COMMON/StMuFwdTrack.h"
#endif

#include "StChain/StMaker.h"
#include "TTree.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "StEvent/StEnumerations.h"

class StMuFwdTrack;
class StMuFwdTrackProjection;
class ForwardTracker;
class StFwdTrack;

/** @brief TClonesArray writer
 * Helper class for writing TClonesArrays to TTree of custom class type
 */
template<class BranchType>
class TClonesArrayWriter {
    public:
	TClonesArrayWriter() {}
	~TClonesArrayWriter() {}

	void createBranch( TTree *tree, const char* name, int buffSize = 256000, int splitLevel = 99){
        _tca = new TClonesArray( BranchType::Class_Name() );
		tree->Branch( name, &this->_tca, buffSize, splitLevel );
	}

	void add( BranchType &branch ){
		if ( nullptr == this->_tca ) return;
		BranchType *new_branch = new ((*this->_tca)[this->_n]) BranchType( );
        *new_branch = branch;
		this->_n++;
	}

	void add( BranchType *branch ){
		if ( nullptr == this->_tca || nullptr == branch) return;
		BranchType *new_branch = new ((*this->_tca)[this->_n]) BranchType( );
        *new_branch = *branch;
		this->_n++;
	}

	void reset(){
		this->_n = 0;
		if( nullptr != this->_tca )
			this->_tca->Clear();
	}

	UInt_t N() const { return _n; }
	BranchType *at( UInt_t i ){
		if ( nullptr == _tca )
			return nullptr;
		return (BranchType*)_tca->At( i );
	}

    protected:
	TClonesArray        * _tca = nullptr;
	UInt_t                _n   = 0;
};

class FwdTreeHeader : public TObject {
    public:
    FwdTreeHeader() : TObject() {
        run = 0;
        event = 0;
        tofmult = 0;
        vpdVz = -999;
        pv.SetXYZ(0, 0, 0);
    }

    void set( int r, int e, int t, TVector3 &p ){
        run = r;
        event = e;
        tofmult = t;
        pv = p;
    }

    void clear() {
        run = 0;
        event = 0;
        tofmult = 0;
        TVector3 pv(-999, -999, -999);
        vpdVz = -999;
    }

    TVector3 pv;
    int run, event, tofmult;
    float vpdVz;

    ClassDef(FwdTreeHeader, 1)
};

class StFcsDb;
class StFcsCluster;
class StFcsHit;
class StMuFcsCluster;
class StMuFcsHit;
class StMuFttCluster;
class StMuFttPoint;
class StMuFwdTrackSeedPoint;

class FwdTreeMonteCarloTrack : public TObject {
    public:

    FwdTreeMonteCarloTrack() : TObject() {
        id = 0;
        q = 0;
        status = 0;
        mom.SetXYZ(0, 0, 0);
    }

    int id, q, status;
    TVector3 mom;

    ClassDef(FwdTreeMonteCarloTrack, 1);
};

/** @brief
* This class is a container for the data that will be written to the output tree.
*/
struct FwdTreeData {

    /** @brief Primary event vertex*/
    FwdTreeHeader header;
    TClonesArrayWriter<StMuFwdTrackSeedPoint> fttSeeds;
    TClonesArrayWriter<StMuFttPoint> fttPoints;
    TClonesArrayWriter<StMuFttCluster> fttClusters;
    TClonesArrayWriter<StMuFwdTrackSeedPoint> fstSeeds;

    TClonesArrayWriter<StMuFcsCluster> wcal;
    TClonesArrayWriter<StMuFcsHit> wcalHits;
    TClonesArrayWriter<StMuFcsCluster> hcal;
    TClonesArrayWriter<StMuFcsHit> hcalHits;
    TClonesArrayWriter<StMuFwdTrack> reco;

    int nSeedTracks;
    TClonesArrayWriter<StMuFwdTrackSeedPoint> seeds;


    void clear();
};


class StMuDstMaker;
class StMuDst;
class StMuFwdTrackCollection;
class StMuFcsCollection;
class StFwdTrackMaker;
class StEvent;

class StFwdQAMaker : public StMaker {

    ClassDef(StFwdQAMaker, 0);

  public:
    StFwdQAMaker();
    ~StFwdQAMaker(){/* nada */};

    int Init();
    int Finish();
    int Make();
    void Clear(const Option_t *opts = "");

    void FillFttClusters();
    void FillFcsStEvent();
    void FillFcsStMuDst();
    void FillTracks();

  protected:
    TFile *mTreeFile = nullptr;
    TTree *mTree     = nullptr;
    FwdTreeData mTreeData;

    StEvent *mStEvent = nullptr;
    StMuDstMaker *mMuDstMaker = nullptr;
    StMuDst *mMuDst = nullptr;
    StMuFwdTrackCollection * mMuForwardTrackCollection = nullptr;
    StMuFcsCollection *mMuFcsCollection = nullptr;
    StFwdTrackMaker *mFwdTrackMaker = nullptr;
    StFcsDb *mFcsDb = nullptr;

};


#endif
