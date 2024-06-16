#ifndef ST_FWD_TRACK_TTREE_H
#define ST_FWD_TRACK_TTREE_H

#include "TClonesArray.h"


/** @brief TClonesArray writer
 * Helper class for writing TClonesArrays to TTree of custom class type
 */
template<class BranchType>
class TClonesArrayWriter {
    public:
	TClonesArrayWriter() {}
	~TClonesArrayWriter() {}

	void createBranch( TTree *tree, const char* name, int buffSize = 256000, int splitLevel = 99){
		_tca = new TClonesArray( BranchType().classname() );
		tree->Branch( name, &this->_tca, buffSize, splitLevel );
	}

	void add( BranchType &branch ){
		if ( nullptr == this->_tca ) return;
		BranchType *new_branch = new ((*this->_tca)[this->_n]) BranchType( );
		new_branch->copy( &branch );
		this->_n++;
	}

	void add( BranchType *branch ){
		if ( nullptr == this->_tca || nullptr == branch) return;
		BranchType *new_branch = new ((*this->_tca)[this->_n]) BranchType( );
		new_branch->copy( branch );
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

class FwdTreeHit : public TObject {
    public: 
    TString classname() { return "FwdTreeHit"; }
    FwdTreeHit() : TObject() {
        pos.SetXYZ(-1, -1, -1);
        id = 0;
        vol = 0;
        det = 0;
    }
    FwdTreeHit(float x, float y, float z, int v, int d) : TObject() {
        pos.SetXYZ(x, y, z);
        id = 0;
        vol = v;
        det = d;
    }
    int id, vol, det;
    TVector3 pos;

    void copy( FwdTreeHit *hit ){
        id = hit->id;
        vol = hit->vol;
        det = hit->det;
        pos = hit->pos;
    }

    void copy( FwdTreeHit &hit ){
        id = hit.id;
        vol = hit.vol;
        det = hit.det;
        pos = hit.pos;
    }

    ClassDef(FwdTreeHit, 1)
};


class FwdTreeRecoTrack : public TObject {
    public:
    FwdTreeRecoTrack() : TObject() {
        id = 0;
        q = 0;
        status = 0;
        mom.SetXYZ(0, 0, 0);
    }

    virtual void Clear(){
        TObject::Clear();
        seeds.clear();
    }

    /**
    * Copies the values of the given FwdTreeRecoTrack object to the current object.
    *
    * @param hit The FwdTreeRecoTrack object to copy from.
    */
    void copy( FwdTreeRecoTrack *hit ){
        id = hit->id;
        q = hit->q;
        status = hit->status;
        mom = hit->mom;
        seeds = hit->seeds;
    }

    /**
    * Copies the values of the given FwdTreeRecoTrack object to the current object.
    *
    * @param hit The FwdTreeRecoTrack object to copy from.
    */
    void copy( FwdTreeRecoTrack &hit ){
        id = hit.id;
        q = hit.q;
        status = hit.status;
        mom = hit.mom;
        seeds = hit.seeds;
    }

    TString classname() { return "FwdTreeRecoTrack"; }
    int id, q, status;
    TVector3 mom;
    vector<TVector3> seeds;
    ClassDef(FwdTreeRecoTrack, 2);
};


const size_t MAX_TREE_ELEMENTS = 4000;
struct FwdTreeData {

    /** @brief Primary event vertex*/
    TVector3 pv;
    TClonesArrayWriter<FwdTreeHit> ftt;
    TClonesArrayWriter<FwdTreeHit> fst;
    TClonesArrayWriter<FwdTreeRecoTrack> reco;

    // fttX, fttY, fttZ;
    // vector<int> fttVolumeId;
    // Note: Below are only avalaible for hits if MC
    // vector<float> fttPt;
    // vector<int> fttTrackId, fttVertexId;

    /** @brief Fst hit related info*/
    int fstN;
    vector<float> fstX, fstY, fstZ;
    vector<int> fstTrackId;

    /** @brief Fcs hit related info*/
    int fcsN;
    vector<float> fcsX, fcsY, fcsZ, fcsE;
    vector<int> fcsDet;

    /** @brief EPD hit related info */
    vector<float> epdX, epdY, epdZ, epdE;

    /** @brief RC track related info*/
    int rcN;
    vector<float> rcPt, rcEta, rcPhi, rcQuality, rcEpdX, rcEpdY;
    vector<int> rcStatus, rcTrackId, rcNumFST, rcCharge, rcNumFTT, rcNumPV;

    /** @brief MC Track related info*/
    int mcN;
    vector<float> mcPt, mcEta, mcPhi;
    vector<int> mcVertexId, mcCharge;
    vector<int> mcNumFtt, mcNumFst;

    /** @brief MC Vertex related info*/
    int vmcN;
    vector<float> vmcX, vmcY, vmcZ;

    /** @brief Track Projection related info*/
    int tprojN;
    vector<float> tprojX, tprojY, tprojZ;
    vector<float> tprojPx, tprojPy, tprojPz;
    vector<int> tprojIdD, tprojIdT, tprojStatus;

    /** @brief Track Seed Point related info*/
    int seedN;
    vector<float> seedX, seedY, seedZ;
    vector<int> seedIdT; // this and the one below are not garanteed to be the same

    vector<float> fitSeedX, fitSeedY, fitSeedZ;
    vector<int> fitSeedIdT;
    vector<int> fitStatus, fitFailedPoints;

    /** @brief RAVE RC Vertex related info*/
    int vrcN;
    vector<float> vrcX, vrcY, vrcZ;

    /** @brief Track-to-hit delta related info*/
    int thdN;
    vector<float> thdX, thdY, thdP, thdR, thaX, thaY, thaZ;

    /** @brief Seed finding Criteria related info*/
    bool saveCrit = false;
    std::map<string, std::vector<float>> Crits;
    std::map<string, std::vector<int>> CritTrackIds;

    void clear();
};


#endif