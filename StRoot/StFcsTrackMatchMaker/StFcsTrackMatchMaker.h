// \class StFcsTrackMatchMaker
// \author Akio Ogawa
//
//   This is FCS-FowardTrack matching maker
// 

#ifndef STAR_StFcsTrackMatchMaker_HH
#define STAR_StFcsTrackMatchMaker_HH

#include "StChain/StMaker.h"


class StFwdTrackCollection;
class StFwdTrack;
class StFcsCollection;
class StFcsDb;
class StEpdGeom;

class StFcsTrackMatchMaker : public StMaker{
public: 
    StFcsTrackMatchMaker(const char* name="FcsTrkMatch");
    ~StFcsTrackMatchMaker();
    int Init();
    int Make();
    int Finish();

    void setFileName(char* file){mFilename=file;} 
    void setMaxDistance(float v) {mMaxDistance=v;}
    void setMinEnergy(float v) {mMinEnergy=v;}

private:
    StFwdTrackCollection* mFwdTrkColl=0;
    StFcsCollection* mFcsColl=0;
    StFcsDb* mFcsDb=0;
    StEpdGeom* mEpdgeo=0;

    TFile* mFile=0;
    char* mFilename=0;

    float mMaxDistance=10.0;
    float mMinEnergy=0.1;
    
    TH1F* mHdx[2];
    TH1F* mHdy[2];
    TH1F* mHdr[2];
    TH1F* mNtrk[2];
    TH1F* mNclu[2];
    
    ClassDef(StFcsTrackMatchMaker,1)
};

#endif
