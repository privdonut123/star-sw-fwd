/***************************************************************************
 *
 * $Id: StFttClusterMaker.cxx,v 1.4 2019/03/08 18:45:40 fseck Exp $
 *
 * Author: Florian Seck, April 2018
 ***************************************************************************
 *
 * Description: StFttClusterMaker - class to fill the StEvent from DAQ reader:
 * unpack raw data & save StETofHeader & StETofDigis in StETofCollection 
 *
 ***************************************************************************/
#include <vector>
#include <map>
#include <array>
#include <algorithm>    // std::is_sorted


#include "StEvent.h"

#include "StFttClusterMaker.h"


#include "StEvent/StFttRawHit.h"
#include "StEvent/StEvent.h"
#include "StEvent/StFttCollection.h"

#include "StFttDbMaker/StFttDb.h"


//_____________________________________________________________
StFttClusterMaker::StFttClusterMaker( const char* name )
: StMaker( name ),
  mEvent( 0 ),          /// pointer to StEvent
  mRunYear( 0 ),        /// year in which the data was taken (switch at 1st Oct)
  mDebug( false ),       /// print out of all full messages for debugging
  mFttDb( nullptr )
{
    LOG_DEBUG << "StFttClusterMaker::ctor"  << endm;
}

//_____________________________________________________________
StFttClusterMaker::~StFttClusterMaker()
{  /* no op */

}

//_____________________________________________________________
Int_t
StFttClusterMaker::Init()
{
    LOG_INFO << "StFttClusterMaker::Init" << endm;

    return kStOk;
}

//_____________________________________________________________
Int_t
StFttClusterMaker::InitRun( Int_t runnumber )
{ 
    mRunYear = ( runnumber + 727000 ) / 1000000 + 1999;

    LOG_INFO << "runnumber: " << runnumber << "  --> year: " << mRunYear << endm;

    return kStOk;
}

//_____________________________________________________________
Int_t
StFttClusterMaker::FinishRun( Int_t runnumber )
{ 
    return kStOk;
}

//_____________________________________________________________
Int_t
StFttClusterMaker::Finish()
{ 
    return kStOk;
}


//_____________________________________________________________
Int_t
StFttClusterMaker::Make()
{ 
    LOG_INFO << "StFttClusterMaker::Make()" << endm;


    mEvent = (StEvent*)GetInputDS("StEvent");
    if(mEvent) {
        LOG_DEBUG<<"Found StEvent"<<endm;
    } else {
        return kStOk;
    }
    mFttCollection=mEvent->fttCollection();
    if(!mFttCollection) {
        return kStOk;
    }

    mFttDb = static_cast<StFttDb*>(GetDataSet("fttDb"));

    assert( mFttDb );
    ApplyHardwareMap();

    // net we need to sort the hits into 1D projections
    // process 1 quadrant at a time,
    // process horizontal, vertical or diagonal strips one at a time

    // key == ROB
    std::map< UChar_t, std::vector<StFttRawHit *> > hStripsPerRob;
    std::map< UChar_t, std::vector<StFttRawHit *> > vStripsPerRob;
    std::map< UChar_t, std::vector<StFttRawHit *> > dStripsPerRob;

    bool hasAny = false;
    for ( StFttRawHit* hit : mFttCollection->rawHits() ) {
        UChar_t rob = mFttDb->rob( hit );
        StFttDb::StripOrientation so = mFttDb->orientation( hit );
        if ( StFttDb::StripOrientation::Horizontal == so ){
            hStripsPerRob[ rob ].push_back(hit);
            hasAny = true;
        }
        if ( StFttDb::StripOrientation::Vertical   == so ){
            vStripsPerRob[ rob ].push_back(hit);
            hasAny = true;
        }
        if ( StFttDb::StripOrientation::Diagonal   == so ){
            dStripsPerRob[ rob ].push_back(hit);
            hasAny = true;
        }
    }

    if ( hasAny ){
    for ( UChar_t iRob = 1; iRob < StFttDb::nRob+1; iRob++ ){
        LOG_INFO << "ROB=" << (int)iRob << " has " << hStripsPerRob[iRob].size() << " horizontal, "
            << vStripsPerRob[iRob].size() << " vertical, "
            << dStripsPerRob[iRob].size() << " diagonal, "
            << " strips hit" << endm;
    }
    }


    return kStOk;
}



void StFttClusterMaker::ApplyHardwareMap(){
    for ( StFttRawHit* rawHit : mFttCollection->rawHits() ) {
        mFttDb->hardwareMap( rawHit );
    }
}