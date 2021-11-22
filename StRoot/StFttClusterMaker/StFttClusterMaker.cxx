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
#include "StEnumerations.h"

#include "StFttClusterMaker.h"


#include "StEvent/StFttRawHit.h"
#include "StEvent/StFttCluster.h"
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

    InjectTestData();

    // net we need to sort the hits into 1D projections
    // process 1 quadrant at a time,
    // process horizontal, vertical or diagonal strips one at a time

    // key == ROB
    std::map< UChar_t, std::vector<StFttRawHit *> > hStripsPerRob;
    std::map< UChar_t, std::vector<StFttRawHit *> > vStripsPerRob;
    std::map< UChar_t, std::vector<StFttRawHit *> > dStripsPerRob;

    size_t nStripsHit = 0;
    for ( StFttRawHit* hit : mFttCollection->rawHits() ) {
        UChar_t rob = mFttDb->rob( hit );
        UChar_t so = mFttDb->orientation( hit );

        if ( kFttHorizontal == so ){
            hStripsPerRob[ rob ].push_back(hit);
            // LOG_INFO << "HORIZONTAL @ ROB = " << (int) rob << endm;
            nStripsHit++;
        }
        if ( kFttVertical   == so ){
            vStripsPerRob[ rob ].push_back(hit);
            // LOG_INFO << "VERTICAL @ ROB = " << (int) rob << endm;
            nStripsHit++;
        }
        if ( kFttDiagonal   == so ){
            dStripsPerRob[ rob ].push_back(hit);
            // LOG_INFO << "DIAGONAL @ ROB = " << (int) rob << endm;
            nStripsHit++;
        }
    } // loop on hit

    if ( nStripsHit > 0 ){ // could make more strict?
        for ( UChar_t iRob = 1; iRob < StFttDb::nRob+1; iRob++ ){
            LOG_INFO << "ROB=" << (int)iRob << " has " << hStripsPerRob[iRob].size() << " horizontal, "
                << vStripsPerRob[iRob].size() << " vertical, "
                << dStripsPerRob[iRob].size() << " diagonal, "
                << " strips hit" << endm;

            auto hClusters = FindClusters( hStripsPerRob[iRob], (UChar_t)kFttHorizontal );
            // Add them to StEvent  
            for ( StFttCluster * clu : hClusters ){
                mFttCollection->addCluster( clu );
            }
            auto vClusters = FindClusters( hStripsPerRob[iRob], (UChar_t)kFttVertical );
            // Add them to StEvent  
            for ( StFttCluster * clu : vClusters ){
                mFttCollection->addCluster( clu );
            }
            auto dClusters = FindClusters( hStripsPerRob[iRob], (UChar_t)kFttDiagonal );
            // Add them to StEvent  
            for ( StFttCluster * clu : dClusters ){
                mFttCollection->addCluster( clu );
            }
        } // loop on iRob
    } // nStripsHit

    return kStOk;
}

void StFttClusterMaker::InjectTestData(){
    mFttCollection->rawHits().clear();

    // TODO: inject clean strip hits to test cluster finder

}


StFttRawHit * StFttClusterMaker::FindMaxAdc( std::vector<StFttRawHit *> hits, size_t &pos ){
    auto itMax = max_element(std::begin(hits), std::end(hits));
    pos = (itMax - hits.begin());
    return *itMax;
}

std::vector<StFttCluster*> StFttClusterMaker::FindClusters( std::vector< StFttRawHit * > hits, UChar_t stripOrientattion ){
    std::vector<StFttCluster*> clusters;

    // Sort the hits by their strip position
    sort(hits.begin(), hits.end(), 
        [](const StFttRawHit* a, const StFttRawHit* b) -> bool
    { 
        return a->strip() > b->strip(); 
    });




    return clusters;
}



void StFttClusterMaker::ApplyHardwareMap(){
    for ( StFttRawHit* rawHit : mFttCollection->rawHits() ) {
        mFttDb->hardwareMap( rawHit );
    }
}