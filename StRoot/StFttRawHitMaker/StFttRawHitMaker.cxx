/***************************************************************************
 *
 * $Id: StFttRawHitMaker.cxx,v 1.4 2019/03/08 18:45:40 fseck Exp $
 *
 * Author: Florian Seck, April 2018
 ***************************************************************************
 *
 * Description: StFttRawHitMaker - class to fill the StEvent from DAQ reader:
 * unpack raw data & save StETofHeader & StETofDigis in StETofCollection 
 *
 ***************************************************************************
 *
 * $Log: StFttRawHitMaker.cxx,v $
 * Revision 1.4  2019/03/08 18:45:40  fseck
 * save middle value of tot bin as raw tot of the digi
 *
 * Revision 1.3  2019/02/19 20:32:09  fseck
 * update for unpacking year 2019+ daq files
 *
 * Revision 1.2  2018/07/27 13:58:12  fseck
 * small change to compile also in 64bit mode
 *
 * Revision 1.1  2018/07/25 14:39:40  jeromel
 * Peer reviewed Raghav+Jerome - code from Florian Seck
 *
 *
 ***************************************************************************/
#include <vector>
#include <map>
#include <array>
#include <algorithm>    // std::is_sorted

#include "StRTSBaseMaker.h"
#include "StDAQMaker/StDAQReader.h"
#include "StRtsTable.h"

#include "StEvent.h"

#include "StFttRawHitMaker.h"


#include "StEvent/StFttRawHit.h"
#include "StEvent/StEvent.h"
#include "StEvent/StFttCollection.h"


//_____________________________________________________________
StFttRawHitMaker::StFttRawHitMaker( const char* name )
: StRTSBaseMaker( "stgc", name ),
  // mEvent( 0 ),          /// pointer to StEvent
  mRunYear( 0 ),        /// year in which the data was taken (switch at 1st Oct)
  mDebug( false )       /// print out of all full messages for debugging
{
    LOG_DEBUG << "StFttRawHitMaker::ctor"  << endm;
}

//_____________________________________________________________
StFttRawHitMaker::~StFttRawHitMaker()
{  /* no op */

}

//_____________________________________________________________
Int_t
StFttRawHitMaker::Init()
{
    LOG_INFO << "StFttRawHitMaker::Init" << endm;

    return kStOk;
}

//_____________________________________________________________
Int_t
StFttRawHitMaker::InitRun( Int_t runnumber )
{ 
    mRunYear = ( runnumber + 727000 ) / 1000000 + 1999;

    LOG_INFO << "runnumber: " << runnumber << "  --> year: " << mRunYear << endm;

    return kStOk;
}

//_____________________________________________________________
Int_t
StFttRawHitMaker::FinishRun( Int_t runnumber )
{ 
    return kStOk;
}

//_____________________________________________________________
Int_t
StFttRawHitMaker::Finish()
{ 
    return kStOk;
}


//_____________________________________________________________
Int_t
StFttRawHitMaker::Make()
{ 
    LOG_INFO << "StFttRawHitMaker::Make()" << endm;


    mEvent = (StEvent*)GetInputDS("StEvent");
    if(mEvent) {
        LOG_DEBUG<<"Found StEvent"<<endm;
    } else {
        mEvent = new StEvent();
        AddData(mEvent);
        LOG_INFO <<"Added StEvent"<<endm;
    }
    mFttCollection=mEvent->fttCollection();
    if(!mFttCollection) {
        mFttCollection = new StFttCollection();
        mEvent->setFttCollection(mFttCollection);
        LOG_INFO <<"Added StFttCollection"<<endm;
    } else {
        mFttCollection = mEvent->fttCollection();
        LOG_DEBUG <<"Found StFttCollection"<<endm;
    }

    mRawVMM.clear();

    StRtsTable* daqdta;

    while (daqdta = GetNext( "vmm" ) ) {

        if ( daqdta == nullptr ) {
            LOG_WARN << "StFttRawHitMaker::Make() - NO STGC DATA found in event" << endm;
            return kStOk;
        }


        // do unpacking of the raw data
        int inputSizeBytes = daqdta->GetSize();
        LOG_DEBUG << "InputSize (bytes): " << inputSizeBytes << endm;
        LOG_DEBUG << "Sector: " << daqdta->Sector() << endm;
        LOG_DEBUG << "Pad: " << daqdta->Pad() << endm;
        LOG_DEBUG << "Row: " << daqdta->Row() << endm;
        LOG_DEBUG << "Rdo: " << daqdta->Rdo() << endm;

        int rdo = daqdta->Rdo();
        int sec = daqdta->Sector();

        if( mDebug ) {
            LOG_DEBUG << "InputSize (bytes): " << inputSizeBytes << endm;
        }

        LOG_DEBUG << "ROWS: " << daqdta->GetNRows() << endm;

        
        struct stgc_vmm_t *the_vmm = (stgc_vmm_t *)  ( *daqdta->begin() );
        int len = 1;
        LOG_DEBUG << "size: " << sizeof(the_vmm) << ", size: " << sizeof(stgc_vmm_t) << endm;
        
        // if ( daqdta->Sector() == 2 && daqdta->Rdo() == 2 )
        //     PrintTheVMM( the_vmm );

        for (auto it = daqdta->begin(); it != daqdta->end(); ++it) {
            
            struct stgc_vmm_t *vmm = (stgc_vmm_t *)  ( *it );
            u_char feb = vmm[0].feb_vmm >> 2 ;  // feb [0..5]
            u_char vm = vmm[0].feb_vmm & 3 ;    // VMM [0..3]


            StFttRawHit *hit = new StFttRawHit( sec, rdo, feb, vm, vmm[0].ch, vmm[0].adc, vmm[0].bcid, vmm[0].tb );
            mFttCollection->addRawHit( hit );
            // printf("  FEB %d:%d, ch %02d: ADC %d, BCID %d, TB %d\n",feb,vm,vmm[0].ch,vmm[0].adc,vmm[0].bcid, vmm[0].tb) ;
            // mRawVMM.push_back( *vmm );
        }
    }

    return kStOk;
}
