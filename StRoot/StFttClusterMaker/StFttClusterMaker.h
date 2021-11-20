/***************************************************************************
 *
 * $Id: StFttClusterMaker.h,v 1.3 2019/02/19 20:32:05 fseck Exp $
 *
 * Author: Florian Seck, April 2018
 ***************************************************************************
 *
 * Description: StFttClusterMaker - class to fill the StEvent from DAQ reader:
 * unpack raw data & save StETofHeader & StETofDigis in StETofCollection 
 *
 ***************************************************************************/
#ifndef STFTTCLUSTERMAKER_H
#define STFTTCLUSTERMAKER_H
#include "StMaker.h"

class StFttDb;
class StEvent;
class StFttCollection;

class StFttClusterMaker: public StMaker {

public:
    StFttClusterMaker( const char* name = "stgcCluster" );

    ~StFttClusterMaker();


    Int_t  Init();
    Int_t  InitRun( Int_t );
    Int_t  FinishRun( Int_t );
    Int_t  Finish();
    Int_t  Make();

private:
    void ApplyHardwareMap();
    StEvent*             mEvent;
    StFttCollection*     mFttCollection;
    Int_t                mRunYear;
    Bool_t               mDebug;
    StFttDb*             mFttDb;


    ClassDef( StFttClusterMaker, 1 )
};

#endif // STFTTCLUSTERMAKER_H