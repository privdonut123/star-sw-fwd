#ifndef STAR_StFttQAMaker_H
#define STAR_StFttQAMaker_H


/***************************************************************************                                          
 *                                                                                                                    
 * $Id: StFttQAMaker.h,v 0.1 2017/02/21 17:50:32 tlusty Exp $                                                    
 * StFttQAMaker - class to fille the StEvent from DAQ reader                                                        
 *--------------------------------------------------------------------------                                          
 *                                                                                                                    
 ***************************************************************************/
#include "StMaker.h"

// ROOT
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"

// STL
#include <vector>

class StEvent;
class StFttCollection;

const Int_t mMax = 10000;
struct STGCVMMData
{
    // event information
    Int_t    EVT;
    Int_t    N;

    //channel information
    Int_t    SEC[mMax];
    Int_t    RDO[mMax];
    Int_t    FEB[mMax];
    Int_t    FEBVMM[mMax];
    Int_t    VMM[mMax];
    Int_t    CH[mMax];
    Int_t    BCID[mMax];
    Int_t    ADC[mMax];
    Int_t    TB[mMax];
    Int_t    ROW[mMax];
    Int_t    STRIP[mMax];

};


class StFttQAMaker: public StMaker
{
private:

public:

/// Default constructor                                                                                          
    StFttQAMaker(const char *name="fttQA");

    ~StFttQAMaker();


    Int_t  Init();
    Int_t  InitRun(Int_t);
    Int_t  FinishRun(Int_t);
    Int_t  Finish();
    Int_t  Make();

    void MakeRawHitQA();

    void WriteHistograms();
    void BookHistograms();
    void BookTree();

    StEvent*             mEvent;
    StFttCollection*     mFttCollection;

    std::map< string, TH1* >    mH1d;
    std::map< string, TH2* >    mH2d;
    TTree * mVMMTree;
    STGCVMMData mVMMData;

    ClassDef(StFttQAMaker, 1)

};

#endif