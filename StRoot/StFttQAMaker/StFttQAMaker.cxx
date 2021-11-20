#include "StFttQAMaker.h"

#include "StFttRawHitMaker/StFttRawHitMaker.h"

#include "StEvent/StFttRawHit.h"
#include "StEvent/StEvent.h"
#include "StEvent/StFttCollection.h"

//_____________________________________________________________                                                       
StFttQAMaker::StFttQAMaker(const char *name):StMaker("etofCluster",name)
{                                            
    LOG_DEBUG << "StFttQAMaker::ctor"  << endm;
}
//_____________________________________________________________                                                       
StFttQAMaker::~StFttQAMaker()
{ 
}

//_____________________________________________________________                                                       
Int_t StFttQAMaker::Init()
{
    LOG_INFO << "StFttQAMaker::Init" << endm;

    BookHistograms();
    BookTree();
    return kStOk;
}
//_____________________________________________________________                                                       
Int_t StFttQAMaker::InitRun(Int_t runnumber)
{ 
    return kStOk;
}

//_____________________________________________________________                                                       
Int_t StFttQAMaker::FinishRun(Int_t runnumber)
{ 
    return kStOk;
}

//-------------------------------------------------------------                                                       
Int_t StFttQAMaker::Finish()
{ 
    LOG_INFO << "StFttQAMaker::Finish()" << endm;
    TFile *histFile = new TFile( "fttQA.root", "RECREATE" );
    histFile->cd();

    WriteHistograms();

    mVMMTree->Write();

    LOG_INFO << "StFttQAMaker::Finish() - writing *.fttQA.root ..." << endm;

    histFile->Close();
    return kStOk;
}

//_____________________________________________________________                                                       
Int_t StFttQAMaker::Make()
{ 
    LOG_INFO << "StFttQAMaker::Make()" << endm;

    mEvent = (StEvent*)GetInputDS("StEvent");
    if(mEvent) {
        LOG_DEBUG<<"Found StEvent"<<endm;
    } else {
        return kStOk;
    }
    mFttCollection=mEvent->fttCollection();
    if(!mFttCollection) {
        return kStOk;
    } else {
        LOG_DEBUG <<"Found StFttCollection"<<endm;
    }
    
    MakeRawHitQA();
    
    return kStOk;
}

//_____________________________________________________________  
void
StFttQAMaker::MakeRawHitQA(){
    

    mVMMData.N = 0;
    for ( auto rawHit : mFttCollection->rawHits() ) {
        
        mH1d[ "ADC" ]->Fill( rawHit->adc() );
        mH1d[ "BCID" ]->Fill( rawHit->bcid() );
        // mH1d[ "FEBVMM" ]->Fill( rawHit->feb_vmm );
        mH1d[ "FEB" ]->Fill( rawHit->feb() );
        mH1d[ "VMM" ]->Fill( rawHit->vmm() );
        mH1d[ "CH" ]->Fill( rawHit->channel() );
        mH1d[ "TB" ]->Fill( rawHit->tb() );

        mVMMData.SEC[mVMMData.N]    = rawHit->sector();
        mVMMData.RDO[mVMMData.N]    = rawHit->rdo();
        mVMMData.ADC[mVMMData.N]    = rawHit->adc();
        mVMMData.FEB[mVMMData.N]    = rawHit->feb();
        // mVMMData.FEBVMM[mVMMData.N] = rawHit->feb_vmm;
        mVMMData.VMM[mVMMData.N]    = rawHit->vmm();
        mVMMData.CH[mVMMData.N]     = rawHit->channel();
        mVMMData.BCID[mVMMData.N]   = rawHit->bcid();
        mVMMData.TB[mVMMData.N]     = rawHit->tb();

        mVMMData.N++;
        // LOG_INFO << "n = " << mVMMData.N << endm;
    }

    mVMMTree->Fill();
}

//_____________________________________________________________  
void
StFttQAMaker::WriteHistograms()
{
    LOG_DEBUG << "StFttQAMaker::WriteHistograms()" << endm;
    for( const auto& kv : mH1d ) {
        if( kv.second->GetEntries() > 0 ) kv.second->Write();
    }
    for( const auto& kv : mH2d ) {
        if( kv.second->GetEntries() > 0 ) kv.second->Write();
    }
}
//_____________________________________________________________  
void
StFttQAMaker::BookHistograms()
{

    mH1d[ "ADC" ]    = new TH1F( "ADC", "ADC", 1025, 0, 1025 );
    mH1d[ "BCID" ]   = new TH1F( "BCID", "BCID", 1024, 0, 4096 );
    mH1d[ "FEB" ]    = new TH1F( "FEB", "FEB", 10, 0, 10 );
    mH1d[ "FEBVMM" ] = new TH1F( "FEBVMM", "FEBVMM", 4096, 0, 4096 );
    mH1d[ "VMM" ]    = new TH1F( "VMM", "VMM", 10, 0, 10 );
    mH1d[ "CH" ]     = new TH1F( "CH", "CH", 100, 0, 100 );
    mH1d[ "TB" ]     = new TH1F( "TB", "TB", 1000, -32000, 32000 );


    for( auto& kv : mH1d ) {
        kv.second->SetDirectory( 0 );
    }

    for( auto& kv : mH2d ) {
        kv.second->SetDirectory( 0 );
    }
}

//_____________________________________________________________  
void
StFttQAMaker::BookTree()
{
    cout << " Booking the event tree " << endl;
    mVMMTree = new TTree("ftt","ftt");
    mVMMTree->SetAutoSave(100000000); // 100 MB

    // Event information
    mVMMTree->Branch("EVT"    , &mVMMData.EVT      , "EVT/I");
    mVMMTree->Branch("N"      , &mVMMData.N        , "N/I");

    // Channel information
    mVMMTree->Branch("SEC"    , mVMMData.SEC       , "SEC[N]/I");
    mVMMTree->Branch("RDO"    , mVMMData.RDO       , "RDO[N]/I");
    mVMMTree->Branch("FEB"    , mVMMData.FEB       , "FEB[N]/I");
    mVMMTree->Branch("FEBVMM" , mVMMData.FEBVMM    , "FEBVMM[N]/I");
    mVMMTree->Branch("VMM"    , mVMMData.VMM       , "VMM[N]/I");
    mVMMTree->Branch("CH"     , mVMMData.CH        , "CH[N]/I");
    mVMMTree->Branch("BCID"   , mVMMData.BCID      , "BCID[N]/I");
    mVMMTree->Branch("ADC"    , mVMMData.ADC       , "ADC[N]/I");
    mVMMTree->Branch("TB"     , mVMMData.TB        , "TB[N]/I");
    mVMMTree->Branch("ROW"    , mVMMData.ROW       , "ROW[N]/I");
    mVMMTree->Branch("STRIP"  , mVMMData.STRIP     , "STRIP[N]/I");
}