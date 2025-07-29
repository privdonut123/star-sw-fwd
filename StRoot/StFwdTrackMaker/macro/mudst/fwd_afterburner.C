//usr/bin/env root4star -l root -l -q  $0; exit $?
// that is a valid shebang to run script as executable, but with only one arg

// For fast fwd tracking run with Db=false, fcs=false, FwdQa=false
bool runDb = true;
bool runFttChain = false;
bool runFcsChain = true;
bool runFwdChain = true;
bool refillMuDst = false;
bool runFwdQa = false;
bool runFitQa = false;
bool runPico = true;

// For EPD QA only
// bool runDb = false;
// bool runFttChain = false;
// bool runFcsChain = true;
// bool runFwdChain = false;
// bool refillMuDst = false;
// bool runFwdQa = false;
// bool runFitQa = true;

// Minimal
bool runDb = true;
bool runFttChain = false;
bool runFcsChain = true;
bool runFwdChain = true;
bool refillMuDst = false;
bool runFwdQa = false;
bool runFitQa = false;
bool runPico = true;

#include "StMemStat.h"


void loadLibs();
void fwd_afterburner( const Char_t * fileList = "/star/data19/reco/forwardCrossSection_2022/ReversedFullField/P25ia/2022/055/23055059/st_physics_23055059_raw_2000001.MuDst.root", int firstEvent = 0, int nEvents = 1000 ){
	cout << "FileList: " << fileList << endl;
	
	cout << "firstEvent: " << firstEvent << endl;
	cout << "nEvents: " << nEvents << endl;

	// First load some shared libraries we need
	loadLibs();

	// create the chain
	StChain *chain  = new StChain("StChain");

	const char* inMuDstFile = fileList;
	// create the StMuDstMaker
	StMuDstMaker *muDstMaker = new StMuDstMaker(  	0,
													0,
													"",
													inMuDstFile,
													"MuDst.root",
													1
												);
	TChain& muDstChain = *muDstMaker.chain();
    printf( "MuDst file has %d events available in tree\n", muDstChain.GetEntries());
	
	/*******************************************************************************************/
	// Initialize the database
		if (runDb){
			cout << endl << "============  Data Base =========" << endl;
			St_db_Maker *dbMk = new St_db_Maker("db","MySQL:StarDb","$STAR/StarDb","StarDb");
			dbMk->SetDateTime(20220225, 0);
			// things will run fine without a timestamp set, but FCS DB will give bad values ...
		}
	/*******************************************************************************************/
	

	/*******************************************************************************************/
	// Create a TEventList to specify the events to process
		TEventList* eventList = new TEventList("selectedEvents");

		// Add event indices from firstEvent to firstEvent + nEvents - 1
		for (int i = firstEvent; i < firstEvent + nEvents; i++) {
			eventList->Enter(i);
		}

		// Set the event list in StMuDstMaker
		muDstMaker->SetEventList(eventList);
	/*******************************************************************************************/

	/*******************************************************************************************/
	// Create the StMuDst2StEventMaker
    StMuDst2StEventMaker * mu2ev = new StMuDst2StEventMaker();
	/*******************************************************************************************/

	/*******************************************************************************************/
	// Setup Fcs Database if needed
	if ( (runFcsChain && runDb) || runFitQa){
		StFcsDbMaker * fcsDb = new StFcsDbMaker();
		chain->AddMaker(fcsDb);
		// fcsDb->SetDebug();
	}
	/*******************************************************************************************/
	

	/*******************************************************************************************/
	// FTT chain
	if (runFttChain){
		StFttDbMaker * fttDbMk = new StFttDbMaker();
		StFttHitCalibMaker * ftthcm = new StFttHitCalibMaker();
		StFttClusterMaker * fttclu = new StFttClusterMaker();
		fttclu->SetTimeCut(1, -40, 40);
		StFttClusterPointMaker *fttCP = new StFttClusterPointMaker();
		// StFttPointMaker * fttpoint = new StFttPointMaker();
	}
	/*******************************************************************************************/

	/*******************************************************************************************/
    // FCS Chain
	if (runFcsChain){
		gSystem->Load("libStFcsWaveformFitMaker.so");
		gSystem->Load("libStFcsClusterMaker.so");
		
		StFcsWaveformFitMaker *fcsWFF = new StFcsWaveformFitMaker();
		fcsWFF->setEnergySelect(0);
		StFcsClusterMaker *fcsclu = new StFcsClusterMaker();
	}
	/*******************************************************************************************/

	/*******************************************************************************************/
	// FwdTrackMaker Chain
	StFwdTrackMaker *fwdTrack = NULL;
	if (runFwdChain){
		// FwdTrackMaker
		fwdTrack = new StFwdTrackMaker();
		fwdTrack->SetDebug(1);
		fwdTrack->setGeoCache( "fGeom.root" );
		fwdTrack->setSeedFindingWithFst();
		fwdTrack->setTrackRefit( false );

		// Fitter Options
		fwdTrack->setFitDebugLvl( 0 );
		fwdTrack->setFitMinIterations( 1 );
		fwdTrack->setFitMaxIterations( 3 );
		
		fwdTrack->setDeltaPval( 1e-3 );
		fwdTrack->setRelChi2Change( 1e-3 );

		// fwdTrack->setSeedFindingOff();
		// fwdTrack->setTrackFittingOff();
		fwdTrack->setFttHitSource( 3 /* = IGNORE */);
	}



		if (runFcsChain){
			// FwdTrack and FcsCluster assciation
			gSystem->Load("StFcsTrackMatchMaker");
			StFcsTrackMatchMaker *match = new StFcsTrackMatchMaker();
			match->setMaxDistance(6,10);
			match->setFileName("fcstrk.root");
		}

		
		
		if (runFwdQa){
			StFwdQAMaker *fwdQA = new StFwdQAMaker();
			fwdQA->SetDebug(2);
			TString fwdqaname( gSystem->BaseName(inMuDstFile) );
			fwdqaname.ReplaceAll(".MuDst.root", ".FwdTree.root");
			cout << fwdqaname.Data() << endl;
			fwdQA->setTreeFilename(fwdqaname);

			gSystem->Load("StFwdUtils.so");
			StFwdAnalysisMaker * fwdAna = new StFwdAnalysisMaker();
			fwdAna->setMuDstInput();
		}


	// The PicoDst
	if (runPico){
		gSystem->Load("libStPicoEvent");
		gSystem->Load("libStPicoDstMaker");
		StPicoDstMaker *picoMk = (StMaker*) (new StPicoDstMaker(StPicoDstMaker::IoWrite, inMuDstFile, "picoDst"));
		cout << "picoMk = " << picoMk << endl;
		picoMk->setVtxMode(StPicoDstMaker::Vtxless);
	}

	if ( runFitQa && runFwdChain){
		StFwdFitQAMaker *fwdFitQA = new StFwdFitQAMaker();
		fwdFitQA->SetDebug();
		TString fitqaoutname(gSystem->BaseName(inMuDstFile));
		fitqaoutname.ReplaceAll(".MuDst.root", ".FwdFitQA.root");
		fwdFitQA->setOutputFilename( fitqaoutname );
	}
	/*******************************************************************************************/


	/*******************************************************************************************/
	// Initialize chain
	chain->SetDebug(1);
	Int_t iInit = chain->Init();
	chain->SetDebug(1);
	cout << "CHAIN INIT DONE? (good==0): " << iInit << endl;
	// ensure that the chain initializes

	if ( iInit ) 
		chain->Fatal(iInit,"on init");

	// print the chain status
	chain->PrintInfo();

	StMemStat stmem;
	stmem.PrintMem("BEFORE Event Loop");
	/*******************************************************************************************/
    // MAIN EVENT LOOP
    /*******************************************************************************************/
	size_t nEntries = muDstChain.GetEntries();
	size_t numProcessed = 0;
	for (int i = 0; i < nEntries; i++) {
		printf("Processing event %d of %d\n", i, nEntries);
		if ( fwdTrack )
			fwdTrack->SetDebug(1);
		if (i > 0) // skip first event to make it consistent
			stmem.Start();
		chain->Clear();
		
        if (kStOK != chain->Make())
            break;

		if (refillMuDst){
			StEvent *mStEvent = static_cast<StEvent *>(muDstMaker->GetInputDS("StEvent"));
			// muDstMaker->fillFwdTrack( mStEvent);
			fwdQA->Make();
		}
		stmem.PrintMem(TString::Format("After Event %d:", i).Data());	
		if (i > 0)
			stmem.Stop();
		// MipMaker->Make();
		// picoMk->Make();
        cout << "EVENT #" << i << " COMPLETED" << endl; 
    }
	stmem.PrintMem("After Event Loop");
	stmem.Summary();
	/*******************************************************************************************/

	// Chain Finish
	// if (nEntries > 1) {
	// 	cout << "FINISH up" << endl;
	// 	chain->Finish();
	// }

	// delete chain;
}



void loadLibs(){	
	// if (gClassTable->GetID("TTable") < 0) {
	// 	gSystem->Load("libStar");
	// 	gSystem->Load("libPhysics");
	// }  
	cout << "LL0" << endl;
	gSystem->Load("libStarClassLibrary.so");
	gSystem->Load("libStarRoot.so");
	cout << "LL1" << endl;
	gROOT->LoadMacro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C");
	loadSharedLibraries();
	cout << "LL2" << endl;
	
	gSystem->Load("StarMagField");
	gSystem->Load("StMagF");
	gSystem->Load("StDetectorDbMaker");
	gSystem->Load("StTpcDb");
	gSystem->Load("StDaqLib");
	gSystem->Load("StDbBroker");
	gSystem->Load("StDbUtilities");
	gSystem->Load("St_db_Maker");

	gSystem->Load("StEvent");
	gSystem->Load("StEventMaker");
	gSystem->Load("StarMagField");
 
	gSystem->Load("libGeom");
	gSystem->Load("St_g2t");
	
	// Added for Run16 And beyond
	gSystem->Load("libGeom.so");
	
	gSystem->Load("St_base.so");
	gSystem->Load("StUtilities.so");
	gSystem->Load("libPhysics.so");
	gSystem->Load("StarAgmlUtil.so");
	gSystem->Load("StarAgmlLib.so");
	gSystem->Load("libStarGeometry.so");
	gSystem->Load("libGeometry.so");
	
	gSystem->Load("xgeometry");
 
	gSystem->Load("St_geant_Maker");


	// needed since I use the StMuTrack
	gSystem->Load("StarClassLibrary");
	gSystem->Load("StStrangeMuDstMaker");
	gSystem->Load("StMuDSTMaker");
	gSystem->Load("StBTofCalibMaker");
	gSystem->Load("StVpdCalibMaker");
	gSystem->Load("StBTofMatchMaker");
	gSystem->Load("StFcsDbMaker");	

	/*******************************************************************************************/
	// loading libraries
	gSystem->Load("StFcsDbMaker");
	gSystem->Load( "StFttDbMaker" );
	gSystem->Load( "StFttHitCalibMaker" );
	gSystem->Load( "StFttClusterMaker" );
	gSystem->Load( "StFttPointMaker" );
    gSystem->Load("libStarGeneratorUtil.so");
    gSystem->Load("libgenfit2");
    gSystem->Load("libKiTrack");
    gSystem->Load("libXMLIO.so");
    gSystem->Load( "StFwdTrackMaker.so" );
	gSystem->Load( "StFwdUtils.so" );
    gSystem->Load("libStEpdUtil.so");
	/*******************************************************************************************/


}
