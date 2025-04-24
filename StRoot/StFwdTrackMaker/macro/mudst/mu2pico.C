//usr/bin/env root4star -l root -l -q  $0; exit $?
// that is a valid shebang to run script as executable, but with only one arg

void loadLibs();
void mu2pico( const Char_t * fileList = "/star/data19/reco/forwardCrossSection_2022/ReversedFullField/P25ia/2022/055/23055058/st_physics_23055058_raw_3500001.MuDst.root", int firstEvent = 0, int nEvents = 10000 ){
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
        cout << endl << "============  Data Base =========" << endl;
        St_db_Maker *dbMk = new St_db_Maker("db","MySQL:StarDb","$STAR/StarDb","StarDb");
        dbMk->SetDateTime(20220225, 0);
        // things will run fine without a timestamp set, but FCS DB will give bad values ...
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
		StFcsDbMaker * fcsDb = new StFcsDbMaker();
		chain->AddMaker(fcsDb);
		// fcsDb->SetDebug();

        gSystem->Load("libStFcsWaveformFitMaker.so");
		gSystem->Load("libStFcsClusterMaker.so");
		
		StFcsWaveformFitMaker *fcsWFF = new StFcsWaveformFitMaker();
		fcsWFF->setEnergySelect(0);

		StFcsClusterMaker *fcsclu = new StFcsClusterMaker();
	/*******************************************************************************************/

    

    // The PicoDst
    gSystem->Load("libStPicoEvent");
    gSystem->Load("libStPicoDstMaker");
    StPicoDstMaker *picoMk = (StMaker*) (new StPicoDstMaker(StPicoDstMaker::IoWrite, inMuDstFile, "picoDst"));
    cout << "picoMk = " << picoMk << endl;
    picoMk->setVtxMode(StPicoDstMaker::Vtxless);
    // chain->SetAttr("TpcVpdVzDiffCut", 600.0);
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

	/*******************************************************************************************/
    // MAIN EVENT LOOP
    /*******************************************************************************************/
	size_t nEntries = muDstChain.GetEntries();
	size_t numProcessed = 0;
	for (int i = 0; i < nEntries; i++) {
		printf("Processing event %d of %d\n", i, nEntries);
        chain->Clear();
        if (kStOK != chain->Make())
            break;

        auto mMuDst = muDstMaker->muDst();
        printf( "BEFORE muprimv: %d \n", mMuDst->primaryVertex());
        StEvent *stEvent = static_cast<StEvent *>(muDstMaker->GetInputDS("StEvent"));
        // muDstMaker->fillVertices(stEvent);
        if (stEvent) {
            printf("StPrimaryVertex::numberOfPrimaryVertices = %d \n", stEvent->numberOfPrimaryVertices());
            printf( "AFTER muprimv: %d \n", mMuDst->primaryVertex());
        }


        cout << "EVENT #" << i << " COMPLETED" << endl; 
    }
	/*******************************************************************************************/

	// Chain Finish
	if (nEntries > 1) {
		cout << "FINISH up" << endl;
		chain->Finish();
	}

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
    gSystem->Load("libStEpdUtil.so");
	/*******************************************************************************************/


}
