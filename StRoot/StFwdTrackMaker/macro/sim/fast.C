//usr/bin/env root4star -l -b -q $0'("'${1:-sim.fzd}'",'${2:-5}')'; exit $?
// that is a valid shebang to run script as executable, but with only one arg


// Run very fast fwd tracking
// generate some input data using genfzd
void loadLibs();
void loadLibs();
TFile *output = 0;

StMemStat stmem;

void fast(       char *inFile =  "sim.fzd",
                int n = 100, // nEvents to run
                bool useFstForSeedFinding = true, // use FTT (default) or FST for track finding
                bool enableTrackRefit = true, // Enable track refit (default off)
                bool realisticSim = true, // enables data-like mode, real track finding and fitting without MC seed
                bool useZeroB = false
            ) {
    stmem.PrintMem("STARTUP");
    // report all of the parameters passed in
    cout << "inFile = " << inFile << endl;
    cout << "n = " << n << endl;
    cout << "useFstForSeedFinding = " << useFstForSeedFinding << endl;
    cout << "enableTrackRefit = " << enableTrackRefit << endl;
    cout << "realisticSim = " << realisticSim << endl;
    cout << "useZeroB = " << useZeroB << endl;
    const char *geom = "";
    TString _geom = geom;
    bool useConstBz = false;

    // to use the geom cache (skip agml build which is faster)
    // set the _geom string to "" and make sure the cache file ("fGeom.root") is present
    // _geom = "";

    // Setup the chain for reading an FZD
    TString _chain;
    
    _chain = Form("fzin %s sdt20211016 fcsDb fwdTrack MakeEvent bigbig evout cmudst tree", _geom.Data() );
    // _chain = Form("fzin %s sdt20211016 fcsDb fwdTrack MakeEvent bigbig evout", _geom.Data() );
    // _chain = Form("fzin %s sdt20211016 fcsDb MakeEvent bigbig ", _geom.Data() );
    

    gSystem->Load( "libStarRoot.so" );
    loadLibs();
    loadLibs();
    gROOT->SetMacroPath(".:/star-sw/StRoot/macros/:./StRoot/macros:./StRoot/macros/graphics:./StRoot/macros/analysis:./StRoot/macros/test:./StRoot/macros/examples:./StRoot/macros/html:./StRoot/macros/qa:./StRoot/macros/calib:./StRoot/macros/mudst:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros/graphics:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros/analysis:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros/test:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros/examples:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros/html:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros/qa:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros/calib:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros/mudst:/afs/rhic.bnl.gov/star/ROOT/36/5.34.38/.sl73_x8664_gcc485/rootdeb/macros:/afs/rhic.bnl.gov/star/ROOT/36/5.34.38/.sl73_x8664_gcc485/rootdeb/tutorials");
    gROOT->LoadMacro("bfc.C");
    stmem.PrintMem("Libs LOADED");
    bfc(-1, _chain, inFile);
    stmem.PrintMem("CHAIN SETUP");

    gSystem->Load( "libStFttSimMaker" );
    
    gSystem->Load( "libStFcsTrackMatchMaker" );
    gSystem->Load( "libMathMore.so" );
    gSystem->Load( "libStarGeneratorUtil" );
    gSystem->Load("StFwdUtils.so");


    gSystem->Load( "libStFttClusterPointMaker" );
    // make an StFttClusterPointMaker
    StFttClusterPointMaker * fttClusterPointMaker = new StFttClusterPointMaker("fttClusterPointMaker");
    fttClusterPointMaker->SetDebug(1);
    fttClusterPointMaker->setUseGeantData( true );

    chain->AddBefore("fwdTrack", fttClusterPointMaker);

    // Configure the Forward Tracker
        StFwdTrackMaker * fwdTrack = (StFwdTrackMaker*) chain->GetMaker( "fwdTrack" );
        if ( fwdTrack ){
            // fwdTrack->SetDebug(0);
            // fwdTrack->SetDebug(0);
            // config file set here for ideal simulation
            if ( _geom == "" ){
                cout << "Using the Geometry cache: fGeom.root" << endl;
                fwdTrack->setGeoCache( "fGeom.root" );
            }

            fwdTrack->setOutputFilename( TString::Format( "%s.output.root", inFile ).Data() );

            // Fitter
            fwdTrack->setFitDebugLvl( 0 );
            fwdTrack->setFitMinIterations( 10 );
            fwdTrack->setFitMaxIterations( 20 );
            
            fwdTrack->setDeltaPval( 1e-1 );
            fwdTrack->setRelChi2Change( 1e-6 );
            
            // fwdTrack->setFttHitSource( 0 /*StFwdHitLoader::GEANT*/ );
            fwdTrack->setFttHitSource( 1 /*StFwdHitLoader::STEVENT*/ );
            // fwdTrack->setFttHitSource( 3 /*StFwdHitLoader::IGNORE*/ );

            fwdTrack->setFstHitSource( 0 /*StFwdHitLoader::GEANT*/ );

            // fwdTrack->setCrit2( "Crit2_DeltaPhi", 0, 2.0 );

            // fwdTrack->setTrackFittingOff();
            
            // fwdTrack->setConfigKeyValue( "TrackFitter:doGlobalTrackFitting", false ); 
            // fwdTrack->setConfigKeyValue( "TrackFitter:findFwdVertices", false ); 
            // fwdTrack->setConfigKeyValue( "TrackFitter:doBeamlineTrackFitting", false ); 
            // fwdTrack->setConfigKeyValue( "TrackFitter:doPrimaryTrackFitting", false ); 
            // fwdTrack->setConfigKeyValue( "TrackFitter:doSecondaryTrackFitting", false );
            fwdTrack->setConfigKeyValue( "TrackFitter:refit", true );


            // float omega = mConfig.get<float>(subsetPath + ".Omega", 0.75);
            // float stableThreshold = mConfig.get<float>(subsetPath + ".StableThreshold", 0.1);
            // float Ti = mConfig.get<float>(subsetPath + ".InitialTemp", 2.1);
            // float Tf = mConfig.get<float>(subsetPath + ".InfTemp", 0.1);

            // fwdTrack->setConfigKeyValue( "TrackFinder.SubsetNN.StableThreshold", 0.01 );
            
            // fwdTrack->setEpdHitSource( 0 /*StFwdHitLoader::GEANT*/ );
            
            bool doFitQA = false;
            if ( doFitQA ){
                StFwdFitQAMaker *fwdFitQA = new StFwdFitQAMaker();
                fwdFitQA->SetDebug();
                TString fitqaoutname(gSystem->BaseName(inFile));
                fitqaoutname.ReplaceAll(".fzd", ".FwdFitQA.root");
                fwdFitQA->setOutputFilename( fitqaoutname );
                chain->AddAfter("fwdTrack", fwdFitQA);
            }
            cout << "fwd tracker setup" << endl;
        }

    // StMuDstMaker * muDstMaker = (StMuDstMaker*)chain->GetMaker( "MuDst" );

    // The PicoDst
    gSystem->Load("libStPicoEvent");
    gSystem->Load("libStPicoDstMaker");
    StPicoDstMaker *picoMk = new StPicoDstMaker(StPicoDstMaker::IoWrite);
    cout << "picoMk = " << picoMk << endl;
    picoMk->setVtxMode(StPicoDstMaker::Vtxless);

    StMemStat stmem;
    stmem.Start();
chain_loop:
	chain->Init();
    stmem.PrintMem("CHAIN INIT COMPLETE");

    //_____________________________________________________________________________
    //
    // MAIN EVENT LOOP
    //_____________________________________________________________________________
    for (int i = 0; i < n; i++) {
        cout << "--------->START EVENT: " << i << endl;
        stmem.PrintMem("BEFORE CHAIN CLEAR");
        chain->Clear();
        stmem.PrintMem("AFTER CHAIN CLEAR COMPLETE");
        stmem.Start();
        if (kStOK != chain->Make())
            break;
    
        cout << "<---------- END EVENT" << endl;
        stmem.Stop();
    } // event loop

    
    stmem.Summary();
    stmem.PrintMem("FINAL MEMORY:");
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