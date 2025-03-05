//usr/bin/env root4star -l -b -q  $0; exit $?
// that is a valid shebang to run script as executable


void daq_track(    int n = 50,
                    const char *inFile = "indaq/st_physics_23072022_raw_4000001.daq",
                    const char *geom = "y2023") {
    cout << "Running daq_track" << endl;
    TString _chain;
    gSystem->Load( "libStarRoot.so" );

    // Simplest chain with fst, fcs, ftt and fwdTracker
    _chain = Form("in, %s, db, StEvent, trgd, btof, fcs, fst, ftt, fwdTrack, fstMuRawHit, EventQA, CMuDst, evout, tree", geom);
    // _chain = Form("in, %s, StEvent, fcs, fst, ftt, fwdTrack, evout, tree", geom);

    // needed in this wonky spack environment / docker container
    gROOT->SetMacroPath(".:/star-sw/StRoot/macros:./StRoot/macros:./StRoot/macros/graphics:./StRoot/macros/analysis:./StRoot/macros/test:./StRoot/macros/examples:./StRoot/macros/html:./StRoot/macros/qa:./StRoot/macros/calib:./StRoot/macros/mudst:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros/graphics:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros/analysis:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros/test:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros/examples:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros/html:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros/qa:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros/calib:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros/mudst:/afs/rhic.bnl.gov/star/ROOT/36/5.34.38/.sl73_x8664_gcc485/rootdeb/macros:/afs/rhic.bnl.gov/star/ROOT/36/5.34.38/.sl73_x8664_gcc485/rootdeb/tutorials");

    gROOT->LoadMacro("bfc.C");

    gROOT->LoadMacro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C");
	loadSharedLibraries();

    printf("Running bfc\n");
    bfc(-1, _chain, inFile);
    printf("Done running bfc -1\n");


    StFttClusterMaker *fttClu = (StFttClusterMaker*) chain->GetMaker("stgcCluster");
    if (fttClu){
        // fttClu->SetDebug(2);
        fttClu->SetTimeCut( 1, -40, 40);
    }

    StMaker * fwdMakerGen = chain->GetMaker("fwdTrack");
    if ( fwdMakerGen ){
        //  Extra configuration  for the Forward Tracking
        StFwdTrackMaker *fwdTrack = (StFwdTrackMaker*) chain->GetMaker("fwdTrack");
        if ( fwdTrack ){ //if it is in the chain
            cout << "Setting up Fwd Tracking in chain" << endl;
            // fwdTrack->SetConfigFile( configFile );
            // fwdTrack->setConfigForData();
            // fwdTrack->setZeroB();
            // fwdTrack->setSeedFindingWithFst();
            // write out wavefront OBJ files
            // fwdTrack->SetVisualize( false );
            fwdTrack->SetDebug(1);
            // fwdTrack->setTrackRefit(true);
            fwdTrack->setGeoCache( "fGeom.root" );
            printf("Loading config for FwdTrackMaker\n");
            fwdTrack->LoadConfiguration();
            // fwdTrack->setDebug();

            
        }
    }

    if (true){
        // The PicoDst
        gSystem->Load("libStPicoEvent");
        gSystem->Load("libStPicoDstMaker");
        StPicoDstMaker *picoMk = new StPicoDstMaker(StPicoDstMaker::IoWrite);
        cout << "picoMk = " << picoMk << endl;
        picoMk->setVtxMode(StPicoDstMaker::Default);

            // Generate FWD QA
        StFwdQAMaker *fwdQAMk = new StFwdQAMaker();
        fwdQAMk->SetDebug(2);
        chain->AddAfter("fwdTrack", fwdQAMk);
    }

    // chain->Print();
    // Initialize the chain
    printf("-------->Initializing chain\n");
    chain->Init();
    printf("-------->Done initializing chain\n");

    //

    //_____________________________________________________________________________
    //
    // MAIN EVENT LOOP
    //_____________________________________________________________________________
    for (int i = 0; i < n; i++) {
        chain->Clear();
        if ( fwdMakerGen )
            fwdMakerGen->SetDebug(1);
        if (kStOK != chain->Make())
            break;
    }
}
