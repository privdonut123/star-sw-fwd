//usr/bin/env root4star -l -b -q  $0; exit $?
// that is a valid shebang to run script as executable


void daq_track(    int n = 50,
                    const char *inFile = "/gpfs01/star/pwg_tasks/FwdCalib/DAQ/zfa/st_physics_23072022_raw_7500009.daq",
                    std::string configFile = "daq/daq_track.xml",
                    const char *geom = "dev2022") {
    TString _chain;
    gSystem->Load( "libStarRoot.so" );

    // Simplest chain with fst, fcs, ftt and fwdTracker
    _chain = Form("in, %s, db, StEvent, fcs, fst, ftt, fwdTrack, fstMuRawHit, CMuDst, evout, tree", geom);
    
    // needed in this wonky spack environment 
    gROOT->SetMacroPath(".:/star-sw/StRoot/macros:./StRoot/macros:./StRoot/macros/graphics:./StRoot/macros/analysis:./StRoot/macros/test:./StRoot/macros/examples:./StRoot/macros/html:./StRoot/macros/qa:./StRoot/macros/calib:./StRoot/macros/mudst:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros/graphics:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros/analysis:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros/test:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros/examples:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros/html:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros/qa:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros/calib:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros/mudst:/afs/rhic.bnl.gov/star/ROOT/36/5.34.38/.sl73_x8664_gcc485/rootdeb/macros:/afs/rhic.bnl.gov/star/ROOT/36/5.34.38/.sl73_x8664_gcc485/rootdeb/tutorials");

    gROOT->LoadMacro("bfc.C");
    bfc(-1, _chain, inFile);


    StMaker * fwdMakerGen = chain->GetMaker("fwdTrack");
    if ( fwdMakerGen ){
        //  Extra configuration  for the Forward Tracking
        StFwdTrackMaker *fwdTrack = (StFwdTrackMaker*) chain->GetMaker("fwdTrack");
        if ( fwdTrack ){ //if it is in the chain
            cout << "Setting up Fwd Tracking in chain" << endl;
            // fwdTrack->SetConfigFile( configFile );
            fwdTrack->setConfigForData();
            fwdTrack->setZeroB();
            // fwdTrack->setSeedFindingWithFst();
            // write debug histograms and ttree?
            fwdTrack->SetGenerateTree( true );
            fwdTrack->SetGenerateHistograms( false );
            // write out wavefront OBJ files
            fwdTrack->SetVisualize( false );
            fwdTrack->SetDebug(2);
            fwdTrack->setGeoCache( "fGeom.root" );
            
            
            // fwdTrack->setDebug();
        }
    }

    // Initialize the chain
    chain->Init();

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
