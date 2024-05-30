//usr/bin/env root4star -l -b -q  $0; exit $?
// that is a valid shebang to run script as executable

void prod( int n = 10000, const char *inFile = "input.daq" ) {
    
    TString _chain;
    gSystem->Load( "libStarRoot.so" );

    // Simplest chain with fst, fcs, ftt and fwdTracker
    _chain = "in, dev2022, db, tpcDB, MakeEvent, StEvent, MuDST, fst, ftt,  CMuDst, evout, tree, fwdTrack";
    
    // needed in this wonky spack environment 
    gROOT->SetMacroPath(".:/star-sw/StRoot/macros:./StRoot/macros:./StRoot/macros/graphics:./StRoot/macros/analysis:./StRoot/macros/test:./StRoot/macros/examples:./StRoot/macros/html:./StRoot/macros/qa:./StRoot/macros/calib:./StRoot/macros/mudst:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros/graphics:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros/analysis:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros/test:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros/examples:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros/html:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros/qa:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros/calib:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros/mudst:/afs/rhic.bnl.gov/star/ROOT/36/5.34.38/.sl73_x8664_gcc485/rootdeb/macros:/afs/rhic.bnl.gov/star/ROOT/36/5.34.38/.sl73_x8664_gcc485/rootdeb/tutorials");

    gROOT->LoadMacro("bfc.C");
    bfc(-1, _chain, inFile, "OUTPUT.root");
    //bfc(-1, _chain, inFile);

    //StFttHitCalibMaker *fttHitCalib = (StFttHitCalibMaker*) chain->GetMaker("fttHitCalib");
    //fttHitCalib->SetMode( StFttHitCalibMaker::CalibMode::Calibration );
    StFwdTrackMaker *fwdTrack = (StFwdTrackMaker*) chain->GetMaker("fwdTrack");
    if ( fwdTrack ){ //if it is in the chain
        //fwdTrack->SetConfigFile( configFile );
        // write debug histograms and ttree?
        //fwdTrack->SetGenerateTree( true );
        fwdTrack->SetGenerateHistograms( true );
        // write out wavefront OBJ files
        fwdTrack->SetVisualize( false );
        chain->AddMaker(fwdTrack);
    }


    // Initialize the chain
    chain->Init();

    //_____________________________________________________________________________
    //
    // MAIN EVENT LOOP
    //_____________________________________________________________________________
    for (int i = 0; i < n; i++) {
        chain->Clear();
        if (kStOK != chain->Make())
            break;
    }
}
