//usr/bin/env root4star -l -b -q  $0; exit $?
// that is a valid shebang to run script as executable

void daq_track(    int n = 1,
                    const char *inFile = "input.daq",
                    std::string configFile = "daq/daq_track.xml",
                    const char *geom = "dev2022") {
    TString _geom = geom;

    bool SiIneff = false;

    TString _chain;

    gSystem->Load( "libStarRoot.so" );

    // Chain with fst, fcs and ftt
    _chain = Form("in, %s, db, StEvent, MuDST, fcs, fst, ftt", _geom.Data());
    
    // needed in this wonky spack environment 
    gROOT->SetMacroPath(".:/star-sw/StRoot/macros:./StRoot/macros:./StRoot/macros/graphics:./StRoot/macros/analysis:./StRoot/macros/test:./StRoot/macros/examples:./StRoot/macros/html:./StRoot/macros/qa:./StRoot/macros/calib:./StRoot/macros/mudst:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros/graphics:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros/analysis:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros/test:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros/examples:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros/html:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros/qa:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros/calib:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros/mudst:/afs/rhic.bnl.gov/star/ROOT/36/5.34.38/.sl73_x8664_gcc485/rootdeb/macros:/afs/rhic.bnl.gov/star/ROOT/36/5.34.38/.sl73_x8664_gcc485/rootdeb/tutorials");

    gROOT->LoadMacro("bfc.C");
    bfc(-1, _chain, inFile);

    // StarMagField::setConstBz(true);


    gSystem->Load("libMathMore.so");
    gSystem->Load("libTMVA.so");
    gSystem->Load("libXMLIO.so"); // needed by FwdTrackerConfig
    gSystem->Load("libStarGeneratorUtil.so"); // needed for StarRandom
    gSystem->Load("libStFstSimMaker.so");
    gSystem->Load("libStFttSimMaker.so");

    gSystem->Load("libgenfit2.so"); // needed for GenFit
    gSystem->Load("libKiTrack.so"); // needed for KiTrack
    gSystem->Load("libStEventUtilities.so");
    gSystem->Load("libStFwdTrackMaker.so");

    gSystem->Load("libStEpdUtil.so");
    gSystem->Load("libStTpcDb.so");

    //  Tracking
    if ( true ) {
        StFwdTrackMaker *gmk = new StFwdTrackMaker();
        gmk->SetConfigFile( configFile );
        
        // write debug histograms and ttree?
        gmk->SetGenerateTree( true );
        gmk->SetGenerateHistograms( true );
        // write out wavefront OBJ files
        gmk->SetVisualize( true );
        
        chain->AddMaker(gmk);
    }

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
