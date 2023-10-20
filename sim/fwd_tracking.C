//usr/bin/env root4star -l -b -q  $0; exit $?
// that is a valid shebang to run script as executable, but with only one arg


// Run very fast fwd tracking with only track seed finding
// generate some input data using genfzd 

TFile *output = 0;

void fwd_tracking(      int n = 500,
                const char *inFile =  "simu/seed.fzd",
                std::string configFile = "simu/seed.xml",
                const char *geom = "dev2022") {
    TString _geom = geom;
    bool SiIneff = false;
    bool useConstBz = false;

    bool useFCS = false;

    // Setup the chain for reading an FZD
    TString _chain;
    
    
    _chain = Form("fzin %s sdt20211016 StEvent ReverseField geant gstar agml usexgeom bigbig fttDb", _geom.Data());

    gSystem->Load( "libStarRoot.so" );
    gROOT->SetMacroPath(".:/star-sw/StRoot/macros/:./StRoot/macros:./StRoot/macros/graphics:./StRoot/macros/analysis:./StRoot/macros/test:./StRoot/macros/examples:./StRoot/macros/html:./StRoot/macros/qa:./StRoot/macros/calib:./StRoot/macros/mudst:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros/graphics:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros/analysis:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros/test:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros/examples:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros/html:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros/qa:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros/calib:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros/mudst:/afs/rhic.bnl.gov/star/ROOT/36/5.34.38/.sl73_x8664_gcc485/rootdeb/macros:/afs/rhic.bnl.gov/star/ROOT/36/5.34.38/.sl73_x8664_gcc485/rootdeb/tutorials");
    gROOT->LoadMacro("bfc.C");
    bfc(-1, _chain, inFile);

    gSystem->Load( "libStFttSimMaker" );

    if ( useConstBz )
        StarMagField::setConstBz(true);
    //fcssim->setDebug(1);
    //fcssim->setLeakyHcal(0);

    // Configure FTT FastSim
        StFttSlowSimMaker *fttSlowSim = new StFttSlowSimMaker();
        // (StFttSlowSimMaker*) chain->GetMaker( "fttSim" );
        cout << "Adding StFttSlowSimMaker to chain" << endl;
        chain->AddMaker(fttSlowSim);

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



// GEXE $STAR_LIB/libStarAgmlUtil.so
// GEXE $STAR_LIB/xgeometry.so