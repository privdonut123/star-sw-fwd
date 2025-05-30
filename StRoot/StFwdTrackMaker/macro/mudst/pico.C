//usr/bin/env root4star -l -b -q $0; exit $?
/////////////////////////////////////////////////////////////////////////////
// $Id: genDst.C,v 1.9 2020/10/10 07:16:56 genevb Exp $
// Author: G. Van Buren (BNL)
//
// Description:
// Process a MuDst for...
// ...creating a PicoDst
// ...re-running vertex-finding to re-create MuDsts
//
// Options are space-separated or comma-separated,
// and case-insensitive. They can be attributed
// whose values are provided after a ':'.
//
// Example options for creating PicoDsts:
// picoDst
// btofMatch
// btofStartless
// mtdMatch
// y2017a
//
// Example lists of options:
// "picoDst"
// "DbV20200125,picoDst,mtdMatch,y2014a"
//
// Example options for vertex-finding:
// beamline, beamline1D, beamline3D (otherwise no beamline)
// useBTOFmatchOnly
// VFstore:100
//
// Example lists of options:
// "VFPPVnoCTB,beamline1D,VFstore:100"
// "VFPPVnoCTB,beamline3D"
//
/////////////////////////////////////////////////////////////////////////////

void afterburner();

void genDst(unsigned int Last,
            const char* options,
            char* infile,
            char* outfile=0);

void genDst(unsigned int First,
            unsigned int Last,
            const char* options,
            char* infile,
            char* outfile=0);

void loadLibs()
{
  gROOT->Macro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C");
  gSystem->Load("StDbBroker");
  gSystem->Load("St_db_Maker");
  gSystem->Load("StEEmcUtil");
}

void loadLibsVF()
{
  gSystem->Load("libMinuit");
  gSystem->Load("Sti");
  gSystem->Load("StBTofUtil");
  gSystem->Load("StGenericVertexMaker");
}

void loadLibsPico()
{

  // EMCs and FMS need DB+converters
  gSystem->Load("StEmcRawMaker");
  gSystem->Load("StEmcADCtoEMaker");
  gSystem->Load("StPreEclMaker");
  gSystem->Load("StEpcMaker");
  gSystem->Load("StEEmcDbMaker");
  gSystem->Load("StFmsUtil");
  gSystem->Load("StFmsDbMaker");
  gSystem->Load("StTriggerUtilities");

  // The PicoDst
  gSystem->Load("libStPicoEvent");
  gSystem->Load("libStPicoDstMaker");
}

void loadLibsAgML()
{
  // load support libraries.  util will fail to load for agml 1.0, but you can ignore the error
  gSystem->Load("libStarAgmlUtil");
  gSystem->Load("libStarAgmlLib");

  // load geometry modules and master steering codes...
  gSystem->Load("libGeometry");
  gSystem->Load("libStarGeometry");
}

void loadLibsMtd()
{
  gSystem->Load("StDetectorDbMaker");
  gSystem->Load("StarMagField");
  gSystem->Load("StMagF");
  gSystem->Load("StMtdUtil");
  gSystem->Load("StMtdMatchMaker");
  gSystem->Load("StMtdCalibMaker");
}

void loadLibsBTof()
{
  gSystem->Load("StBTofUtil");
  gSystem->Load("StVpdCalibMaker");
  gSystem->Load("StBTofCalibMaker");
  gSystem->Load("StBTofMatchMaker");
}

void loadLibsETof()
{
  gSystem->Load("StETofUtil");
  gSystem->Load("StETofCalibMaker");
  gSystem->Load("StETofHitMaker");
  gSystem->Load("StETofMatchMaker");
}

void loadLibsFwd()
{
  gSystem->Load("libStDetectorDbMaker.so");
  gSystem->Load("StEvent");
  gSystem->Load("StEventMaker");

  // Fwd
  gSystem->Load("StFwdUtils");
  gSystem->Load("libStarGeneratorUtil.so");
  gSystem->Load("libXMLIO");
  gSystem->Load("libgenfit2");
  gSystem->Load("libKiTrack");
  // gSystem->Load("StarGeneratorUtil")  ;
  gSystem->Load("libMathMore");
  gSystem->Load("StEventUtilities");
  gSystem->Load("StEpdUtil");
  gSystem->Load("StFwdTrackMaker");
  gSystem->Load("StFwdUtils");
  //Ftt
  gSystem->Load("StFttDbMaker");
  gSystem->Load("StFttRawHitMaker");
  gSystem->Load("StFttHitCalibMaker");
  gSystem->Load("StFttClusterMaker");
  gSystem->Load("StFttPointMaker");
  // Fcs
  gSystem->Load("libMinuit.so");
  gSystem->Load("StFcsDbMaker");
  gSystem->Load("StFcsRawHitMaker");
  gSystem->Load("StFcsWaveformFitMaker");
  gSystem->Load("StFcsClusterMaker");
  gSystem->Load("StFcsPointMaker");
  gSystem->Load("StFcsTrackMatchMaker");

  // Fst
  gSystem->Load("StFstUtil");
  gSystem->Load("StFstDbMaker");
  gSystem->Load("StFstRawHitMaker");
  gSystem->Load("StFstClusterMaker");
  gSystem->Load("StFstHitMaker");
}

void procGeoTag(TObjArray* optionTokens)
{
  if (TClass::GetClass("AgBlock")) return; // arbitrarily chosen AgML class
  loadLibsAgML();

  const char* tag = 0;
  for (int tk=0; tk < optionTokens->GetEntries(); tk++) {
    TString& tok = ((TObjString*) (optionTokens->At(tk)))->String();
    if (tok.BeginsWith("y20")){
      tag = tok.Data();
      optionTokens->RemoveAt(tk);
      optionTokens->Compress();
      break;
    }
  }

  // Let agml know we want the ROOT geometry
  AgModule::SetStacker( new StarTGeoStacker );

  // now pass the geometry "tag" and build.  If the class StarGeometry exists, we have
  // AgML 2.0 and can run using the new steering.  Otherwise, old steering codes...
  if (tag) {
    if ( TClass::GetClass("StarGeometry") ) {      StarGeometry::Construct( tag ); }
    else                                    {      ( new Geometry() )->ConstructGeometry(tag); }
  } else {
    gMessMgr->Warning() << "No geometry tag passed! (e.g. y2017a)" << endm;
  }
}

bool findAndRemoveOption(const char* optionName, TObjArray* optionTokens)
{
  TString optName = optionName;
  optName.ToLower();
  TObject* obj = optionTokens->FindObject(optName.Data());
  if (obj) {
    optionTokens->Remove(obj);
    optionTokens->Compress();
    return true;
  }
  return false;
}

void genDst(unsigned int First,
            unsigned int Last,
            const char* options,
            char* infile,
            char* outfile)
{
  loadLibs();

  StChain fullChain("genDst");
  fullChain.SetDebug(1);

  StMuDstMaker muDstMaker(0, 0, "", infile, "st:MuDst.root", 1e9); // set up maker in read mode
  //                      0, 0                        this means read mode
  //                           dir                    read all files in this directory
  //                               file               bla.lis read all file in this list, if (file!="") dir is ignored
  //                                    filter        apply filter to filenames, multiple filters are separated by ':'
  //                                          10      maximum number of file to read


  TChain& muDstChain = *muDstMaker.chain();
  unsigned int nEntries      = muDstChain.GetEntries();
  unsigned int LastToRead = Last > 0 ? min(Last, nEntries) : nEntries;
  gMessMgr->Info() << nEntries << " events in chain, " << LastToRead-First+1 << " will be read." << endm;

  // St_db_Maker* db = new St_db_Maker("db", "StarDb", "MySQL:StarDb", "$STAR/StarDb");

  // Initialize some values and pointers
  StMaker* processMaker = 0;
  TFile* outFile = 0;
  TTree* muDstTreeOut = 0;

  // Basic decisions based on options
  TString CasedOptions = options;
  TString Options = options;
  Options.ToLower();
  TString optDelim = " ,";
  TObjArray* optionTokens = Options.Tokenize(optDelim);
  optionTokens->SetOwner(kTRUE);

  // Determine database flavors
  TString flavors = "ofl"; // default flavor for offline

  // simulation flavors
  if (findAndRemoveOption("Simu",optionTokens) && ! findAndRemoveOption("NoSimuDb",optionTokens))
    flavors.Prepend("sim+");

  // filestream flavors
  TObject* firstFile = muDstChain.GetListOfFiles()->At(0);
  if (firstFile) {
    TString firstFileName = firstFile->GetTitle();
    firstFileName = firstFileName(firstFileName.Last('/')+1,firstFileName.Length());
    if (firstFileName.BeginsWith("st_")) {
      TString fileStream = firstFileName(3,firstFileName.Index('_',3)-3);
      if (fileStream.Length()>0) flavors.Prepend(fileStream += '+');
    }
  }

  // gMessMgr->Info() << "Using DB flavors: " << flavors << endm;
  // db->SetFlavor(flavors.Data());

  if (findAndRemoveOption("picodst",optionTokens)) {
    // _________________________________________________________________
    // Processing with generation of PicoDsts

    loadLibsPico();

    // Specify active branches but first disable all branches
    muDstMaker.SetStatus("*", 0);
    muDstMaker.SetStatus("MuEvent", 1);
    muDstMaker.SetStatus("PrimaryVertices", 1);
    muDstMaker.SetStatus("PrimaryTracks", 1);
    muDstMaker.SetStatus("GlobalTracks", 1);
    muDstMaker.SetStatus("CovGlobTrack", 1);
    muDstMaker.SetStatus("BTof*", 1);
    muDstMaker.SetStatus("Emc*", 1);
    muDstMaker.SetStatus("MTD*", 1);
    muDstMaker.SetStatus("ETof*", 1);
    muDstMaker.SetStatus("Epd*", 1);
    muDstMaker.SetStatus("Fms*", 1);
    muDstMaker.SetStatus("MCAll", 1);
    muDstMaker.SetStatus("Fwd*", 1);
    muDstMaker.SetStatus("Fcs*", 1);
    muDstMaker.SetStatus("Ftt*", 1);
    muDstMaker.SetStatus("Fst*", 1);

    // EMCs
    // StEEmcDbMaker* eemcDb = new StEEmcDbMaker;
    // StEmcADCtoEMaker* adc2e = new StEmcADCtoEMaker();
    // adc2e->saveAllStEvent(true);
    // StPreEclMaker* pre_ecl = new StPreEclMaker();
    // StEpcMaker* epc = new StEpcMaker();

    // FMS
    // StFmsDbMaker* fmsDb = new StFmsDbMaker("fmsDb");

    // Trigger simulator
    // StTriggerSimuMaker* trigSimu = new StTriggerSimuMaker;
    // trigSimu->setMC(false);
    // trigSimu->useBemc();
    // trigSimu->useEemc();
    // trigSimu->useOfflineDB();
    // trigSimu->bemc->setConfig(StBemcTriggerSimu::kOffline);

    loadLibsFwd();

    // If you need the geometry you can load it with this:
    // gROOT->LoadMacro("/star-sw/StarVMC/Geometry/macros/loadStarGeometry.C");
    // loadStarGeometry( "y2023" );

    // This makes the Fcs and Ftt hits show up in StEvent,
    // and produces an StEvent for StFwdTracks to go into
    StMuDst2StEventMaker * mu2ev = new StMuDst2StEventMaker();

    // Re-run the FCS chain
    StFcsDbMaker * fcsDbMk = new StFcsDbMaker();
    // StFcsWaveformFitMaker * fcsWFF = new StFcsWaveformFitMaker();
    // StFcsClusterMaker * fcsClu = new StFcsClusterMaker();
    // StFcsPointMaker * fcsPoint = new StFcsPointMaker();

    // FTT chain
    StFttDbMaker * fttDbMk = new StFttDbMaker();
    StFttHitCalibMaker * ftthcm = new StFttHitCalibMaker();
    StFttClusterMaker * fttclu = new StFttClusterMaker();
    fttclu->SetTimeCut(1, -40, 40);
    StFttPointMaker * fttpoint = new StFttPointMaker();


    StFwdTrackMaker * fwdTrack = new StFwdTrackMaker();
    fwdTrack->setConfigForData( );
    fwdTrack->setGeoCache( "fGeom.root" );
    fwdTrack->setSeedFindingWithFst();
    fwdTrack->setIncludePrimaryVertexInFit(true);
    fwdTrack->setMaxFailedHitsInFit(2);
    fwdTrack->SetDebug(2);

    StFcsTrackMatchMaker *match = new StFcsTrackMatchMaker();
    match->setMaxDistance(6,10);
    match->setFileName("fcstrk.root");
    match->SetDebug();

    StFwdAnalysisMaker *fwdAna = new StFwdAnalysisMaker();
    fwdAna->SetDebug();
    fwdAna->setLocalOutputFile( "StFwdAnalysisMaker.root");

    StFwdFitQAMaker *fwdFitQA = new StFwdFitQAMaker();
    fwdFitQA->SetDebug();

    StFwdQAMaker *fwdQAMk = new StFwdQAMaker();
    fwdQAMk->SetDebug(2);
    // chain->AddAfter("fwdTrack", fwdQAMk);

    if (findAndRemoveOption("btofmatch",optionTokens)) {

      procGeoTag(optionTokens);
      loadLibsBTof();

      // instantiate both VPD and BTOF CalibMakers and MatchMaker and point them to the MuDST
      StBTofMatchMaker* btofMatch = new StBTofMatchMaker();
      btofMatch->setMuDstIn();
      StVpdCalibMaker *vpdCalib = new StVpdCalibMaker();
      vpdCalib->setMuDstIn();
      StBTofCalibMaker *btofCalib = new StBTofCalibMaker();
      btofCalib->setMuDstIn();

      if (findAndRemoveOption("btofstartless",optionTokens)) {
        //Disable the VPD as start detector, BTOF calib maker will switch to the "start-less" algorithm.
        vpdCalib->setUseVpdStart(kFALSE);
      }

    }

    if (findAndRemoveOption("etofmatch",optionTokens)) {

      procGeoTag(optionTokens);
      loadLibsETof();

      // instantiate eTOF CalibMakers, HitMaker and MatchMaker
      StETofCalibMaker* etofCalib = new StETofCalibMaker();
      StETofHitMaker*   etofHit   = new StETofHitMaker();
      StETofMatchMaker* etofMatch = new StETofMatchMaker();

    }

    if (findAndRemoveOption("mtdmatch",optionTokens)) {

      procGeoTag(optionTokens);
      loadLibsMtd();

      StMagFMaker* magfMk = new StMagFMaker;
      StMtdMatchMaker* mtdMatchMaker = new StMtdMatchMaker();
      StMtdCalibMaker* mtdCalibMaker = new StMtdCalibMaker("mtdcalib");

    }

    processMaker = (StMaker*) (new StPicoDstMaker(StPicoDstMaker::IoWrite, infile, "picoDst"));
    StPicoDstMaker *picoMk = (StPicoDstMaker*) (processMaker);
    picoMk->setVtxMode(StPicoDstMaker::Vpd);
    gMessMgr->Info() << "PicoSetup Complete" << endm;
  } else if (Options.Contains("fwd")) {
    // _________________________________________________________________
    // Processing with new vertex-finding

    // loadLibsVF();

    // Specify inactive branches but first enable all branches
    muDstMaker.SetStatus("*",1);
    muDstMaker.SetStatus("PrimaryTracks",0);
    muDstMaker.SetStatus("PrimaryVertices",0);

    // Create new branch
    // TClonesArray* verticesRefitted = new TClonesArray("StMuPrimaryVertex", 1000);

    // Specify output
    if (outfile) {
      outFile = new TFile(outfile, "RECREATE");
    } else {
      // Use the same filename for output as was given by input
      TString fileStr = infile;
      Ssiz_t dir = fileStr.Last('/');
      if (dir<0) {
        gMessMgr->Error() << "No specification for output when input is in local directory!" << endm;
        return;
      }
      fileStr.Remove(0,dir+1);
      outFile = new TFile(fileStr.Data(), "RECREATE");
    }
    muDstTreeOut = muDstChain.CloneTree(0);
    muDstTreeOut->Branch("PrimaryVertices", &verticesRefitted, 65536, 99);

    // processMaker = (StMaker*) (new StGenericVertexMaker());
    // processMaker->ToWhiteConst("vtxArray",verticesRefitted);
    // processMaker->SetAttr("useMuDst",1);

  } else {

    gMessMgr->Info() << "No processing specified - just reading a MuDst?" << endm;
    // User code may be inserted here

  }

  // Set additional options (except DbV) as maker attributes
  if (processMaker) {
    for (int tk=0; tk < optionTokens->GetEntries(); tk++) {
      TString& Tag = ((TObjString*) (optionTokens->At(tk)))->String();

      // copy of DbV code from StBFChain.cxx
      if (Tag.BeginsWith("dbv")) {
        int FDate=0,FTime=0;
        if (Tag.Length() == 11)  (void) sscanf(Tag.Data(),"dbv%8d",&FDate);
        if (Tag.Length() == 18)  (void) sscanf(Tag.Data(),"dbv%8d.%6d",&FDate,&FTime);
        if (FDate) {
          // db->SetMaxEntryTime(FDate,FTime);
          // gMessMgr->Info() << "\tSet DataBase max entry time " << FDate << "/" << FTime
          //                  << " for St_db_Maker(\"" << db->GetName() <<"\")" << endm;
        }
        continue;
      }

      // assign attributes
      StMaker* attrMaker = processMaker;
      Ssiz_t delim = Tag.First(':');
      // look for "::" to set attributes for a different maker
      if (delim > 0 && Tag[delim+1] == ':') {
        TString altMakerName = Tag(0,delim);
        // GetMaker...() functions are case sensitive, so find original case
        Ssiz_t casedMakerNameIdx = CasedOptions.Index(altMakerName,0,TString::ECaseCompare::kIgnoreCase);
        if (casedMakerNameIdx >= 0) altMakerName = CasedOptions(casedMakerNameIdx,delim);
        StMaker* altMaker = fullChain.GetMaker(altMakerName.Data());
        if (!altMaker) altMaker = fullChain.GetMakerInheritsFrom(altMakerName.Data());
        if (!altMaker) {
          gMessMgr->Warning() << "No maker found with name or class " << altMakerName.Data() << endm;
          continue;
        }
        attrMaker = altMaker;
        Tag.Remove(0,delim+2);
        delim = Tag.First(':');
      }
      if (delim < 0) {
        attrMaker->SetAttr(Tag.Data(),1);
      } else {
        TString key(Tag(0,delim));
        TString& val = Tag.Remove(0,delim+1);
        if (val.IsDigit()) { attrMaker->SetAttr(key.Data(),val.Atoi()); }
        else if (val.IsFloat()) { attrMaker->SetAttr(key.Data(),val.Atof()); }
        else { attrMaker->SetAttr(key.Data(),val.Data()); }
      }
    }
    processMaker->PrintAttr();
  }

  {
    TDatime t;
    gMessMgr->QAInfo() << Form("Run is started at Date/Time %i/%i",t.GetDate(),t.GetTime()) << endm;
  }
  gMessMgr->QAInfo() << Form("Run on %s in %s",gSystem->HostName(),gSystem->WorkingDirectory()) << endm;
  gMessMgr->QAInfo() << Form("with %s", fullChain.GetCVS()) << endm;

  gMessMgr->Info() << "Chain setup Complete" << endm;

  // Main loop over events
  int iInit = fullChain.Init();
  if (iInit >= kStEOF) {fullChain.FatalErr(iInit,"on init"); return;}
  if (Last == 0) return;
  int eventCount = 0;
  // Skip, if any
  if (First > 1) fullChain.Skip(First - 1);
  for (unsigned int iEvent = First; iEvent <= LastToRead; iEvent++)
  {
      // make an StEvent so that Fwd can save tracks, etc.
      // StEvent * stEv = new StEvent();
      // fullChain.AddData( stEv );
      fullChain.SetDebug(2);
      fwdTrack->SetDebug(2);

     int iMake = fullChain.Make();
     if (iMake) {fullChain.FatalErr(iMake,"on make"); return;}

     if (muDstTreeOut) muDstTreeOut->Fill();

     int iClear = fullChain.Clear();
     if (iClear) {fullChain.FatalErr(iClear,"on clear"); return;}
     eventCount++;
  }
  fullChain.Finish();

  //
  // ATTENTION - please DO NOT change the format of the next 2 lines,
  //   they are used by our DataManagement parsers to detect a generation
  //   was succesful and thereafter Catalog the produced files.
  //   Thank you.
  //
  gMessMgr->QAInfo() << "NumberOfEvents= " << eventCount << endm;
  gMessMgr->QAInfo() << "Run completed " << endm;

  if (outFile) {
    outFile->Write();
    outFile->Close();
    delete outFile;
  }

  // delete db;
  delete processMaker;
  delete optionTokens;
}

//__________________________________________________________
void genDst(unsigned int Last,
            const char* options,
            char* infile,
            char* outfile)
{
  cout << TString::Format("genDst( %u, '%s', '%s', '%s' )", Last, options, infile, outfile ) << endl;
  genDst(1,Last,options,infile,outfile);
}

void pico(){
  genDst(5000, "y2023a picodst PicoVtxMode:PicoVtxDefault", "/work/st_fwd_22355048_raw_1000012.MuDst.root", "PROD.root");
}

void pico( TString f, int n = 500){
  genDst(n, "y2023a picodst PicoVtxMode:PicoVtxDefault", f.Data(), "PROD.root" /*Not used?*/);
}

/////////////////////////////////////////////////////////////////////////////
//
// $Log: genDst.C,v $
// Revision 1.10  2022/05/12 07:16:56  weidenkaff
// Added eTOF support
//
// $Log: genDst.C,v $
// Revision 1.9  2020/10/10 07:16:56  genevb
// Specify makers to set attributes
//
// Revision 1.8  2020/08/19 15:27:33  genevb
// Add DB flavors
//
// Revision 1.7  2020/01/25 05:10:00  genevb
// Include DbV, more like BFC
//
// Revision 1.6  2019/09/18 17:54:48  genevb
// Acivate additional branches by default
//
// Revision 1.5  2019/03/21 18:53:34  jeromel
// Added ATTENTION message
//
// Revision 1.4  2019/01/15 17:24:29  genevb
// Added FMS
//
// Revision 1.3  2018/03/16 18:41:14  genevb
// Add BTof-matching
//
// Revision 1.2  2017/12/15 18:36:53  genevb
// Remove explicit function of StPicoDstMaker...params should be passed by attribute
//
// Revision 1.1  2017/12/05 16:47:58  genevb
// Introduce genDst.C for creating new Dsts from MuDsts
//
//
/////////////////////////////////////////////////////////////////////////////
