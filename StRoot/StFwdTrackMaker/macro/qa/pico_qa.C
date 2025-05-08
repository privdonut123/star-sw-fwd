
string branches[] = {
    "FwdTracks.mPVal",
    "FwdTracks.mChi2",
    "FwdTracks.mNumberOfFitPoints",
    "FwdTracks.mNumberOfSeedPoints",
    "FwdTracks.mVtxIndex",
    "FwdTracks.mEcalMatchIndex",
    "FwdTracks.mHcalMatchIndex",
    "FwdTracks.mECalX",
    "FwdTracks.mECalY",
    "FwdTracks.mHCalX",
    "FwdTracks.mHCalY",
    "FwdTracks.mDCAXY",
    "FwdTracks.mDCAZ",
    "FwdTracks.mStatus",
    "FwdTracks.mMomentumX",
    "FwdTracks.mMomentumY",
    "FwdTracks.mMomentumZ"
};

TChain *chain = nullptr;
TCanvas *c1 = nullptr;

void a_plot( TString b, TString cut, TString prefix = "" ){

    TString bNamed = b+" >>hb";
    chain->Draw( bNamed, cut, "hist" );
    TH1F *hb = (TH1F*)gDirectory->Get("hb");
    hb->SetTitle( Form("%s %s", b.Data(), cut.Data()) );
    hb->GetXaxis()->SetTitle( b.Data() );
    hb->GetYaxis()->SetTitle( Form("Entries / %f", hb->GetXaxis()->GetBinWidth(1)) );
    hb->SetLineColor( kBlack );
    hb->SetLineWidth(2);
    hb->SetFillColorAlpha( kRed, 0.6 );
    hb->Draw();
    // if the yrange is larger than 1e2, set log scale
    if ( hb->GetMaximum() - hb->GetMinimum() > 1e2 ) c1->SetLogy();
    // c1->SetLogy();
    c1->Print( Form( "plots/%s%s.png", prefix.Data(), b.Data()) );
    c1->SetLogy(0);

    hb->Delete();
    gDirectory->Delete("hb");

}


void pico_qa(){
    chain = new TChain("PicoDst");
    chain->Add("*picoDst.root");
    c1 = new TCanvas("c1", "c1", 1600, 1200);
    

    // no cuts
    for ( auto b : branches ){
        a_plot( b, "" );
    }

    for ( auto b : branches ){
        a_plot( b, "FwdTracks.mPVal > 0.1", "pval_gt_0p1_" );
    }

    for ( auto b : branches ){
        a_plot( b, "FwdTracks.mVtxIndex == 0", "prim_" );
    }

    for ( auto b : branches ){
        a_plot( b, "FwdTracks.mVtxIndex != 0", "glob_" );
    }

    for ( auto b : branches ){
        a_plot( b, "FwdTracks.mVtxIndex == 0 && FwdTracks.mPVal>0.1", "good_prim_" );
    }

    for ( auto b : branches ){
        a_plot( b, "FwdTracks.mVtxIndex != 0 && FwdTracks.mPVal>0.1", "good_glob_" );
    }

    

}