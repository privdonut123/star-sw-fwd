

    

void dndeta_plot(){
    gStyle->SetOptStat(0);
    TFile *fPythia = new TFile("dndch1_pythia.root");
    TFile *fData = new TFile("dndch1_data.root");

    TH1 * hPythia = (TH1*)fPythia->Get("hFwdMultRecoMC")->Clone("hPythia");
    hPythia->SetTitle("N_Ch in 2.5 < #eta < 4.0; Multiplicity; arb. units");
    hPythia->SetLineColor(kRed);
    hPythia->SetFillColorAlpha(kRed, 0.3);
    hPythia->SetLineWidth(4);
    TH1 * hData = (TH1*)fData->Get("hFwdMultReco")->Clone("hData");
    hData->SetTitle("Data N_Ch; Multiplicity; dN/d#Delta#eta");
    hData->SetLineColor(kBlue);
    hData->SetLineWidth(4);
    hData->SetMarkerStyle(20);
    hData->SetMarkerColor(kBlue);
    hData->SetMarkerSize(1.2);

    TCanvas *c = new TCanvas("c", "c", 800*3, 600*3);
    
    hPythia->SetLineColor(kRed);
    hData->SetLineColor(kBlue);
    hPythia->GetXaxis()->SetRangeUser(1, 21);
    hPythia->DrawNormalized();
    
    hData->DrawNormalized("same pe");
    c->BuildLegend();
    c->SaveAs("nch_plot_linear.png");
    // gPad->SetLogy();
    // c->SaveAs("nch_plot_logy.png");
}