


void eff_from_FwdFitQA(const char *filename = "/Users/brandenburg.89/star/ssw/data/muon_minus/FwdFitQA.root", const char *output = "eff_muon_minus.root") {
    TFile *f = TFile::Open(filename);
    if (!f || f->IsZombie()) {
    std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }

    TFile *fOut = new TFile(output, "RECREATE");
    //   MC Histograms without requiring acceptance cuts
    TH1 *hMcEta = (TH1*)f->Get("McEta");
    TH1 *hMcPt = (TH1*)f->Get("McPt");
    TH1 *hMcPhi = (TH1*)f->Get("McPhi");

    TH1 *McEta3FST = (TH1*)f->Get("McEta3FST");
    TH1 *McPt3FST = (TH1*)f->Get("McPt3FST");
    TH1 *McPhi3FST = (TH1*)f->Get("McPhi3FST");

    TH1 *hRcMatchedMcEta = (TH1*)f->Get("RcMatchedMcEta");
    TH1 *hRcMatchedMcPt = (TH1*)f->Get("RcMatchedMcPt");
    TH1 *hRcMatchedMcPhi = (TH1*)f->Get("RcMatchedMcPhi");

    TH1 *RcMatched3FSTMcEtaPrimary = (TH1*)f->Get("RcMatched3FSTMcEtaPrimary");
    TH1 *RcMatched3FSTMcPtPrimary = (TH1*)f->Get("RcMatched3FSTMcPtPrimary");
    TH1 *RcMatched3FSTMcPhiPrimary = (TH1*)f->Get("RcMatched3FSTMcPhiPrimary");

    TH1 *RcMatched3FSTMcEtaGlobal = (TH1*)f->Get("RcMatched3FSTMcEtaGlobal");
    TH1 *RcMatched3FSTMcPtGlobal = (TH1*)f->Get("RcMatched3FSTMcPtGlobal");
    TH1 *RcMatched3FSTMcPhiGlobal = (TH1*)f->Get("RcMatched3FSTMcPhiGlobal");


    TH1 *hEffEta = (TH1*)RcMatched3FSTMcEtaPrimary->Clone("EffEta");
    hEffEta->Divide(McEta3FST);
    // hEffEta->Scale(1.0 / 3.0);
    hEffEta->SetTitle("Efficiency vs #eta (3FST, Primary); #eta; Efficiency");
    hEffEta->SetLineColor(kBlue);
    hEffEta->SetLineWidth(4);
    hEffEta->Draw("");

    TH1 *hEffPt = (TH1*)RcMatched3FSTMcPtPrimary->Clone("EffPt");
    hEffPt->Divide(McPt3FST);
    hEffPt->SetTitle("Efficiency vs p_{T} (3FST, Primary); p_{T} (GeV/c); Efficiency");
    hEffPt->SetLineColor(kRed);
    hEffPt->SetLineWidth(4);
    hEffPt->Draw("");

    TH1 *hEffPhi = (TH1*)RcMatched3FSTMcPhiPrimary->Clone("EffPhi");
    hEffPhi->Divide(McPhi3FST);
    hEffPhi->SetTitle("Efficiency vs #phi (3FST, Primary); #phi (rad); Efficiency");
    hEffPhi->SetLineColor(kBlack);
    hEffPhi->SetLineWidth(4);
    // hEffPhi->Draw("");


    TH1 *hEffEtaGlobal = (TH1*)RcMatched3FSTMcEtaGlobal->Clone("EffEtaGlobal");
    hEffEtaGlobal->Divide(McEta3FST);
    hEffEtaGlobal->Scale(0.5); // since it includes Global and beamline
    hEffEtaGlobal->SetTitle("Efficiency vs #eta (3FST, Global); #eta; Efficiency");
    hEffEtaGlobal->SetLineColor(kBlue);
    hEffEtaGlobal->SetLineWidth(4);
    hEffEtaGlobal->Draw("");

    TH1 *hEffPtGlobal = (TH1*)RcMatched3FSTMcPtGlobal->Clone("EffPtGlobal");
    hEffPtGlobal->Divide(McPt3FST);
    hEffPtGlobal->Scale(0.5); // since it includes Global and beamline
    hEffPtGlobal->SetTitle("Efficiency vs p_{T} (3FST, Global); p_{T} (GeV/c); Efficiency");
    hEffPtGlobal->SetLineColor(kRed);
    hEffPtGlobal->SetLineWidth(4);
    hEffPtGlobal->Draw("hist");

    TH1 *hEffPhiGlobal = (TH1*)RcMatched3FSTMcPhiGlobal->Clone("EffPhiGlobal");
    hEffPhiGlobal->Divide(McPhi3FST);
    hEffPhiGlobal->Scale(0.5); // since it includes Global and beamline
    hEffPhiGlobal->SetTitle("Efficiency vs #phi (3FST, Global); #phi (rad); Efficiency");
    hEffPhiGlobal->SetLineColor(kBlack);
    hEffPhiGlobal->SetLineWidth(4);
    // hEffPhi->Draw("");

    TCanvas *cEff = new TCanvas("cEff", "Efficiency Plots", 800, 600);
    cEff->cd();

    hEffEtaGlobal->Draw("hist");
    hEffEta->SetLineColor(kRed);
    hEffEta->Draw("same");
    TLegend *leg = new TLegend(0.1, 0.7, 0.3, 0.9);
    leg->AddEntry(hEffEtaGlobal, "Global Efficiency", "l");
    leg->AddEntry(hEffEta, "Primary Efficiency", "l");
    leg->Draw();

    // hEffEta->SetLineColor(kBlack);
    // hEffEta->Draw("same hist");

    cEff->Print("efficiency_eta.pdf");


    hEffPtGlobal->SetLineColor(kBlue);
    hEffPtGlobal->Draw("hist");
    hEffPt->SetLineColor(kRed);
    hEffPt->Draw("hist same");
    leg = new TLegend(0.6, 0.2, 0.9, 0.45);
    leg->AddEntry(hEffPtGlobal, "Global Efficiency", "l");
    leg->AddEntry(hEffPt, "Primary Efficiency", "l");
    leg->Draw();

    cEff->Print("efficiency_Pt.pdf");


    hEffPhiGlobal->SetLineColor(kBlue);
    hEffPhiGlobal->Draw("hist");
    hEffPhi->SetLineColor(kRed);
    hEffPhi->Draw("same");
    leg = new TLegend(0.6, 0.2, 0.9, 0.45);
    leg->AddEntry(hEffPhiGlobal, "Global Efficiency", "l");
    leg->AddEntry(hEffPhi, "Primary Efficiency", "l");
    leg->Draw();

    cEff->Print("efficiency_Phi.pdf");
    
}