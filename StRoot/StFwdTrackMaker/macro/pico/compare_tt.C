

TH1 * getResolution(TH2 * h2, TString nn){
    printf("Fitting %s\n", h2->GetName());
    ((TH2*)h2->Clone(nn))->FitSlicesY();
    TH1 * h1 = (TH1*)gDirectory->Get( Form("%s_2", nn.Data()) );
    assert(h1);
    return h1;
}



void compare_tt( TString label = "", TString note = "" ){
    // open files for each track type
    TFile *fGlobal = TFile::Open( TString::Format("efficiency_trackType0_%s.root", label.Data()));
    TFile *fBLC = TFile::Open( TString::Format("efficiency_trackType1_%s.root", label.Data()));
    TFile *fPrimary = TFile::Open( TString::Format("efficiency_trackType2_%s.root", label.Data()));
    

    // now draw each one with a different color
    TCanvas *c = new TCanvas("c", "Efficiency Comparison", 3840, 2160/2.0);
    c->Divide(3, 1);
    gStyle->SetOptStat(0);
    {
        c->cd(1);
        TH1 *hGlobalEta = (TH1*)fGlobal->Get("Eta_MatchedPrimary_McPrimary;1");
        TH1 *hBLCEta = (TH1*)fBLC->Get("Eta_MatchedPrimary_McPrimary;1");
        TH1 *hPrimaryEta = (TH1*)fPrimary->Get("Eta_MatchedPrimary_McPrimary;1");
        
        hGlobalEta->SetLineColor(kBlack);
        hBLCEta->SetLineColor(kRed);
        hPrimaryEta->SetLineColor(kBlue);

        hGlobalEta->SetTitle(label + "; #eta; Efficiency (Matched McPrimary / McPrimary)");
        hGlobalEta->GetXaxis()->SetRangeUser(2.0, 4.5);
        hGlobalEta->Draw();
        hBLCEta->Draw("same");
        hPrimaryEta->Draw("same");

        TLine *line = new TLine(1.0, 0.0, 5.0, 0.0);
        line->SetLineColor(kBlack);
        line->SetLineStyle(2);
        line->DrawLine(1.0, 1.0, 5.0, 1.0);

        TLegend *leg1 = new TLegend(0.6, 0.2, 0.9, 0.5);
        leg1->AddEntry(hGlobalEta, "Global", "l");
        leg1->AddEntry(hBLCEta, "Beamline", "l");
        leg1->AddEntry(hPrimaryEta, "Primary", "l");
        leg1->Draw();
    }

    {
        auto pad = c->cd(2);
        TH1 *hGlobalPt = (TH1*)fGlobal->Get("Pt_MatchedPrimary_McPrimary;1");
        TH1 *hBLCPt = (TH1*)fBLC->Get("Pt_MatchedPrimary_McPrimary;1");
        TH1 *hPrimaryPt = (TH1*)fPrimary->Get("Pt_MatchedPrimary_McPrimary;1");
        
        hGlobalPt->SetLineColor(kBlack);
        hBLCPt->SetLineColor(kRed);
        hPrimaryPt->SetLineColor(kBlue);

        hGlobalPt->SetTitle(TString::Format("%s; Pt; Efficiency (Matched McPrimary / McPrimary)", label.Data()));
        // hGlobalPt->GetXaxis()->SetRangeUser(1.0, 5.0);
        hGlobalPt->GetYaxis()->SetRangeUser(0.0, 1.2);
        hGlobalPt->Draw();
        hBLCPt->Draw("same");
        hPrimaryPt->Draw("same");

        TLine *line = new TLine(0.0, 0.0, 5.0, 0.0);
        line->SetLineColor(kBlack);
        line->SetLineStyle(2);
        line->DrawLine(0.0, 1.0, 5.0, 1.0);

        TLegend *leg1 = new TLegend(0.6, 0.2, 0.9, 0.5);
        leg1->AddEntry(hGlobalPt, "Global", "l");
        leg1->AddEntry(hBLCPt, "Beamline", "l");
        leg1->AddEntry(hPrimaryPt, "Primary", "l");
        leg1->Draw();
    }

    {
        c->cd(3);
        TH1 *hGlobalPhi = (TH1*)fGlobal->Get("Phi_MatchedPrimary_McPrimary;1");
        TH1 *hBLCPhi = (TH1*)fBLC->Get("Phi_MatchedPrimary_McPrimary;1");
        TH1 *hPrimaryPhi = (TH1*)fPrimary->Get("Phi_MatchedPrimary_McPrimary;1");
        
        hGlobalPhi->SetLineColor(kBlack);
        hBLCPhi->SetLineColor(kRed);
        hPrimaryPhi->SetLineColor(kBlue);

        hGlobalPhi->SetTitle(label + "; #phi; Efficiency (Matched McPrimary / McPrimary)");
        // hGlobalPhi->GetXaxis()->SetRangeUser(1.0, 5.0);

        hGlobalPhi->Draw();
        hBLCPhi->Draw("same");
        hPrimaryPhi->Draw("same");

        TLine *line = new TLine(0.0, 0.0, 5.0, 0.0);
        line->SetLineColor(kBlack);
        line->SetLineStyle(2);
        line->DrawLine(-3.24159, 1.0, 3.24159, 1.0);

        TLegend *leg1 = new TLegend(0.6, 0.2, 0.9, 0.5);
        leg1->AddEntry(hGlobalPhi, "Global", "l");
        leg1->AddEntry(hBLCPhi, "Beamline", "l");
        leg1->AddEntry(hPrimaryPhi, "Primary", "l");
        leg1->Draw();
    }

    c->Print( TString::Format("plots/%s_efficiency_comparison.png", label.Data()));

    {
        c->cd(1);
        TH2 *hGlobalEtaCharge = (TH2*)fGlobal->Get("hMatchedMcPrimaryEtaCharge");
        TH2 *hBLCEtaCharge = (TH2*)fBLC->Get("hMatchedMcPrimaryEtaCharge");
        TH2 *hPrimaryEtaCharge = (TH2*)fPrimary->Get("hMatchedMcPrimaryEtaCharge");
        
        hGlobalEtaCharge->GetXaxis()->SetRangeUser(1.0, 5.0);
        hGlobalEtaCharge->SetLineColor(kBlack);
        hBLCEtaCharge->SetLineColor(kRed);
        hPrimaryEtaCharge->SetLineColor(kBlue);

        TH1 * hGlobalEtaChargeProfile =  hGlobalEtaCharge->ProfileX("hGlobalEtaChargeProfile");
        TH1 * hBLCEtaChargeProfile = hBLCEtaCharge->ProfileX("hBLCEtaChargeProfile");
        TH1 * hPrimaryEtaChargeProfile = hPrimaryEtaCharge->ProfileX("hPrimaryEtaChargeProfile");

        hGlobalEtaChargeProfile->SetLineWidth(5);
        hBLCEtaChargeProfile->SetLineWidth(5);
        hPrimaryEtaChargeProfile->SetLineWidth(5);

        hGlobalEtaChargeProfile->SetTitle( label + "; #eta; Charge MisId");
        hGlobalEtaChargeProfile->Draw();
        hBLCEtaChargeProfile->Draw("same");
        hPrimaryEtaChargeProfile->Draw("same");
        
        TLegend *leg2 = new TLegend(0.1, 0.6, 0.4, 0.9);
        leg2->AddEntry(hGlobalEtaChargeProfile, "Global", "l");
        leg2->AddEntry(hBLCEtaChargeProfile, "Beamline", "l");
        leg2->AddEntry(hPrimaryEtaChargeProfile, "Primary", "l");
        leg2->Draw();
    }

    {
        c->cd(2);
        TH2 *hGlobalPtCharge = (TH2*)fGlobal->Get("hMatchedMcPrimaryPtCharge");
        TH2 *hBLCPtCharge = (TH2*)fBLC->Get("hMatchedMcPrimaryPtCharge");
        TH2 *hPrimaryPtCharge = (TH2*)fPrimary->Get("hMatchedMcPrimaryPtCharge");
        
        hGlobalPtCharge->GetXaxis()->SetRangeUser(0.0, 5.0);
        hGlobalPtCharge->GetYaxis()->SetRangeUser(0.0, 1.1);
        hGlobalPtCharge->SetLineColor(kBlack);
        hBLCPtCharge->SetLineColor(kRed);
        hPrimaryPtCharge->SetLineColor(kBlue);

        TH1 * hGlobalPtChargeProfile =  hGlobalPtCharge->ProfileX("hGlobalPtChargeProfile");
        TH1 * hBLCPtChargeProfile = hBLCPtCharge->ProfileX("hBLCPtChargeProfile");
        TH1 * hPrimaryPtChargeProfile = hPrimaryPtCharge->ProfileX("hPrimaryPtChargeProfile");

        hGlobalPtChargeProfile->SetLineWidth(5);
        hBLCPtChargeProfile->SetLineWidth(5);
        hPrimaryPtChargeProfile->SetLineWidth(5);

        hGlobalPtChargeProfile->SetTitle( label + "; P_{T}; Charge MisId");
        hGlobalPtChargeProfile->GetYaxis()->SetRangeUser(0.0, 1.1);
        hGlobalPtChargeProfile->Draw();
        hBLCPtChargeProfile->Draw("same");
        hPrimaryPtChargeProfile->Draw("same");
        
        TLegend *leg2 = new TLegend(0.1, 0.6, 0.4, 0.9);
        leg2->AddEntry(hGlobalPtChargeProfile, "Global", "l");
        leg2->AddEntry(hBLCPtChargeProfile, "Beamline", "l");
        leg2->AddEntry(hPrimaryPtChargeProfile, "Primary", "l");
        leg2->Draw();
    }

    {
        c->cd(3);
        TH2 *hGlobalPhiCharge = (TH2*)fGlobal->Get("hMatchedMcPrimaryPhiCharge");
        TH2 *hBLCPhiCharge = (TH2*)fBLC->Get("hMatchedMcPrimaryPhiCharge");
        TH2 *hPrimaryPhiCharge = (TH2*)fPrimary->Get("hMatchedMcPrimaryPhiCharge");
        
        // hGlobalPhiCharge->GetXaxis()->SetRangeUser(0.0, 5.0);
        // hGlobalPhiCharge->GetYaxis()->SetRangeUser(0.0, 1.1);
        hGlobalPhiCharge->SetLineColor(kBlack);
        hBLCPhiCharge->SetLineColor(kRed);
        hPrimaryPhiCharge->SetLineColor(kBlue);

        TH1 * hGlobalPhiChargeProfile =  hGlobalPhiCharge->ProfileX("hGlobalPhiChargeProfile");
        TH1 * hBLCPhiChargeProfile = hBLCPhiCharge->ProfileX("hBLCPhiChargeProfile");
        TH1 * hPrimaryPhiChargeProfile = hPrimaryPhiCharge->ProfileX("hPrimaryPhiChargeProfile");

        hGlobalPhiChargeProfile->SetLineWidth(5);
        hBLCPhiChargeProfile->SetLineWidth(5);
        hPrimaryPhiChargeProfile->SetLineWidth(5);

        hGlobalPhiChargeProfile->SetTitle( label + "; #phi; Charge MisId");
        hGlobalPhiChargeProfile->GetYaxis()->SetRangeUser(0.0, 1.1);
        hGlobalPhiChargeProfile->Draw();
        hBLCPhiChargeProfile->Draw("same");
        hPrimaryPhiChargeProfile->Draw("same");
        
        TLegend *leg2 = new TLegend(0.1, 0.6, 0.4, 0.9);
        leg2->AddEntry(hGlobalPhiChargeProfile, "Global", "l");
        leg2->AddEntry(hBLCPhiChargeProfile, "Beamline", "l");
        leg2->AddEntry(hPrimaryPhiChargeProfile, "Primary", "l");
        leg2->Draw();
    }


    c->Print( TString::Format("plots/%s_charge_comparison_eta.png", label.Data()));


    c = new TCanvas("c", "Efficiency Comparison", 3840, 2160/2.0);
    c->Divide(4, 1);
    {
        c->cd(1);
        TH1 *hGlobalCurveRes = (TH1*)fGlobal->Get("hCurveResolution");
        TH1 *hBLCCurveRes = (TH1*)fBLC->Get("hCurveResolution");
        TH1 *hPrimaryCurveRes = (TH1*)fPrimary->Get("hCurveResolution");
        
        hGlobalCurveRes->SetLineColor(kBlack);
        hBLCCurveRes->SetLineColor(kRed);
        hPrimaryCurveRes->SetLineColor(kBlue);

        hGlobalCurveRes->SetLineWidth(5);
        hBLCCurveRes->SetLineWidth(5);
        hPrimaryCurveRes->SetLineWidth(5);

        hGlobalCurveRes->GetYaxis()->SetRangeUser(0.0, hPrimaryCurveRes->GetMaximum() * 1.2);
        hGlobalCurveRes->SetTitle(label + "; Curve Resolution");
        hGlobalCurveRes->Draw();
        hBLCCurveRes->Draw("same");
        hPrimaryCurveRes->Draw("same");

        TLegend *leg3 = new TLegend(0.1, 0.1, 0.4, 0.3);
        leg3->AddEntry(hGlobalCurveRes, "Global", "l");
        leg3->AddEntry(hBLCCurveRes, "Beamline", "l");
        leg3->AddEntry(hPrimaryCurveRes, "Primary", "l");
        leg3->Draw();
    }

    {
        c->cd(2);
        TH2 *hGlobalCurveResEta = (TH2*)fGlobal->Get("hCurveResolutionEta");
        TH2 *hBLCCurveResEta = (TH2*)fBLC->Get("hCurveResolutionEta");
        TH2 *hPrimaryCurveResEta = (TH2*)fPrimary->Get("hCurveResolutionEta");
        
        TH1 * hGlobalCurveResProfileEta =  getResolution(hGlobalCurveResEta, "hGlobalCurveResEta");
        TH1 * hBLCCurveResProfileEta = getResolution(hBLCCurveResEta, "hBLCCurveResEta");
        TH1 * hPrimaryCurveResProfileEta = getResolution(hPrimaryCurveResEta, "hPrimaryCurveResEta");

        hGlobalCurveResProfileEta->SetLineWidth(5);
        hBLCCurveResProfileEta->SetLineWidth(5);
        hPrimaryCurveResProfileEta->SetLineWidth(5);

        hGlobalCurveResProfileEta->SetLineColor(kBlack);
        hBLCCurveResProfileEta->SetLineColor(kRed);
        hPrimaryCurveResProfileEta->SetLineColor(kBlue);

        hGlobalCurveResProfileEta->SetTitle(label + "; #eta; Curve Resolution");
        hGlobalCurveResProfileEta->GetYaxis()->SetRangeUser(0.0, 1.0);
        hGlobalCurveResProfileEta->GetXaxis()->SetRangeUser(2.0, 4.5);
        hGlobalCurveResProfileEta->Draw();
        hBLCCurveResProfileEta->Draw("same");
        hPrimaryCurveResProfileEta->Draw("same");
        
        TLegend *leg3 = new TLegend(0.6, 0.1, 0.9, 0.3);
        leg3->AddEntry(hGlobalCurveResProfileEta, "Global", "l");
        leg3->AddEntry(hBLCCurveResProfileEta, "Beamline", "l");
        leg3->AddEntry(hPrimaryCurveResProfileEta, "Primary", "l");
        leg3->Draw();
    }
    {
        c->cd(3);
        TH2 *hGlobalCurveResPt = (TH2*)fGlobal->Get("hCurveResolutionPt");
        TH2 *hBLCCurveResPt = (TH2*)fBLC->Get("hCurveResolutionPt");
        TH2 *hPrimaryCurveResPt = (TH2*)fPrimary->Get("hCurveResolutionPt");
        
        TH1 * hGlobalCurveResProfilePt =  getResolution(hGlobalCurveResPt, "hGlobalCurveResPt");
        TH1 * hBLCCurveResProfilePt = getResolution(hBLCCurveResPt, "hBLCCurveResPt");
        TH1 * hPrimaryCurveResProfilePt = getResolution(hPrimaryCurveResPt, "hPrimaryCurveResPt");

        hGlobalCurveResProfilePt->SetLineWidth(5);
        hBLCCurveResProfilePt->SetLineWidth(5);
        hPrimaryCurveResProfilePt->SetLineWidth(5);

        hGlobalCurveResProfilePt->SetLineColor(kBlack);
        hBLCCurveResProfilePt->SetLineColor(kRed);
        hPrimaryCurveResProfilePt->SetLineColor(kBlue);

        hGlobalCurveResProfilePt->SetTitle(label + "; P_{T}; Curve Resolution");
        hGlobalCurveResProfilePt->GetYaxis()->SetRangeUser(0.0, 1.0);
        hGlobalCurveResProfilePt->GetXaxis()->SetRangeUser(0.0, 5.0);
        hGlobalCurveResProfilePt->Draw();
        hBLCCurveResProfilePt->Draw("same");
        hPrimaryCurveResProfilePt->Draw("same");
        
        TLegend *leg3 = new TLegend(0.6, 0.1, 0.9, 0.3);
        leg3->AddEntry(hGlobalCurveResProfilePt, "Global", "l");
        leg3->AddEntry(hBLCCurveResProfilePt, "Beamline", "l");
        leg3->AddEntry(hPrimaryCurveResProfilePt, "Primary", "l");
        leg3->Draw();
    }
    {
        c->cd(4);
        TH2 *hGlobalCurveResPhi = (TH2*)fGlobal->Get("hCurveResolutionPhi");
        TH2 *hBLCCurveResPhi = (TH2*)fBLC->Get("hCurveResolutionPhi");
        TH2 *hPrimaryCurveResPhi = (TH2*)fPrimary->Get("hCurveResolutionPhi");
        
        TH1 * hGlobalCurveResProfilePhi =  getResolution(hGlobalCurveResPhi, "hGlobalCurveResPhi");
        TH1 * hBLCCurveResProfilePhi = getResolution(hBLCCurveResPhi, "hBLCCurveResPhi");
        TH1 * hPrimaryCurveResProfilePhi = getResolution(hPrimaryCurveResPhi, "hPrimaryCurveResPhi");

        hGlobalCurveResProfilePhi->SetLineWidth(5);
        hBLCCurveResProfilePhi->SetLineWidth(5);
        hPrimaryCurveResProfilePhi->SetLineWidth(5);

        hGlobalCurveResProfilePhi->SetLineColor(kBlack);
        hBLCCurveResProfilePhi->SetLineColor(kRed);
        hPrimaryCurveResProfilePhi->SetLineColor(kBlue);

        hGlobalCurveResProfilePhi->SetTitle(label + "; #phi; Curve Resolution");
        hGlobalCurveResProfilePhi->GetYaxis()->SetRangeUser(0.0, 1.0);
        hGlobalCurveResProfilePhi->GetXaxis()->SetRangeUser(0.0, 5.0);
        hGlobalCurveResProfilePhi->Draw();
        hBLCCurveResProfilePhi->Draw("same");
        hPrimaryCurveResProfilePhi->Draw("same");
        
        TLegend *leg3 = new TLegend(0.6, 0.1, 0.9, 0.3);
        leg3->AddEntry(hGlobalCurveResProfilePhi, "Global", "l");
        leg3->AddEntry(hBLCCurveResProfilePhi, "Beamline", "l");
        leg3->AddEntry(hPrimaryCurveResProfilePhi, "Primary", "l");
        leg3->Draw();
    }

    c->Print( TString::Format("plots/%s_curve_resolution.png", label.Data()));


    c = new TCanvas("cMult", "Efficiency Comparison", 3840/2.0, 2160/2.0);
    // c->Divide(1, 1);
    {
        c->cd();
        TH1 *hGlobalMultEff = (TH1*)fGlobal->Get("Mult_MatchedPrimary_McPrimary");
        TH1 *hBLCMultEff = (TH1*)fBLC->Get("Mult_MatchedPrimary_McPrimary");
        TH1 *hPrimaryMultEff = (TH1*)fPrimary->Get("Mult_MatchedPrimary_McPrimary");
        
        hGlobalMultEff->SetLineColor(kBlack);
        hBLCMultEff->SetLineColor(kRed);
        hPrimaryMultEff->SetLineColor(kBlue);

        hGlobalMultEff->SetLineWidth(5);
        hBLCMultEff->SetLineWidth(5);
        hPrimaryMultEff->SetLineWidth(5);

        hGlobalMultEff->GetXaxis()->SetRangeUser(0.0, 30.0);
        hGlobalMultEff->SetTitle(label + "; Multiplicity; Efficiency (Matched McPrimary / McPrimary)");
        hGlobalMultEff->Draw();
        hBLCMultEff->Draw("same");
        hPrimaryMultEff->Draw("same");

        TLegend *leg3 = new TLegend(0.6, 0.1, 0.9, 0.3);
        leg3->AddEntry(hGlobalMultEff, "Global", "l");
        leg3->AddEntry(hBLCMultEff, "Beamline", "l");
        leg3->AddEntry(hPrimaryMultEff, "Primary", "l");
        leg3->Draw();
    }
    c->Print( TString::Format("plots/%s_eff_vs_mult.png", label.Data()));
}