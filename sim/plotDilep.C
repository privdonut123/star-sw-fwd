static const int mNCut=7;
const char* nameCut[mNCut] = {"All","ETOT","HTOT","Cone","SigmaMax","TrackMatch","ChargeSign"};

void plotDilep(int plt=0, char* data="JPsi", int run=0, int cut=6, int log=1){
  TString Data(data);
  char file[100];
  if(Data.Contains("data")){
    if(run==0){
      sprintf(file,"hist/st_fwd.dilep.root");
    }else{
      sprintf(file,"hist/st_fwd_23024030_raw_7500071.dilep.root");
    }
  }else{
    if(run==0){
      sprintf(file,"hist/pythia_%s_vz0.dilep.root",data);      
    }else{
      sprintf(file,"hist/pythia_%s_vz0_run%d.dilep.root",data,run);
    }
  }
  printf("Reading %s\n",file);
  TFile *F = new TFile(file,"old");
  
  c1 = new TCanvas("c1","FCS DiLepton",50,0,1500,1200);
  gStyle->SetLabelSize(0.1,"xy");
  gStyle->SetPalette(1);
  gStyle->SetStatW(0.4);
    
  TH1F* h1;
  TH2F* h2;
  char hname[100];
  if(plt==0 || plt==1) {
    c1->Clear();
    c1->Divide(2,2);
    
    c1->cd(1); h2=(TH2F*)F->Get(Form("XFPT_%s",nameCut[cut]));   h2->Draw("colz");
    c1->cd(2); h2=(TH2F*)F->Get(Form("ET12_%s",nameCut[cut]));   h2->Draw("colz");
    c1->cd(3); h2=(TH2F*)F->Get(Form("XY_%s",nameCut[cut]));     h2->Draw("colz");
    c1->cd(4); h2=(TH2F*)F->Get(Form("PTET_%s",nameCut[cut]));   h2->Draw("colz");
    c1->cd(0); t = new TText(0.4,0.95,Form("Cut=%s",nameCut[cut])); t->Draw();

    c1->SaveAs(Form("dilep.cut%d.%s.%d.png",cut,data,run));
  }

  if(plt==0 || plt==2) {
    c1->Clear();
    c1->Divide(2,3);
    
    int c=0;
    c1->cd(1)->SetLogy(log); h1=(TH1F*)F->Get(Form("RETot_%s",nameCut[c]));   h1->Draw();
    c1->cd(2)->SetLogy(log); h1=(TH1F*)F->Get(Form("RHTot_%s",nameCut[c]));   h1->Draw();
    c1->cd(3)->SetLogy(log); h1=(TH1F*)F->Get(Form("RCone_%s",nameCut[c]));   h1->Draw();
    c1->cd(4)->SetLogy(log); h1=(TH1F*)F->Get(Form("Sigmax_%s",nameCut[c]));  h1->Draw();
    c1->cd(5)->SetLogy(log); h1=(TH1F*)F->Get(Form("PToverET_%s",nameCut[c]));  h1->SetMinimum(0.2); h1->Draw();
    c1->cd(6)->SetLogy(log); h1=(TH1F*)F->Get(Form("ChargeSum_%s",nameCut[c]));  h1->SetMinimum(0.2); h1->Draw();

    c=4;
    c1->cd(1)->SetLogy(log); h1=(TH1F*)F->Get(Form("RETot_%s",nameCut[c]));   h1->SetLineColor(kBlue); h1->Draw("same");
    c1->cd(2)->SetLogy(log); h1=(TH1F*)F->Get(Form("RHTot_%s",nameCut[c]));   h1->SetLineColor(kBlue); h1->Draw("same");
    c1->cd(3)->SetLogy(log); h1=(TH1F*)F->Get(Form("RCone_%s",nameCut[c]));   h1->SetLineColor(kBlue); h1->Draw("same");
    c1->cd(4)->SetLogy(log); h1=(TH1F*)F->Get(Form("Sigmax_%s",nameCut[c]));  h1->SetLineColor(kBlue); h1->Draw("same");
    c1->cd(5)->SetLogy(log); h1=(TH1F*)F->Get(Form("PToverET_%s",nameCut[c]));  h1->SetLineColor(kBlue); h1->Draw("same");
    c1->cd(6)->SetLogy(log); h1=(TH1F*)F->Get(Form("ChargeSum_%s",nameCut[c])); h1->SetLineColor(kBlue);  h1->Draw("same");

    c=5;
    c1->cd(1)->SetLogy(log); h1=(TH1F*)F->Get(Form("RETot_%s",nameCut[c]));   h1->SetLineColor(kMagenta); h1->Draw("same");
    c1->cd(2)->SetLogy(log); h1=(TH1F*)F->Get(Form("RHTot_%s",nameCut[c]));   h1->SetLineColor(kMagenta); h1->Draw("same");
    c1->cd(3)->SetLogy(log); h1=(TH1F*)F->Get(Form("RCone_%s",nameCut[c]));   h1->SetLineColor(kMagenta); h1->Draw("same");
    c1->cd(4)->SetLogy(log); h1=(TH1F*)F->Get(Form("Sigmax_%s",nameCut[c]));  h1->SetLineColor(kMagenta); h1->Draw("same");
    c1->cd(5)->SetLogy(log); h1=(TH1F*)F->Get(Form("PToverET_%s",nameCut[c]));  h1->SetLineColor(kMagenta); h1->Draw("same");
    c1->cd(6)->SetLogy(log); h1=(TH1F*)F->Get(Form("ChargeSum_%s",nameCut[c])); h1->SetLineColor(kMagenta);  h1->Draw("same");

    c=6;
    c1->cd(1)->SetLogy(log); h1=(TH1F*)F->Get(Form("RETot_%s",nameCut[c]));   h1->SetLineColor(kRed); h1->Draw("same");
    c1->cd(2)->SetLogy(log); h1=(TH1F*)F->Get(Form("RHTot_%s",nameCut[c]));   h1->SetLineColor(kRed); h1->Draw("same");
    c1->cd(3)->SetLogy(log); h1=(TH1F*)F->Get(Form("RCone_%s",nameCut[c]));   h1->SetLineColor(kRed); h1->Draw("same");
    c1->cd(4)->SetLogy(log); h1=(TH1F*)F->Get(Form("Sigmax_%s",nameCut[c]));  h1->SetLineColor(kRed); h1->Draw("same");
    c1->cd(5)->SetLogy(log); h1=(TH1F*)F->Get(Form("PToverET_%s",nameCut[c]));  h1->SetLineColor(kRed); h1->Draw("same");
    c1->cd(6)->SetLogy(log); h1=(TH1F*)F->Get(Form("ChargeSum_%s",nameCut[c])); h1->SetLineColor(kRed);  h1->Draw("same");

    c1->cd(1);
    t=new TText(0.15, 0.8,"No cut");        t->SetTextColor(kBlack);   t->SetNDC(); t->Draw();
    t=new TText(0.15, 0.7,"FCS cuts");      t->SetTextColor(kBlue);    t->SetNDC(); t->Draw();
    t=new TText(0.15, 0.6,"Track Matched"); t->SetTextColor(kMagenta); t->SetNDC(); t->Draw();
    t=new TText(0.15, 0.5,"Charge Sign");   t->SetTextColor(kRed);     t->SetNDC(); t->Draw();
    c1->SaveAs(Form("dilep1.%s.%d.png",data,run));
  }
  
  if(plt==0 || plt==3) {
    c1->Clear();
    c1->Divide(2,3);    

    int c=0;
    c1->cd(1)->SetLogy(log); h1=(TH1F*)F->Get(Form("ET_%s",nameCut[c]));  h1->SetMinimum(0.2);   h1->Draw();
    c1->cd(2)->SetLogy(log); h1=(TH1F*)F->Get(Form("EZ_%s",nameCut[c]));  h1->SetMinimum(0.2);    h1->Draw();
    c1->cd(3)->SetLogy(log); h1=(TH1F*)F->Get(Form("M_%s",nameCut[c]));   h1->SetMinimum(0.2);    h1->Draw();
    c1->cd(4)->SetLogy(log); h1=(TH1F*)F->Get(Form("Z_%s",nameCut[c]));   h1->SetMinimum(0.2);    h1->Draw();
    c1->cd(5)->SetLogy(log); h1=(TH1F*)F->Get(Form("CosT_%s",nameCut[c]));h1->SetMinimum(0.2); h1->Draw();
    c1->cd(6)->SetLogy(log); h1=(TH1F*)F->Get(Form("Phi_%s",nameCut[c])); h1->SetMinimum(0.2); h1->Draw();

    c=4;
    c1->cd(1)->SetLogy(log); h1=(TH1F*)F->Get(Form("ET_%s",nameCut[c]));   h1->SetLineColor(kBlue); h1->Draw("same");
    c1->cd(2)->SetLogy(log); h1=(TH1F*)F->Get(Form("EZ_%s",nameCut[c]));   h1->SetLineColor(kBlue); h1->Draw("same");
    c1->cd(3)->SetLogy(log); h1=(TH1F*)F->Get(Form("M_%s",nameCut[c]));    h1->SetLineColor(kBlue); h1->Draw("same");
    c1->cd(4)->SetLogy(log); h1=(TH1F*)F->Get(Form("Z_%s",nameCut[c]));    h1->SetLineColor(kBlue); h1->Draw("same");
    c1->cd(5)->SetLogy(log); h1=(TH1F*)F->Get(Form("CosT_%s",nameCut[c])); h1->SetLineColor(kBlue); h1->Draw("same");
    c1->cd(6)->SetLogy(log); h1=(TH1F*)F->Get(Form("Phi_%s",nameCut[c]));  h1->SetLineColor(kBlue); h1->Draw("same");

    c=5;
    c1->cd(1)->SetLogy(log); h1=(TH1F*)F->Get(Form("ET_%s",nameCut[c]));   h1->SetLineColor(kMagenta); h1->Draw("same");
    c1->cd(2)->SetLogy(log); h1=(TH1F*)F->Get(Form("EZ_%s",nameCut[c]));   h1->SetLineColor(kMagenta); h1->Draw("same");
    c1->cd(3)->SetLogy(log); h1=(TH1F*)F->Get(Form("M_%s",nameCut[c]));    h1->SetLineColor(kMagenta); h1->Draw("same");
    c1->cd(4)->SetLogy(log); h1=(TH1F*)F->Get(Form("Z_%s",nameCut[c]));    h1->SetLineColor(kMagenta); h1->Draw("same");
    c1->cd(5)->SetLogy(log); h1=(TH1F*)F->Get(Form("CosT_%s",nameCut[c])); h1->SetLineColor(kMagenta); h1->Draw("same");
    c1->cd(6)->SetLogy(log); h1=(TH1F*)F->Get(Form("Phi_%s",nameCut[c]));  h1->SetLineColor(kMagenta); h1->Draw("same");

    c=6;
    c1->cd(1)->SetLogy(log); h1=(TH1F*)F->Get(Form("ET_%s",nameCut[c]));   h1->SetLineColor(kRed); h1->Draw("same");
    c1->cd(2)->SetLogy(log); h1=(TH1F*)F->Get(Form("EZ_%s",nameCut[c]));   h1->SetLineColor(kRed); h1->Draw("same");
    c1->cd(3)->SetLogy(log); h1=(TH1F*)F->Get(Form("M_%s",nameCut[c]));    h1->SetLineColor(kRed); h1->Draw("same");
    c1->cd(4)->SetLogy(log); h1=(TH1F*)F->Get(Form("Z_%s",nameCut[c]));    h1->SetLineColor(kRed); h1->Draw("same");
    c1->cd(5)->SetLogy(log); h1=(TH1F*)F->Get(Form("CosT_%s",nameCut[c])); h1->SetLineColor(kRed); h1->Draw("same");
    c1->cd(6)->SetLogy(log); h1=(TH1F*)F->Get(Form("Phi_%s",nameCut[c]));  h1->SetLineColor(kRed); h1->Draw("same");
    
    c1->SaveAs(Form("dilep2.%s.%d.png",data,run));
  }
}
