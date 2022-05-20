void plot(int plt=0, int pid=5, int job=1){
  const char* cpid[10] = {"?","Gamma","e+","e-","N","mu+","mu-","pi0","pi+","pi-"};
  char file[100];
  sprintf(file,"out/fast_track_pid%d/job%d_fcstrk.root",pid,job);
  printf("Reading %s\n",file);
  TFile *F = new TFile(file,"old");
  
  c1 = new TCanvas("c1","FCS-TRK",50,0,1500,1200);
  gStyle->SetLabelSize(0.1,"xy");
  gStyle->SetPalette(1);
  gStyle->SetStatW(0.4);
  
  TH1F* h1;
  TH2F* h2;
  char hname[100];

  if(plt==0 || plt==1) {
    c1->Clear();
    c1->Divide(3,3);    
    c1->cd(1); h1=(TH1F*)F->Get("dx_EcalTrk");   h1->Draw();
    c1->cd(2); h1=(TH1F*)F->Get("dy_EcalTrk");   h1->Draw();
    c1->cd(3); h1=(TH1F*)F->Get("dr_EcalTrk");   h1->Draw();
    c1->cd(4); h1=(TH1F*)F->Get("NTrk_Ecal");    h1->Draw();
    c1->cd(5); h1=(TH1F*)F->Get("NEcalClu_Trk"); h1->Draw();
    c1->cd(6); h1=(TH1F*)F->Get("ETovPT_E");     h1->Draw();
    c1->cd(7); h2=(TH2F*)F->Get("ETPT_E");       h2->Draw("colz");
    c1->cd(8); h1=(TH1F*)F->Get("Charge_E");     h1->Draw();
    c1->SaveAs(Form("ecaltrk_%s.png",cpid[pid]));
  }

  if(plt==0 || plt==2) {
    c1->Clear();
    c1->Divide(3,3);    
    c1->cd(1); h1=(TH1F*)F->Get("dx_HcalTrk");   h1->Draw();
    c1->cd(2); h1=(TH1F*)F->Get("dy_HcalTrk");   h1->Draw();
    c1->cd(3); h1=(TH1F*)F->Get("dr_HcalTrk");   h1->Draw();
    c1->cd(4); h1=(TH1F*)F->Get("NTrk_Hcal");    h1->Draw();
    c1->cd(5); h1=(TH1F*)F->Get("NHcalClu_Trk"); h1->Draw();
    c1->cd(6); h1=(TH1F*)F->Get("ETovPT_EH");    h1->Draw();
    c1->cd(7); h2=(TH2F*)F->Get("ETPT_EH");      h2->Draw("colz");
    c1->cd(8); h1=(TH1F*)F->Get("Charge_EH");    h1->Draw();
    c1->SaveAs(Form("hcaltrk_%s.png",cpid[pid]));
  }
}
