const char* CPID[10]={"?","gamma","e+","e-","nu","mu+","mu-","pi0","pi+","pi-"};
int PID=0;

void plot(int plt=0, int pid, int run=0){
  PID=pid;
  plot(plt,CPID[pid],run);
}

void plot(int plt=0, char* data="JPsi", int run=0){
  
  TString Data(data);
  char file[100];
  if(Data.Contains("data")){
    if(run==0){
      sprintf(file,"hist/st_fwd.match.root");
    }else{
      sprintf(file,"hist/st_fwd_23024030_raw_7500071.match.root");
    }
  }else if(PID>0){
    sprintf(file,"out/fast_track_pid%d/job%d_fcstrk.root",PID,run); 
  }else{
    if(run==0){
      sprintf(file,"hist/pythia_%s_vz0.match.root",data);      
    }else{
      sprintf(file,"hist/pythia_%s_vz0_run%d.match.root",data,run);
    }
  }
  printf("Reading %s\n",file);
  TFile *F = new TFile(file,"old");
  
  c1 = new TCanvas("c1","FCS-TRK",50,0,1500,1200);
  gStyle->SetLabelSize(0.1,"xy");
  gStyle->SetPalette(1);
  gStyle->SetStatW(0.4);
  
  TH1F* h1;
  TH2F* h2;
  char hname[100];
  TText* t;

  if(plt==0 || plt==1) {
    c1->Clear();
    c1->Divide(3,2);    
    c1->cd(1); h1=(TH1F*)F->Get("NTrk");     h1->Draw();
    c1->cd(2); h1=(TH1F*)F->Get("NGoodTrk"); h1->Draw();
    c1->cd(3); h1=(TH1F*)F->Get("Charge");   h1->Draw();
    c1->cd(4); h2=(TH2F*)F->Get("XY");       h2->Draw("colz");
    c1->cd(5); h1=(TH1F*)F->Get("NEcalClu"); h1->Draw();
    c1->cd(6); h1=(TH1F*)F->Get("NHcalClu"); h1->Draw();
    c1->SaveAs(Form("www/trackmatch/trk_%s.png",data));
  }

  if(plt==0 || plt==2) {
    c1->Clear();
    c1->Divide(3,3);    
    c1->cd(1); h1=(TH1F*)F->Get("dx_EcalTrk");   h1->Draw();
    c1->cd(2); h1=(TH1F*)F->Get("dy_EcalTrk");   h1->Draw();
    c1->cd(3); h1=(TH1F*)F->Get("dr_EcalTrk");   h1->Draw();
    c1->cd(4); h1=(TH1F*)F->Get("NTrk_Ecal");    h1->Draw();
    c1->cd(5); h1=(TH1F*)F->Get("NEcalClu_Trk"); h1->Draw();
    c1->cd(6); h1=(TH1F*)F->Get("ETovPT_E");     h1->Draw();
    c1->cd(7)->SetLogz(); h2=(TH2F*)F->Get("ETPT_E"); h2->Draw("colz");
    c1->cd(8); h1=(TH1F*)F->Get("Charge_E");     h1->Draw();
    c1->cd(9); h2=(TH2F*)F->Get("XY_E");         h2->Draw("colz");
    c1->SaveAs(Form("www/trackmatch/ecaltrk_%s.png",data));
  }

  if(plt==0 || plt==3) {
    c1->Clear();
    c1->Divide(3,3);    
    c1->cd(1); h1=(TH1F*)F->Get("dx_HcalTrk");   h1->Draw();
    c1->cd(2); h1=(TH1F*)F->Get("dy_HcalTrk");   h1->Draw();
    c1->cd(3); h1=(TH1F*)F->Get("dr_HcalTrk");   h1->Draw();
    c1->cd(4); h1=(TH1F*)F->Get("NTrk_Hcal");    h1->Draw();
    c1->cd(5); h1=(TH1F*)F->Get("NHcalClu_Trk"); h1->Draw();
    c1->cd(6); h1=(TH1F*)F->Get("ETovPT_EH");    h1->Draw();
    c1->cd(7)->SetLogz();; h2=(TH2F*)F->Get("ETPT_EH");      h2->Draw("colz");
    c1->cd(8); h1=(TH1F*)F->Get("Charge_E");     h1->Draw();
    c1->cd(9); h2=(TH2F*)F->Get("XY_H");         h2->Draw("colz");
    c1->SaveAs(Form("www/trackmatch/hcaltrk_%s.png",data));
  }
}
