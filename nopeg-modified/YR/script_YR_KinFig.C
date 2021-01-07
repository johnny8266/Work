void script(const char* name, int part, int config)
{
  //Please label the plot with indicating the particle ID. Normalize your counts to an integrated luminosity of 10 fb-1
 TString title;//(#intL=10 fb^{-1} )"); // Title of your plot
  switch(part){
case 1: title="Electron Momentum vs Theta"; break;//(#intL=10 fb^{-1} )"); // Title of your plot
case 2: title="Photon Momentum vs Theta";   break;//(#intL=10 fb^{-1} )"); // Title of your plot
case 3: title="Helium-4 Momentum vs Theta"; break;//(#intL=10 fb^{-1} )"); // Title of your plot
}

  switch(config){
case 1: title+=" ( 5x 41)"; break;
case 2: title+=" (10x110)"; break;
case 3: title+=" (18x110)"; break;
}

  TFile *f=new TFile(name);
  TTree *tree=(TTree*)f->Get("TOPEG");

  //Please keep the following binning so that we can compare across different WG
  TH2D *h=new TH2D("h",title.Data(),36,0.,TMath::Pi(),100.,0.,50.);

  //Fill your 2D histogram with Theta (in deg) and Momenta (in GeV) for a given particle ID
  switch(part){
case 1: tree->Draw("sqrt(part_px*part_px+part_py*part_py+part_pz*part_pz):acos(part_pz/sqrt(part_px*part_px+part_py*part_py+part_pz*part_pz))>>h","part_id==11");         break;
case 2: tree->Draw("sqrt(part_px*part_px+part_py*part_py+part_pz*part_pz):acos(part_pz/sqrt(part_px*part_px+part_py*part_py+part_pz*part_pz))>>h","part_id==22");         break;
case 3: tree->Draw("sqrt(part_px*part_px+part_py*part_py+part_pz*part_pz):acos(part_pz/sqrt(part_px*part_px+part_py*part_py+part_pz*part_pz))>>h","part_id==1000020040"); break;
}

 
  // --------- Plotting ---------------- //
  //     (code by R. Seidl) 

  gStyle->SetOptStat(0);
  gStyle->SetTitleX(.5);
  gStyle->SetTitleY(.87);

  TH2F *hdummy=new TH2F("hdummy",title.Data(),10,-18,5., 10, 0., 7.);

  Double_t xtit = h->GetYaxis()->GetTitleSize();
  xtit=xtit*1.6;
  
  Double_t maxy = h->GetMaximum()*1.1;
  Double_t miny = h->GetMinimum(0.);
  
  hdummy->GetZaxis()->SetRangeUser(miny,maxy);
  hdummy->GetYaxis()->SetTitleSize(xtit);
  hdummy->GetXaxis()->SetTitleSize(xtit);
  
  hdummy->GetYaxis()->SetLabelSize(xtit);
  hdummy->GetXaxis()->SetLabelSize(xtit);
  
  // hdummy->GetYaxis()->SetTitle("Momentum (GeV)  ");
  hdummy->GetXaxis()->SetTitle("Momentum (GeV)  ");
  
  hdummy->GetXaxis()->SetTitleOffset(1.2);     
  hdummy->GetYaxis()->SetTitleOffset(1.3);
  hdummy->SetLineColor(1);
  hdummy->SetLineStyle(1);
  hdummy->SetLineWidth(2);
  
  hdummy->DrawCopy("AXIS");
  delete hdummy;

  gPad->SetRightMargin(0.15);
  gPad->SetBottomMargin(0.15);
  gPad->SetLogz();
 
  h->DrawCopy("POL COLZ SAME");

  // ---- Output -------

  TFile *fout=new TFile("output_plot.root","recreate");
  h->Write(title.Data());
  gPad->Write(title.Data());
  title+=".gif";
  gPad->Print(title.Data()); 

}
