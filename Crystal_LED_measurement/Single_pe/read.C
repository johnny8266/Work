void read()
{
  TCanvas *c1 = new TCanvas("c1", "c1", 800, 800);
  TH1D* h1;
  
  TFile* gfile = TFile::Open("./pe_distribution_clear/Really_clear_pe_5minutes_per_intensity_930_1018mV_8mv_interval.root");
  TDirectory* dir = gFile->GetDirectory("Energy");

  dir->GetObject("_R_EnergyCH0@DT5730_1204", h1);
  h1->SetStats(0);
  h1->RebinX(3);
  h1->GetXaxis()->SetRangeUser(0, 500);
  h1->Draw();

}
