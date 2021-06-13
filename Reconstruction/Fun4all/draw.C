

void draw()
{
  //  TFile *f = new TFile("inner_result.root");
  //  TFile *f = new TFile("outer_result.root");
  TFile *f = new TFile("result.root");

  TGraphErrors *Pion_reject_glass_3pars = (TGraphErrors*)f->Get("Pion_reject_glass_3pars");
  TGraphErrors *Pion_reject_glass_2pars = (TGraphErrors*)f->Get("Pion_reject_glass_2pars");
  TGraphErrors *Pion_reject_crystal_3pars = (TGraphErrors*)f->Get("Pion_reject_crystal_3pars");
  TGraphErrors *Pion_reject_crystal_2pars = (TGraphErrors*)f->Get("Pion_reject_crystal_2pars");
  //  TH1F *P_E_distribution_1, *P_E_distribution_2;
  // TH2F *P_E_distribution_pi = (TH2F*)f->Get("P_E_distribution_pi");
  // TH2F *P_E_distribution_gamma = (TH2F*)f->Get("P_E_distribution_gamma");
  

  TCanvas *c1 = new TCanvas("c1", "c1", 800, 800);
  TLegend *legend_2_par;
  legend_2_par = new TLegend(0.55, 0.35, 0.75, 0.5);
  legend_2_par->SetBorderSize(0);
  gPad->SetLogy();
  gPad->SetLogx();
  Pion_reject_glass_2pars->SetTitle("Rejection factor [2 pars case]");
  Pion_reject_glass_2pars->SetLineColor(1);
  Pion_reject_glass_2pars->SetMarkerStyle(24);
  Pion_reject_glass_2pars->SetMarkerSize(1);
  Pion_reject_glass_2pars->SetMarkerColor(1);
  Pion_reject_glass_2pars->Draw("");
  Pion_reject_crystal_2pars->SetLineColor(2);
  Pion_reject_crystal_2pars->SetMarkerStyle(25);
  Pion_reject_crystal_2pars->SetMarkerSize(1);
  Pion_reject_crystal_2pars->SetMarkerColor(2);
  Pion_reject_crystal_2pars->Draw("same");
  legend_2_par->AddEntry(Pion_reject_glass_2pars, "glass", "p");
  legend_2_par->AddEntry(Pion_reject_crystal_2pars, "crystal", "p");
  legend_2_par->Draw("same");
  
  TCanvas *c2 = new TCanvas("c2", "c2", 800, 800);
  TLegend *legend_3_par;
  legend_3_par = new TLegend(0.55, 0.35, 0.75, 0.5);
  legend_3_par->SetBorderSize(0);
  gPad->SetLogy();
  gPad->SetLogx();
  Pion_reject_glass_3pars->SetTitle("Rejection factor [3 pars case]");
  Pion_reject_glass_3pars->SetLineColor(1);
  Pion_reject_glass_3pars->SetMarkerStyle(24);
  Pion_reject_glass_3pars->SetMarkerSize(1);
  Pion_reject_glass_3pars->SetMarkerColor(1);
  Pion_reject_glass_3pars->Draw("");
  Pion_reject_crystal_3pars->SetLineColor(2);
  Pion_reject_crystal_3pars->SetMarkerStyle(25);
  Pion_reject_crystal_3pars->SetMarkerSize(1);
  Pion_reject_crystal_3pars->SetMarkerColor(2);
  Pion_reject_crystal_3pars->Draw("same");
  legend_3_par->AddEntry(Pion_reject_glass_3pars, "glass", "p");
  legend_3_par->AddEntry(Pion_reject_crystal_3pars, "crystal", "p");
  legend_3_par->Draw("same");
    
}
