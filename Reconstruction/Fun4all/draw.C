double_t E_resolu_fit_quardratic_sum_3(double *x, double *par)
{
  return sqrt( par[0] * par[0] + (par[1] * par[1]) / x[0] + (par[2] * par[2]) / (x[0] * x[0]) );
}

double_t E_resolu_fit_quardratic_sum_2(double *x, double *par)
{
  return sqrt( par[0] * par[0] + (par[1] * par[1]) / x[0] );
}



void draw()
{
  TFile *f = new TFile("result.root");
  TFile *g = new TFile("E_Reso_result.root");
  


  // ================================
  // Draw the position resolution
  // ================================



  /*  
  // ================================
  // Draw the E_resolution resolution 
  // ================================
  TGraphErrors *Eresolution_glass = (TGraphErrors*)g->Get("Eresolution_glass");
  TGraphErrors *Eresolution_crystal = (TGraphErrors*)g->Get("Eresolution_crystal");

  TCanvas *c2 = new TCanvas("c2", "c2", 900, 900);
  TLegend *legend_e_reso;
  legend_e_reso = new TLegend(0.3, 0.45, 0.85, 0.85);
  legend_e_reso->SetBorderSize(0);
  
  auto *fun_recons_e_res_2_glass = new TF1("E_rfit_quardratic_sum_2_glass", E_resolu_fit_quardratic_sum_2, 0., 19., 2);
  fun_recons_e_res_2_glass->SetParLimits(0, 0.1, 10.);
  fun_recons_e_res_2_glass->SetParLimits(1, 0.01, 10.);
  
  auto *fun_recons_e_res_3_glass = new TF1("E_rfit_quardratic_sum_3_glass", E_resolu_fit_quardratic_sum_3, 0., 19., 3);
  fun_recons_e_res_3_glass->SetParLimits(0, 0.1, 10.);
  fun_recons_e_res_3_glass->SetParLimits(1, 0.01, 10.);
  fun_recons_e_res_3_glass->SetParLimits(2, 0.01, 10.);

  auto *fun_recons_e_res_2_crystal = new TF1("E_rfit_quardratic_sum_2_crystal", E_resolu_fit_quardratic_sum_2, 0., 19., 2);
  fun_recons_e_res_2_crystal->SetParLimits(0, 0.1, 10.);
  fun_recons_e_res_2_crystal->SetParLimits(1, 0.01, 10.);
  
  auto *fun_recons_e_res_3_crystal = new TF1("E_rfit_quardratic_sum_3_crystal", E_resolu_fit_quardratic_sum_3, 0., 19., 3);
  fun_recons_e_res_3_crystal->SetParLimits(0, 0.1, 10.);
  fun_recons_e_res_3_crystal->SetParLimits(1, 0.01, 10.);
  fun_recons_e_res_3_crystal->SetParLimits(2, 0.01, 10.);



  Eresolution_glass->SetTitle("E resolution of inner[crystal] and outer[sci-glass]");
  Eresolution_glass->SetMarkerStyle(32);
  Eresolution_glass->SetMarkerColor(2);
  Eresolution_glass->SetMarkerSize(1.5);
  Eresolution_glass->GetXaxis()->SetTitle("E [GeV]");
  Eresolution_glass->GetYaxis()->SetTitle(" #frac{#sigma}{E} [%]");

  Eresolution_glass->Draw("AP");
  fun_recons_e_res_2_glass->SetLineColor(2);
  fun_recons_e_res_2_glass->SetLineStyle(4);
  Eresolution_glass->Fit("E_rfit_quardratic_sum_2_glass", "", "R", 0., 19.);
  fun_recons_e_res_2_glass->Draw("same");
  fun_recons_e_res_3_glass->SetLineColor(2);
  fun_recons_e_res_3_glass->SetLineStyle(1);
  Eresolution_glass->Fit("E_rfit_quardratic_sum_3_glass", "", "R", 0., 19.);
  fun_recons_e_res_3_glass->Draw("same");
  cout << fun_recons_e_res_2_glass->GetChisquare() << endl;

  double alpha_2, beta_2;
  double alpha_3, beta_3, gamma_3;

  alpha_2 = fun_recons_e_res_2_glass->GetParameter(0);  beta_2 = fun_recons_e_res_2_glass->GetParameter(1);
  alpha_3 = fun_recons_e_res_3_glass->GetParameter(0);  beta_3 = fun_recons_e_res_3_glass->GetParameter(1);  gamma_3 = fun_recons_e_res_3_glass->GetParameter(2);
  
  legend_e_reso->AddEntry("E_rfit_quardratic_sum_2_glass",
			  Form("#frac{#sigma}{E} = %.2f #oplus #frac{%.2f}{#sqrt{E}}        [sci-glass]", alpha_2, beta_2),
			  "l");
  legend_e_reso->AddEntry((TObject*)0, "", "");
  legend_e_reso->AddEntry("E_rfit_quardratic_sum_3_glass",
			  Form("#frac{#sigma}{E} = %.2f #oplus #frac{%.2f}{#sqrt{E}} #oplus #frac{%.2f}{E}", alpha_3, beta_3, gamma_3),
			  "l");
  legend_e_reso->AddEntry((TObject*)0, "", "");
  legend_e_reso->AddEntry((TObject*)0, "", "");
  //  legend_e_reso->Draw("same");

  

  Eresolution_crystal->SetMarkerStyle(24);
  Eresolution_crystal->SetMarkerColor(4);
  Eresolution_crystal->SetMarkerSize(1.5);
  Eresolution_crystal->GetXaxis()->SetTitle("E [GeV]");
  Eresolution_crystal->GetYaxis()->SetTitle("E_sig / E [%]");
  Eresolution_crystal->Draw("Psame");
  
  fun_recons_e_res_2_crystal->SetLineColor(4);
  fun_recons_e_res_2_crystal->SetLineStyle(4);
  Eresolution_crystal->Fit("E_rfit_quardratic_sum_2_crystal", "", "R", 0., 19.);
  fun_recons_e_res_2_crystal->Draw("same");
  
  fun_recons_e_res_3_crystal->SetLineColor(4);
  fun_recons_e_res_3_crystal->SetLineStyle(1);
  Eresolution_crystal->Fit("E_rfit_quardratic_sum_3_crystal", "", "R", 0., 19.);
  fun_recons_e_res_3_crystal->Draw("same");
  
  cout << fun_recons_e_res_2_crystal->GetChisquare() << endl;
  alpha_2 = fun_recons_e_res_2_crystal->GetParameter(0);  beta_2 = fun_recons_e_res_2_crystal->GetParameter(1);
  alpha_3 = fun_recons_e_res_3_crystal->GetParameter(0);  beta_3 = fun_recons_e_res_3_crystal->GetParameter(1);  gamma_3 = fun_recons_e_res_3_crystal->GetParameter(2);

  legend_e_reso->AddEntry("E_rfit_quardratic_sum_2_crystal",
			  Form("#frac{#sigma}{E} = %.2f #oplus #frac{%.2f}{#sqrt{E}}        [crystal]", alpha_2, beta_2),
			  "l");
  legend_e_reso->AddEntry((TObject*)0, "", "");
  legend_e_reso->AddEntry("E_rfit_quardratic_sum_3_crystal",
			  Form("#frac{#sigma}{E} = %.2f #oplus #frac{%.2f}{#sqrt{E}} #oplus #frac{%.2f}{E}", alpha_3, beta_3, gamma_3),
			  "l");
  legend_e_reso->AddEntry((TObject*)0, "", "");
  legend_e_reso->AddEntry((TObject*)0, "", "");
  legend_e_reso->AddEntry(Eresolution_glass, "sci-glass", "p");
  legend_e_reso->AddEntry(Eresolution_crystal, "crystal", "p");
  legend_e_reso->Draw("same");
  */
  






  
  // ================================
  // Draw the pion rejection
  // ================================
  
  TGraphErrors *Pion_reject_glass_3pars = (TGraphErrors*)f->Get("Pion_reject_glass_3pars");
  TGraphErrors *Pion_reject_glass_2pars = (TGraphErrors*)f->Get("Pion_reject_glass_2pars");
  TGraphErrors *Pion_reject_crystal_3pars = (TGraphErrors*)f->Get("Pion_reject_crystal_3pars");
  TGraphErrors *Pion_reject_crystal_2pars = (TGraphErrors*)f->Get("Pion_reject_crystal_2pars");  

  
  TCanvas *c1 = new TCanvas("c1", "c1", 900, 900);
  TLegend *legend_2_par;
  legend_2_par = new TLegend(0.5, 0.3, 0.75, 0.5);
  legend_2_par->SetBorderSize(0);
  gPad->SetLogy();
  gPad->SetLogx();
  Pion_reject_glass_2pars->SetTitle("Pion rejection factor inner[crystal] and outer[sci-glass]");
  Pion_reject_glass_2pars->SetLineColor(2);
  Pion_reject_glass_2pars->SetLineStyle(4);
  Pion_reject_glass_2pars->SetMarkerStyle(24);
  Pion_reject_glass_2pars->SetMarkerSize(1);
  Pion_reject_glass_2pars->SetMarkerColor(2);
  Pion_reject_glass_2pars->Draw("");
  Pion_reject_crystal_2pars->SetLineColor(4);
  Pion_reject_crystal_2pars->SetLineStyle(4);
  Pion_reject_crystal_2pars->SetMarkerStyle(25);
  Pion_reject_crystal_2pars->SetMarkerSize(1);
  Pion_reject_crystal_2pars->SetMarkerColor(4);
  Pion_reject_crystal_2pars->Draw("same");
  legend_2_par->AddEntry(Pion_reject_glass_2pars, "glass  [ #frac{#sigma}{E} = #alpha #oplus #frac{#beta}{#sqrt{E}} ]", "pl");
  legend_2_par->AddEntry((TObject*)0, "", "");
  legend_2_par->AddEntry(Pion_reject_crystal_2pars, "crystal  [ #frac{#sigma}{E} = #alpha #oplus #frac{#beta}{#sqrt{E}} ]", "pl");
  //  legend_2_par->AddEntry(Pion_reject_glass_2pars, "glass  #frac{#sigma}{E} fit 2 pars", "pl");
  //  legend_2_par->AddEntry(Pion_reject_crystal_2pars, "crystal  #frac{#sigma}{E} fit 2 pars", "pl");
  legend_2_par->AddEntry((TObject*)0, "", "");
  legend_2_par->AddEntry((TObject*)0, "", "");
  //  legend_2_par->Draw("same");
  
  
  //  TCanvas *c2 = new TCanvas("c2", "c2", 800, 800);
  TLegend *legend_3_par;
  legend_3_par = new TLegend(0.55, 0.35, 0.75, 0.5);
  legend_3_par->SetBorderSize(0);
  gPad->SetLogy();
  gPad->SetLogx();
  Pion_reject_glass_3pars->SetTitle("Pion rejection factor inner[crystal] and outer[sci-glass]");
  Pion_reject_glass_3pars->SetLineColor(2);
  Pion_reject_glass_3pars->SetLineStyle(1);
  Pion_reject_glass_3pars->SetMarkerStyle(24);
  Pion_reject_glass_3pars->SetMarkerSize(1);
  Pion_reject_glass_3pars->SetMarkerColor(2);
  Pion_reject_glass_3pars->Draw("same");
  Pion_reject_crystal_3pars->SetLineColor(4);
  Pion_reject_crystal_3pars->SetLineStyle(1);
  Pion_reject_crystal_3pars->SetMarkerStyle(25);
  Pion_reject_crystal_3pars->SetMarkerSize(1);
  Pion_reject_crystal_3pars->SetMarkerColor(4);
  Pion_reject_crystal_3pars->Draw("same");
  // legend_2_par->AddEntry(Pion_reject_glass_3pars, "glass  #frac{#sigma}{E} fit 3 pars", "pl");
  // legend_2_par->AddEntry(Pion_reject_crystal_3pars, "crystal  #frac{#sigma}{E} fit 3 pars", "pl");
  legend_2_par->AddEntry(Pion_reject_glass_3pars, "glass  [ #frac{#sigma}{E} = #alpha #oplus #frac{#beta}{#sqrt{E}} #oplus #frac{#gamma}{E}]", "pl");
  legend_2_par->AddEntry((TObject*)0, "", "");
  legend_2_par->AddEntry(Pion_reject_crystal_3pars, "crystal  [ #frac{#sigma}{E} = #alpha #oplus #frac{#beta}{#sqrt{E}} #oplus #frac{#gamma}{E}]", "pl");
  legend_2_par->Draw("same");
  


  
}
