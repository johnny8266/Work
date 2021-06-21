double GausM(double *x, double *par)
{
  return par[0] * exp(-0.5 * TMath::Power(((x[0] - par[1]) / par[2]), 2));
}

double_t E_resolu_fit_quardratic_sum_3(double *x, double *par)
{
  return sqrt( par[0] * par[0] + (par[1] * par[1]) / x[0] + (par[2] * par[2]) / (x[0] * x[0]) );
}

double_t E_resolu_fit_quardratic_sum_2(double *x, double *par)
{
  return sqrt( par[0] * par[0] + (par[1] * par[1]) / x[0] );
}

double E_shift_fit(double *x, double *par)
{
  return par[0] * x[0] + par[1];
}



void E_resolution()
{

  const int N_E_bin = 12;
  string File_num[N_E_bin] = {"0.5", "1", "2", "3", "4", "5", "7", "10", "13", "18", "30", "50"};
  double e_per_bin[N_E_bin] = {0.5, 1., 2., 3., 4., 5., 7., 10., 13., 18., 30., 50.};
  double e_statement[13] = {0.5, 1., 2., 3., 4., 5., 7., 10., 13., 18., 30., 50., 52.};
  int flag = 1; // 0: glass, 1: crystal
  
  TH1F *h1_e_diff[N_E_bin];
  for(int i = 0 ; i < N_E_bin ; i++)
    {
      /*
      if( i < 2 )
	h1_e_diff[i] = new TH1F(Form("h1_e_diff_%d", i), Form("h1_e_diff_%d", i), 80, 0., 0.4);
      else if( (i >= 2) && (i < 5) )
	h1_e_diff[i] = new TH1F(Form("h1_e_diff_%d", i), Form("h1_e_diff_%d", i), 100, 0., 1.);
      else if( (i >= 5) && (i < 7) )
	h1_e_diff[i] = new TH1F(Form("h1_e_diff_%d", i), Form("h1_e_diff_%d", i), 200, 0., 2.);
      else if( (i >= 7) && (i < 10) )
      	h1_e_diff[i] = new TH1F(Form("h1_e_diff_%d", i), Form("h1_e_diff_%d", i), 80, 0.5, 4.5);
      else if( (i >= 10) && (i < 11) )
	h1_e_diff[i] = new TH1F(Form("h1_e_diff_%d", i), Form("h1_e_diff_%d", i), 120, 1., 7.);
      else if( (i >= 11) && (i < 12) )
	h1_e_diff[i] = new TH1F(Form("h1_e_diff_%d", i), Form("h1_e_diff_%d", i), 240, 3., 15.);
      */
      if( i < 2 )
	h1_e_diff[i] = new TH1F(Form("h1_e_diff_%d", i), Form("h1_e_diff_%d", i), 100, 0., 1.);
      else if( (i >= 2) && (i < 5) )
	h1_e_diff[i] = new TH1F(Form("h1_e_diff_%d", i), Form("h1_e_diff_%d", i), 400, 0., 2.);
      else if( (i >= 5) && (i < 7) )
	h1_e_diff[i] = new TH1F(Form("h1_e_diff_%d", i), Form("h1_e_diff_%d", i), 400, 0., 4.);
      else if( (i >= 7) && (i < 10) )
      	h1_e_diff[i] = new TH1F(Form("h1_e_diff_%d", i), Form("h1_e_diff_%d", i), 1000, 0., 10.);
      else if( (i >= 10) && (i < 11) )
	h1_e_diff[i] = new TH1F(Form("h1_e_diff_%d", i), Form("h1_e_diff_%d", i), 1400, 1., 15.);
      else if( (i >= 11) && (i < 12) )
	h1_e_diff[i] = new TH1F(Form("h1_e_diff_%d", i), Form("h1_e_diff_%d", i), 2200, 3., 25.);
    }


  
  for( int i = 0 ; i < N_E_bin ; i++ )
    {
      string root_file_name;

      if( flag == 1 )
	{
	  if( i == 0)
	    root_file_name = "/vol0/pwang-l/Singularity/my_det/sub_crystal/data/g4eemc_crystal_eval_mono_" + File_num[i] + "_GeV.root";
	  else
	    root_file_name = "/vol0/pwang-l/Singularity/my_det/sub_crystal/data/g4eemc_crystal_eval_mono_" + File_num[i] + ".0_GeV.root";
	}
      else
	{
	  if( i == 0)
	    root_file_name = "/vol0/pwang-l/Singularity/my_det/sub_glass/data/g4eemc_glass_eval_mono_" + File_num[i] + "_GeV.root";
	  else
	    root_file_name = "/vol0/pwang-l/Singularity/my_det/sub_glass/data/g4eemc_glass_eval_mono_" + File_num[i] + ".0_GeV.root";	  
	}

      

      const char *rfn = root_file_name.c_str();
      TFile *f = new TFile(rfn);

      TTree *ntp_cluster = (TTree*)f->Get("ntp_cluster");
      Int_t n_cluster = (Int_t)ntp_cluster->GetEntries();

      //declare the variables used in the following ntuples
      //
      TLorentzVector v_g;
      Float_t ge, gpt, geta, gphi, gvx, gvy, gvz;
      Float_t event, clusterID, x, y, z, e, ntowers, e_base = 0., x_base = 0., y_base = 0., z_base = 0., n_tow_base = 0.;
      vector<float> x_reco, y_reco, z_reco, e_reco, N_towers;
      vector<int> i_th_cluster;
      int count_single_clu = 0, count = 0;

 
      //Load ntuple of ntp_cluster
      //
      ntp_cluster->SetBranchAddress("event", &event);
      ntp_cluster->SetBranchAddress("clusterID", &clusterID);
      ntp_cluster->SetBranchAddress("x", &x);
      ntp_cluster->SetBranchAddress("y", &y);
      ntp_cluster->SetBranchAddress("z", &z);
      ntp_cluster->SetBranchAddress("e", &e);
      ntp_cluster->SetBranchAddress("ntowers", &ntowers);
      ntp_cluster->SetBranchAddress("gpt", &gpt);
      ntp_cluster->SetBranchAddress("geta", &geta);
      ntp_cluster->SetBranchAddress("gphi", &gphi);
      ntp_cluster->SetBranchAddress("ge", &ge);

      
      for(int j = 0 ; j < n_cluster ; j++)
	{
	  ntp_cluster->GetEntry(j);

	  int j_eve = (int)event;
      
	  if( j_eve == count )
	    {
	      if( e > e_base )
		{
		  e_base = e;
		  x_base = x;
		  y_base = y;
		  z_base = z;
		  n_tow_base = ntowers;
		  //	      cout << j_eve << " " << count << " " << e << endl;
		}
	    }
	  else
	    {
	      double e_diff = ge - e_base;
	      for(int k = 0 ; k < N_E_bin ; k++)
		{
		  if( (ge > (e_statement[k] - 0.1)) && ( ge < e_statement[k+1] ) )
		    h1_e_diff[k]->Fill(e_diff);
		}

	      count++;
	      count_single_clu++;
	      j = j - 1;

	      e_base = 0.;  x_base = 0.;  y_base = 0.; z_base = 0.; n_tow_base = 0.;
	    }
	}
      cout << count_single_clu << endl;
      
    }

 
  

  double fit_E_mean[N_E_bin] = {};
  double fit_E_sigma[N_E_bin] = {}, fit_E_RMS[N_E_bin] = {}, fit_E_FWHM[N_E_bin] = {};;
  double x_err_sigma[N_E_bin] = {0.}, y_err_sigma[N_E_bin] = {0.};
  double x_err_RMS[N_E_bin] = {0.}, y_err_RMS[N_E_bin] = {0.};
  double x_err_FWHM[N_E_bin] = {0.}, y_err_FWHM[N_E_bin] = {0.};
  
  TCanvas *c3 = new TCanvas("c3", "c3", 1280, 960);
  c3->Divide(4,3);
  auto *fun_e_res = new TF1("fun_e_res", GausM, -0.5, 2., 3);
  for(int i = 0 ; i < N_E_bin ; i++)
    {
      c3->cd(i+1);

      TLegend *legend_e_bin;
      legend_e_bin = new TLegend(0.55, 0.6, 0.85, 0.85);
      legend_e_bin->SetBorderSize(0);
      
      double up_lim = h1_e_diff[i]->GetEntries();
      double plot_mean = h1_e_diff[i]->GetMean();
      double plot_std = h1_e_diff[i]->GetStdDev();
      fun_e_res->SetParLimits(0, 5., up_lim);
      fun_e_res->SetParLimits(1, 0.05, 7.5);
      fun_e_res->SetParLimits(2, 0.01, 3.);
      //      fun_e_res->SetParLimits(3, -500., 100.);
      if( i < 3 )
        h1_e_diff[i]->Fit("fun_e_res", "", "R", ( plot_mean - 2.6 * plot_std), (plot_mean + 0.1 * plot_std) );
      else if( (i >= 3) && (i < 5) )
        h1_e_diff[i]->Fit("fun_e_res", "", "R", ( plot_mean - 1.5 * plot_std), (plot_mean + 0. * plot_std) );
      else if( (i >= 5) && (i < 7) )
        h1_e_diff[i]->Fit("fun_e_res", "", "R", ( plot_mean - 1.4 * plot_std), (plot_mean + 0.1 * plot_std) );
      else if( (i >= 7) && (i < 9) )
	h1_e_diff[i]->Fit("fun_e_res", "", "R", ( plot_mean - 1.3 * plot_std), (plot_mean + 0. * plot_std) );
      else if( (i >= 9) && (i < N_E_bin) )
	h1_e_diff[i]->Fit("fun_e_res", "", "R", ( plot_mean - 1.6 * plot_std), (plot_mean + 0. * plot_std) );

      h1_e_diff[i]->SetTitle(Form("Mono_E %.1f GeV reco - primary", e_statement[i]));
      h1_e_diff[i]->SetStats(0);
      h1_e_diff[i]->GetXaxis()->SetTitle("e diff [GeV]");
      h1_e_diff[i]->Draw();

      legend_e_bin->AddEntry((TObject*)0, Form(" #sigma:  %.3f", (fun_e_res->GetParameter(2)) ), "");
      //      legend_e_bin->AddEntry((TObject*)0, Form(" #sigma error:  %.5f", (fun_e_res->GetParError(2)) ), "");
      //      legend_e_bin->AddEntry((TObject*)0, Form(" #chi^{2}:  %.3f", (fun_e_res->GetChisquare() / fun_e_res->GetNDF()) ), "");
      //      legend_e_bin->Draw("same");

      double e_resolu = (fun_e_res->GetParameter(2)) / e_per_bin[i] * 100.;  // x100% 
      fit_E_mean[i] = fun_e_res->GetParameter(1);
      fit_E_sigma[i] = e_resolu;
      y_err_sigma[i] = (fun_e_res->GetParError(2)) * 100.;
      x_err_sigma[i] = 0.;
      
      fit_E_RMS[i] = (h1_e_diff[i]->GetRMS()) / e_per_bin[i] * 100.;

      int bin1 = h1_e_diff[i]->FindFirstBinAbove(h1_e_diff[i]->GetMaximum()/2);
      int bin2 = h1_e_diff[i]->FindLastBinAbove(h1_e_diff[i]->GetMaximum()/2);
      double fwhm_E = h1_e_diff[i]->GetBinCenter(bin2) - h1_e_diff[i]->GetBinCenter(bin1);
      fit_E_FWHM[i] = fwhm_E / e_per_bin[i] * 100.;

      legend_e_bin->AddEntry((TObject*)0, Form(" RMS:  %.3f", (h1_e_diff[i]->GetRMS()) ), "");
      legend_e_bin->AddEntry((TObject*)0, Form(" FWHM:  %.3f", fwhm_E ), "");
      legend_e_bin->Draw("same");
      //      cout << (fun_e_res->GetParameter(2)) << " " << (h1_e_diff[i]->GetRMS()) << " " << fwhm_E << endl;

    }


  TCanvas *c4 = new TCanvas("c4", "c4", 600, 600);
  TGraphErrors* E_reso_RMS = new TGraphErrors(N_E_bin, e_per_bin, fit_E_RMS, x_err_RMS, y_err_RMS);
  E_reso_RMS->GetHistogram()->SetMinimum(0.);
  E_reso_RMS->GetHistogram()->SetMaximum(15.);
  E_reso_RMS->SetMarkerSize(1.);
  E_reso_RMS->SetMarkerStyle(21);
  E_reso_RMS->SetMarkerColor(1);
  E_reso_RMS->SetTitle("E resolution from RMS");
  E_reso_RMS->GetXaxis()->SetTitle("E [GeV]");
  E_reso_RMS->GetYaxis()->SetTitle("#frac{RMS}{E}");
  E_reso_RMS->Draw("AP");
  
  
  /*
  cout << endl << endl << "Start fit E resolution......." << endl << endl;

  //  TFile *g = new TFile("E_Reso_result.root", "UPDATE");
  
  TGraphErrors* E_reso_sigma = new TGraphErrors(N_E_bin, e_per_bin, fit_E_sigma, x_err_sigma, y_err_sigma);
  TGraphErrors* E_reso_RMS = new TGraphErrors(N_E_bin, e_per_bin, fit_E_RMS, x_err_RMS, y_err_RMS);
  TGraphErrors* E_reso_FWHM = new TGraphErrors(N_E_bin, e_per_bin, fit_E_FWHM, x_err_FWHM, y_err_FWHM);
  TGraph* E_shift = new TGraph(N_E_bin, e_per_bin, fit_E_mean);
  E_shift->GetHistogram()->SetMinimum(0.);
  E_reso_sigma->GetHistogram()->SetMinimum(0.);
  E_reso_sigma->GetHistogram()->SetMaximum(10.);
  E_reso_RMS->GetHistogram()->SetMinimum(0.);
  E_reso_RMS->GetHistogram()->SetMaximum(15.);
  E_reso_FWHM->GetHistogram()->SetMinimum(0.);
  E_reso_FWHM->GetHistogram()->SetMaximum(20.);
    
  TCanvas *c4 = new TCanvas("c4", "c4", 1800, 600);
  c4->Divide(3,1);
  c4->cd(1);
  TLegend *legend_e_sigma;
  legend_e_sigma = new TLegend(0.3, 0.55, 0.8, 0.85);
  legend_e_sigma->SetBorderSize(0);
  
  auto *eres_sigma_fit_2 = new TF1("eres_sigma_fit_2", E_resolu_fit_quardratic_sum_2, 0., 51., 2);
  eres_sigma_fit_2->SetParLimits(0, 0.1, 10.);
  eres_sigma_fit_2->SetParLimits(1, 0.01, 10.);
  eres_sigma_fit_2->SetLineColor(2);
  auto *eres_sigma_fit_3 = new TF1("eres_sigma_fit_3", E_resolu_fit_quardratic_sum_3, 0., 51., 3);
  eres_sigma_fit_3->SetParLimits(0, 0.1, 10.);
  eres_sigma_fit_3->SetParLimits(1, 0.01, 10.);
  eres_sigma_fit_3->SetParLimits(2, 0.01, 10.);
  eres_sigma_fit_3->SetLineColor(4);
  
  E_reso_sigma->SetMarkerStyle(21);
  E_reso_sigma->SetMarkerSize(1);
  if( flag == 1 )
    {
      E_reso_sigma->SetTitle("crystal pe/GeV = 15000, noise = 0, 0");
      E_reso_sigma->SetName("Eresolution_crystal");
    }
  else
    {
      E_reso_sigma->SetTitle("glass pe/GeV = 5000, noise = 0, 0");
      E_reso_sigma->SetName("Eresolution_glass");
    }

  E_reso_sigma->GetXaxis()->SetTitle("E [GeV]");
  E_reso_sigma->GetYaxis()->SetTitle("#frac{#sigma}{E}");
  E_reso_sigma->Draw("AP");
  E_reso_sigma->Fit("eres_sigma_fit_2", "", "R", 0., 51.);
  eres_sigma_fit_2->Draw("same");
  E_reso_sigma->Fit("eres_sigma_fit_3", "", "R", 0., 51.);
  eres_sigma_fit_3->Draw("same");
  //  cout << eres_sigma_fit_2->GetChisquare() / eres_sigma_fit_2->GetNDF() << endl;

  double alpha_2 = eres_sigma_fit_2->GetParameter(0),
         beta_2 = eres_sigma_fit_2->GetParameter(1);
  double alpha_3 = eres_sigma_fit_3->GetParameter(0),
         beta_3 = eres_sigma_fit_3->GetParameter(1),
         gamma_3 = eres_sigma_fit_3->GetParameter(2);
  
  legend_e_sigma->AddEntry((TObject*)0, "", "");
  legend_e_sigma->AddEntry("E_resolu_fit_quardratic_sum_2",
			  Form("#frac{#sigma}{E} = %.2f #oplus #frac{%.2f}{#sqrt{E}}", alpha_2, beta_2),
			  "l");
  legend_e_sigma->AddEntry((TObject*)0, "", "");
  legend_e_sigma->AddEntry("E_resolu_fit_quardratic_sum_3",
			  Form("#frac{#sigma}{E} = %.2f #oplus #frac{%.2f}{#sqrt{E}} #oplus #frac{%.2f}{E}", alpha_3, beta_3, gamma_3),
			  "l");
  legend_e_sigma->Draw("same");

  
  c4->cd(2);
  TLegend *legend_e_RMS;
  legend_e_RMS = new TLegend(0.3, 0.55, 0.8, 0.85);
  legend_e_RMS->SetBorderSize(0);
  E_reso_RMS->SetMarkerSize(1.);
  E_reso_RMS->SetMarkerStyle(21);
  E_reso_RMS->SetMarkerColor(1);
  E_reso_RMS->SetTitle("E resolution from RMS");
  E_reso_RMS->GetXaxis()->SetTitle("E [GeV]");
  E_reso_RMS->GetYaxis()->SetTitle("#frac{RMS}{E}");
  E_reso_RMS->Draw("AP");
  auto *eres_RMS_fit_2 = new TF1("eres_RMS_fit_2", E_resolu_fit_quardratic_sum_2, 0., 51., 2);
  eres_RMS_fit_2->SetParLimits(0, 0.1, 10.);
  eres_RMS_fit_2->SetParLimits(1, 0.01, 10.);
  eres_RMS_fit_2->SetLineColor(2);  
  auto *eres_RMS_fit_3 = new TF1("eres_RMS_fit_3", E_resolu_fit_quardratic_sum_3, 0., 51., 3);
  eres_RMS_fit_3->SetParLimits(0, 0.1, 10.);
  eres_RMS_fit_3->SetParLimits(1, 0.01, 10.);
  eres_RMS_fit_3->SetParLimits(2, 0.01, 10.);
  eres_RMS_fit_3->SetLineColor(4);
  E_reso_RMS->Fit("eres_RMS_fit_2", "", "R", 0., 51.);
  eres_RMS_fit_2->Draw("same");
  E_reso_RMS->Fit("eres_RMS_fit_3", "", "R", 0., 51.);
  eres_RMS_fit_3->Draw("same");
  alpha_2 = eres_RMS_fit_2->GetParameter(0);  beta_2 = eres_RMS_fit_2->GetParameter(1);
  alpha_3 = eres_RMS_fit_3->GetParameter(0);  beta_3 = eres_RMS_fit_3->GetParameter(1);  gamma_3 = eres_RMS_fit_3->GetParameter(2);
  legend_e_RMS->AddEntry((TObject*)0, "", "");
  legend_e_RMS->AddEntry("E_resolu_fit_quardratic_sum_2",
			  Form("#frac{#sigma}{E} = %.2f #oplus #frac{%.2f}{#sqrt{E}}", alpha_2, beta_2),
			  "l");
  legend_e_RMS->AddEntry((TObject*)0, "", "");
  legend_e_RMS->AddEntry("E_resolu_fit_quardratic_sum_3",
			  Form("#frac{#sigma}{E} = %.2f #oplus #frac{%.2f}{#sqrt{E}} #oplus #frac{%.2f}{E}", alpha_3, beta_3, gamma_3),
			  "l");
  legend_e_RMS->Draw("same");

  c4->cd(3);
  TLegend *legend_e_FWHM;
  legend_e_FWHM = new TLegend(0.3, 0.55, 0.8, 0.85);
  legend_e_FWHM->SetBorderSize(0);
  E_reso_FWHM->SetMarkerSize(1.);
  E_reso_FWHM->SetMarkerStyle(21);
  E_reso_FWHM->SetMarkerColor(1);
  E_reso_FWHM->SetTitle("E resolution from FWHM");
  E_reso_FWHM->GetXaxis()->SetTitle("E [GeV]");
  E_reso_FWHM->GetYaxis()->SetTitle("#frac{FWHM}{E}");
  E_reso_FWHM->Draw("AP");
  auto *eres_FWHM_fit_2 = new TF1("eres_FWHM_fit_2", E_resolu_fit_quardratic_sum_2, 0., 51., 2);
  eres_FWHM_fit_2->SetParLimits(0, 0.1, 10.);
  eres_FWHM_fit_2->SetParLimits(1, 0.01, 10.);
  eres_FWHM_fit_2->SetLineColor(2);  
  auto *eres_FWHM_fit_3 = new TF1("eres_FWHM_fit_3", E_resolu_fit_quardratic_sum_3, 0., 51., 3);
  eres_FWHM_fit_3->SetParLimits(0, 0.1, 10.);
  eres_FWHM_fit_3->SetParLimits(1, 0.01, 10.);
  eres_FWHM_fit_3->SetParLimits(2, 0.01, 10.);
  eres_FWHM_fit_3->SetLineColor(4);
  E_reso_FWHM->Fit("eres_FWHM_fit_2", "", "R", 0., 51.);
  eres_FWHM_fit_2->Draw("same");
  E_reso_FWHM->Fit("eres_FWHM_fit_3", "", "R", 0., 51.);
  eres_FWHM_fit_3->Draw("same");
  alpha_2 = eres_FWHM_fit_2->GetParameter(0);  beta_2 = eres_FWHM_fit_2->GetParameter(1);
  alpha_3 = eres_FWHM_fit_3->GetParameter(0);  beta_3 = eres_FWHM_fit_3->GetParameter(1);  gamma_3 = eres_FWHM_fit_3->GetParameter(2);
  legend_e_FWHM->AddEntry((TObject*)0, "", "");
  legend_e_FWHM->AddEntry("E_resolu_fit_quardratic_sum_2",
			  Form("#frac{#sigma}{E} = %.2f #oplus #frac{%.2f}{#sqrt{E}}", alpha_2, beta_2),
			  "l");
  legend_e_FWHM->AddEntry((TObject*)0, "", "");
  legend_e_FWHM->AddEntry("E_resolu_fit_quardratic_sum_3",
			  Form("#frac{#sigma}{E} = %.2f #oplus #frac{%.2f}{#sqrt{E}} #oplus #frac{%.2f}{E}", alpha_3, beta_3, gamma_3),
			  "l");
  legend_e_FWHM->Draw("same");

  //  E_reso_sigma->Write();
  //  g->Close();

  */

  /* 
  c4->cd(2);
  auto *fun_E_shift = new TF1("E_shift_fit", E_shift_fit, 0., 20., 2);
  fun_E_shift->SetParLimits(0, 0.01, 10.);
  fun_E_shift->SetParLimits(1, 0.01, 10.);
  E_shift->Fit("E_shift_fit", "", "R", 0., 20.);
  //  E_shift->SetStats(0);
  E_shift->SetMarkerStyle(3);
  E_shift->SetMarkerSize(3);
  E_shift->Draw();
  */
  
  
  
}
