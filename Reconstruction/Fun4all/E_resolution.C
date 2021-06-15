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

  const int N_E_bin = 10;
  string File_num[N_E_bin] = {"0.5", "1", "2", "3", "4", "5", "7", "10", "13", "18"};
  double e_per_bin[N_E_bin] = {0.5, 1., 2., 3., 4., 5., 7., 10., 13., 18.};
  double e_statement[11] = {0.5, 1., 2., 3., 4., 5., 7., 10., 13., 18., 20.};
  int flag = 0; // 0: glass, 1: crystal
  
  TH1F *h1_e_diff[N_E_bin];
  for(int i = 0 ; i < N_E_bin ; i++)
    {
      if( i < 2 )
	h1_e_diff[i] = new TH1F(Form("h1_e_diff_%d", i), Form("h1_e_diff_%d", i), 80, 0., 0.4);
      else if( (i >= 2) && (i < 5) )
	h1_e_diff[i] = new TH1F(Form("h1_e_diff_%d", i), Form("h1_e_diff_%d", i), 100, 0., 1.);
      else if( (i >= 5) && (i < 7) )
	h1_e_diff[i] = new TH1F(Form("h1_e_diff_%d", i), Form("h1_e_diff_%d", i), 200, 0., 2.);
      else if( (i >= 7) && (i < N_E_bin) )
      	h1_e_diff[i] = new TH1F(Form("h1_e_diff_%d", i), Form("h1_e_diff_%d", i), 80, 0.5, 4.5);
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
  double fit_E_resolu[N_E_bin] = {};
  double x_err[N_E_bin] = {0.}, y_err[N_E_bin] = {};
  
  TCanvas *c3 = new TCanvas("c3", "c3", 1800, 720);
  c3->Divide(5,2);
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
      fun_e_res->SetParLimits(1, 0.05, 2.5);
      fun_e_res->SetParLimits(2, 0.01, 0.5);
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
      legend_e_bin->AddEntry((TObject*)0, Form(" #sigma error:  %.5f", (fun_e_res->GetParError(2)) ), "");
      legend_e_bin->AddEntry((TObject*)0, Form(" #chi^{2}:  %.3f", (fun_e_res->GetChisquare() / fun_e_res->GetNDF()) ), "");
      legend_e_bin->Draw("same");
      //      cout << "Chisquare: " << (fun_e_res->GetChisquare() / fun_e_res->GetNDF()) << endl;

      double e_resolu = (fun_e_res->GetParameter(2)) / e_per_bin[i] * 100.;  // x100% 

      fit_E_mean[i] = fun_e_res->GetParameter(1);
      fit_E_resolu[i] = e_resolu;
      y_err[i] = (fun_e_res->GetParError(2)) * 100.;
      x_err[i] = 0.;
    }







  TFile *g = new TFile("E_Reso_result.root", "UPDATE");
  
  TGraphErrors* E_reso = new TGraphErrors(N_E_bin, e_per_bin, fit_E_resolu, x_err, y_err);
  TGraph* E_shift = new TGraph(N_E_bin, e_per_bin, fit_E_mean);
  E_shift->GetHistogram()->SetMinimum(0.);
  E_reso->GetHistogram()->SetMinimum(0.);
  E_reso->GetHistogram()->SetMaximum(7.);
    
  TCanvas *c4 = new TCanvas("c4", "c4", 800, 800);
  // c4->Divide(2,1);
  // c4->cd(1);
  TLegend *legend_e_reso;
  legend_e_reso = new TLegend(0.3, 0.55, 0.8, 0.85);
  legend_e_reso->SetBorderSize(0);
  
  auto *fun_recons_e_res_2 = new TF1("E_resolu_fit_quardratic_sum_2", E_resolu_fit_quardratic_sum_2, 0., 19., 2);
  fun_recons_e_res_2->SetParLimits(0, 0.1, 10.);
  fun_recons_e_res_2->SetParLimits(1, 0.01, 10.);
  //  fun_recons_e_res_2->SetLineWidth(10);
  fun_recons_e_res_2->SetLineColor(2);
  
  auto *fun_recons_e_res_3 = new TF1("E_resolu_fit_quardratic_sum_3", E_resolu_fit_quardratic_sum_3, 0., 19., 3);
  fun_recons_e_res_3->SetParLimits(0, 0.1, 10.);
  fun_recons_e_res_3->SetParLimits(1, 0.01, 10.);
  fun_recons_e_res_3->SetParLimits(2, 0.01, 10.);
  fun_recons_e_res_3->SetLineColor(4);
  
  E_reso->SetMarkerStyle(21);
  E_reso->SetMarkerSize(1);
  if( flag == 1 )
    {
      E_reso->SetTitle("crystal pe/GeV = 15000, noise = 0, 0");
      E_reso->SetName("Eresolution_crystal");
    }
  else
    {
      E_reso->SetTitle("glass pe/GeV = 5000, noise = 0, 0");
      E_reso->SetName("Eresolution_glass");
    }

  E_reso->GetXaxis()->SetTitle("E [GeV]");
  E_reso->GetYaxis()->SetTitle("E_sig / E [%]");
  E_reso->Draw("AP");
  E_reso->Fit("E_resolu_fit_quardratic_sum_2", "", "R", 0., 19.);
  fun_recons_e_res_2->Draw("same");
  E_reso->Fit("E_resolu_fit_quardratic_sum_3", "", "R", 0., 19.);
  fun_recons_e_res_3->Draw("same");
  //  cout << fun_recons_e_res_2->GetChisquare() / fun_recons_e_res_2->GetNDF() << endl;
  cout << fun_recons_e_res_2->GetChisquare() << endl;

  double alpha_2 = fun_recons_e_res_2->GetParameter(0),
         beta_2 = fun_recons_e_res_2->GetParameter(1);
  double alpha_3 = fun_recons_e_res_3->GetParameter(0),
         beta_3 = fun_recons_e_res_3->GetParameter(1),
         gamma_3 = fun_recons_e_res_3->GetParameter(2);
  
  //  legend_e_reso->AddEntry((TObject*)0, Form(" #alpha:  %.3f", (fun_recons_e_res_2->GetParameter(0)) ), "");
  //  legend_e_reso->AddEntry((TObject*)0, Form(" #beta:  %.3f", (fun_recons_e_res_2->GetParameter(1)) ), "");
  legend_e_reso->AddEntry((TObject*)0, "", "");
  legend_e_reso->AddEntry("E_resolu_fit_quardratic_sum_2",
			  Form("#frac{#sigma}{E} = %.2f #oplus #frac{%.2f}{#sqrt{E}}", alpha_2, beta_2),
			  "l");
  // legend_e_reso->AddEntry((TObject*)0, Form(" #alpha:  %.3f", (fun_recons_e_res_3->GetParameter(0)) ), "");
  // legend_e_reso->AddEntry((TObject*)0, Form(" #beta:  %.3f", (fun_recons_e_res_3->GetParameter(1)) ), "");
  // legend_e_reso->AddEntry((TObject*)0, Form(" #gamma:  %.3f", (fun_recons_e_res_3->GetParameter(2)) ), "");
  legend_e_reso->AddEntry((TObject*)0, "", "");
  legend_e_reso->AddEntry("E_resolu_fit_quardratic_sum_3",
			  Form("#frac{#sigma}{E} = %.2f #oplus #frac{%.2f}{#sqrt{E}} #oplus #frac{%.2f}{E}", alpha_3, beta_3, gamma_3),
			  "l");
  legend_e_reso->Draw("same");
  E_reso->Write();
  g->Close();



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
