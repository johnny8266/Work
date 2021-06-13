double GausM(double *x, double *par)
{
  return par[0] * exp(-0.5 * TMath::Power(((x[0] - par[1]) / par[2]), 2));
}

double_t E_resolu_fit_quardratic_sum(double *x, double *par)
{
  return sqrt( par[0] * par[0] + (par[1] * par[1]) / x[0] + (par[2] * par[2]) / (x[0] * x[0]) );
  //  return sqrt( par[0] * par[0] + (par[1] * par[1]) / x[0] );

  //  return sqrt( (par[0] * par[0]) / x[0] + (par[1] * par[1]) / (x[0] * x[0]) );
}

double E_shift_fit(double *x, double *par)
{
  return par[0] * x[0] + par[1];
}


void read()
{

  //  TFile *f = new TFile("/vol0/pwang-l/Work/Data/fun4all_simulation/e_pion/outer/gamma_g4eemc_eval_10k_pe100_lower_noise.root");
  TFile *f = new TFile("/vol0/pwang-l/Work/Data/fun4all_simulation/e_pion/inner/gamma_g4eemc_eval_25k_pe15000_noise_80.root");

  TTree *ntp_gpoint = (TTree*)f->Get("ntp_gpoint");
  TTree *ntp_gshower = (TTree*)f->Get("ntp_gshower");
  TTree *ntp_tower = (TTree*)f->Get("ntp_tower");
  TTree *ntp_cluster = (TTree*)f->Get("ntp_cluster");

  Int_t n_gpoint = (Int_t)ntp_gpoint->GetEntries();
  Int_t n_gshower = (Int_t)ntp_gshower->GetEntries();
  Int_t n_tower = (Int_t)ntp_tower->GetEntries();
  Int_t n_cluster = (Int_t)ntp_cluster->GetEntries();

  

  //declare the variables used in the following ntuples
  //
  TLorentzVector v_g;
  Float_t ge, gpt, geta, gphi, gvx, gvy, gvz;
  Float_t event, clusterID, x, y, z, e, ntowers, e_base = 0., x_base = 0., y_base = 0., z_base = 0., n_tow_base = 0.;
  vector<float> x_reco, y_reco, z_reco, e_reco, N_towers;
  vector<int> i_th_cluster;
  int count_single_clu = 0, count = 0;

  
  /*
  //Load ntuple of ntp_gpoint
  //
  ntp_gpoint->SetBranchAddress("event",&event);
  ntp_gpoint->SetBranchAddress("gvx",&gvx);
  ntp_gpoint->SetBranchAddress("gvy",&gvy);
  ntp_gpoint->SetBranchAddress("gvz",&gvz);

  for(int i = 0 ; i < n_gpoint ; i++)
    {
      ntp_gpoint->GetEntry(i);
      cout << "ntp_gpoint: " << event << " " << gvx << " " << gvy << " " << gvz << endl;
    }
  */


  //Load ntuple of ntp_cluster
  //
  ntp_cluster->SetBranchAddress("event", &event);
  ntp_cluster->SetBranchAddress("clusterID", &clusterID);
  ntp_cluster->SetBranchAddress("x", &x);
  ntp_cluster->SetBranchAddress("y", &y);
  ntp_cluster->SetBranchAddress("z", &z);
  ntp_cluster->SetBranchAddress("e", &e);
  ntp_cluster->SetBranchAddress("ntowers", &ntowers);
  
  for(int i = 0 ; i < n_cluster ; i++)
    {
      ntp_cluster->GetEntry(i);

      int i_eve = (int)event;
      
      if( i_eve == count )
	{
	  if( e > e_base )
	    {
	      e_base = e;
	      x_base = x;
	      y_base = y;
	      z_base = z;
	      n_tow_base = ntowers;
	      //	      cout << i_eve << " " << count << " " << e << endl;
	    }
	}
      else
	{
	  //	  cout << i_eve << endl;
	  i_th_cluster.push_back(count);
	  x_reco.push_back(x_base);
	  y_reco.push_back(y_base);
	  z_reco.push_back(z_base);
	  e_reco.push_back(e_base);
	  N_towers.push_back(n_tow_base);
	  count++;
	  count_single_clu++;
	  i = i - 1;

	  e_base = 0.;  x_base = 0.;  y_base = 0.; z_base = 0.; n_tow_base = 0.;
	}
    }
  cout << count_single_clu << endl;
  


  //Load ntuple of ntp_gshower
  //
  ntp_gshower->SetBranchAddress("event",&event);
  ntp_gshower->SetBranchAddress("gvx",&gvx);
  ntp_gshower->SetBranchAddress("gvy",&gvy);
  ntp_gshower->SetBranchAddress("gvz",&gvz);
  ntp_gshower->SetBranchAddress("gpt",&gpt);
  ntp_gshower->SetBranchAddress("geta",&geta);
  ntp_gshower->SetBranchAddress("gphi",&gphi);
  ntp_gshower->SetBranchAddress("ge",&ge);

  auto h1_x_diff = new TH1F("h1_x_diff", "h1_x_diff", 60, -3., 3.);
  auto h1_y_diff = new TH1F("h1_y_diff", "h1_y_diff", 60, -3., 3.);
  TH1F *h1_e_diff[9];
  for(int i = 0 ; i < 9 ; i++)
    {
      if( i == 0 )
	h1_e_diff[i] = new TH1F(Form("h1_e_diff_%d", i), Form("h1_e_diff_%d", i), 160, -0.5, 3.5);
      else
	h1_e_diff[i] = new TH1F(Form("h1_e_diff_%d", i), Form("h1_e_diff_%d", i), 160, -0.5, 3.5);
    }


  auto E_reso = new TH1F("E_reso", "E_reso", 11, 0., 11.);
  auto E_shift = new TH1F("E_shift", "E_shift", 11, 0., 11.);

  auto h1_pjtx_pjty = new TH2F("h1_pjtx_pjty", "h1_pjtx_pjty", 60, -150., 150., 60, -150., 150.);
  auto h1_x_diff_pjtx = new TH2F("h1_x_diff_pjtx", "h1_x_diff_pjtx", 60, -3., 3., 45, -135., 135.);
  auto h1_y_diff_pjty = new TH2F("h1_y_diff_pjty", "h1_y_diff_pjty", 60, -3., 3., 45, -135., 135.);
  auto h1_e_diff_pjtx = new TH2F("h1_e_diff_pjtx", "h1_e_diff_pjtx", 50, -1., 4., 45, -90., 90.);
  auto h1_e_diff_pjty = new TH2F("h1_e_diff_pjty", "h1_e_diff_pjty", 50, -1., 4., 45, -90., 90.);

  auto btwn_inner_outer = new TH3F("btwn_inner_outer", "btwn_inner_outer", 140, -140, 140, 140, -140, 140, 30, -190, -220);
  
  
  for(int i = 0 ; i < count_single_clu ; i++)
    {
      int i_eve = i_th_cluster[i];
      ntp_gshower->GetEntry(i_eve);
      
      v_g.SetPtEtaPhiE(gpt, geta, gphi, ge);
      //      cout << gpt << " " << geta << " " << gphi << " " << ge << endl;

      
      double gx = v_g.Px(), gy = v_g.Py(), gz = v_g.Pz(), g_E = v_g.E();
      double x_pjt = gx / gz * (z_reco[i]), y_pjt = gy / gz * (z_reco[i]);
      double x_D = x_pjt - x_reco[i], y_D = y_pjt - y_reco[i];
      double dis_reco = std::sqrt(x_D * x_D + y_D * y_D);
      //     cout << gx << " " << gy << " " << gz << " " << g_E << endl << endl << endl;

      //      h1_pjtx_pjty->Fill(x_pjt, y_pjt);
      h1_pjtx_pjty->Fill(x_reco[i], y_reco[i]);
      
      //      if(e_reco[i] > 8.)
      if(N_towers[i] > 4.)
	{
	  if( e_reco[i] > (ge * 0.1) )
	    {
	      h1_x_diff->Fill((x_pjt - x_reco[i]));
	      h1_y_diff->Fill((y_pjt - y_reco[i]));

	      h1_x_diff_pjtx->Fill((x_pjt - x_reco[i]), x_pjt);
	      h1_y_diff_pjty->Fill((y_pjt - y_reco[i]), y_pjt);

	      btwn_inner_outer->Fill(x_reco[i], y_reco[i], z_reco[i]);
	    }

	  if( dis_reco < 2.)
	    {
	      int e_bin = ge - 1;

	      
	      if(e_bin == 9)
		{
		  cout << e_bin << endl;	  
		  e_bin = 8;
		}

	      h1_e_diff[e_bin]->Fill((ge - e_reco[i]));
	      //      if(e_bin == 8)
		//		cout << ge << ",  " << e_reco[i] << endl;
	    }
	}
     

      h1_e_diff_pjtx->Fill((ge - e_reco[i]), x_pjt);
      h1_e_diff_pjty->Fill((ge - e_reco[i]), y_pjt);
    }

  
  //  h1_pjtx_pjty->Draw("colorz");
  //  TCanvas *c1 = new TCanvas("c1", "c1", 900, 900);  
  //  btwn_inner_outer->Draw("box");

  
  // Star Draw the plots and fit the results
  
  TCanvas *c1 = new TCanvas("c1", "c1", 900, 900);
  c1->Divide(2,2);
  TLegend *legend_posx, *legend_posy;
  legend_posx = new TLegend(0.55, 0.6, 0.9, 0.8);
  legend_posy = new TLegend(0.55, 0.6, 0.9, 0.8);
  legend_posx->SetBorderSize(0);
  legend_posy->SetBorderSize(0);
  auto *fun_pos_res = new TF1("fun_pos_res", GausM, -2.5, 2.5, 4);
  fun_pos_res->SetParLimits(0, 3000, 5000);
  fun_pos_res->SetParLimits(1, -0.5, 0.5);
  fun_pos_res->SetParLimits(2, 0.1, 1.);
  fun_pos_res->SetParLimits(3, -1000., 50.);
  
  c1->cd(1);
  //  h1_x_diff->Fit("fun_pos_res", "", "R", -2.5, 2.5);
  h1_x_diff->SetStats(0);
  h1_x_diff->SetTitle("sci-glass outer [4, 4, 40] X resco - 2X0");
  h1_x_diff->GetXaxis()->SetTitle("projection - reconstruction[cm]");
  h1_x_diff->Draw();
  //  legend_posx->AddEntry((TObject*)0, Form("Posx resolution:  %.3f cm", (fun_pos_res->GetParameter(2)) ), "");
  legend_posx->Draw("same");
  
  c1->cd(2);
  //  h1_y_diff->Fit("fun_pos_res", "", "R", -2.5, 2.5);
  h1_y_diff->SetStats(0);
  h1_y_diff->SetTitle("sci-glass outer [4, 4, 40] Y resco - 2X0");
  h1_y_diff->GetXaxis()->SetTitle("projection - reconstruction[cm]");
  h1_y_diff->Draw();
  //  legend_posy->AddEntry((TObject*)0, Form("Posx resolution:  %.3f cm", (fun_pos_res->GetParameter(2)) ), "");
  legend_posy->Draw("same");


  //  TCanvas *c2 = new TCanvas("c2", "c2", 1200, 600);
  //  c2->Divide(2,1);
  c1->cd(3);
  h1_x_diff_pjtx->SetTitle("pjt_x  v.s.  difference [pjt_x, reco_x]");
  h1_x_diff_pjtx->GetXaxis()->SetTitle("projection - reconstruction[cm]");
  h1_x_diff_pjtx->SetStats(0);
  h1_x_diff_pjtx->Draw("colorz");
  c1->cd(4);
  h1_y_diff_pjty->SetTitle("pjt_y  v.s.  difference [pjt_y, reco_y]");
  h1_y_diff_pjty->GetXaxis()->SetTitle("projection - reconstruction[cm]");
  h1_y_diff_pjty->SetStats(0);
  h1_y_diff_pjty->Draw("colorz");


    
  TCanvas *c3 = new TCanvas("c3", "c3", 900, 900);
  c3->Divide(3,3);
  auto *fun_e_res = new TF1("fun_e_res", GausM, -0.5, 2., 3);
  for(int i = 0 ; i < 9 ; i++)
    {
      c3->cd(i+1);

      double up_lim = h1_e_diff[i]->GetEntries();
      double plot_mean = h1_e_diff[i]->GetMean();
      double plot_std = h1_e_diff[i]->GetStdDev();
      fun_e_res->SetParLimits(0, 5., up_lim);
      fun_e_res->SetParLimits(1, 0.05, 1.5);
      fun_e_res->SetParLimits(2, 0.05, 1.5);
      //      fun_e_res->SetParLimits(3, -500., 100.);
      if( i < 1 )
        h1_e_diff[i]->Fit("fun_e_res", "", "R", ( plot_mean - 3. * plot_std), (plot_mean + 1.6 * plot_std) );
      else
        h1_e_diff[i]->Fit("fun_e_res", "", "R", ( plot_mean - 1.5 * plot_std), (plot_mean + 0.15 * plot_std) );
      //      h1_e_diff[i]->SetStats(0);

      h1_e_diff[i]->SetTitle(Form("Energy [%d ~ %d] GeV reco - primary", i+1, i+2));
      h1_e_diff[i]->GetXaxis()->SetTitle("e diff [GeV]");
      //      h1_e_diff[i]->SetStats(0);
      h1_e_diff[i]->Draw();

      cout << "Chisquare: " << fun_e_res->GetChisquare() << endl;
      cout << fun_e_res->GetNDF() << endl;
      cout << fun_e_res->GetParameter(2) << endl;
      
      double e_resolu = (fun_e_res->GetParameter(2)) / (i + 1.5) * 100.;  // x100% 
      
      //      if( (i >= 1) )
      //        E_reso->SetBinContent( i+1, e_resolu );
      //      if( i != 8 )
      E_reso->SetBinContent( i + 2, e_resolu );
      E_shift->SetBinContent( i + 2, (fun_e_res->GetParameter(1)) );
      //      E_reso->SetBinError( i + 2, (fun_e_res->GetParError(2)) );

      
    }

  
  TCanvas *c4 = new TCanvas("c4", "c4", 800, 800);
  //  c4->Divide(2,1);
  //  c4->cd(1);
  TLegend *legend_e_reso;
  legend_e_reso = new TLegend(0.45, 0.2, 0.9, 0.4);
  legend_e_reso->SetBorderSize(0);
  auto *fun_recons_e_res = new TF1("E_resolu_fit_quardratic_sum", E_resolu_fit_quardratic_sum, 1.5, 9.5, 3);
  fun_recons_e_res->SetParLimits(0, 0.1, 10.);
  fun_recons_e_res->SetParLimits(1, 0.01, 10.);
  fun_recons_e_res->SetParLimits(2, 0.01, 30.);
  E_reso->SetStats(0);
  E_reso->SetMarkerStyle(2);
  E_reso->SetMarkerSize(1);
  E_reso->SetTitle("glass pe/GeV = 15000, noise = 8, 16");
  E_reso->GetXaxis()->SetTitle("E [GeV]");
  E_reso->GetYaxis()->SetTitle("E_sig / E [%]");
  E_reso->Fit("E_resolu_fit_quardratic_sum", "", "R", 1.3, 9.5);
  E_reso->Draw("p");
  cout << fun_recons_e_res->GetChisquare() / fun_recons_e_res->GetNDF() << endl;
  legend_e_reso->AddEntry((TObject*)0, Form(" #alpha:  %.3f", (fun_recons_e_res->GetParameter(0)) ), "");
  legend_e_reso->AddEntry((TObject*)0, Form(" #beta:  %.3f", (fun_recons_e_res->GetParameter(1)) ), "");
  legend_e_reso->AddEntry((TObject*)0, Form(" #gamma:  %.3f", (fun_recons_e_res->GetParameter(2)) ), "");
  legend_e_reso->AddEntry((TObject*)0, "", "");
  //  legend_e_reso->AddEntry((TObject*)0, "#frac{#sigma}{E} = #alpha #oplus #frac{#beta}{#sqrt{E}}", "");
  legend_e_reso->AddEntry((TObject*)0, "#frac{#sigma}{E} = #alpha #oplus #frac{#beta}{#sqrt{E}} #oplus #frac{#gamma}{E}", "");
  legend_e_reso->Draw("same");

  /*			  
  c4->cd(2);
  auto *fun_E_shift = new TF1("E_shift_fit", E_shift_fit, 1., 10., 2);
  fun_E_shift->SetParLimits(0, 0.01, 10.);
  fun_E_shift->SetParLimits(1, 0.01, 10.);
  E_shift->Fit("E_shift_fit", "", "R", 1., 10.);
  E_shift->SetStats(0);
  E_shift->SetMarkerStyle(3);
  E_shift->SetMarkerSize(3);
  E_shift->Draw("P");
  */
  
  
  /*
  legend[0] = new TLegend(0.5, 0.6, 0.75, 0.8);
  legend[0]->SetBorderSize(0);
  legend[0]->AddEntry((TObject*)0, Form("Par A:  %.3f", (fun_recons_e_res->GetParameter(0)) ), "");
  legend[0]->AddEntry((TObject*)0, Form("Par B:  %.3f", (fun_recons_e_res->GetParameter(1)) ), "");
  legend[0]->AddEntry((TObject*)0, Form("Par C:  %.3f", (fun_recons_e_res->GetParameter(2)) ), "");
  legend[0]->Draw("same");
  */

  
}
