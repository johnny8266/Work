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



void pi0_2_gamma()
{

  const int N_E_bin = 7;
  string File_num[N_E_bin] = {"2.0", "4.0", "7.0", "10.0", "14.0", "20.0", "30.0"};
  double e_per_bin[N_E_bin] = {2.0, 4.0, 7.0, 10.0, 14.0, 20.0, 30.0};
  // string File_num[N_E_bin] = {"2.0", "4.0", "6.0", "8.0", "10.0", "14.0", "20.0"};
  // double e_per_bin[N_E_bin] = {2.0, 4.0, 6.0, 8.0, 10.0, 14.0, 20.0};
  // string File_num[N_E_bin] = {"10.0"};
  // double e_per_bin[N_E_bin] = {10.0};
  double e_statement[13] = {0.5, 1., 2., 3., 4., 5., 7., 10., 13., 18., 30., 50., 52.};
  double angle_statement[9] = {0., 0.3, 0.6, 0.9, 1.2, 1.6, 2., 3., 4.};

  double E_cut[N_E_bin] = {1.52833, 3.16879, 5.46432, 7.9867, 10.8647, 15.3859, 22.7036};
  
  int flag = 1; // 0: glass, 1: crystal
  TLorentzVector LV;
  double pri_Px = 0., pri_Py = 0., pri_Pz = 0.;
  double pjx = 0., pjy = 0., pjz = 0., diff_x = 0., diff_y = 0.;
  double x_graph[N_E_bin], y_graph[N_E_bin], x_err[N_E_bin], y_err[N_E_bin], possible_pi0_count = 0.;
  

  // Declare the plots
  //
  TH1F *h1_e_diff[N_E_bin];
  TH1F *h1_posx_diff[N_E_bin], *h1_posy_diff[N_E_bin];
  TH1F *recons_1_angle[N_E_bin];
  TH1F *recons_2_angle[N_E_bin];
  TH1F *recons_1_pos[N_E_bin];
  TH1F *recons_2_pos[N_E_bin];
  TH1F *E_reco_1[N_E_bin];
  TH1F *E_reco_2[N_E_bin];
  TH1F *total_E_reco[N_E_bin];
  TH1F *zpos_reco_1[N_E_bin];
  TH1F *zpos_reco_2[N_E_bin];
  TH1I *number_of_clus[N_E_bin];
  
  TH1F *angle_range[N_E_bin][8];
  const int N_bin_plot = 220;
  for(int i = 0 ; i < N_E_bin ; i++)
    for(int j = 0 ; j < 8 ; j++)
      angle_range[i][j] = new TH1F(Form("angle_range_%d_%d", i, j), Form("angle_range_%d_%d", i, j), N_bin_plot, 0., (e_per_bin[i]+10.0) * 1.1 );

  
  TH2F *h2_posx_diff[N_E_bin], *h2_posy_diff[N_E_bin];
  TH2F *cluster_1 = new TH2F("cluster_1", "cluster_1", 30, 0., 60., 30, 0., 60.);
  TH2F *cluster_2 = new TH2F("cluster_2", "cluster_2", 30, 0., 60., 30, 0., 60.);



  // Start read the files and fill data
  //
  for(int i = 0 ; i < N_E_bin ; i++)
    {
      recons_1_angle[i] = new TH1F(Form("recons_1_angle_%d", i), Form("recons_1_angle_%d", i), 120, 0., 6.);
      recons_1_angle[i]->SetLineColor(2);
      recons_2_angle[i] = new TH1F(Form("recons_2_angle_%d", i), Form("recons_2_angle_%d", i), 120, 0., 6.);
      recons_2_angle[i]->SetLineColor(4);
      
      E_reco_1[i] = new TH1F(Form("E_reco_1_%d", i), Form("E_reco_1_%d", i), 100, 0., 20.);
      E_reco_2[i] = new TH1F(Form("E_reco_2_%d", i), Form("E_reco_2_%d", i), 100, 0., 20.);

      //      total_E_reco[i] = new TH1F(Form("total_E_reco_%d", i), Form("total_E_reco_%d", i), 240, e_per_bin[i]*0.5, e_per_bin[i]*1.1);
      total_E_reco[i] = new TH1F(Form("total_E_reco_%d", i), Form("total_E_reco_%d", i), 220, 0., e_per_bin[i]*1.1);

      number_of_clus[i] = new TH1I(Form("number_of_clus_%d", i), Form("number_of_clus_%d", i), 5, 0, 5);
      
      zpos_reco_1[i] = new TH1F(Form("zpos_reco_1_%d", i), Form("zpos_reco_1_%d", i), 100, -200., -190.);
      zpos_reco_1[i]->GetXaxis()->SetRangeUser(-200., -190.);
      zpos_reco_2[i] = new TH1F(Form("zpos_reco_2_%d", i), Form("zpos_reco_2_%d", i), 100, -200., -190.);
      zpos_reco_2[i]->GetXaxis()->SetRangeUser(-200., -190.);
      
      h1_posx_diff[i] = new TH1F(Form("h1_posx_diff_%d", i), Form("h1_posx_diff_%d", i), 120, -3., 3.);
      h1_posy_diff[i] = new TH1F(Form("h1_posy_diff_%d", i), Form("h1_posy_diff_%d", i), 120, -3., 3.);
      
      h2_posx_diff[i] = new TH2F(Form("h2_posx_diff_%d", i), Form("h2_posx_diff_%d", i), 70, -70., 70., 70, -70., 70.);
      h2_posy_diff[i] = new TH2F(Form("h2_posy_diff_%d", i), Form("h2_posy_diff_%d", i), 100, -100., 100., 120, -3., 3.);

    }

  


  // Read the data and calculate the reconstruction
  //
  for( int i = 0 ; i < N_E_bin ; i++ )
    {
      string root_file_name;

      if( flag == 1 )
	{
	  root_file_name = "/vol0/pwang-l/Singularity/my_det/sub_crystal/data/pi0/g4eemc_crystal_eval_mono_" + File_num[i] + "_GeV.root";
	}
      else
	{
	  root_file_name = "/vol0/pwang-l/Singularity/my_det/sub_glass/data/pos/g4eemc_glass_eval_mono_" + File_num[i] + "_GeV_close_event_study.root";
	}

      

      const char *rfn = root_file_name.c_str();
      TFile *f = new TFile(rfn);

      TTree *ntp_gshower = (TTree*)f->Get("ntp_gshower");
      Int_t n_gshower = (Int_t)ntp_gshower->GetEntries();      
      TTree *ntp_cluster = (TTree*)f->Get("ntp_cluster");
      Int_t n_cluster = (Int_t)ntp_cluster->GetEntries();


      
      
      //declare the variables used in the following ntuples
      //
      Float_t ge, gpt, geta, gphi, gvx, gvy, gvz, gparticleID;
      Float_t event, clusterID, x, y, z, e, ntowers;
      Float_t total_e = 0., e_base = 0., x_base = 0., y_base = 0., z_base = 0., n_tow_base = 0.;
      Float_t e_1_1, e_2_1, e_2_2, e_3_1, e_3_2, e_3_3;
      Float_t e_4_1, e_4_2, e_4_3, e_4_4;
      Float_t total_pi0 = n_gshower, denominator = 0., numerator = 0.;
      int count_single_clu = 0, count_double_clu = 0, count = 0, ID_sum = 0, good_event = 0, good_event_1 = 0, good_event_2 = 0;
      int N_cluster = 0;
      //      cout << total_pi0 << " " << n_cluster << endl;      
      
      
      ntp_cluster->SetBranchAddress("event", &event);
      ntp_cluster->SetBranchAddress("clusterID", &clusterID);
      ntp_cluster->SetBranchAddress("x", &x);
      ntp_cluster->SetBranchAddress("y", &y);
      ntp_cluster->SetBranchAddress("z", &z);
      ntp_cluster->SetBranchAddress("e", &e);
      ntp_cluster->SetBranchAddress("ntowers", &ntowers);
      ntp_cluster->SetBranchAddress("gparticleID", &gparticleID);

      
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
		}
	      total_e = total_e + e;
	      N_cluster++;
	    }
	  else
	    {
	      //	      total_E_reco[i]->Fill(total_e);
	      
	      if( (total_e > (e_per_bin[i] * 0.4)) && (total_e < (e_per_bin[i] * 0.7)) )
		number_of_clus[i]->Fill(N_cluster);
		//		h2_posx_diff[i]->Fill(x_base, y_base);

	      float dis = std::sqrt(x_base * x_base + y_base * y_base);
	      if( (dis < 15.) || (dis > 48.) )
	       	total_e = 0.;	      
	      //		  cout << x_base << ", " << y_base << endl;
	      

	      
	      //	      if( total_e > (e_per_bin[i] * 0.7) )
	      if( total_e > E_cut[i] )
	      //	      if( total_e > 0. )
		{
		  denominator++;
		  //		  total_E_reco[i]->Fill(total_e);	       
	      	      
		  ntp_cluster->GetEntry(j - 1);
		  j_eve = (int)event;
		  int Nth_clus = clusterID;

	      
		  if( Nth_clus == 0 )
		    {
		      //		  cout << "Only 1 cluster reconstructed: " << j_eve << endl;		  
		      // recons_1_angle[i]->Fill(angle_2p[j_eve]);
		      // E_reco_1[i]->Fill(e);
		      // zpos_reco_1[i]->Fill(z);
		      count_single_clu++;
		    }
		  else if( Nth_clus == 1 )
		    {
		      double xi[2] = {}, yi[2] = {}, D_2_clu = 0.;
		      
		      ntp_cluster->GetEntry(j - 1);
		      e_2_1 = e;
		      xi[0] = x;  yi[0] = y; 
		      // E_reco_2[i]->Fill(e_2_1);
		      // zpos_reco_2[i]->Fill(z);
		  
		      ntp_cluster->GetEntry(j - 2);
		      e_2_2 = e;
		      xi[1] = x;  yi[1] = y; 
		      // E_reco_2[i]->Fill(e_2_2);
		      // zpos_reco_2[i]->Fill(z);
		      //		  total_E_reco[i]->Fill((e_2_1+e_2_2));
		      if( (e_2_1 > 0.1) && (e_2_2 > 0.1) )
			{
			  total_E_reco[i]->Fill(e_2_1);
			  total_E_reco[i]->Fill(e_2_2);
			  // if( i == 2 )
			  //   if( (e_2_1 > e_per_bin[i] * 0.5) && (e_2_1 < e_per_bin[i] * 0.6) )
			  //     cout << e_2_1 << " " << e_2_2 << " " << std::sqrt( (xi[1] - xi[0]) * (xi[1] - xi[0]) + (yi[1] - yi[0]) * (yi[1] - yi[0]) ) << endl;
			  numerator++;
			}			

		    
		      ID_sum = 0;
		      count_double_clu++;

		    }
		  else if( Nth_clus == 2 )
		    {
		      int flag = 0;
		  
		      ntp_cluster->GetEntry(j - 1);
		      e_3_1 = e;
		      if( e_3_1 > 0.1 ) flag++;
		  
		      ntp_cluster->GetEntry(j - 2);
		      e_3_2 = e;
		      if( e_3_2 > 0.1 ) flag++;
		  
		      ntp_cluster->GetEntry(j - 3);
		      e_3_3 = e;
		      if( e_3_3 > 0.1 ) flag++;

		      if( flag == 2 ) 
			numerator++;
		      // if( total_e > (e_per_bin[i] * 0.7) )
		      //   cout << e_3_1 << " " << e_3_2 << " " << e_3_3 << " " << (e_3_1 + e_3_2 + e_3_3) << " || " << total_e << " " << N_cluster << endl;
		      //		    good_event_2++;
		    }
		  else if( Nth_clus == 3 )
		    {
		      int flag = 0;
		  
		      ntp_cluster->GetEntry(j - 1);
		      e_4_1 = e;
		      if( e_4_1 > 0.1 ) flag++;
		  
		      ntp_cluster->GetEntry(j - 2);
		      e_4_2 = e;
		      if( e_4_2 > 0.1 ) flag++;
		  
		      ntp_cluster->GetEntry(j - 3);
		      e_4_3 = e;
		      if( e_4_3 > 0.1 ) flag++;

		      ntp_cluster->GetEntry(j - 4);
		      e_4_4 = e;
		      if( e_4_4 > 0.1 ) flag++;

		      if( flag == 2 ) 
			numerator++;
		      // if( total_e > (e_per_bin[i] * 0.7) )
		      //   cout << e_3_1 << " " << e_3_2 << " " << e_3_3 << " " << (e_3_1 + e_3_2 + e_3_3) << " || " << total_e << " " << N_cluster << endl;
		      //		    good_event_2++;
		    }
	      
		}
	      
	      // if( total_e > (e_per_bin[i] * 0.7) )
	      // 	good_event++;
	      
	      count++;
	      j = j - 1;

	      e_base = 0.;  x_base = 0.;  y_base = 0.; z_base = 0.; n_tow_base = 0.;
	      total_e = 0.;  N_cluster = 0;
	    }
	}// Finish the cluster ntuple
      
      x_graph[i] = e_per_bin[i];
      x_err[i] = 0.;
      y_graph[i] = numerator / denominator * 100.;	    
      //      y_graph[i] = (double)count_double_clu / total_pi0 * 100.;
      y_err[i] = y_graph[i] * std::sqrt(1. / (double)count_double_clu + 1. / total_pi0);
      
      cout << endl << endl << "Total events: " << n_gshower << endl;      
      cout << "Single cluster count: " << count_single_clu << " || good events: " << good_event_1 << endl;
      cout << "Double cluster count: " << count_double_clu << " || good events: " << good_event_2 << endl;
      cout << "numerator: " << numerator << " || denominator: " << denominator << endl;
      //      cout << "Good event: " << good_event << endl;
      cout << "Primary pi0 " << x_graph[i] << " GeV || E cut: " << E_cut[i] << " GeV." << endl;
      cout << "2 photon min angle: " << 0.135 / x_graph[i] / M_PI * 180. << endl;
      cout << "Pi0 reco effi: " << y_graph[i] << "% || error: " << y_err[i] << endl << endl;

      count_single_clu = 0;  count_double_clu = 0;
      good_event = 0;  good_event_1 = 0;  good_event_2 = 0;
      
    }// Finish the files read      



  
  TCanvas *c1 = new TCanvas("c1", "c1", 800, 800);
  TGraphErrors *Pi0_effi = new TGraphErrors(N_E_bin, x_graph, y_graph, x_err, y_err);
  Pi0_effi->GetYaxis()->SetTitle("#pi^{0} effi [%]");
  Pi0_effi->GetXaxis()->SetTitle("#pi^{0} E[GeV]");
  Pi0_effi->SetTitle("#pi^{0} reconstruct efficiency");
  Pi0_effi->GetHistogram()->SetMaximum(100.);
  Pi0_effi->GetHistogram()->SetMinimum(0.);
  Pi0_effi->SetMarkerSize(1);
  Pi0_effi->SetMarkerStyle(24);
  Pi0_effi->SetMarkerColor(1);
  Pi0_effi->Draw("AP");
  



  TCanvas *c2 = new TCanvas("c2", "c2", 1600, 800);
  c2->Divide(4,2);
  for(int i = 0 ; i < N_E_bin ; i++)
    {
      c2->cd(i+1);


      /*
      auto *fit_gauss = new TF1("fit_gauss", GausM, e_per_bin[i]*0.82, e_per_bin[i]*0.97, 3);
      if( i == 0 )
	{
	  fit_gauss->SetParLimits(0, 600., 2000.);
	  fit_gauss->SetParLimits(1, e_per_bin[i]*0.8, e_per_bin[i]*0.9);
	  fit_gauss->SetParLimits(2, 0.02, 3.);
	  total_E_reco[i]->Fit("fit_gauss", "", "R", e_per_bin[i]*0.82, e_per_bin[i]*0.95);	  
	}
      else if( (i >= 1) && (i < 3) )
	{
	  fit_gauss->SetParLimits(0, 600., 1000.);
	  fit_gauss->SetParLimits(1, e_per_bin[i]*0.86, e_per_bin[i]*0.9);
	  fit_gauss->SetParLimits(2, 0.02, 3.);
	  total_E_reco[i]->Fit("fit_gauss", "", "R", e_per_bin[i]*0.84, e_per_bin[i]*0.915);	  
	}
      else if( (i >= 3) && (i < 5) )
	{
	  fit_gauss->SetParLimits(0, 600., 1000.);
	  fit_gauss->SetParLimits(1, e_per_bin[i]*0.8, e_per_bin[i]*0.98);
	  fit_gauss->SetParLimits(2, 0.02, 3.);
	  total_E_reco[i]->Fit("fit_gauss", "", "R", e_per_bin[i]*0.86, e_per_bin[i]*0.925);	  
	}

      else if( (i >= 5) && (i < N_E_bin) )
	{
	  fit_gauss->SetParLimits(0, 400., 900.);
	  fit_gauss->SetParLimits(1, e_per_bin[i]*0.8, e_per_bin[i]*0.98);
	  fit_gauss->SetParLimits(2, 0.1, 3.);
	  total_E_reco[i]->Fit("fit_gauss", "", "R", e_per_bin[i]*0.85, e_per_bin[i]*0.93);	  
	}

      cout << endl << endl << (fit_gauss->GetParameter(1)) - 3. * (fit_gauss->GetParameter(2)) << endl << endl;

      double gauss_cut = (fit_gauss->GetParameter(1)) - 3. * (fit_gauss->GetParameter(2));
      
      TLine *l1 = new TLine(gauss_cut, 0, gauss_cut, (fit_gauss->GetParameter(0)));
      l1->SetLineColor(2);
      l1->SetLineStyle(9);
      */
	
      //      gPad->SetLogy();
      total_E_reco[i]->SetTitle(Form("#pi^{0} %.1f GeV and Min angle: %.3f#circ", e_per_bin[i], (0.135 / x_graph[i] / M_PI * 180.)));
      //      total_E_reco[i]->GetXaxis()->SetTitle("Reconstructed E[GeV]");
      total_E_reco[i]->GetXaxis()->SetTitle("Deposited E of g1 and g2[GeV]");
      //      total_E_reco[i]->SetStats(0);
      total_E_reco[i]->Draw();
      //      l1->Draw("same");
      

      //      h2_posx_diff[i]->Draw("colorz");

      //      number_of_clus[i]->Draw();
    }


  
}
