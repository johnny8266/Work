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



void close_event_diffE()
{

  const int N_E_bin = 3;
  string File_num[N_E_bin] = {"1.0", "5.0", "10.0"};
  double e_per_bin[N_E_bin] = {1.0, 5.0, 10.0};
  // string File_num[N_E_bin] = {"1.0"};
  // double e_per_bin[N_E_bin] = {1.0};
  double e_statement[13] = {0.5, 1., 2., 3., 4., 5., 7., 10., 13., 18., 30., 50., 52.};
  double angle_statement[9] = {0., 0.3, 0.6, 0.9, 1.2, 1.6, 2., 3., 4.};
  
  int flag = 1; // 0: glass, 1: crystal
  TLorentzVector LV;
  double pri_Px = 0., pri_Py = 0., pri_Pz = 0.;
  double pjx = 0., pjy = 0., pjz = 0., diff_x = 0., diff_y = 0.;
  

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
  TH1F *zpos_reco_1[N_E_bin];
  TH1F *zpos_reco_2[N_E_bin];

  TH1F *angle_range[N_E_bin][8];
  const int N_bin_plot = 220;
  for(int i = 0 ; i < N_E_bin ; i++)
    for(int j = 0 ; j < 8 ; j++)
      angle_range[i][j] = new TH1F(Form("angle_range_%d_%d", i, j), Form("angle_range_%d_%d", i, j), N_bin_plot, 0., (e_per_bin[i]+10.0) * 1.1 );

  
  TH2F *h2_posx_diff[N_E_bin], *h2_posy_diff[N_E_bin];
  TH2F *cluster_1 = new TH2F("cluster_1", "cluster_1", 30, 0., 60., 30, 0., 60.);
  TH2F *cluster_2 = new TH2F("cluster_2", "cluster_2", 30, 0., 60., 30, 0., 60.);

  TProfile *hprof_x[N_E_bin], *hprof_y[N_E_bin];

  TGraph *clusters_effi[N_E_bin];



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

      zpos_reco_1[i] = new TH1F(Form("zpos_reco_1_%d", i), Form("zpos_reco_1_%d", i), 100, -200., -190.);
      zpos_reco_1[i]->GetXaxis()->SetRangeUser(-200., -190.);
      zpos_reco_2[i] = new TH1F(Form("zpos_reco_2_%d", i), Form("zpos_reco_2_%d", i), 100, -200., -190.);
      zpos_reco_2[i]->GetXaxis()->SetRangeUser(-200., -190.);
      
      h1_posx_diff[i] = new TH1F(Form("h1_posx_diff_%d", i), Form("h1_posx_diff_%d", i), 120, -3., 3.);
      h1_posy_diff[i] = new TH1F(Form("h1_posy_diff_%d", i), Form("h1_posy_diff_%d", i), 120, -3., 3.);
      h2_posx_diff[i] = new TH2F(Form("h2_posx_diff_%d", i), Form("h2_posx_diff_%d", i), 100, -100., 100., 120, -3., 3.);
      h2_posy_diff[i] = new TH2F(Form("h2_posy_diff_%d", i), Form("h2_posy_diff_%d", i), 100, -100., 100., 120, -3., 3.);
      hprof_x[i] = new TProfile(Form("Tprofile_x_%d", i), Form("Tprofile_x_%d", i), 60, -60., 60., -3., 3.);
      hprof_y[i] = new TProfile(Form("Tprofile_y_%d", i), Form("Tprofile_y_%d", i), 60, -60., 60., -3., 3.);

    }

  


  // Read the data and calculate the reconstruction
  //
  for( int i = 0 ; i < N_E_bin ; i++ )
    {
      string root_file_name;

      
      if( flag == 1 )
	{
	  root_file_name = "/vol0/pwang-l/Singularity/my_det/sub_crystal/data/pos/E_portion/g4eemc_crystal_eval_E_portion_10.0_" + File_num[i] + "_GeV.root";
	}
      else
	{
	  root_file_name = "/vol0/pwang-l/Singularity/my_det/sub_glass/data/pos/g4eemc_glass_eval_mono_" + File_num[i] + "_GeV_close_event_study.root";
	}

      

      const char *rfn = root_file_name.c_str();
      TFile *f = new TFile(rfn);

      TTree *ntp_shower = (TTree*)f->Get("ntp_gshower");
      Int_t n_shower = (Int_t)ntp_shower->GetEntries();
      TTree *ntp_cluster = (TTree*)f->Get("ntp_cluster");
      Int_t n_cluster = (Int_t)ntp_cluster->GetEntries();
      // TTree *ntp_tower = (TTree*)f->Get("ntp_tower");
      // Int_t n_tower = (Int_t)ntp_tower->GetEntries();

      
      
      //declare the variables used in the following ntuples
      //
      TLorentzVector v_g, pr_1, pr_2;
      TVector3 p1, p2;
      Float_t ge, gpt, geta, gphi, gvx, gvy, gvz, gparticleID;
      Float_t event, clusterID, x, y, z, e, ntowers, e_base = 0., x_base = 0., y_base = 0., z_base = 0., n_tow_base = 0.;
      Float_t e_1_1, e_2_1, e_2_2;
      vector<float> x_reco, y_reco, z_reco, e_reco, N_towers, angle_2p;
      vector<float> p1x, p1y, p1z, p2x, p2y, p2z;
      vector<int> i_th_cluster;
      int count_single_clu = 0, count_double_clu = 0, count = 0, ID_sum = 0;




      

      ntp_shower->SetBranchAddress("event", &event);
      ntp_shower->SetBranchAddress("gpt", &gpt);
      ntp_shower->SetBranchAddress("geta", &geta);
      ntp_shower->SetBranchAddress("gphi", &gphi);
      ntp_shower->SetBranchAddress("ge", &ge);

      for(int j = 0 ; j < n_shower ; j = j + 2)
	{
	  ntp_shower->GetEntry(j);
	  pr_1.SetPtEtaPhiE(gpt, geta, gphi, ge);
	  p1.SetXYZ(pr_1.Px(), pr_1.Py(), pr_1.Pz());
	  p1x.push_back(pr_1.Px());  p1y.push_back(pr_1.Py());  p1z.push_back(pr_1.Pz());

	  
	  ntp_shower->GetEntry(j+1);
	  pr_2.SetPtEtaPhiE(gpt, geta, gphi, ge);
	  p2.SetXYZ(pr_2.Px(), pr_2.Py(), pr_2.Pz());
	  p2x.push_back(pr_2.Px());  p2y.push_back(pr_2.Py());  p2z.push_back(pr_2.Pz());
	  //	  cout << "Angle of 2 primary particle of event " << event << " th: " << p1.Angle(p2) / TMath::Pi() * 180. << endl;

	  angle_2p.push_back(p1.Angle(p2) / TMath::Pi() * 180.);

	  int j_eve = (int)event;

	  if( j_eve == 1054 )
	    {
	      cluster_1->SetTitle(Form("p1[%.2f, %.2f, %.2f]      p2[%.2f, %.2f, %.2f]", pr_1.Px(), pr_1.Py(), pr_1.Pz(), pr_2.Px(), pr_2.Py(), pr_2.Pz()));	    
	    }


	  if( j_eve == 1053 )
	    {
	      cluster_2->SetTitle(Form("p1[%.2f, %.2f, %.2f]      p2[%.2f, %.2f, %.2f]", pr_1.Px(), pr_1.Py(), pr_1.Pz(), pr_2.Px(), pr_2.Py(), pr_2.Pz()));	    
	    }
	  
	}




      /*
      ntp_tower->SetBranchAddress("event", &event);
      ntp_tower->SetBranchAddress("x", &x);
      ntp_tower->SetBranchAddress("y", &y);
      ntp_tower->SetBranchAddress("z", &z);
      ntp_tower->SetBranchAddress("e", &e);
      
      for(int j = 0 ; j < n_tower ; j++) 
	{
	  ntp_tower->GetEntry(j);

	  int j_eve = (int)event;

	  //	  if( j_eve == 1054 )
	    //	    cout << x << " " << y << " " << e << endl;
	  cluster_1->Fill(x, y, e);

	  //	  if( j_eve == 1053 )
	  //	    cluster_2->Fill(x, y, e);
	}
      */



      
      
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
	    }
	  else
	    {
	      ntp_cluster->GetEntry(j - 1);
	      j_eve = (int)event;
	      
	      //	      cout << clusterID << endl;
	      
	      if( clusterID == 0 )
		{
		  //		  cout << "Only 1 cluster reconstructed: " << j_eve << endl;		  
		  recons_1_angle[i]->Fill(angle_2p[j_eve]);
		  E_reco_1[i]->Fill(e);
		  zpos_reco_1[i]->Fill(z);
		  count_single_clu++;
		}

	      if( clusterID == 1 )
		{
		  //		  cout << "2 clusters reconstructed: " << j_eve << endl;
		  //		  recons_2_angle[i]->Fill(angle_2p[j_eve]);
		    
		  e_2_1 = e;
		  E_reco_2[i]->Fill(e_2_1);
		  zpos_reco_2[i]->Fill(z);
		  ID_sum = ID_sum + gparticleID;
		  float D_1_1 = std::sqrt( std::pow((z / p1z[j_eve] * p1x[j_eve] - x), 2) + std::pow((z / p1z[j_eve] * p1y[j_eve] - y), 2) );
		  float D_1_2 = std::sqrt( std::pow((z / p2z[j_eve] * p2x[j_eve] - x), 2) + std::pow((z / p2z[j_eve] * p2y[j_eve] - y), 2) );

		  //		  if( (D_1_1 > 10.) && (D_1_2 > 10.) )
		  //		    cout << j_eve << " " << z / p1z[j_eve] * p1x[j_eve] << " " << z / p1z[j_eve] * p1y[j_eve] << " || " << x << " " << y << endl;
		  /*
		  if( gparticleID < 1.5 )
		    {
		      if( D_1_1 > D_1_2 )
			cout << "1 Wrong pair !!!  " << D_1_1 << " " << D_1_2 << " " << gparticleID << " " << angle_2p[j_eve] << endl;		      
		    }
		  else
		    {
		      if( D_1_1 < D_1_2 )
			cout << "1 Wrong pair !!!" << endl;		      		      
		    }
		  */

		  
		  ntp_cluster->GetEntry(j - 2);
		  e_2_2 = e;
		  E_reco_2[i]->Fill(e_2_2);
		  zpos_reco_2[i]->Fill(z);
		  ID_sum = ID_sum + gparticleID;
		  float D_2_1 = std::sqrt( std::pow((z / p1z[j_eve] * p1x[j_eve] - x), 2) + std::pow((z / p1z[j_eve] * p1y[j_eve] - y), 2) );
		  float D_2_2 = std::sqrt( std::pow((z / p2z[j_eve] * p2x[j_eve] - x), 2) + std::pow((z / p2z[j_eve] * p2y[j_eve] - y), 2) );

		  //		  if( (D_2_1 > 10.) && (D_2_2 > 10.) )
		  //		  if( (D_1_1 > 10.) && (D_1_2 > 10.) )
		  //		    cout << j_eve << " " << z / p1z[j_eve] * p1x[j_eve] << " " << z / p1z[j_eve] * p1y[j_eve] << " || " << x << " " << y << endl << endl;

		  /*
		  if( gparticleID < 1.5 )
		    {
		      if( D_2_1 > D_2_2 )
			cout << "2 Wrong pair !!!  " << D_2_1 << " " << D_2_2 << " " << gparticleID << " " << angle_2p[j_eve] << endl;		      
		    }
		  else
		    {
		      if( D_2_1 < D_2_2 )
			cout << "2 Wrong pair !!!" << endl;		      		      
		    }
		  */

		  for(int k = 0 ; k < 8 ; k++)
		    {
		      if( (angle_2p[j_eve] >= angle_statement[k]) && (angle_2p[j_eve] < angle_statement[k+1]) )
			{
			  angle_range[i][k]->Fill(e_2_1);
			  angle_range[i][k]->Fill(e_2_2);			  
			}		      
		    }
		  
		  
		  // if( (ID_sum != 3) && (angle_2p[j_eve] > 1.5) )
		  //   cout << e_2_1 << " " << e_2_2 << " " << ID_sum << " " << angle_2p[j_eve] << endl;

		  ID_sum = 0;
		  count_double_clu++;

		  
		  //		  E_reco_2[i]->Fill((e_2_1 + e_2_2));
		  //		  if( (e_2_1 < e_per_bin[i]) && (e_2_2 < e_per_bin[i]) )
		  if( (e_2_1 < e_per_bin[i]) || (e_2_2 < e_per_bin[i]) )		 
		    {
		      //		      cout << e_2_1 << " " << e_2_2 << endl;
		      //		      count_double_clu++;
		      recons_2_angle[i]->Fill(angle_2p[j_eve]);
		    }


		}

	     	      
	      count++;
	      //	      count_double_clu++;
	      j = j - 1;

	      e_base = 0.;  x_base = 0.;  y_base = 0.; z_base = 0.; n_tow_base = 0.;
	    }
	}// Finish the cluster ntuple


      
      
      cout << endl << endl << "Single cluster count: " << count_single_clu << endl;
      cout << "Double cluster count: " << count_double_clu << endl << endl;
      count_single_clu = 0;  count_double_clu = 0;

      
    }// Finish the files read      





  
  //  TF1 fit_gauss = new TF1();
  TGraphErrors *angle_effi[N_E_bin];
  TLegend *legend_bin[3];
  TCanvas *c1 = new TCanvas("c1", "c1", 1800, 900);
  c1->Divide(4,2);
  for(int i = 0 ; i < N_E_bin ; i++)
    {
      auto *fun_e_res = new TF1("fun_e_res", GausM, e_per_bin[i]*0.6, e_per_bin[i]*1.1, 3);
      double x[8], y[8], x_err[8], y_err[8];
      for(int j = 0 ; j < 8 ; j++)
	{
	  c1->cd(j+1);

	  int mean_bin = e_per_bin[i] / (e_per_bin[i] * 1.1 / N_bin_plot);
	  int width_bin = mean_bin * 0.5;
	  int fit_flag = angle_range[i][j]->Integral((mean_bin - width_bin), (mean_bin + width_bin));

	  //	  cout << "INT: " << mean_bin << " " << width_bin << " " << fit_flag << endl;
	  angle_range[i][j]->SetTitle(Form("E[%.1f, 10.0]    Angle: %.1f ~ %.1f", e_per_bin[i], angle_statement[j], angle_statement[j+1]));
	  angle_range[i][j]->Draw();
	  c1->Update();

	  double y_axis_max = gPad->GetUymax();
	  double all_2_clus_events = angle_range[i][j]->GetEntries();
	  double selected_2_clus_events = 0.;
	  int mean_separa = (e_per_bin[i] * 1.4) / ((e_per_bin[i]+10.) * 1.1 / N_bin_plot);
	  int low_separa = (e_per_bin[i] * 0.35) / ((e_per_bin[i]+10.) * 1.1 / N_bin_plot);
	  TLine *l_1 = new TLine((e_per_bin[i] * 0.3), 0, (e_per_bin[i] * 0.3), y_axis_max);
	  TLine *l_2 = new TLine((e_per_bin[i] * 1.2), 0, (e_per_bin[i] * 1.2), y_axis_max);
	  l_1->SetLineColor(2);
	  l_1->SetLineWidth(2);
	  l_1->SetLineStyle(9);
	  l_1->Draw("same");
	  l_2->SetLineColor(2);
	  l_2->SetLineWidth(2);
	  l_2->SetLineStyle(9);
	  l_2->Draw("same");
	  gPad->SetLogy();


	  
	  //	  cout << mean_separa << " " << low_separa << endl;
	  int low_e_count = 0, high_e_count = 0;
	  
	  if( j < 5 )
	    {
	      low_e_count = angle_range[i][j]->Integral(low_separa, mean_separa);
	      if(i < 2)
		selected_2_clus_events = low_e_count * 2.;
	      else
		selected_2_clus_events = low_e_count;

	      cout << j << ": " << mean_separa << " || " << low_e_count << endl;
	    }
	  else 
	    {
	      low_e_count = angle_range[i][j]->Integral(low_separa, mean_separa);
	      if(i < 2)
		selected_2_clus_events = low_e_count * 2.;
	      else
		selected_2_clus_events = low_e_count;

	      high_e_count = angle_range[i][j]->Integral(mean_separa, N_bin_plot);

	      legend_bin[j-5] = new TLegend(0.25, 0.6, 0.55, 0.85);
	      legend_bin[j-5]->SetBorderSize(0);
	      legend_bin[j-5]->AddEntry((TObject*)0, Form("Lower E count: %d", low_e_count), "");
	      legend_bin[j-5]->AddEntry((TObject*)0, Form("Higher E count: %d", high_e_count), "");
	      legend_bin[j-5]->Draw("same");
	      //	      selected_2_clus_events = low_e_count + high_e_count;
	      //	      cout << j << ": " << mean_separa << " || " << low_e_count << ", " << high_e_count << endl;
	    }

	  x[j] = (angle_statement[j] + angle_statement[j+1]) / 2.;
	  x_err[j] = 0.;
	  y[j] = selected_2_clus_events / all_2_clus_events;
	  y_err[j] = y[j] * std::sqrt( 1. / selected_2_clus_events + 1. / all_2_clus_events );
	  cout << "y error: " << y_err[j] << endl;
	}
      angle_effi[i] = new TGraphErrors(8, x, y, x_err, y_err);
    }

  
  TLegend *legend_e_reso;
  legend_e_reso = new TLegend(0.5, 0.45, 0.85, 0.65);
  legend_e_reso->SetBorderSize(0);
  TCanvas *c3 = new TCanvas("c3", "c3", 800, 800);
  //  c3->Divide(3,1);
  for(int i = 0 ; i < N_E_bin ; i++)
    {
      //      c3->cd(i+1);
      //      angle_effi[i]->SetTitle(Form("Efficiency of %.1f GeV gamma", e_per_bin[i]));
      angle_effi[i]->SetTitle("Efficiency of different energy gamma separation");
      angle_effi[i]->SetName(Form("angle_effi_%d", i));
      angle_effi[i]->SetMarkerSize(1.);
      angle_effi[i]->SetMarkerStyle(24+i);
      angle_effi[i]->SetMarkerColor(i+1);

      //      angle_effi[i]->Draw("PSame");
      
      if( i == 0 )
	angle_effi[i]->Draw("AP");
      else
	angle_effi[i]->Draw("P same");

      legend_e_reso->AddEntry(Form("angle_effi_%d", i), Form("%.1f GeV and 10.0 GeV", e_per_bin[i]), "p");
    }
  legend_e_reso->Draw("same");
  

  
  /*
  //Drawing the plots
  //
  TCanvas *c1 = new TCanvas("c1", "c1", 1800, 600);
  c1->Divide(3,1);
  for(int i = 0 ; i < N_E_bin ; i++)
    {
      c1->cd( (i+1) );
      
      for(int j = 0 ; j < 2 ; j++)      
	{
	  //	  c1->cd( (i+1) + 3 * j );

	  if(j == 0)
	    recons_2_angle[i]->Draw();
	  else
	    recons_1_angle[i]->Draw("same");
	}
    }

  
  TCanvas *c2 = new TCanvas("c2", "c2", 1350, 900);
  c2->Divide(3,2);
  for(int i = 0 ; i < N_E_bin ; i++)
    {
      for(int j = 0 ; j < 2 ; j++)      
	{
	  c2->cd( (i+1) + 3 * j );

	  if(j == 0)
	    //	    zpos_reco_2[i]->Draw();
	    E_reco_2[i]->Draw();
	  else
	    //	    zpos_reco_1[i]->Draw();
	    E_reco_1[i]->Draw();
	}
    }
  */
  
}
