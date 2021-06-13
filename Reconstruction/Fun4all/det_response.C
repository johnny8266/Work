

void det_response()
{

  //  TFile *f = new TFile("/vol0/pwang-l/Singularity/my_det/emcal/data/e_pion/outer/gamma_g4eemc_eval_10k.root");
  //  TFile *f = new TFile("/vol0/pwang-l/Singularity/my_det/emcal/data/e_pion/outer/pi_g4eemc_eval_50k.root");
  TFile *f = new TFile("/vol0/pwang-l/Singularity/my_det/emcal/g4eemc_eval.root");
  
  TTree *ntp_gpoint = (TTree*)f->Get("ntp_gpoint");
  TTree *ntp_gshower = (TTree*)f->Get("ntp_gshower");
  TTree *ntp_tower = (TTree*)f->Get("ntp_tower");
  TTree *ntp_cluster = (TTree*)f->Get("ntp_cluster");

  Int_t n_gpoint = (Int_t)ntp_gpoint->GetEntries();
  Int_t n_gshower = (Int_t)ntp_gshower->GetEntries();
  Int_t n_tower = (Int_t)ntp_tower->GetEntries();
  Int_t n_cluster = (Int_t)ntp_cluster->GetEntries();

  Int_t colors[] = {0, 1, 2, 3, 4, 5, 6}; // #colors >= #levels - 1
  gStyle->SetPalette((sizeof(colors)/sizeof(Int_t)), colors);
  Double_t levels[] = {-1.79e308, 1.17e-38, 0.90, 0.95, 1.00, 1.05, 1.10, 1.79e308};

  const Int_t Number = 3;
  Double_t Red[Number]    = { 1.00, 1.00, 0.00};
  Double_t Green[Number]  = { 0.00, 1.00, 0.00};
  Double_t Blue[Number]   = { 0.80, 0.00, 0.80};
  Double_t Length[Number] = { 0.00, 0.50, 1.00 };
  Int_t nb=50;
  TColor::CreateGradientColorTable(Number,Length,Red,Green,Blue,nb);

  

  TLorentzVector v_g;
  Float_t ge, gpt, geta, gphi, gvx, gvy, gvz;
  Float_t event, clusterID, x, y, z, e, ntowers;
  Float_t gparticleID, gflavor;
  Float_t e_e_base = 0., e_x_base = 0., e_y_base = 0., e_z_base = 0., e_n_tow_base = 0.;
  Float_t pi_e_base = 0., pi_x_base = 0., pi_y_base = 0., pi_z_base = 0., pi_n_tow_base = 0.;
  
  int count = 0;
  vector<float> e_x_reco, e_y_reco, e_z_reco, e_e_reco, e_N_towers;
  vector<float> pi_x_reco, pi_y_reco, pi_z_reco, pi_e_reco, pi_N_towers;
  vector<int> i_th_cluster;
  int count_single_clu = 0;
  
  
  ntp_cluster->SetBranchAddress("event", &event);
  ntp_cluster->SetBranchAddress("clusterID", &clusterID);
  ntp_cluster->SetBranchAddress("gparticleID",&gparticleID);
  ntp_cluster->SetBranchAddress("gflavor",&gflavor);
  ntp_cluster->SetBranchAddress("x", &x);
  ntp_cluster->SetBranchAddress("y", &y);
  ntp_cluster->SetBranchAddress("z", &z);
  ntp_cluster->SetBranchAddress("e", &e);
  ntp_cluster->SetBranchAddress("ntowers", &ntowers);

  //  cout << "     i || clusterID event" << endl;
  for(int i = 0 ; i < n_cluster ; i++)
    //  for(int i = 0 ; i < 10 ; i++)
    {
      ntp_cluster->GetEntry(i);

      //      cout << i << " " << clusterID << " " << event << ": [" << x << ", " << y << ", " << z << "] || " << e << endl;
      cout << i << " " << clusterID << " " << event << " " << e << " || ";
      cout << gparticleID << " " << gflavor << " " << endl;
      
      
      int i_eve = (int)event;
      
      if( i_eve == count )
	{
	  if( gflavor > 0 )  // gamma = 22, this is not int type, comparision is safer than equal
	    {
	      if( e > e_e_base )
		{
		  e_e_base = e;
		  e_x_base = x;
		  e_y_base = y;
		  e_z_base = z;
		  e_n_tow_base = ntowers;
		  //	      cout << i_eve << " " << count << " " << e << endl;
		}
	    }
	  else if( gflavor < 0 )  // pi- = -211, this is not int type, comparision is safer than equal
	    {
	      if( e > pi_e_base )
		{	
		  pi_e_base = e;
		  pi_x_base = x;
		  pi_y_base = y;
		  pi_z_base = z;
		  pi_n_tow_base = ntowers;
		  //	      cout << i_eve << " " << count << " " << e << endl;
		}
	    }
	  
	  
	}
      else
	{
	  //	  cout << i_eve << endl;
	  i_th_cluster.push_back(count);
	  x_reco.push_back(e_x_base);
	  y_reco.push_back(e_y_base);
	  z_reco.push_back(e_z_base);
	  e_reco.push_back(e_e_base);
	  N_towers.push_back(e_n_tow_base);
	  count_single_clu++;
	  
	  i_th_cluster.push_back(count);
	  x_reco.push_back(pi_x_base);
	  y_reco.push_back(pi_y_base);
	  z_reco.push_back(pi_z_base);
	  e_reco.push_back(pi_e_base);
	  N_towers.push_back(pi_n_tow_base);
	  count_single_clu++;
	  
	  count++;
	  i = i - 1;
	  e_e_base = 0.;  e_x_base = 0.;  e_y_base = 0.;  e_z_base = 0.;  e_n_tow_base = 0.;
	  pi_e_base = 0.;  pi_x_base = 0.;  pi_y_base = 0.;  pi_z_base = 0.;  pi_n_tow_base = 0.;
	}

      
    }
  //  cout << endl << endl << count_single_clu << endl << endl;

  //  for(int i = 0 ; i < 10 ; i++)
  //    cout << i_th_cluster[i] << " " << e_reco[i] << " || [" << x_reco[i] << ", " << y_reco[i] << ", " << z_reco[i] << "]" << endl;
  
      
  /*  
  auto hit_distribution = new TH2F("hit_distribution", "hit_distribution", 150, -150., 150., 150, -150., 150.);

  
  ntp_tower->SetBranchAddress("event", &event);
  ntp_tower->SetBranchAddress("x", &x);
  ntp_tower->SetBranchAddress("y", &y);
  ntp_tower->SetBranchAddress("z", &z);
  
  for(int i = 0 ; i < n_tower ; i++)
    {
      ntp_tower->GetEntry(i);

      if(event == 4)
	hit_distribution->Fill(x, y);
    }
  
  //  hit_distribution->Draw("colorz");
  */
  
  auto P_E_distribution = new TH2F("P_E_distribution", "P_E_distribution", 110, 0., 11., 110, 0., 11.);
  auto event_e = new TH2F("event_e", "event_e", 110, 0., 11., 110, 0., 11.);
  
  
  ntp_gshower->SetBranchAddress("event",&event);
  ntp_gshower->SetBranchAddress("gparticleID",&gparticleID);
  ntp_gshower->SetBranchAddress("gflavor",&gflavor);
  ntp_gshower->SetBranchAddress("gvx",&gvx);
  ntp_gshower->SetBranchAddress("gvy",&gvy);
  ntp_gshower->SetBranchAddress("gvz",&gvz);
  ntp_gshower->SetBranchAddress("gpt",&gpt);
  ntp_gshower->SetBranchAddress("geta",&geta);
  ntp_gshower->SetBranchAddress("gphi",&gphi);
  ntp_gshower->SetBranchAddress("ge",&ge);


  //  for(int i = 0 ; i < count_single_clu ; i++)
  for(int i = 0 ; i < n_gshower ; i++)
    {
      //      int i_eve = i_th_cluster[i];
      //      ntp_gshower->GetEntry(i_eve);
      ntp_gshower->GetEntry(i);
      v_g.SetPtEtaPhiE(gpt, geta, gphi, ge);
      double gx = v_g.Px(), gy = v_g.Py(), gz = v_g.Pz(), g_E = v_g.E();

      //      cout << i << " " << event << " " << gflavor << " " << gparticleID << " || ";
      //      cout << gx << " " << gy << " " << gz << endl;
      
      /*
      v_g.SetPtEtaPhiE(gpt, geta, gphi, ge);
      //      cout << gpt << " " << geta << " " << gphi << " " << ge << endl;

      double gx = v_g.Px(), gy = v_g.Py(), gz = v_g.Pz(), g_E = v_g.E();
      //      cout << gx << " " << gy << " " << gz << endl;
      double total_p = std::sqrt(gx * gx + gy * gy + gz * gz);
      double x_pjt = gx / gz * z_reco[i];
      double y_pjt = gy / gz * z_reco[i];
      double dis_reco = std::sqrt( (x_pjt - x_reco[i]) * (x_pjt - x_reco[i]) + (y_pjt - y_reco[i]) * (y_pjt - y_reco[i]) );
      //      cout << x_pjt << " " << x_reco[i] << " " << dis_reco << endl;

      //      if( (N_towers[i] > 4.) && (dis_reco < 3.) )
      if( N_towers[i] > 4. )
	{
	  if( e_reco[i] > (ge * 0.1) )
	    {
	      P_E_distribution->Fill(total_p, e_reco[i]);
	    }
	}

      */
    }



  /*
  // Draw and save the plots

  TFile *g = new TFile("result.root", "UPDATE");

  TCanvas *c1 = new TCanvas("c1", "c1", 900, 900);
  P_E_distribution->SetStats(0);
  P_E_distribution->SetName("P_E_distribution_pi");
  P_E_distribution->SetTitle("gamma energy deposition v.s. primary momentum");
  P_E_distribution->GetXaxis()->SetTitle("primary momentum [GeV]");
  P_E_distribution->GetYaxis()->SetTitle("E deposition in EEMC [GeV]");
  P_E_distribution->Draw("colorz");
  
  P_E_distribution->Write();
  g->Close();
*/

  
}
