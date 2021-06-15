#include <iostream>
using namespace std;


double E_reso_2(double par_a, double par_b, double e)
{
  return sqrt( par_a * par_a + (par_b * par_b) / e ) / 100.;  
}

double E_reso_3(double par_a, double par_b, double par_c, double e)
{
  return sqrt( par_a * par_a + (par_b * par_b) / e + (par_c * par_c) / (e * e) ) / 100.;  
}

double E_correct(double par_a, double par_b, double e)
{
  return par_a * e + par_b;  
}


void rejection_efficiency()
{
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


  

  const int E_bin = 9;
  int flag = 0; // 0: glass, 1: crystal
  int par_N = 2;
  string File_num[E_bin] = {"0.5", "1", "2", "3", "5", "7", "10", "13", "18"};
  double e_per_bin[E_bin] = {0.5, 1., 2., 3., 5., 7., 10., 13., 18.};
  double pion_rejection_factor[E_bin] = {};
  double e_corr[2] = {}, e_reso_pars_2[2] = {}, e_reso_pars_3[3] = {};
  double err_x[E_bin] = {0.}, err_y[E_bin] = {}; 

  if( flag == 1 )
    {
      e_corr[0] = 0.0884224;  e_corr[1] = 0.0950983;
      e_reso_pars_2[0] = 1.38664;  e_reso_pars_2[1] = 2.51748;
      e_reso_pars_3[0] = 1.51918;  e_reso_pars_3[1] = 2.24762;  e_reso_pars_3[2] = 0.760284;
    }
  else
    {
      e_corr[0] = 0.115516;  e_corr[1] = 0.0308727;
      e_reso_pars_2[0] = 1.02691;  e_reso_pars_2[1] = 3.01323;
      e_reso_pars_3[0] = 1.77746,  e_reso_pars_3[1] = 1.08107;  e_reso_pars_3[2] = 1.9217;
    }
  

  TLorentzVector v_g;
  TH1F *h1[E_bin];
  TCanvas *c1 = new TCanvas("c1", "c1", 960, 960);
  c1->Divide(3,3);

  
  for(int i = 0 ; i < E_bin ; i ++)
    {
      c1->cd(i+1);
      gPad->SetLogy();
      string root_file_name;

      if( flag == 1 )
	{
	  if(i == 0)
	    root_file_name = "/vol0/pwang-l/Singularity/my_det/sub_crystal/data/pion/g4eemc_crystal_eval_mono_" + File_num[i] + "_GeV.root";
	  else
	    root_file_name = "/vol0/pwang-l/Singularity/my_det/sub_crystal/data/pion/g4eemc_crystal_eval_mono_" + File_num[i] + ".0_GeV.root";	  
	}
      else
	{
	  if(i == 0)
	    root_file_name = "/vol0/pwang-l/Singularity/my_det/sub_glass/data/pion/g4eemc_glass_eval_mono_" + File_num[i] + "_GeV.root";
	  else
	    root_file_name = "/vol0/pwang-l/Singularity/my_det/sub_glass/data/pion/g4eemc_glass_eval_mono_" + File_num[i] + ".0_GeV.root";	  
	}

      const char *rfn = root_file_name.c_str();
      TFile *f = new TFile(rfn);
  
      TTree *ntp_gshower = (TTree*)f->Get("ntp_gshower");
      TTree *ntp_cluster = (TTree*)f->Get("ntp_cluster");

      Int_t n_gshower = (Int_t)ntp_gshower->GetEntries();
      Int_t n_cluster = (Int_t)ntp_cluster->GetEntries();
      
      h1[i] = new TH1F(Form("h1_%d", i), Form("h1_%d", i), 110, 0., 1.1);

      
      cout << File_num[i] << " GeV have: " << n_gshower << " pions." << endl;
      
      TLorentzVector v_g;
      Float_t ge, gpt, geta, gphi, gvx, gvy, gvz;
      Float_t event, clusterID, x, y, z, e, ntowers;
      Float_t gparticleID, gflavor;
      Float_t e_p_cut;
      Float_t e_base = 0., x_base = 0., y_base = 0., z_base = 0., n_tow_base = 0.;
      Float_t ge_base, gpt_base, geta_base, gphi_base;
  
      int count = 0, reject_count = 0;
      vector<float> x_reco, y_reco, z_reco, e_reco, N_towers;
      vector<int> i_th_cluster;
      int count_single_clu = 0;


      ntp_cluster->SetBranchAddress("event", &event);
      ntp_cluster->SetBranchAddress("gflavor",&gflavor);
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

		  ge_base = ge;
		  gpt_base = gpt;
		  geta_base = geta;
		  gphi_base = gphi;
		  //	      cout << i_eve << " " << count << " " << e << endl;
		}
	    }
	  else
	    {
	      v_g.SetPtEtaPhiE(gpt_base, geta_base, gphi_base, ge_base);
	      double total_P = v_g.P();
	      if(total_P > 0.1)
		{
		  double e_after_correct = E_correct(e_corr[0], e_corr[1], e_base) + e_base;
		  double e_over_p = e_after_correct / total_P;
		  if( par_N == 2 )
		    e_p_cut = 1. - 2. * E_reso_2(e_reso_pars_2[0], e_reso_pars_2[1], e_per_bin[i]);
		  else if( par_N == 3 )
		    e_p_cut = 1. - 2. * E_reso_3(e_reso_pars_3[0], e_reso_pars_3[1], e_reso_pars_3[2], e_per_bin[i]);
	      
		  if( isnan(e_over_p) == false || (e_over_p) < 10. )
		    {
		      h1[i]->Fill(e_over_p);
		  
		      if( e_over_p > e_p_cut)
			{
			  //			  cout << e_base << " " << e_after_correct << " " << total_P << " || " << e_over_p << "  " << e_p_cut << endl;
			  reject_count++;		      
			}
		    }		  
		}
	      
	      count++;
	      j = j - 1;
	      e_base = 0.;  x_base = 0.;  y_base = 0.;  z_base = 0.;  n_tow_base = 0.;
	      gpt_base = 0.; geta_base = 0.; gphi_base = 0.; ge_base = 0.;
	    }
	}

      //      cout << count << endl;
      cout << "cut for this energy: " << e_p_cut << endl;
      cout << reject_count << " pions pass the cut..." << endl;
      
      if( reject_count != 0 )
	{
	  double drc = reject_count;
	  cout << "rejection factor: " << n_gshower / reject_count << endl;
	  pion_rejection_factor[i] = n_gshower / reject_count;
	  err_y[i] = n_gshower / reject_count * sqrt(drc) / drc;
	  cout << "Error: " << n_gshower / reject_count * sqrt(drc) / drc << endl << endl; 
	}
      else
	{
	  cout << "rejection factor might greater than 10000..... " << endl << endl;
	  //	  pion_reject->SetBinContent(e_per_bin[i], 10000 / reject_count);
	}
      
      reject_count = 0;
      
      
      string plot_name = "pion " + File_num[i] + " GeV and E/P matching cut: " + to_string(e_p_cut);
      const char *pn = plot_name.c_str();
      
      h1[i]->SetTitle(pn);
      h1[i]->SetStats(0);
      h1[i]->GetXaxis()->SetTitle("E over P");
      h1[i]->Draw();
      TLine *line_cut = new TLine(e_p_cut, 0, e_p_cut, 10000);
      line_cut->SetLineColor(1);
      line_cut->SetLineStyle(9);
      line_cut->Draw("same");
    } // Loop all the different energy bin



  
  TGraphErrors *pion_reject = new TGraphErrors(E_bin, e_per_bin, pion_rejection_factor, err_x, err_y);

  TFile *g = new TFile("result.root", "UPDATE");
  
  TCanvas *c2 = new TCanvas("c2", "c2", 800, 800);
  gPad->SetLogy();
  gPad->SetLogx();
  pion_reject->GetXaxis()->SetLimits(0.2, 30.);
  pion_reject->SetMaximum(30000.);
  pion_reject->SetMinimum(1.);
  if( flag == 1 )
    {
      pion_reject->SetTitle(Form("Pion rejection [crystal] cut(E: %d pars)", par_N));
      pion_reject->SetName(Form("Pion_reject_crystal_%dpars", par_N));
      pion_reject->SetMarkerStyle(25);
      pion_reject->SetMarkerSize(1);
    }
  else
    {
      pion_reject->SetTitle(Form("Pion rejection [glass] cut(E: %d pars)", par_N));
      pion_reject->SetName(Form("Pion_reject_glass_%dpars", par_N));
      pion_reject->SetMarkerStyle(24);
      pion_reject->SetMarkerSize(1);
    }
  pion_reject->GetYaxis()->SetTitle("Rejection factor");
  pion_reject->GetXaxis()->SetTitle("E [GeV]");
  pion_reject->Draw("AP");
  pion_reject->Write();

  g->Close();

}
