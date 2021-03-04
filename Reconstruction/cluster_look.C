#include <stdio.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TF1.h>
#include <string.h>
#include "TLatex.h"
#include <TTree.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <TGraphErrors.h>
#include <TFile.h>
#include "TLatex.h"

void cluster_look()
{

  std::string path = "./data/";
  std::string fileName = path + "outCluster_90k.root";
  //  std::string fileName = path + "outCluster_multiclus_8k.root";
    
  cout << "************ " << fileName << " ************" << endl << endl;

  
  TFile *f1 = new TFile(fileName.c_str(), "Read");
  TTree *outTree = (TTree *) f1->Get("outTree");

  TTreeReader fReader("outTree", f1);

  TTreeReaderValue<int>        N_cluster = {fReader, "N_cluster"};
  TTreeReaderArray<double>     Cl_seed_energy = {fReader, "Cl_seed_energy"};
  TTreeReaderArray<double>     Cl_seed_x = {fReader, "Cl_seed_x"};
  TTreeReaderArray<double>     Cl_seed_y = {fReader, "Cl_seed_y"};
  TTreeReaderArray<double>     Cl_seed_z = {fReader, "Cl_seed_z"};
  TTreeReaderArray<double>     Cl_x = {fReader, "Cl_x"};
  TTreeReaderArray<double>     Cl_y = {fReader, "Cl_y"};
  TTreeReaderArray<double>     Cl_radius = {fReader, "Cl_radius"};
  TTreeReaderArray<double>     Cl_theta = {fReader, "Cl_theta"};
  TTreeReaderArray<double>     Cl_phi = {fReader, "Cl_phi"};
  TTreeReaderArray<double>     Cl_Energy_tot_simul = {fReader, "Cl_Energy_tot_simul"};
  TTreeReaderArray<double>     Cl_size_simul = {fReader, "Cl_size_simul"};
  TTreeReaderArray<double>     Cl_par_a = {fReader, "Cl_par_a"};
  TTreeReaderArray<double>     Cl_x_corr = {fReader, "Cl_x_corr"};
  TTreeReaderArray<double>     Cl_y_corr = {fReader, "Cl_y_corr"};
  
  TTreeReaderValue<double>     e_hit_emcal_x = {fReader, "e_hit_emcal_x"};
  TTreeReaderValue<double>     e_hit_emcal_y = {fReader, "e_hit_emcal_y"};
  TTreeReaderValue<double>     e_hit_emcal_z = {fReader, "e_hit_emcal_z"};
  TTreeReaderValue<double>     g_hit_emcal_x = {fReader, "g_hit_emcal_x"};
  TTreeReaderValue<double>     g_hit_emcal_y = {fReader, "g_hit_emcal_y"};
  TTreeReaderValue<double>     g_hit_emcal_z = {fReader, "g_hit_emcal_z"};
  TTreeReaderValue<double>     e_pjt_emcal_x = {fReader, "e_pjt_emcal_x"};
  TTreeReaderValue<double>     e_pjt_emcal_y = {fReader, "e_pjt_emcal_y"};
  TTreeReaderValue<double>     g_pjt_emcal_x = {fReader, "g_pjt_emcal_x"};
  TTreeReaderValue<double>     g_pjt_emcal_y = {fReader, "g_pjt_emcal_y"};
  TTreeReaderValue<double>     e_px = {fReader, "e_px"};
  TTreeReaderValue<double>     e_py = {fReader, "e_py"};
  TTreeReaderValue<double>     e_pz = {fReader, "e_pz"};
  TTreeReaderValue<double>     e_E = {fReader, "e_E"};
  TTreeReaderValue<double>     g_px = {fReader, "g_px"};
  TTreeReaderValue<double>     g_py = {fReader, "g_py"};
  TTreeReaderValue<double>     g_pz = {fReader, "g_pz"};
  TTreeReaderValue<double>     g_E = {fReader, "g_E"};
  TTreeReaderValue<int>        e_flag_emcal = {fReader, "e_flag_emcal"};
  TTreeReaderValue<int>        g_flag_emcal = {fReader, "g_flag_emcal"};
  TTreeReaderValue<int>        N_hit_emcal = {fReader, "N_hit_emcal"};


  auto Ratio_E_seed_E_recons = new TH1F("Ratio_E_seed_E_recons", "Ratio_E_seed_E_recons", 50, 0., 1.);
  auto E_recons = new TH1F("E_recons", "E_recons", 100, 0., 10.);
  auto diff_posx_g_res_pri = new TH1F("diff_posx_g_res_pri", "diff_posx_g_res_pri", 80, -80., 80.);
  auto diff_posy_g_res_pri = new TH1F("diff_posy_g_res_pri", "diff_posy_g_res_pri", 80, -80., 80.);
  
  auto N_clus_3_x_y = new TH2F("N_clus_3_x_y", "N_clus_3_x_y", 150, -1500., 1500., 100, -1000., 1000.);


  size_t events_numer = 0;  
  double g_eD, g_poxD, g_poyD, g_poxD_cor, g_poyD_cor;
  double e_eD, e_poxD, e_poyD, e_poxD_cor, e_poyD_cor;
  double Eseed_vs_Et;
  int count_3_clus = 0, e_g_simultaneous_count = 0, e_count = 0, g_count = 0, no_primary_hit = 0;


  
  while (fReader.Next())
    {
      if(++events_numer > 89000)
	break;
      
      auto count = static_cast<size_t>(*N_cluster.Get());	      
      double ge = *g_E.Get(), g_pjx = *g_pjt_emcal_x.Get(), g_pjy = *g_pjt_emcal_y.Get();

      if(count == 3)
	{
	  count_3_clus++;
	  cout << events_numer << endl;
	  
	  for(int i = 0 ; i < count ; i++)
	    {
	      //	  cout <<  Cl_seed_energy[i] << " " << Cl_x[i] << " " << Cl_y[i] << " ";
	      Eseed_vs_Et = Cl_seed_energy[i] / Cl_Energy_tot_simul[i];

	      e_eD = *e_E.Get() - Cl_Energy_tot_simul[i] / 1000.;
	      e_poxD = *e_pjt_emcal_x.Get() - Cl_x[i];
	      e_poyD = *e_pjt_emcal_y.Get() - Cl_y[i];
	      e_poxD_cor = *e_pjt_emcal_x.Get() - Cl_x_corr[i];
	      e_poyD_cor = *e_pjt_emcal_y.Get() - Cl_y_corr[i];

	      g_eD = *g_E.Get() - Cl_Energy_tot_simul[i] / 1000.;
	      g_poxD = *g_pjt_emcal_x.Get() - Cl_x[i];
	      g_poyD = *g_pjt_emcal_y.Get() - Cl_y[i];
	      g_poxD_cor = *g_pjt_emcal_x.Get() - Cl_x_corr[i];
	      g_poyD_cor = *g_pjt_emcal_y.Get() - Cl_y_corr[i];


	      if( ((*e_flag_emcal.Get()) == 1) && ((*g_flag_emcal.Get()) == 1) )
		{
		  cout << Cl_x[i] << " " << Cl_y[i] << endl;
		}
	      
	      E_recons->Fill( (Cl_Energy_tot_simul[i] / 1000.) );
	      Ratio_E_seed_E_recons->Fill(Eseed_vs_Et);
	      if( ((*e_flag_emcal.Get()) == 1) && ((*g_flag_emcal.Get()) == 0) )
		{
		  N_clus_3_x_y->Fill(Cl_x[i], Cl_y[i]);
		  if( Cl_x[i] < 0. )
		    {
		      diff_posx_g_res_pri->Fill(g_poxD);
		      diff_posy_g_res_pri->Fill(g_poyD);
		    }
		}
	  	  	  
	    }// loop the N cluster
	  cout << endl << endl;

	  if( ((*e_flag_emcal.Get()) == 1) && ((*g_flag_emcal.Get()) == 1) )
	    e_g_simultaneous_count++;

	  if( ((*e_flag_emcal.Get()) == 1) && ((*g_flag_emcal.Get()) == 0) )
	    e_count++;
	    
	  if( ((*e_flag_emcal.Get()) == 0) && ((*g_flag_emcal.Get()) == 1) )
	    g_count++;

	  if( ((*e_flag_emcal.Get()) == 0) && ((*g_flag_emcal.Get()) == 0) )
	    no_primary_hit++;
	    
	}// N cluster == 3
    }


  cout << "Number of events with 3 cluster reconstructed: " << count_3_clus << endl;
  cout << "Number of events with e & g hit emcal simultaneous: " << e_g_simultaneous_count << endl;
  cout << e_count << " " << g_count << " " << no_primary_hit << endl;

  auto c1 = new TCanvas("c1", "c1", 800, 400);
  c1->Divide(2,1);
  c1->cd(1);
  diff_posx_g_res_pri->Draw();
  c1->cd(2);
  diff_posy_g_res_pri->Draw();
  
  auto c2 = new TCanvas("c2", "c2", 800, 800);
  N_clus_3_x_y->Draw("colorz");
  
}
