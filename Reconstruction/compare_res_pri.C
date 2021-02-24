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

using namespace std;

double GausM(double *x, double *par)
{
  return par[0] * exp(-0.5 * TMath::Power(((x[0] - par[1]) / par[2]), 2)) + 0.11;
}


Double_t Eresolution_fit(double *x, double *par)
{
  return sqrt(TMath::Power(par[0] / sqrt(x[0]), 2) + TMath::Power(par[1], 2));
}


void compare_res_pri()
{
  
  std::string path = "./data/";
  std::string fileName = path + "outCluster.root";
  //  std::string fileName = path + "outCluster_multiclus_8k.root";
    
  cout << "************ " << fileName << " ************" << endl;

  
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

  auto diffE_e_res_pri = new TH1F("diffE_e_res_pri", "diffE_e_res_pri", 100, -5., 5.);
  auto diffE_g_res_pri = new TH1F("diffE_g_res_pri", "diffE_g_res_pri", 100, -5., 5.);
  auto diff_posx_g_res_pri = new TH1F("diff_posx_g_res_pri", "diff_posx_g_res_pri", 40, -200., 200.);
  auto diff_posy_g_res_pri = new TH1F("diff_posy_g_res_pri", "diff_posy_g_res_pri", 40, -200., 200.);
  auto diff_posx_cor_g_res_pri = new TH1F("diff_posx_cor_g_res_pri", "diff_posx_cor_g_res_pri", 40, -200., 200.);
  auto diff_posy_cor_g_res_pri = new TH1F("diff_posy_cor_g_res_pri", "diff_posy_cor_g_res_pri", 40, -200., 200.);
  auto hit_cl_size = new TH1F("hit_cl_size", "hit_cl_size", 20, 0., 20.);

  auto diffx_g_E = new TH2F("diffx_g_E", "diffx_g_E", 40, -200., 200., 10, 0., 10.);
  auto diffx_g_x = new TH2F("diffx_g_x", "diffx_g_x", 40, 0., 200., 50, -1500., 0.);
  
  size_t events_numer = 0;  
  
  while (fReader.Next())
    {
      if(++events_numer > 8500)
	break;
      auto count = static_cast<size_t>(*N_cluster.Get());

      for(int i = 0 ; i < count ; i++)
	{
	  //	  cout <<  Cl_seed_energy[i] << " " << Cl_x[i] << " " << Cl_y[i] << " ";
	  if( (Cl_x[i] < 0.) && ((*g_flag_emcal.Get()) == 1) )
	    {
	      double g_eD = *g_E.Get() - Cl_Energy_tot_simul[i] / 1000.;
	      double g_poxD = *g_pjt_emcal_x.Get() - Cl_x[i];
	      double g_poyD = *g_pjt_emcal_y.Get() - Cl_y[i];
	      double g_poxD_cor = *g_pjt_emcal_x.Get() - Cl_x_corr[i];
	      double g_poyD_cor = *g_pjt_emcal_y.Get() - Cl_y_corr[i];
	      
	      double ge = *g_E.Get(), g_pjx = *g_pjt_emcal_x.Get(), g_pjy = *g_pjt_emcal_y.Get();
	      //	      double g_poxD = *g_hit_emcal_x.Get() - Cl_seed_x[i];
	      //	      double g_poyD = *g_hit_emcal_y.Get() - Cl_seed_y[i];
	      
	      diffE_g_res_pri->Fill(g_eD);
	      diff_posx_g_res_pri->Fill(g_poxD);
	      diff_posy_g_res_pri->Fill(g_poyD);
	      diff_posx_cor_g_res_pri->Fill(g_poxD_cor);
	      diff_posy_cor_g_res_pri->Fill(g_poyD_cor);
		      
	      diffx_g_E->Fill(g_poxD, g_eD);
	      diffx_g_x->Fill(g_poyD, g_pjy);
	    }
	}
    }

  auto c1 = new TCanvas("c1", "c1", 1000, 1000);
  c1->Divide(2,2);
  c1->cd(1);
  diffE_g_res_pri->SetStats(0);
  diffE_g_res_pri->GetXaxis()->SetTitle("[GeV]");
  diffE_g_res_pri->Draw();
  c1->cd(2);
  diffx_g_x->SetTitle("project Y pos v.s. delta Y");
  diffx_g_x->SetStats(0);
  diffx_g_x->GetXaxis()->SetTitle("project-recons [mm]");
  diffx_g_x->GetYaxis()->SetTitle("ypos [mm]");
  diffx_g_x->Draw("colorz");
  //  diffx_g_E->Draw("colorz");
  c1->cd(3);
  diff_posx_g_res_pri->SetStats(0);
  diff_posx_g_res_pri->GetXaxis()->SetTitle("[mm]");
  diff_posx_g_res_pri->Draw();
  c1->cd(4);
  diff_posy_g_res_pri->SetStats(0);
  diff_posy_g_res_pri->GetXaxis()->SetTitle("[mm]");
  diff_posy_g_res_pri->Draw();


  auto c2 = new TCanvas("c2", "c2", 1000, 500);
  c2->Divide(2,1);
  c2->cd(1);
  diff_posx_cor_g_res_pri->Draw();
  c2->cd(2);
  diff_posy_cor_g_res_pri->Draw();
}
