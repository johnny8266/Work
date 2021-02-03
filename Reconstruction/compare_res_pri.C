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
  double Cl_seed_energy = 0;
  int Cl_seed_npe = 0;
  double Cl_energy = 0;
  double Cl_seed_x = 0;
  double Cl_seed_y = 0;
  double Cl_seed_z = 0;
  double Cl_x = 0;
  double Cl_y = 0;
  double Cl_radius = 0;
  double Cl_theta = 0;
  double Cl_phi = 0;
  int Cl_size = 0;
  int Cl_size_simul = 0;
  double Cl_Energy_tot_simul = 0;
  double g_px = 0., g_py = 0., g_pz = 0., g_E = 0., e_px = 0., e_py = 0., e_pz = 0., e_E = 0.;
  double e_hit_emcal_x = 0., e_hit_emcal_y = 0., e_hit_emcal_z = 0.;
  double g_hit_emcal_x = 0., g_hit_emcal_y = 0., g_hit_emcal_z = 0.;
  int e_flag_emcal=0, g_flag_emcal=0; 

  std::string path = "./data/";
  //  std::string fileName_out = path + "outEnergy_resolution_test.root";
  //  cout << fileName_out << endl;

  //  std::string fileName = path + "outCluster." + std::to_string(ene) + ".root";
  std::string fileName = path + "outCluster.root";
    
  cout << "************ " << fileName << " ************" << endl;

  
  TFile *f1 = new TFile(fileName.c_str(), "Read");
  TTree *outTree = (TTree *) f1->Get("outTree");
  outTree->SetBranchAddress("Cl_seed_energy", &Cl_seed_energy);
  outTree->SetBranchAddress("Cl_seed_npe", &Cl_seed_npe);
  outTree->SetBranchAddress("Cl_seed_x", &Cl_seed_x);
  outTree->SetBranchAddress("Cl_seed_y", &Cl_seed_y);
  outTree->SetBranchAddress("Cl_energy", &Cl_energy);
  outTree->SetBranchAddress("Cl_x", &Cl_x);
  outTree->SetBranchAddress("Cl_y", &Cl_y);
  outTree->SetBranchAddress("Cl_radius", &Cl_radius);
  outTree->SetBranchAddress("Cl_theta", &Cl_theta);
  outTree->SetBranchAddress("Cl_phi", &Cl_phi);
  outTree->SetBranchAddress("Cl_size", &Cl_size);
  outTree->SetBranchAddress("Cl_size_simul", &Cl_size_simul);
  outTree->SetBranchAddress("Cl_Energy_tot_simul", &Cl_Energy_tot_simul);
  outTree->SetBranchAddress("e_hit_emcal_x", &e_hit_emcal_x);
  outTree->SetBranchAddress("e_hit_emcal_y", &e_hit_emcal_y);
  outTree->SetBranchAddress("e_hit_emcal_z", &e_hit_emcal_z);
  outTree->SetBranchAddress("g_hit_emcal_x", &g_hit_emcal_x);
  outTree->SetBranchAddress("g_hit_emcal_y", &g_hit_emcal_y);
  outTree->SetBranchAddress("g_hit_emcal_z", &g_hit_emcal_z);
  outTree->SetBranchAddress("g_px", &g_px);
  outTree->SetBranchAddress("g_py", &g_py);
  outTree->SetBranchAddress("g_pz", &g_pz);
  outTree->SetBranchAddress("g_E", &g_E);
  outTree->SetBranchAddress("e_px", &e_px);
  outTree->SetBranchAddress("e_py", &e_py);
  outTree->SetBranchAddress("e_pz", &e_pz);
  outTree->SetBranchAddress("e_E", &e_E);
  outTree->SetBranchAddress("e_flag_emcal", &e_flag_emcal);
  outTree->SetBranchAddress("g_flag_emcal", &g_flag_emcal);

  auto diffE_e_res_pri = new TH1F("diffE_e_res_pri", "diffE_e_res_pri", 100, -5., 5.);
  auto diffE_g_res_pri = new TH1F("diffE_g_res_pri", "diffE_g_res_pri", 100, -5., 5.);
  auto hit_cl_size = new TH1F("hit_cl_size", "hit_cl_size", 20, 0., 20.);

  auto diffE_e_res_pri_vs_ori = new TH2F("diffE_e_res_pri_vs_ori", "diffE_e_res_pri_vs_ori", 140, -4., 10., 110, 0., 10.);
  auto diffE_g_res_pri_vs_ori = new TH2F("diffE_g_res_pri_vs_ori", "diffE_g_res_pri_vs_ori", 100, -5., 0., 50, 5., 10.);
  auto diffE_e_res_pri_vs_size = new TH2F("diffE_e_res_pri_vs_size", "diffE_e_res_pri_vs_size", 140, -4., 10., 30, 0., 30.);
  auto good_res_e_vs_dis = new TH2F("good_res_e_vs_dis", "good_res_e_vs_dis", 120, 0., 1200, 50, -5., 0.);
  auto bad_res_e_vs_dis = new TH2F("bad_res_e_vs_dis", "bad_res_e_vs_dis", 80, 400., 1200, 100, 0., 10.);
  
  
  cout << "Number of events: " << outTree->GetEntries() << endl << endl;

  int e_g_count = 0, e_count = 0, g_count = 0, bad_e_recons = 0;
  double diffE_e = 0., diffE_g = 0.;
  
  for(int i = 0 ; i < outTree->GetEntries() ; i++)
    {
      outTree->GetEntry(i);

      if( (e_flag_emcal == 1) && (g_flag_emcal == 1) )
	{
	  e_g_count++;
	  continue;
	}

      if( (e_flag_emcal == 1) && (g_flag_emcal == 0) )
	{
	  double Dis_clu_pri = (e_hit_emcal_x - Cl_x) * (e_hit_emcal_x - Cl_x) + (e_hit_emcal_y - Cl_y) * (e_hit_emcal_y - Cl_y);
	  Dis_clu_pri = TMath::Sqrt(Dis_clu_pri);

	  diffE_e = (Cl_Energy_tot_simul / 1000.) - e_E;
	  diffE_e_res_pri->Fill(diffE_e);
	  diffE_e_res_pri_vs_ori->Fill(diffE_e, e_E);
	  diffE_e_res_pri_vs_size->Fill(diffE_e, Cl_size_simul);
	  if( diffE_e < 0. )
	    good_res_e_vs_dis->Fill(Dis_clu_pri, diffE_e);
	  
	  //	  if( Dis_clu_pri > (3. * 20.5) )
	  if( diffE_e > 0. )
	    {
	      cout << e_E << " " << (Cl_Energy_tot_simul / 1000.) << endl;
	      cout << "[" << e_hit_emcal_x << " " << e_hit_emcal_y << "] || " << "[" << Cl_x << " " << Cl_y << "]" << endl;
	      cout << "Pos difference of pri e- and recon e-: " << Dis_clu_pri << endl << endl;
	      bad_res_e_vs_dis->Fill(Dis_clu_pri, diffE_e);
	      bad_e_recons++;
	    }
	  e_count++;
	}

      if( (e_flag_emcal == 0) && (g_flag_emcal == 1) )
	{
	  diffE_g = (Cl_Energy_tot_simul / 1000.) - g_E;
	  diffE_g_res_pri->Fill(diffE_g);
	  diffE_g_res_pri_vs_ori->Fill(diffE_g, g_E);
	  g_count++;
	}
	
    }


  cout << "Both e & g hit on Emcal: " << e_g_count << endl;
  cout << "e hit on Emcal: " << e_count << "& bad electron reconstruction: " << bad_e_recons << endl;
  cout << "g hit on Emcal: " << g_count << endl;

    
  TCanvas *c0 = new TCanvas("c0", "c0", 1000, 1000);
  c0->Divide(2,2);
  c0->cd(1);
  diffE_e_res_pri->SetTitle("Electron reconstruct - primary [energy]");
  diffE_e_res_pri->SetStats(0);
  diffE_e_res_pri->Draw();
  c0->cd(2);
  diffE_g_res_pri->SetTitle("Photon reconstruct - primary [energy]");
  diffE_g_res_pri->SetStats(0);
  diffE_g_res_pri->Draw();
  c0->cd(3);
  diffE_e_res_pri_vs_ori->SetTitle("ELectron delta E v.s. primary E");
  diffE_e_res_pri_vs_ori->SetStats(0);
  diffE_e_res_pri_vs_ori->Draw("colorz");
  c0->cd(4);
  diffE_g_res_pri_vs_ori->SetTitle("Photon delta E v.s. primary E");
  diffE_g_res_pri_vs_ori->SetStats(0);
  diffE_g_res_pri_vs_ori->Draw("colorz");

  TCanvas *c1 = new TCanvas("c1", "c1", 1000, 500);
  c1->Divide(2,1);
  c1->cd(1);
  good_res_e_vs_dis->Draw("colorz");
  c1->cd(2);
  bad_res_e_vs_dis->Draw("colorz");
  //  diffE_e_res_pri_vs_size->SetTitle("ELectron delta E v.s. cluster size");
  //  diffE_e_res_pri_vs_size->SetStats(0);
  //  diffE_e_res_pri_vs_size->Draw("colorz");
  
}
