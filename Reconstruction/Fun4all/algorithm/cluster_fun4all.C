//=================================================================================
// Reconstruction the photon information with the cluster hits in the Calorimeter
//
// 21.01.2021
//
// Pu-Kai
//
//=================================================================================

#include <ROOT/RDataFrame.hxx>
#include <stdio.h>
#include <vector>
#include <TMath.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TH1D.h>
#include <TH2D.h>
#include "TF1.h"
#include <string.h>
#include "TLatex.h"
#include <TTree.h>
#include <vector>
#include <map>
#include <iostream>
#include <fstream>
#include <TGraphErrors.h>
#include <TFile.h>
#include "TLatex.h"

using namespace std;

struct Hit
{
  double x_crs;
  double y_crs;
  double z_crs;
  double Et_dep;     // ce_emcal ETot deposit
  double E_digi;     // ce_emcal ADC
  double time;
  int npe;
  int c_row, c_col;  // crystal row & col
  int c_ID, c_sec;   // crystal ID: 1000*row+col(section:0), 1000000+1000*row+col(section:1)
};

struct Cluster
{
  int num_cluster;

  double C_seed_energy, C_energy, C_seed_x, C_seed_y, C_seed_z;
  int C_size;
};

Cluster ComputeCluster(vector<Hit> hit)
{
  int Size = hit.size(), Clus_size = 0;
  double ClusSeed_xcrs, ClusSeed_ycrs, ClusSeed_zcrs, ClusSeed_Ene;
  double Clus_Etot, Rmoliere = 4., Ethr = 0.01;
  
  Cluster cluster;


  for(int i = 0 ; i < Size ; i++)
    {	
      if (hit.at(i).Et_dep > Ethr && hit.at(i).Et_dep > ClusSeed_Ene)  // iterate the energy to the larger one.
	{
	  ClusSeed_Ene = hit.at(i).Et_dep;
	  ClusSeed_xcrs = hit.at(i).x_crs;
	  ClusSeed_ycrs = hit.at(i).y_crs;
	  ClusSeed_zcrs = hit.at(i).z_crs;
	  //	  cout << " run " << endl;
	}
    }
  //  cout << "Seed ene: " << ClusSeed_Ene << endl;
  
  
  for (int i = 0; i < Size; i++)
    {
      if (hit.at(i).Et_dep > Ethr)
	{
	  double Dx = hit.at(i).x_crs - ClusSeed_xcrs;
	  double Dy = hit.at(i).y_crs - ClusSeed_ycrs;

	  if (sqrt(Dx * Dx + Dy * Dy) <= 3. * Rmoliere)
	    {
	      //   cout <<hit.at(i).E_digi<< " "<<hit.at(i).x_crs<< " "<<hit.at(i).y_crs<<" "<< Dx << " "<< Dy<< " "<<sqrt(Dx*Dx+Dy*Dy)<<endl;
	      Clus_Etot += hit.at(i).Et_dep;
	      Clus_size++;
	    }
	} //end if ethr
    } //for energy tot


  
  Cluster clust;

  clust.C_seed_energy = ClusSeed_Ene;                 //
  clust.C_energy = Clus_Etot;                         //
  clust.C_seed_x = ClusSeed_xcrs;
  clust.C_seed_y = ClusSeed_ycrs;
  clust.C_seed_z = ClusSeed_zcrs;
  clust.C_size = Clus_size;
  
  return clust;
}


//.......oooooooo000000ooooooooOOOOOOOOOoooooooo000000oooooooo.......


void cluster_fun4all()
{
  
  //===================================
  // Open the g4e output file
  //===================================
  
  //  TFile *file = TFile::Open("/vol0/pwang-l/Singularity/my_det/emcal/data/material_test/g4eemc_eval_200_inner_crystal.root");
  //  TFile *file = TFile::Open("/vol0/pwang-l/Singularity/my_det/emcal/data/material_test/g4eemc_eval_Nathaly_glass_4_4_40cm_thickness.root");
  //  TFile *file = TFile::Open("/vol0/pwang-l/Singularity/my_det/emcal/data/radiation_length_test/g4eemc_eval_500_2_5X0.root");

  TFile *file = TFile::Open("/vol0/pwang-l/Singularity/my_det/emcal/data/g4eemc_eval.root");
  
  TTree *ntp_gpoint = (TTree*)file->Get("ntp_gpoint");
  TTree *ntp_gshower = (TTree*)file->Get("ntp_gshower");
  TTree *ntp_tower = (TTree*)file->Get("ntp_tower");
  TTree *ntp_cluster = (TTree*)file->Get("ntp_cluster");

  Int_t n_gpoint = (Int_t)ntp_gpoint->GetEntries();
  Int_t n_gshower = (Int_t)ntp_gshower->GetEntries();
  Int_t n_tower = (Int_t)ntp_tower->GetEntries();
  Int_t n_cluster = (Int_t)ntp_cluster->GetEntries();
  

  double mass_electron = 0.00051099895;
  double mass_proton = 0.0938271998;
  double hit_energy = 0.;
  int count = 0;
  
  
  //==================================
  // Histgram
  //==================================

  auto *C_size = new TH1I("C_size", "C_size", 50, 0, 50);

  auto *reco_E = new TH1F("reco_E", "reco_E", 55, 0, 11);
  auto *real_E = new TH1F("real_E", "real_E", 55, 0, 11);
    
  //==================================
  // Save the output
  //==================================

  /*
  int ene = 0, N_cluster = 0;
  int g_possible_hit_emcal_count = 0;
  double g_px = 0., g_py = 0., g_pz = 0., g_E = 0., e_px = 0., e_py = 0., e_pz = 0., e_E = 0.;
  double e_hit_emcal_x = 0., e_hit_emcal_y = 0., e_hit_emcal_z = 0.;
  double g_hit_emcal_x = 0., g_hit_emcal_y = 0., g_hit_emcal_z = 0.;
  double e_pjt_emcal_x = 0., e_pjt_emcal_y = 0.;
  double g_pjt_emcal_x = 0., g_pjt_emcal_y = 0.;
  int e_flag_emcal = 0, g_flag_emcal = 0, N_hit_emcal = 0;

  vector<double> Cl_seed_energy, Cl_seed_x, Cl_seed_y, Cl_seed_z, Cl_x, Cl_y;
  vector<double> Cl_radius, Cl_theta, Cl_phi;
  vector<double> Cl_Energy_tot_simul, Cl_Energy_pe, Cl_size_simul;
  vector<double> Cl_par_a, Cl_x_corr, Cl_y_corr;

  vector<double> g_pjt_x, g_pjt_y, g_E_all;
  */
  
  TLorentzVector v_g;
  Float_t ge, gpt, geta, gphi, gvx, gvy, gvz;
  Float_t event, clusterID, towerID, x, y, z, e, ntowers, e_tot = 0.;
  int count_single_clu = 0, n_event = 0, count_per_event = 0, sum = 0, success_reco = 0;

  vector<float> hit_x, hit_y, hit_z, hit_e;
  vector<double> reco_e;
  vector<int> event_n;


  
  /*
  string fileName_out = "../data/outCluster.root";
  TFile *fout = new TFile(fileName_out.c_str(), "Recreate");
  TTree *outTree = new TTree("outTree", "outTree");

  outTree->Branch("N_cluster", &N_cluster, "N_cluster/I");

  outTree->Branch("Cl_seed_energy", &Cl_seed_energy);
  outTree->Branch("Cl_seed_x", &Cl_seed_x);
  outTree->Branch("Cl_seed_y", &Cl_seed_y);
  outTree->Branch("Cl_seed_z", &Cl_seed_z);
  outTree->Branch("Cl_x", &Cl_x);
  outTree->Branch("Cl_y", &Cl_y);
  outTree->Branch("Cl_radius", &Cl_radius);
  outTree->Branch("Cl_theta", &Cl_theta);
  outTree->Branch("Cl_phi", &Cl_phi);
  outTree->Branch("Cl_Energy_tot_simul", &Cl_Energy_tot_simul);
  outTree->Branch("Cl_Energy_pe", &Cl_Energy_pe);
  outTree->Branch("Cl_size_simul", &Cl_size_simul);
  outTree->Branch("Cl_par_a", &Cl_par_a);
  outTree->Branch("Cl_x_corr", &Cl_x_corr);
  outTree->Branch("Cl_y_corr", &Cl_y_corr);

  outTree->Branch("g_pjt_x", &g_pjt_x);
  outTree->Branch("g_pjt_y", &g_pjt_y);
  outTree->Branch("g_E_all", &g_E_all);
  
  outTree->Branch("N_hit_emcal", &N_hit_emcal, "N_hit_emcal/I");
  */

  
  
  
  //==================================
  // Loop the Event
  //==================================


  //Load ntuple of ntp_cluster
  //
  ntp_tower->SetBranchAddress("event", &event);
  ntp_tower->SetBranchAddress("towerID", &towerID);
  ntp_tower->SetBranchAddress("x", &x);
  ntp_tower->SetBranchAddress("y", &y);
  ntp_tower->SetBranchAddress("z", &z);
  ntp_tower->SetBranchAddress("e", &e);

  cout << "N towers are hit: " << n_tower << endl;

  Cluster clust;
  
  for(int i = 0 ; i < n_tower ; i++)
    {

      ntp_tower->GetEntry(i);
      //      cout << i << " || " << clusterID << " " << event << ": [" << x << ", " << y << ", " << z << "]" << endl;
      //      cout << event << endl;
      if(event == n_event)
	{
	  count_per_event++;
	  hit_x.push_back(x);
	  hit_y.push_back(y);
	  hit_z.push_back(z);
	  hit_e.push_back(e);
	  e_tot = e_tot + e;
	}
      //      else if( (event != n_event) || ( i == (n_tower - 1)) )
      else
	{
	  vector<Hit> hhit;
	  
	  for(int j = 0 ; j < count_per_event ; j++)
	    {
	      Hit hit;
	      hit.x_crs = hit_x[j];
	      hit.y_crs = hit_y[j];
	      hit.z_crs = hit_z[j];
	      hit.Et_dep = hit_e[j];
	  
	      hhit.push_back(hit);
	    }

	  clust = ComputeCluster(hhit);

	  double CSE = clust.C_seed_energy, CE = clust.C_energy;
	  double CSx = clust.C_seed_x, CSy = clust.C_seed_y, CSz = clust.C_seed_z;
	  int CS = clust.C_size;

	  if( CE > 0.1 )
	    {
	      reco_e.push_back(CE);
	      event_n.push_back(n_event);
	      success_reco++;
	      C_size->Fill(CS);
	    }

	  cout << "etot: " << e_tot << endl;
	  e_tot = 0.;
	  /*
	  cout << n_event << " has " << count_per_event << " hits....." << endl;
	  cout << CSE << " " << CE << " " << CS << " || " ;
	  cout << CSx << ", " << CSy << ", " << CSz << endl << endl; 
	  */
	  
	  hit_x.clear();  hit_y.clear();  hit_z.clear();  hit_e.clear();
	  sum = sum + count_per_event;
	  count_per_event = 0;
	  i = i - 1;
	  n_event++;

	}
            

      
    } //End of the for loop
  cout << success_reco << endl;

  //  cout << endl << sum << endl;



  ntp_gshower->SetBranchAddress("event",&event);
  ntp_gshower->SetBranchAddress("gvx",&gvx);
  ntp_gshower->SetBranchAddress("gvy",&gvy);
  ntp_gshower->SetBranchAddress("gvz",&gvz);
  ntp_gshower->SetBranchAddress("gpt",&gpt);
  ntp_gshower->SetBranchAddress("geta",&geta);
  ntp_gshower->SetBranchAddress("gphi",&gphi);
  ntp_gshower->SetBranchAddress("ge",&ge);
  
  for(int i = 0 ; i < success_reco ; i++)
    {
      int se = event_n[i];
      //      cout << se << endl;
      ntp_gshower->GetEntry(se);

      reco_E->Fill(reco_e[i]);
      real_E->Fill(ge);
      cout << se << ": " << reco_e[i] << ", " << ge << endl;
    }


  TCanvas *c1 = new TCanvas("c1", "c1", 500, 500);
  C_size->SetTitle("Glass cluster size[n hits]");
  C_size->SetStats(0);
  C_size->Draw();
  
  TCanvas *c2 = new TCanvas("c2", "c2", 1000, 500);
  c2->Divide(2,1);
  c2->cd(1);
  reco_E->SetTitle("glass_Nathaly_4_4_40cm_thickness reco E");
  reco_E->SetStats(0);
  reco_E->Draw();
  c2->cd(2);
  real_E->SetTitle("glass_Nathaly_4_4_40cm_thickness real E");
  real_E->SetStats(0);
  real_E->Draw();
  
}
           



