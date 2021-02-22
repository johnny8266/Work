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

struct Hit {
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
  vector<double> C_seed_energy, C_seed_x, C_seed_y, C_seed_z, C_x, C_y;
  vector<double> C_energy, C_size;  // these 2 variables were not used.
  vector<double> C_radius, C_theta, C_phi;
  vector<double> C_Energy_tot_simul, C_size_simul;
};

Cluster ComputeCluster(vector<Hit> hit)
{
  int Size = hit.size(), seed_total=0, seed_cry=0, seed_gla=0;
  double CRYS_ZPOS = 2110; // ECAL zpos from the center of EIC detector
  double Ethr = 10.; // threshold[MeV]
  double Rmoliere = 20.01; // in mm    ->Rmolier for PbWO is 20 mm
  int Clus_size_simul = 0;
  double Clus_Energy_tot_simul = 0, Clus_Etot = 0.;
  double Clus_x = 0., Clus_y = 0., Clus_xx = 0., Clus_yy = 0.;
  double Clus_sigmaX = 0., Clus_sigmaY = 0.;
  double Clus_Radius = 0, Clus_Theta = 0, Clus_phi = 0; //in deg;

  vector<double> ClusSeed_xcrs, ClusSeed_ycrs, ClusSeed_zcrs, ClusEtotsimu;
  map<int, double> map_crystal;  // ID, row, col, E
  Cluster cluster;

  //    cout << "=============================================" << endl;
  //    cout << "Ce_Emcal hit size of this event: " << Size << endl;


  //========================================================
  // Save all hit into map form and find the multiple seed
  //========================================================
  for(int i = 0 ; i < Size ; i++)
    map_crystal.insert(pair<int, double>(hit.at(i).c_ID, hit.at(i).Et_dep));
  

  //  cout << "This event has " << Size << " hits and...." << endl;
  for(int i = 0 ; i < Size ; i++)
    {
      int count = 0;
      int i_row = hit.at(i).c_row, i_col = hit.at(i).c_col;
      double i_xcrs = hit.at(i).x_crs, i_ycrs = hit.at(i).y_crs, i_zcrs = hit.at(i).z_crs;
      double i_et_dep = hit.at(i).Et_dep;
      double i_et_digi = hit.at(i).E_digi;
      
      if( hit.at(i).c_sec == 0 )    //For crystal region
	{
	  //	  cout << "initial: [" << i_row << ", " << i_col << "]" << endl;	  
	  for(int j = -1 ; j < 2 ; j++)
	    {
	      for(int k = -1 ; k < 2 ; k++)
		{
		  if( (j == 0) && (k == 0) )
		    continue;
		  i_row = hit.at(i).c_row;  i_row = i_row + j;
		  i_col = hit.at(i).c_col;  i_col = i_col + k;
		  //		  cout << i_row << " " << i_col << endl;
		  int i_ID = i_row * 1000 + i_col;
		  auto iter = map_crystal.find(i_ID);
		  if( iter != map_crystal.end() )
		    if( i_et_dep > (iter->second) )
		      count++;
		} //scan col
	    } //scan row
	  if( count == 8 && i_et_dep > 200. )
	    {
	      cluster.C_seed_energy.push_back(i_et_dep);
	      cluster.C_seed_x.push_back(i_xcrs);  cluster.C_seed_y.push_back(i_ycrs);  cluster.C_seed_z.push_back(i_zcrs);

	      ClusSeed_xcrs.push_back(i_xcrs);  ClusSeed_ycrs.push_back(i_ycrs);  ClusSeed_zcrs.push_back(i_zcrs); 
	      //	      cout << hit.at(i).c_col << " " << hit.at(i).c_row << " count: " << count << endl;
	      //	      cout << "crystal seed: " << i_et_digi << "  " << i_et_dep << " [" << i_xcrs << " " << i_ycrs << "] || row, col: [" << hit.at(i).c_col << " " << hit.at(i).c_row << "]" << endl;
	      seed_cry++;
	    }
	} //crystal section
      
      else    //For glass region
	{
	  //	  cout << "initial: [" << i_row << ", " << i_col << "]" << endl;	  
	  for(int j = -1 ; j < 2 ; j++)
	    {
	      for(int k = -1 ; k < 2 ; k++)
		{
		  if( (j == 0) && (k == 0) )
		    continue;
		  i_row = hit.at(i).c_row;  i_row = i_row + j;
		  i_col = hit.at(i).c_col;  i_col = i_col + k;
		  //		  cout << i_row << " " << i_col << endl;
		  //		  cout << hit.at(i)
		  int i_ID = i_row * 1000 + i_col + 1000000;
		  auto iter = map_crystal.find(i_ID);
		  if( iter != map_crystal.end() )
		    if( i_et_dep > (iter->second) )
		      count++;
		} //scan col
	    } //scan row
	  if( count == 8 && i_et_dep > 200. )
	    {
	      cluster.C_seed_energy.push_back(i_et_dep);
	      cluster.C_seed_x.push_back(i_xcrs); cluster.C_seed_y.push_back(i_ycrs); cluster.C_seed_z.push_back(i_zcrs);

	      ClusSeed_xcrs.push_back(i_xcrs);  ClusSeed_ycrs.push_back(i_ycrs);  ClusSeed_zcrs.push_back(i_zcrs);
	      //	      cout << hit.at(i).c_col << " " << hit.at(i).c_row << " count: " << count << endl;
	      //	      cout << "glass seed: " << i_et_digi << "  " << i_et_dep << " [" << i_xcrs << " " << i_ycrs << "] || row, col: [" << hit.at(i).c_col << " " << hit.at(i).c_row << "]" << endl;
	      seed_gla++;
	    }
	} //glass section
    } //loop ce_emcal hits

  seed_total = seed_gla + seed_cry;
  cluster.num_cluster = seed_total;
  
  //  if(seed_total >= 1 )
  //    cout << "Seed counts: " << seed_total << ", glass: " << seed_gla << ", crystal: " << seed_cry << endl;
  
  //
  //build the cluster energy and size
  //
  
  for(int i = 0 ; i < seed_total ; i++)
    {
      Clus_Energy_tot_simul = 0.;  Clus_size_simul = 0;
      
      for(int k = 0 ; k < Size ; k++)
	{
	  double Dx = hit.at(k).x_crs - ClusSeed_xcrs[i];
	  double Dy = hit.at(k).y_crs - ClusSeed_ycrs[i];

	  if (sqrt(Dx * Dx + Dy * Dy) <= 3. * Rmoliere)
	    {
	      Clus_Energy_tot_simul += hit.at(k).Et_dep;
	      Clus_size_simul++;
	    }
	}
      //      cout << "Reconstructed energy: " << Clus_Energy_tot_simul << " || cluster size: " << Clus_size_simul << endl;
      cluster.C_Energy_tot_simul.push_back(Clus_Energy_tot_simul);  ClusEtotsimu.push_back(Clus_Energy_tot_simul);
      cluster.C_size_simul.push_back(Clus_size_simul);
    }


  double w_tot, x, y;
 
  for(int i = 0 ; i < seed_total ; i++ )
    {
      w_tot = 0.;  x = 0.;  y = 0.;
      Clus_Etot = ClusEtotsimu[i];
      Clus_x = 0.; Clus_xx = 0.;
      Clus_y = 0.; Clus_yy = 0.;
      
      for(int k = 0 ; k < Size ; k++)
	{
	  double Dx = hit.at(k).x_crs - ClusSeed_xcrs[i];
	  double Dy = hit.at(k).y_crs - ClusSeed_ycrs[i];

	  if (sqrt(Dx * Dx + Dy * Dy) <= 3. * Rmoliere)
	    {
	      double w1 = std::max(0., (3.45 + std::log(hit.at(k).E_digi / Clus_Etot))); 
	      x += w1 * hit.at(k).x_crs;
	      y += w1 * hit.at(k).y_crs;
	      Clus_xx += w1 * hit.at(k).x_crs * hit.at(k).x_crs;
	      Clus_yy += w1 * hit.at(k).y_crs * hit.at(k).y_crs;
	      w_tot += w1;
	    }
	}
      Clus_x = x / w_tot;  Clus_y = y / w_tot;
      cluster.C_x.push_back(Clus_x);  cluster.C_y.push_back(Clus_y);

      Clus_xx /= w_tot;  Clus_yy /= w_tot;
      
      double sigmax2 = Clus_xx - std::pow(Clus_x, 2.);  // <x^2> - <x>^2
      if (sigmax2 < 0) sigmax2 = 0;
      Clus_sigmaX = std::sqrt(sigmax2);

      double sigmay2 = Clus_yy - std::pow(Clus_y, 2.);
      if (sigmay2 < 0) sigmay2 = 0;
      Clus_sigmaY = std::sqrt(sigmay2);

      //Cluster radius
      double radius2 = (sigmax2 + sigmay2);
      if (radius2 < 0) radius2 = 0;
      Clus_Radius = std::sqrt(radius2);
      cluster.C_radius.push_back(Clus_Radius);

      //Cluster theta
      Clus_Theta = (std::atan((std::sqrt(std::pow(Clus_x, 2.) + std::pow(Clus_y, 2.))) / (CRYS_ZPOS + ClusSeed_zcrs[i]))) * (180. / M_PI);
      cluster.C_theta.push_back(Clus_Theta);

      //Cluster phi
      Clus_phi = std::atan2(Clus_x, Clus_y) * (180. / M_PI);
      cluster.C_phi.push_back(Clus_phi);

      //      cout << "Reconstruction pos: [" << Clus_x << " " << Clus_y << "]  ";
      //      cout << "Sigmax: " << Clus_sigmaX << " " << "Sigmay: " << Clus_sigmaY << endl;
      //      cout << "Radius: " << Clus_Radius << " " << "theta: " << Clus_Theta << " " << "phi: " << Clus_phi << endl; 
    }
  //  cout << endl << endl;

  return cluster;
}


//.......oooooooo000000ooooooooOOOOOOOOOoooooooo000000oooooooo.......


void multi_cluster_recons()
{
  
  //===================================
  // Open the g4e output file
  //===================================
  
  //  TFile *file = TFile::Open("../Data/g4e_simulation/g4e_output_10k_events_crossing_angle.root");
  TFile *file = TFile::Open("../Data/g4e_simulation/g4e_output_foam_imposed_10k_events.root");
  TTree *events = (TTree *) file->Get("events");

  TTreeReader fReader("events", file);

  // Readers to access the data (delete the ones you do not need).
  TTreeReaderValue<ULong64_t>     event_id = {fReader, "event_id"};

  // Hits collections
  TTreeReaderValue<ULong64_t>     hit_count = {fReader, "hit_count"};
  TTreeReaderArray<unsigned long> hit_id = {fReader, "hit_id"};
  TTreeReaderArray<unsigned long> hit_trk_id = {fReader, "hit_trk_id"};
  TTreeReaderArray<unsigned long> hit_ptr_id = {fReader, "hit_ptr_id"};
  TTreeReaderArray<unsigned long> hit_parent_trk_id = {fReader, "hit_parent_trk_id"};  // PDG code of that particle
  TTreeReaderArray<string>        hit_vol_name = {fReader, "hit_vol_name"};
  TTreeReaderArray<double>        hit_x = {fReader, "hit_x"};
  TTreeReaderArray<double>        hit_y = {fReader, "hit_y"};
  TTreeReaderArray<double>        hit_z = {fReader, "hit_z"};
  TTreeReaderArray<double>        hit_e_loss = {fReader, "hit_e_loss"};

  // Tracks data
  TTreeReaderValue<ULong64_t>     trk_count = {fReader, "trk_count"};
  TTreeReaderArray<unsigned long> trk_id = {fReader, "trk_id"};
  TTreeReaderArray<long>          trk_pdg = {fReader, "trk_pdg"};
  TTreeReaderArray<unsigned long> trk_parent_id = {fReader, "trk_parent_id"};
  TTreeReaderArray<long>          trk_create_proc = {fReader, "trk_create_proc"};
  TTreeReaderArray<unsigned long> trk_level = {fReader, "trk_level"};
  TTreeReaderArray<double>        trk_vtx_x = {fReader, "trk_vtx_x"};
  TTreeReaderArray<double>        trk_vtx_y = {fReader, "trk_vtx_y"};
  TTreeReaderArray<double>        trk_vtx_z = {fReader, "trk_vtx_z"};
  TTreeReaderArray<double>        trk_vtx_dir_x = {fReader, "trk_vtx_dir_x"};
  TTreeReaderArray<double>        trk_vtx_dir_y = {fReader, "trk_vtx_dir_y"};
  TTreeReaderArray<double>        trk_vtx_dir_z = {fReader, "trk_vtx_dir_z"};
  TTreeReaderArray<double>        trk_mom = {fReader, "trk_mom"};

  // 'Copy' of generated particle data
  TTreeReaderValue<ULong64_t>     gen_prt_count = {fReader, "gen_prt_count"};
  TTreeReaderArray<unsigned long> gen_prt_id = {fReader, "gen_prt_id"};
  TTreeReaderArray<unsigned long> gen_prt_vtx_id = {fReader, "gen_prt_vtx_id"};
  TTreeReaderArray<unsigned long> gen_prt_pdg = {fReader, "gen_prt_pdg"};
  TTreeReaderArray<unsigned long> gen_prt_trk_id = {fReader, "gen_prt_trk_id"};
  TTreeReaderArray<double>        gen_prt_charge = {fReader, "gen_prt_charge"};
  TTreeReaderArray<double>        gen_prt_dir_x = {fReader, "gen_prt_dir_x"};
  TTreeReaderArray<double>        gen_prt_dir_y = {fReader, "gen_prt_dir_y"};
  TTreeReaderArray<double>        gen_prt_dir_z = {fReader, "gen_prt_dir_z"};
  TTreeReaderArray<double>        gen_prt_tot_mom = {fReader, "gen_prt_tot_mom"};
  TTreeReaderArray<double>        gen_prt_tot_e = {fReader, "gen_prt_tot_e"};
  TTreeReaderArray<double>        gen_prt_time = {fReader, "gen_prt_time"};
  TTreeReaderValue<ULong64_t>     gen_vtx_count = {fReader, "gen_vtx_count"};
  TTreeReaderArray<unsigned long> gen_vtx_id = {fReader, "gen_vtx_id"};
  TTreeReaderArray<unsigned long> gen_vtx_part_count = {fReader, "gen_vtx_part_count"};
  TTreeReaderArray<double>        gen_vtx_x = {fReader, "gen_vtx_x"};
  TTreeReaderArray<double>        gen_vtx_y = {fReader, "gen_vtx_y"};
  TTreeReaderArray<double>        gen_vtx_z = {fReader, "gen_vtx_z"};

  // The information of the emcal
  TTreeReaderArray<int>           ce_emcal_npe = {fReader, "ce_emcal_npe"};
  TTreeReaderArray<int>           ce_emcal_id = {fReader, "ce_emcal_id"};
  TTreeReaderArray<int>           ce_emcal_row = {fReader, "ce_emcal_row"};
  TTreeReaderArray<int>           ce_emcal_col = {fReader, "ce_emcal_col"};
  TTreeReaderArray<int>           ce_emcal_section = {fReader, "ce_emcal_section"};
  TTreeReaderArray<double>        ce_emcal_etot_dep = {fReader, "ce_emcal_etot_dep"};
  TTreeReaderArray<double>        ce_emcal_adc = {fReader, "ce_emcal_adc"};
  TTreeReaderArray<double>        ce_emcal_tdc = {fReader, "ce_emcal_tdc"};
  TTreeReaderArray<double>        ce_emcal_xcrs = {fReader, "ce_emcal_xcrs"};  // x, y, z are translation value relative to the mother frame
  TTreeReaderArray<double>        ce_emcal_ycrs = {fReader, "ce_emcal_ycrs"};
  TTreeReaderArray<double>        ce_emcal_zcrs = {fReader, "ce_emcal_zcrs"};

  double mass_electron = 0.00051099895;
  double mass_proton = 0.0938271998;
  double hit_energy = 0.;
  int count = 0;
  
  
  //==================================
  // Histgram
  //==================================

  auto eloss_of_bad_res = new TH1F("eloss_of_bad_res", "eloss_of_bad_res", 100, 0., 0.1);
  auto emcal_e_of_bad_res = new TH1F("emcal_e_of_bad_res", "emcal_e_of_bad_res", 100, 0., 10000.);
  
  auto emcal_hit_xy_photon = new TH2F("emcal_hit_xy_photon", "emcal_hit_xy_photon", 300, -1500., 1500., 300, -1500., 1500.);
  auto emcal_hit_xy_electron = new TH2F("emcal_hit_xy_electron", "emcal_hit_xy_electron", 300, -1500., 1500., 300, -1500., 1500.);
  auto e_res_pri_xy_pos = new TH2F("e_res_pri_xy_pos", "e_res_pri_xy_pos", 300, -1500., 1500., 300, -1500., 1500.);
  auto e_res_bad_xy_pos = new TH2F("e_res_bad_xy_pos", "e_res_bad_xy_pos", 300, -1500., 1500., 300, -1500., 1500.);
  auto hit_row_col_cry = new TH2I("hit_row_col_cry", "hit_row_col_cry", 100, 0, 100, 100, 0, 100);
  auto hit_row_col_gla = new TH2I("hit_row_col_gla", "hit_row_col_gla", 100, 0, 100, 100, 0, 100);
  auto hit_pos_crystal = new TH2F("hit_pos_crystal", "hit_pos_crystal", 140, -1400, 1400, 140, -1400, 1400);
  auto hit_pos_glass = new TH2F("hit_pos_glass", "hit_pos_glass", 140, -1400, 1400, 140, -1400, 1400);

  
  //==================================
  // Save the output
  //==================================

  int ene = 0, N_cluster = 0;
  double g_px = 0., g_py = 0., g_pz = 0., g_E = 0., e_px = 0., e_py = 0., e_pz = 0., e_E = 0.;
  double e_hit_emcal_x = 0., e_hit_emcal_y = 0., e_hit_emcal_z = 0.;
  double g_hit_emcal_x = 0., g_hit_emcal_y = 0., g_hit_emcal_z = 0.;
  double e_pjt_emcal_x = 0., e_pjt_emcal_y = 0.;
  double g_pjt_emcal_x = 0., g_pjt_emcal_y = 0.;
  int e_flag_emcal = 0, g_flag_emcal = 0, N_hit_emcal = 0;

  vector<double> Cl_seed_energy, Cl_seed_x, Cl_seed_y, Cl_seed_z, Cl_x, Cl_y;
  vector<double> Cl_radius, Cl_theta, Cl_phi;
  vector<double> Cl_Energy_tot_simul, Cl_size_simul;

  
  string fileName_out = "./data/outCluster.root";
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
  outTree->Branch("Cl_size_simul", &Cl_size_simul);
  
  outTree->Branch("e_hit_emcal_x", &e_hit_emcal_x, "e_hit_emcal_x/D");
  outTree->Branch("e_hit_emcal_y", &e_hit_emcal_y, "e_hit_emcal_y/D");
  outTree->Branch("e_hit_emcal_z", &e_hit_emcal_z, "e_hit_emcal_z/D");
  outTree->Branch("g_hit_emcal_x", &g_hit_emcal_x, "g_hit_emcal_x/D");
  outTree->Branch("g_hit_emcal_y", &g_hit_emcal_y, "g_hit_emcal_y/D");
  outTree->Branch("g_hit_emcal_z", &g_hit_emcal_z, "g_hit_emcal_z/D");
  outTree->Branch("e_pjt_emcal_x", &e_pjt_emcal_x, "e_pjt_emcal_x/D");
  outTree->Branch("e_pjt_emcal_y", &e_pjt_emcal_y, "e_pjt_emcal_y/D");
  outTree->Branch("g_pjt_emcal_x", &g_pjt_emcal_x, "g_pjt_emcal_x/D");
  outTree->Branch("g_pjt_emcal_y", &g_pjt_emcal_y, "g_pjt_emcal_y/D");
  outTree->Branch("g_px", &g_px, "g_px/D");
  outTree->Branch("g_py", &g_py, "g_py/D");
  outTree->Branch("g_pz", &g_pz, "g_pz/D");
  outTree->Branch("g_E", &g_E, "g_E/D");
  outTree->Branch("e_px", &e_px, "e_px/D");
  outTree->Branch("e_py", &e_py, "e_py/D");
  outTree->Branch("e_pz", &e_pz, "e_pz/D");
  outTree->Branch("e_E", &e_E, "e_E/D");
  outTree->Branch("e_flag_emcal", &e_flag_emcal, "e_flag_emcal/I");
  outTree->Branch("g_flag_emcal", &g_flag_emcal, "g_flag_emcal/I");
  outTree->Branch("N_hit_emcal", &N_hit_emcal, "N_hit_emcal/I");


  
  //==================================
  // Loop the Event
  //==================================

  size_t events_numer = 0;  
  while (fReader.Next())
    {
      if(++events_numer != 89)
      	continue;
      
      if(++events_numer > 200)
	break;
	//	continue;
      
      if(events_numer%100 == 0)
	cout << "Read " << events_numer << " th events..." << endl;
      
      std::unordered_set<uint64_t> track_ids_in_ecap_emcal;  // Get tracks information that have hits in ion EMCAL
      std::unordered_set<uint64_t> track_ids_in_ffi_RPOTS;  // Get tracks information that have hits in Roman Pots(FarForward ion direction area)
      
    
      // Read basic values
      auto hits_count = static_cast<size_t>(*hit_count.Get());         
      auto tracks_count = static_cast<size_t>(*trk_count.Get());       
      //      cout << endl << "This event has: " << tracks_count << " tracks || " << hits_count << " hits." << endl << endl;

      g_hit_emcal_x = 0.;  g_hit_emcal_y = 0.;  g_hit_emcal_z = 0.;
      e_hit_emcal_x = 0.;  e_hit_emcal_y = 0.;  e_hit_emcal_z = 0.;

      // =============================
      // Save the track hit on Emcal
      // =============================
      for(size_t i = 0 ; i < hits_count ; i++)
	{
	  // This is is of a track that made this hit
	  uint64_t hit_track_id = static_cast<uint64_t>(hit_trk_id[i]);
	  uint64_t hit_parent_track_id = static_cast<uint64_t>(hit_parent_trk_id[i]);
	  
	  // This is a volume name of a hit
	  std::string vol_name = static_cast<std::string>(hit_vol_name[i]);
		  
	  double x = hit_x[i], y = hit_y[i], z = hit_z[i];
	  double e_loss = hit_e_loss[i];
	 
        
	  // Check that the name starts with "ce_EMCAL"
	  if(vol_name.rfind("ce_EMCAL", 0) == 0)
	    {
	      //	      if( hit_parent_track_id == 22 && hit_track_id == 2 )
	      //	      cout << hit_parent_track_id << " " << hit_track_id << " " << hit_vol_name[i] << " " << "[" << x << ", " << y << ", " << z << "] " << endl;
	      //	      cout << hit_e_loss[i] * 1000. << "[MeV]" << endl;
	      //	      hit_energy+=hit_e_loss[i];
		      
	      
	      track_ids_in_ecap_emcal.insert(hit_track_id);

	      if( hit_track_id == 2 )
		{
		  emcal_hit_xy_photon->Fill(x, y);
		  if( g_flag_emcal == 0 )
		    {
		      //		      cout << x << " " << y << " " << z << endl;
		      g_hit_emcal_x = x;
		      g_hit_emcal_y = y;
		      g_hit_emcal_z = z;
		      g_flag_emcal = 1;
		    }
		}
	      if( hit_track_id == 1 )
		{
		  emcal_hit_xy_electron->Fill(x, y);
		  if( e_flag_emcal == 0 )
		    {
		      //		      cout << x << " " << y << " " << z << endl;
		      e_hit_emcal_x = x;
		      e_hit_emcal_y = y;
		      e_hit_emcal_z = z;
		      e_flag_emcal = 1;
		    }
		}
	      count++;
	    }
	}

      //      cout << "e flag: " << e_flag_emcal << " || g flag: " << g_flag_emcal << endl;
      //      cout << "hit counts in emcal: " << count << " || energy loss: " << hit_energy * 1000. << endl << endl;
      
      
      // =============================
      // Calculate the hit cluster information
      // =============================

      vector<Hit> hhit;
      double hit_e_check = 0.;
      
      N_hit_emcal = ce_emcal_etot_dep.GetSize();
      //      cout << endl << N_hit_emcal << ": " << endl;
      
      for(int i = 0 ; i < N_hit_emcal ; i++)
	{
	  Hit hit;
	  hit.x_crs = ce_emcal_xcrs[i];
	  hit.y_crs = ce_emcal_ycrs[i];
	  hit.z_crs = ce_emcal_zcrs[i];
	  hit.Et_dep = ce_emcal_etot_dep[i];
	  hit.E_digi = ce_emcal_adc[i];
	  hit.time = ce_emcal_tdc[i];
	  hit.npe = ce_emcal_npe[i];
	  hit.c_row = ce_emcal_row[i];
	  hit.c_col = ce_emcal_col[i];
	  hit.c_ID = ce_emcal_id[i];
	  hit.c_sec = ce_emcal_section[i];
	  
	  hhit.push_back(hit);
	  /*
	  cout << ce_emcal_id[i] << " ";
          cout << ce_emcal_row[i] << " " << ce_emcal_col[i] << " " << ce_emcal_section[i] << " ";
	  cout << ce_emcal_etot_dep[i] << endl;
          cout << "[" << ce_emcal_xcrs[i] << " " << ce_emcal_ycrs[i] << " " << ce_emcal_zcrs[i] << "]" << endl;
	  */
	  //	  if( ce_emcal_section[i] == 0 )
	  //	    hit_row_col_cry->Fill(ce_emcal_row[i], ce_emcal_col[i]);


	  if(ce_emcal_section[i] == 0)
	    {
	      hit_pos_crystal->Fill(ce_emcal_xcrs[i], ce_emcal_ycrs[i], ce_emcal_etot_dep[i]);
	      //	      hit_pos_crystal->Fill(ce_emcal_xcrs[i], ce_emcal_ycrs[i]);
	      hit_row_col_cry->Fill(ce_emcal_row[i], ce_emcal_col[i], ce_emcal_etot_dep[i]);
	      //	      cout << ce_emcal_id[i] << " " << ce_emcal_row[i] << " " << ce_emcal_col[i] << endl;
	      //	      cout << ce_emcal_id[i] << " [" << ce_emcal_xcrs[i] << " " << ce_emcal_ycrs[i] << " " << ce_emcal_zcrs[i] << "]" << endl;
	    }
	  else
	    {
	      hit_pos_glass->Fill(ce_emcal_xcrs[i], ce_emcal_ycrs[i], ce_emcal_etot_dep[i]);
	      //	      hit_pos_glass->Fill(ce_emcal_xcrs[i], ce_emcal_ycrs[i]);
	      hit_row_col_gla->Fill(ce_emcal_row[i], ce_emcal_col[i], ce_emcal_etot_dep[i]);
	      //	      cout << ce_emcal_id[i] << " " << ce_emcal_row[i] << " " << ce_emcal_col[i] << endl;
	      //	      cout << ce_emcal_id[i] << " [" << ce_emcal_xcrs[i] << " " << ce_emcal_ycrs[i] << " " << ce_emcal_zcrs[i] << "]" << endl;
	    }
	}

      g_px = gen_prt_dir_x[1] * gen_prt_tot_mom[1];
      g_py = gen_prt_dir_y[1] * gen_prt_tot_mom[1];
      g_pz = gen_prt_dir_z[1] * gen_prt_tot_mom[1];
      g_E = gen_prt_tot_e[1];

      e_px = gen_prt_dir_x[0] * gen_prt_tot_mom[0];
      e_py = gen_prt_dir_y[0] * gen_prt_tot_mom[0];
      e_pz = gen_prt_dir_z[0] * gen_prt_tot_mom[0];
      e_E = gen_prt_tot_e[0];
      
      Cluster cluster;

      /*
      cout << "e- hit: {" << e_flag_emcal << ", " << e_E << "} || g hit: {" << g_flag_emcal << ", " << g_E << "}" << endl;
      //      cout << "e- direction: " << gen_prt_dir_x[0] << " " << gen_prt_dir_y[0] << " " << gen_prt_dir_z[0] << endl;
      cout << "e- project pos: [" << 2240 * gen_prt_dir_x[0] / TMath::Abs(gen_prt_dir_z[0]) << ", " << 2240 * gen_prt_dir_y[0] / TMath::Abs(gen_prt_dir_z[0]) << "] || ";
      cout << "g project pos: [" << 2240 * gen_prt_dir_x[1] / TMath::Abs(gen_prt_dir_z[1]) << ", " << 2240 * gen_prt_dir_y[1] / TMath::Abs(gen_prt_dir_z[1]) << "]" << endl << endl; 
      */
      e_pjt_emcal_x = 2110. * gen_prt_dir_x[0] / TMath::Abs(gen_prt_dir_z[0]);
      e_pjt_emcal_y = 2110. * gen_prt_dir_y[0] / TMath::Abs(gen_prt_dir_z[0]);
      g_pjt_emcal_x = 2110. * gen_prt_dir_x[1] / TMath::Abs(gen_prt_dir_z[1]);
      g_pjt_emcal_y = 2110. * gen_prt_dir_y[1] / TMath::Abs(gen_prt_dir_z[1]);
      
      cluster = ComputeCluster(hhit);

      N_cluster = cluster.num_cluster;
      //      cout <<  N_cluster << endl;
      if(N_cluster > 0)
	{
	  for(int i = 0 ; i < N_cluster ; i++)
	    {
	      //	  cout << "N cluster: " << N_cluster << endl << endl;
	  
	      double cce = cluster.C_seed_energy[i];
	      double ccx = cluster.C_seed_x[i], ccy = cluster.C_seed_y[i], ccz = cluster.C_seed_z[i];
	      double cx = cluster.C_x[i], cy = cluster.C_y[i];
	      double cR = cluster.C_radius[i], cT = cluster.C_theta[i], cP = cluster.C_phi[i]; 
	      double cets = cluster.C_Energy_tot_simul[i], css = cluster.C_size_simul[i];
	  
	      Cl_seed_energy.push_back(cce);
	      Cl_seed_x.push_back(ccx);
	      Cl_seed_y.push_back(ccy);
	      Cl_seed_z.push_back(ccz);
	      Cl_x.push_back(cx);
	      Cl_y.push_back(cy);
	      Cl_radius.push_back(cR);
	      Cl_theta.push_back(cT);
	      Cl_phi.push_back(cP);
	      Cl_Energy_tot_simul.push_back(cets);
	      Cl_size_simul.push_back(css);
	  

	      //	  cout << cluster.C_seed_x[i] << endl;
	    }
	}
      else
	{
	  Cl_seed_energy.clear();  Cl_seed_x.clear();  Cl_seed_y.clear();  Cl_seed_z.clear();
	  Cl_x.clear();  Cl_y.clear();
	  Cl_radius.clear();  Cl_theta.clear();  Cl_phi.clear();
	  Cl_Energy_tot_simul.clear();  Cl_size_simul.clear();
	  cout << events_numer << " th have " << N_hit_emcal << ".........." << endl;
		
	  continue;
	}
      
      outTree->Fill();

      Cl_seed_energy.clear();  Cl_seed_x.clear();  Cl_seed_y.clear();  Cl_seed_z.clear();
      Cl_x.clear();  Cl_y.clear();
      Cl_radius.clear();  Cl_theta.clear();  Cl_phi.clear();
      Cl_Energy_tot_simul.clear();  Cl_size_simul.clear();


      
      g_flag_emcal = 0;  e_flag_emcal = 0;  hit_energy = 0.;
      
    }// end for the read g4e_output
  fout->Write();
  fout->Close();

  
  
  // ======================================
  // Draw the Results
  // ======================================

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

  /*
  auto *c1 = new TCanvas("c1", "c1", 1600, 800);
  c1->Divide(2, 1);
  c1->cd(1);
  emcal_hit_xy_electron->SetTitle("Primary electron hits on the EMcal");
  emcal_hit_xy_electron->SetStats(0);
  emcal_hit_xy_electron->SetContour(nb);
  emcal_hit_xy_electron->Draw("colorz");
  c1->cd(2);
  emcal_hit_xy_photon->SetTitle("Primary photon hits on the EMcal");
  emcal_hit_xy_photon->SetStats(0);
  emcal_hit_xy_photon->SetContour(nb);
  emcal_hit_xy_photon->Draw("colorz");
  
  auto *c2 = new TCanvas("c2", "c2", 1600, 800);
  c2->Divide(2, 1);
  c2->cd(1);
  e_res_pri_xy_pos->SetTitle("Hits distribution on Emcal for bad reconstruction event");
  e_res_pri_xy_pos->SetStats(0);
  e_res_pri_xy_pos->SetContour(nb);
  e_res_pri_xy_pos->Draw("colorz");
  c2->cd(2);
  e_res_bad_xy_pos->SetStats(0);
  e_res_bad_xy_pos->SetContour(nb);
  e_res_bad_xy_pos->Draw("colorz");
  */
  auto *c3 = new TCanvas("c3", "c3", 1600, 800);
  c3->Divide(2, 1);
  c3->cd(1);
  //  hit_row_col_cry->Draw("colorz");
  hit_pos_crystal->SetStats(0);
  hit_pos_crystal->Draw("colorz");
  //  emcal_e_of_bad_res->Draw();
  c3->cd(2);
  //  hit_row_col_gla->Draw("colorz");
  hit_pos_glass->SetStats(0);
  hit_pos_glass->Draw("colorz");
  //  eloss_of_bad_res->Draw();
  
}
