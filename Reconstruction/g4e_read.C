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

struct Cluster {

    double C_seed_energy;
    int C_seed_npe;
    double C_seed_x;
    double C_seed_y;
    double C_seed_z;
    double C_energy;
    double C_x;
    double C_y;
    double C_radius;
    double C_theta;
    double C_phi;
    double C_size;
    double C_Energy_tot_simul;
    double C_size_simul;

};

Cluster ComputeCluster(vector<Hit> hit)
{
  int Size = hit.size();
  double CRYS_ZPOS = 2110; // ECAL zpos from the center of EIC detector
  double Ethr = 10.; // threshold[MeV]
  double Rmoliere = 20.01; // in mm    ->Rmolier for PbWO is 20 mm
  
  double ClusSeed_Ene = 0;
  int ClusSeed_npe = 0;
  double ClusSeed_xcrs = 0;
  double ClusSeed_ycrs = 0;
  double ClusSeed_zcrs = 0;
  double Clus_Etot = 0;
  double Clus_xx = 0;
  double Clus_yy = 0;
  double Clus_x = 0;
  double Clus_y = 0;
  int Clus_size = 0;
  double Clus_sigmaX = 0;
  double Clus_sigmaY = 0;
  double Clus_Radius = 0;
  double Clus_Theta = 0; //in deg;
  double Clus_phi = 0; //in deg;
  int Clus_size_simul = 0;
  double Clus_Energy_tot_simul = 0;

  map<int, double> map_crystal;  // ID, row, col, E
  
  //    cout << "=============================================" << endl;
  //    cout << "Ce_Emcal hit size of this event: " << Size << endl;


  //================================
  // Save all hit into map form
  //================================
  for(int i = 0 ; i < Size ; i++)
    if( hit.at(i).c_sec == 0 )
      {
	//	cout << i << endl;
	map_crystal.insert(pair<int, double>(hit.at(i).c_ID, hit.at(i).E_digi));
      }

  cout << "This event has " << Size << " hits and...." << endl;
  for(int i = 0 ; i < Size ; i++)
    {
      int count = 0;
      if( hit.at(i).c_sec == 0 )    // For crystal region
	{
	  int i_row = hit.at(i).c_row;
	  int i_col = hit.at(i).c_col;
	  double i_edigi = hit.at(i).E_digi;
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
		    {
		      if( i_edigi > (iter->second) )
			count++;
		      //		      else
		      //			continue;
		    }
		}
	    }
	  //	  if( count > 5 )
	    //	    cout << hit.at(i).c_col << " " << hit.at(i).c_row << " count: " << count << endl;
	}
    }
  cout << endl << endl;

  
  //================================
  // Loop all hits per event
  // Find Cluster Seed;
  //================================
  for(int i = 0 ; i < Size ; i++)
    {

	
      if (hit.at(i).E_digi > Ethr && hit.at(i).E_digi > ClusSeed_Ene)  // iterate the energy to the larger one.
	{
	  ClusSeed_Ene = hit.at(i).E_digi;
	  ClusSeed_xcrs = hit.at(i).x_crs;
	  ClusSeed_ycrs = hit.at(i).y_crs;
	  ClusSeed_zcrs = hit.at(i).z_crs;
	  ClusSeed_npe = hit.at(i).npe;
	}
    }

  // this seed has the largest energy and print its information
  //    cout << "Most powerful hit: " << ClusSeed_xcrs << " "<< ClusSeed_ycrs << " " << ClusSeed_Ene << endl << endl;
    

  //ENERGY TOT simul starting fro Et_dep
  for (int i = 0; i < Size; i++) {
    if (hit.at(i).Et_dep > Ethr) {
      double Dx = hit.at(i).x_crs - ClusSeed_xcrs;
      double Dy = hit.at(i).y_crs - ClusSeed_ycrs;
      //	    cout << "SEED " << i << " : " << Dx << " " << Dy << " " << ClusSeed_xcrs << " " << ClusSeed_ycrs << endl;
	    
      if (sqrt(Dx * Dx + Dy * Dy) <= 3. * Rmoliere)  // find the hits close to the powerful hit
	{
	  Clus_Energy_tot_simul += hit.at(i).Et_dep;
	  Clus_size_simul++;
	}
    }
  }


  //Cluster Energy tot
  for (int i = 0; i < Size; i++) {

    if (hit.at(i).E_digi > Ethr) {
      double Dx = hit.at(i).x_crs - ClusSeed_xcrs;
      double Dy = hit.at(i).y_crs - ClusSeed_ycrs;

      if (sqrt(Dx * Dx + Dy * Dy) <= 3. * Rmoliere) {

	//   cout <<hit.at(i).E_digi<< " "<<hit.at(i).x_crs<< " "<<hit.at(i).y_crs<<" "<< Dx << " "<< Dy<< " "<<sqrt(Dx*Dx+Dy*Dy)<<endl;
	Clus_Etot += hit.at(i).E_digi;
	Clus_size++;
      }
    } //end if ethr
  } //for energy tot

    //    cout << "Clus Etot: " << Clus_Etot << " || Clus size: " << Clus_size << endl << endl;


    // Cluster Center

  double w_tot = 0;
  double x, y;
  x = 0;
  y = 0;

  for (int i = 0; i < Size; i++)
    {
      // Consider hits near the seed! if hit energy is so small, give it weight 0
      if (hit.at(i).E_digi > Ethr)
	{
	  double Dx = hit.at(i).x_crs - ClusSeed_xcrs;
	  double Dy = hit.at(i).y_crs - ClusSeed_ycrs;

	  if (sqrt(Dx * Dx + Dy * Dy) <= 3. * Rmoliere)
	    {
	      double w1 = std::max(0., (3.45 + std::log(hit.at(i).E_digi / Clus_Etot))); 
	      x += w1 * hit.at(i).x_crs;
	      y += w1 * hit.at(i).y_crs;
	      Clus_xx += w1 * hit.at(i).x_crs * hit.at(i).x_crs;
	      Clus_yy += w1 * hit.at(i).y_crs * hit.at(i).y_crs;
	      w_tot += w1;
	    }
	}
    }
  Clus_x = x / w_tot;
  Clus_y = y / w_tot;
  Clus_xx /= w_tot;
  Clus_yy /= w_tot;

    
  // Cluster sigma

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

  //Cluster theta
  Clus_Theta = (std::atan((std::sqrt(std::pow(Clus_x, 2.) + std::pow(Clus_y, 2.))) / (CRYS_ZPOS + ClusSeed_zcrs))) *
    (180. / M_PI);

  //Cluster phi
  Clus_phi = std::atan2(Clus_x, Clus_y) * (180. / M_PI); //

  Cluster cluster;

  cluster.C_seed_energy = ClusSeed_Ene;                 //
  cluster.C_energy = Clus_Etot;                         //
  cluster.C_seed_x = ClusSeed_xcrs;
  cluster.C_seed_y = ClusSeed_ycrs;
  cluster.C_seed_z = ClusSeed_zcrs;
  cluster.C_x = Clus_x;
  cluster.C_y = Clus_y;
  cluster.C_radius = Clus_Radius;
  cluster.C_theta = Clus_Theta;
  cluster.C_phi = Clus_phi;
  cluster.C_size = Clus_size;
  cluster.C_Energy_tot_simul = Clus_Energy_tot_simul;   //
  cluster.C_size_simul = Clus_size_simul;               
  cluster.C_seed_npe = ClusSeed_npe;

  return cluster;
}




void g4e_read()
{
  
  //===================================
  // Open the g4e output file
  //===================================
  
  //  TFile *file = TFile::Open("../Data/g4e_simulation/g4e_output_10k_events_crossing_angle.root");
  TFile *file = TFile::Open("../Data/g4e_simulation/g4e_output_foam_imposed_1k_events.root");
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
  auto hit_pos_crystal = new TH2F("hit_pos_crystal", "hit_pos_crystal", 300, -1500, 1500, 300, -1500, 1500);
  auto hit_pos_glass = new TH2F("hit_pos_glass", "hit_pos_glass", 300, -1500, 1500, 300, -1500, 1500);

  
  //==================================
  // Save the output
  //==================================

  double Cl_seed_energy = 0, Cl_Energy_tot_simul = 0, Cl_energy = 0;
  double Cl_seed_x = 0, Cl_seed_y = 0, Cl_seed_z = 0;
  double Cl_x = 0, Cl_y = 0, Cl_radius = 0, Cl_theta = 0, Cl_phi = 0;
  int Cl_size = 0, Cl_seed_npe = 0, Cl_size_simul = 0, ene = 0;
  double g_px = 0., g_py = 0., g_pz = 0., g_E = 0., e_px = 0., e_py = 0., e_pz = 0., e_E = 0.;
  double e_hit_emcal_x = 0., e_hit_emcal_y = 0., e_hit_emcal_z = 0.;
  double g_hit_emcal_x = 0., g_hit_emcal_y = 0., g_hit_emcal_z = 0.;
  int e_flag_emcal = 0, g_flag_emcal = 0, Etot_size = 0;


  
  string fileName_out = "./data/outCluster.root";
  TFile *fout = new TFile(fileName_out.c_str(), "Recreate");
  TTree *outTree = new TTree("outTree", "outTree");

  outTree->Branch("Cl_seed_energy", &Cl_seed_energy, "Cl_seed_energy/D");
  outTree->Branch("Cl_seed_npe", &Cl_seed_npe, "Cl_seed_npe/I");
  outTree->Branch("Cl_seed_x", &Cl_seed_x, "Cl_seed_x/D");
  outTree->Branch("Cl_seed_y", &Cl_seed_y, "Cl_seed_y/D");
  outTree->Branch("Cl_energy", &Cl_energy, "Cl_energy/D");
  outTree->Branch("Cl_x", &Cl_x, "Cl_x/D");
  outTree->Branch("Cl_y", &Cl_y, "Cl_y/D");
  outTree->Branch("Cl_radius", &Cl_radius, "Cl_radius/D");
  outTree->Branch("Cl_theta", &Cl_theta, "Cl_theta/D");
  outTree->Branch("Cl_phi", &Cl_phi, "Cl_phi/D");
  outTree->Branch("Cl_size", &Cl_size, "Cl_size/I");
  outTree->Branch("Cl_phi", &Cl_phi, "Cl_phi/D");
  outTree->Branch("Cl_size", &Cl_size, "Cl_size/I");
  outTree->Branch("Cl_size_simul", &Cl_size_simul, "Cl_size_simul/I");
  outTree->Branch("Cl_Energy_tot_simul", &Cl_Energy_tot_simul, "Cl_Energy_tot_simul/D");
  outTree->Branch("e_hit_emcal_x", &e_hit_emcal_x, "e_hit_emcal_x/D");
  outTree->Branch("e_hit_emcal_y", &e_hit_emcal_y, "e_hit_emcal_y/D");
  outTree->Branch("e_hit_emcal_z", &e_hit_emcal_z, "e_hit_emcal_z/D");
  outTree->Branch("g_hit_emcal_x", &g_hit_emcal_x, "g_hit_emcal_x/D");
  outTree->Branch("g_hit_emcal_y", &g_hit_emcal_y, "g_hit_emcal_y/D");
  outTree->Branch("g_hit_emcal_z", &g_hit_emcal_z, "g_hit_emcal_z/D");
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
  outTree->Branch("Etot_size", &Etot_size, "Etot_size/I");


  
  //==================================
  // Loop the Event
  //==================================

  size_t events_numer = 0;  
  while (fReader.Next())
    {
      if(++events_numer > 5)
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
	  
	  //	  if( hit_track_id == 1 )
	  //	    cout << hit_parent_track_id << " " << hit_vol_name[i] << endl;

        
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
		      e_hit_emcal_x = x;
		      e_hit_emcal_y = y;
		      e_hit_emcal_z = z;
		      e_flag_emcal = 1;
		    }
		}
	      //		  cout << hit_track_id << " " << hit_parent_track_id << endl;
	      count++;
	    }
	}

      //      cout << "e flag: " << e_flag_emcal << " || g flag: " << g_flag_emcal << endl;
      //      cout << "hit counts in emcal: " << count << " || energy loss: " << hit_energy * 1000. << endl << endl;

      // iterate over the hit emcal tracks
      for(size_t i=0; i < hits_count; i++)
	{
	  if (trk_pdg[i] != 11) continue;       // Take only electrons for now
	  if (trk_parent_id[i] != 0) continue;  // Take only particles from a generator
                
	  // Check track has hits in ce_EMCAL
	  if (!track_ids_in_ecap_emcal.count(trk_id[i])) continue;
        
	  // Construct TLorenz vector
	  double px = trk_vtx_dir_x[i] * trk_mom[i];
	  double py = trk_vtx_dir_y[i] * trk_mom[i];
	  double pz = trk_vtx_dir_z[i] * trk_mom[i];

	  TLorentzVector lv;
	  lv.SetXYZM(px, py, pz, mass_electron);
	  //	  h1_el_e_tot->Fill(lv.Energy());
	}

      
      
      // =============================
      // Calculate the hit cluster information
      // =============================

      vector<Hit> hhit;
      double hit_e_check = 0.;
      
      Etot_size = ce_emcal_etot_dep.GetSize();
      //      cout << endl << Etot_size << ": " << endl;
      
      for(int i = 0 ; i < Etot_size ; i++)
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
	      hit_pos_crystal->Fill(ce_emcal_xcrs[i], ce_emcal_ycrs[i]);
	      hit_row_col_cry->Fill(ce_emcal_row[i], ce_emcal_col[i]);
				    //	      cout << ce_emcal_id[i] << " " << ce_emcal_row[i] << " " << ce_emcal_col[i] << endl;
	      cout << ce_emcal_id[i] << " [" << ce_emcal_xcrs[i] << " " << ce_emcal_ycrs[i] << " " << ce_emcal_zcrs[i] << "]" << endl;
	    }
	  else
	    {
	      hit_pos_glass->Fill(ce_emcal_xcrs[i], ce_emcal_ycrs[i]);
	      hit_row_col_gla->Fill(ce_emcal_row[i], ce_emcal_col[i]);
				    //	      cout << ce_emcal_id[i] << " " << ce_emcal_row[i] << " " << ce_emcal_col[i] << endl;
	      cout << ce_emcal_id[i] << " [" << ce_emcal_xcrs[i] << " " << ce_emcal_ycrs[i] << " " << ce_emcal_zcrs[i] << "]" << endl;
	    }
		    
	  //	  else
	    //	    hit_row_col_cry->Fill(ce_emcal_row[i], ce_emcal_col[i], 10);
	  //	  hit_e_check+=ce_emcal_etot_dep[i];
	  //	  cout << ce_emcal_xcrs[i] << " " << ce_emcal_ycrs[i] << " " << ce_emcal_zcrs[i] << "     " << ce_emcal_etot_dep[i] << endl;
	}
      //      cout << "double check: " << hit_e_check << endl;

      Cluster cluster;
      cluster = ComputeCluster(hhit);

      Cl_seed_energy = cluster.C_seed_energy;               // Max. energy of hits crystal(digi)
      Cl_seed_npe = cluster.C_seed_npe;                     // Max. npe of hits crystals
      Cl_energy = cluster.C_energy;                         // Total energy of the cluster(digi)
      Cl_Energy_tot_simul = cluster.C_Energy_tot_simul;     // Total energy of the cluster(deposit)

      Cl_seed_x = cluster.C_seed_x;                         // Position of hit with max. energy
      Cl_seed_y = cluster.C_seed_y;
      Cl_seed_z = cluster.C_seed_z;
      
      Cl_x = cluster.C_x;                                   // Cluster position after weighted
      Cl_y = cluster.C_y;
      
      Cl_radius = cluster.C_radius;
      Cl_theta = cluster.C_theta;
      Cl_phi = cluster.C_phi;

      Cl_size = cluster.C_size;                             // Cluster size(number of hits in this cluster)
      Cl_size_simul = cluster.C_size_simul;                 // Same as the above one


      g_px = gen_prt_dir_x[1] * gen_prt_tot_mom[1];
      g_py = gen_prt_dir_y[1] * gen_prt_tot_mom[1];
      g_pz = gen_prt_dir_z[1] * gen_prt_tot_mom[1];
      g_E = gen_prt_tot_e[1];

      e_px = gen_prt_dir_x[0] * gen_prt_tot_mom[0];
      e_py = gen_prt_dir_y[0] * gen_prt_tot_mom[0];
      e_pz = gen_prt_dir_z[0] * gen_prt_tot_mom[0];
      e_E = gen_prt_tot_e[0];
      //      cout << g_px << " " << g_py << " " << g_pz << " " << g_E << endl; 

      outTree->Fill();


      if( (e_flag_emcal == 1) && (g_flag_emcal == 0) )
	{
	  //	  cout << events_numer << endl;
	  double Dis_re_pri = (e_hit_emcal_x - Cl_x) * (e_hit_emcal_x - Cl_x) + (e_hit_emcal_y - Cl_y) * (e_hit_emcal_y - Cl_y);
	  //	  double Dis_re_pri = (g_hit_emcal_x - Cl_x) * (g_hit_emcal_x - Cl_x) + (g_hit_emcal_y - Cl_y) * (g_hit_emcal_y - Cl_y);
	  Dis_re_pri = TMath::Sqrt(Dis_re_pri);
	  if( ((Cl_Energy_tot_simul / 1000.) - e_E) > 0. )
	    {
	      if( (Dis_re_pri > 515.) && (Dis_re_pri < 520.) )
		{
		  cout << "Cluster energy: " << Cl_Energy_tot_simul << " || Primary e- energy: " << e_E << endl; 
		  cout << "Seed pos: " << Cl_seed_x << " " << Cl_seed_y << endl;
		  cout << "Cluster center pos: " << Cl_x << " " << Cl_y << endl;
		  cout << "Electron hit pos: " << e_hit_emcal_x << " " << e_hit_emcal_y << endl;
		  cout << "Photon primary direction: " << gen_prt_dir_x[1] << " " << gen_prt_dir_y[1] << " " << gen_prt_dir_z[1] << endl; 
		  //	      cout << "Photon hit pos: " << g_hit_emcal_x << " " << g_hit_emcal_y << endl;
		  cout << Dis_re_pri << " !!!" << endl;

		  
		  for(int k = 0 ; k < Etot_size ; k++)
		    {
		      e_res_pri_xy_pos->Fill(ce_emcal_xcrs[k], ce_emcal_ycrs[k], ce_emcal_etot_dep[k]);
		      emcal_e_of_bad_res->Fill(ce_emcal_etot_dep[k]);
		    }

		  for(size_t i = 0 ; i < hits_count ; i++)
		    {
		      uint64_t hit_track_id = static_cast<uint64_t>(hit_trk_id[i]);
		      uint64_t hit_parent_track_id = static_cast<uint64_t>(hit_parent_trk_id[i]);
		      std::string vol_name = static_cast<std::string>(hit_vol_name[i]);
		      
		      double x = hit_x[i], y = hit_y[i], z = hit_z[i];
		      double e_loss = hit_e_loss[i];
		      if(vol_name.rfind("ce_EMCAL", 0) == 0)
			{
			  e_res_bad_xy_pos->Fill(x, y, e_loss);
			  eloss_of_bad_res->Fill(e_loss * 1000.);
			  cout << hit_parent_track_id << "   " << hit_track_id << "   " << hit_vol_name[i] << "   " << "[" << x << ", " << y << ", " << z << "]   ||   " << e_loss << endl;
			  hit_energy = hit_energy + hit_e_loss[i];
			}
		    }
		  cout << "Emcal total hit energy: " << hit_energy << endl;
		}
	    }
	} // End if loop 


      
      g_flag_emcal = 0;  e_flag_emcal = 0;  hit_energy = 0.;
      
    }// end for the read g4e_output

  fout->cd();
  outTree->Write();
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
  hit_row_col_cry->Draw("colorz");
  //  hit_pos_crystal->Draw("colorz");
  //  hit_row_col_cry->Draw("colorz");
  //  emcal_e_of_bad_res->Draw();
  c3->cd(2);
  hit_row_col_gla->Draw("colorz");
  //  hit_pos_glass->Draw("colorz");
  //  eloss_of_bad_res->Draw();
  
}
