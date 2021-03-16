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
  //  TFile *file = TFile::Open("../../Data/g4e_simulation/g4e_output_foam_imposed_0_to_20k.root");
  TFile *file = TFile::Open("./g4e_output.root");
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
  
  auto RPOT_2_hit_xyz_proton = new TH3F("RPOT_2_hit_xyz_proton", "RPOT_2_hit_xyz_proton", 40, 550., 950., 40, -200., 200., 50, 26000., 26250.);
  auto RPOT_3_hit_xyz_proton = new TH3F("RPOT_3_hit_xyz_proton", "RPOT_3_hit_xyz_proton", 40, 750., 1150., 20, -100., 100., 50, 28080., 28150.);
  

  
  //==================================
  // Save the output
  //==================================
  double g_px = 0., g_py = 0., g_pz = 0., g_E = 0., e_px = 0., e_py = 0., e_pz = 0., e_E = 0.;
  double e_hit_emcal_x = 0., e_hit_emcal_y = 0., e_hit_emcal_z = 0.;
  double g_hit_emcal_x = 0., g_hit_emcal_y = 0., g_hit_emcal_z = 0.;
  int e_flag_emcal = 0, g_flag_emcal = 0, Etot_size = 0;

  /*
  string fileName_out = "../data/outCluster.root";
  TFile *fout = new TFile(fileName_out.c_str(), "Recreate");
  TTree *outTree = new TTree("outTree", "outTree");

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
  */

  
  //==================================
  // Loop the Event
  //==================================

  cout << "start" << endl;
  size_t events_numer = 0;  
  while (fReader.Next())
    {
      if(++events_numer > 2000)
	break;
      
      if(events_numer%1 == 0)
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
	  if(vol_name.rfind("ffi_RPOT", 0) == 0)
	    {
	      //	      if( hit_parent_track_id == 22 && hit_track_id == 2 )
	      cout << hit_parent_track_id << " " << hit_track_id << " " << hit_vol_name[i] << " " << "[" << x << ", " << y << ", " << z << "] " << endl;
	      //	      cout << hit_e_loss[i] * 1000. << "[MeV]" << endl;
	      //	      hit_energy+=hit_e_loss[i];
		      	      
	      track_ids_in_ecap_emcal.insert(hit_track_id);

	      //	      cout << x << ", " << y << ", " << z << endl;
	      
	      if( hit_track_id == 3 )
		{
		  if(z < 27000.)
		    RPOT_2_hit_xyz_proton->Fill(x, y, z);
		  else if(z > 27000.)
		    RPOT_3_hit_xyz_proton->Fill(x, y, z);
		}

	      /*
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
	      */
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
      
      
      g_flag_emcal = 0;  e_flag_emcal = 0;  hit_energy = 0.;
      
    }// end for the read g4e_output

  // fout->cd();
  // outTree->Write();
  // fout->Close();

  
  
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
  
  auto *c1 = new TCanvas("c1", "c1", 1600, 800);
  c1->Divide(2,1);
  c1->cd(1);
  RPOT_2_hit_xyz_proton->SetStats(0);
  RPOT_2_hit_xyz_proton->SetContour(nb);
  RPOT_2_hit_xyz_proton->Draw("box");
  c1->cd(2);
  RPOT_3_hit_xyz_proton->SetStats(0);
  RPOT_3_hit_xyz_proton->SetContour(nb);
  RPOT_3_hit_xyz_proton->Draw("box");
  //  RPOT_2_hit_xyz_proton->Draw("colorz");
  /*
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
  */
  
}
