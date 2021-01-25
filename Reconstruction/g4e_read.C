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
  double Et_dep;   // ce_emcal ETot deposit
  double E_digi;   // ce_emcal ADC
  double time;
  int npe;

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

Cluster ComputeCluster(vector<Hit> hit) {
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

    //    cout << "=============================================" << endl;
    //    cout << "Ce_Emcal hit size of this event: " << Size << endl;

    //  Loop all hits per event
    //  Cluster Seed;
    for (int i = 0; i < Size; i++) {
      if (hit.at(i).E_digi > Ethr && hit.at(i).E_digi > ClusSeed_Ene)  // iterate the energy to the larger one.
	{
	  ClusSeed_Ene = hit.at(i).E_digi;
	  ClusSeed_xcrs = hit.at(i).x_crs;
	  ClusSeed_ycrs = hit.at(i).y_crs;
	  ClusSeed_zcrs = hit.at(i).z_crs;
	  ClusSeed_npe = hit.at(i).npe;
	  
	  //	  cout << "SEED " << i << " : " << ClusSeed_xcrs << " "<< ClusSeed_ycrs << " " << ClusSeed_Ene << endl;
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
	    
            if (sqrt(Dx * Dx + Dy * Dy) <= 3 * Rmoliere)  // find the hits close to the powerful hit
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

            if (sqrt(Dx * Dx + Dy * Dy) <= 3 * Rmoliere) {

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

    for (int i = 0; i < Size; i++) { // Consider all ce_emcal hits! if hit energy so small, give it weight 0
      double w1 = std::max(0., (3.45 + std::log(hit.at(i).E_digi / Clus_Etot))); 
      x += w1 * hit.at(i).x_crs;
      y += w1 * hit.at(i).y_crs;
      Clus_xx += w1 * hit.at(i).x_crs * hit.at(i).x_crs;
      Clus_yy += w1 * hit.at(i).y_crs * hit.at(i).y_crs;
      w_tot += w1;
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
  
  TFile *file = TFile::Open("../Data/g4e_simulation/g4e_output_10k_events_crossing_angle.root");
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
  TTreeReaderArray<int>           ce_emcal_Npe = {fReader, "ce_emcal_Npe"};   
  TTreeReaderArray<double>        ce_emcal_Etot_dep = {fReader, "ce_emcal_Etot_dep"};
  TTreeReaderArray<double>        ce_emcal_ADC = {fReader, "ce_emcal_ADC"};
  TTreeReaderArray<double>        ce_emcal_TDC = {fReader, "ce_emcal_TDC"};
  TTreeReaderArray<double>        ce_emcal_xcrs = {fReader, "ce_emcal_xcrs"};
  TTreeReaderArray<double>        ce_emcal_ycrs = {fReader, "ce_emcal_ycrs"};
  TTreeReaderArray<double>        ce_emcal_zcrs = {fReader, "ce_emcal_zcrs"};

  double mass_electron = 0.00051099895;
  double mass_proton = 0.0938271998;

  
  //==================================
  // Histgram
  //==================================


  //==================================
  // Save the output
  //==================================

  double Cl_seed_energy = 0, Cl_Energy_tot_simul = 0, Cl_energy = 0;
  double Cl_seed_x = 0, Cl_seed_y = 0, Cl_seed_z = 0;
  double Cl_x = 0, Cl_y = 0, Cl_radius = 0, Cl_theta = 0, Cl_phi = 0;
  int Cl_size = 0, Cl_seed_npe = 0, Cl_size_simul = 0, ene = 0;
  double g_px = 0., g_py = 0., g_pz = 0., g_E = 0.;
  
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
  outTree->Branch("g_px", &g_px, "g_px/D");
  outTree->Branch("g_py", &g_py, "g_py/D");
  outTree->Branch("g_pz", &g_pz, "g_pz/D");
  outTree->Branch("g_E", &g_E, "g_E/D");
  
  //==================================
  // Loop the Event
  //==================================

  size_t events_numer = 0;  
  while (fReader.Next())
    {
      if(++events_numer > 7000)
	break;
      
      if(events_numer%100 == 0)
	cout << "Read " << events_numer << " th events..." << endl;
      
      std::unordered_set<uint64_t> track_ids_in_ecap_emcal;  // Get tracks information that have hits in ion EMCAL
      std::unordered_set<uint64_t> track_ids_in_ffi_RPOTS;  // Get tracks information that have hits in Roman Pots(FarForward ion direction area)
      
    
      // Read basic values
      auto hits_count = static_cast<size_t>(*hit_count.Get());         
      auto tracks_count = static_cast<size_t>(*trk_count.Get());       
      //      cout << endl << "This event has: " << tracks_count << " tracks || " << hits_count << " hits." << endl << endl;
      
	      
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
	      track_ids_in_ecap_emcal.insert(hit_track_id);
	    }
		  
	  if(vol_name.rfind("ffi", 0) == 0)
	    {
	      track_ids_in_ffi_RPOTS.insert(hit_track_id);
	    }
	    
	}


      
      // =============================
      // Calculate the hit cluster information
      // =============================
      auto Etot_size = 0;
      //      cout << Etot_size << endl;
      vector<Hit> hhit;
      
      for(int j = 0 ; j < 3 ; j++)
	{
	  if(gen_prt_charge[j] == 0.)
	    {
	      Etot_size = ce_emcal_Etot_dep.GetSize();
	      
	      for(int i = 0 ; i < Etot_size ; i++)
		{
		  Hit hit;
		  hit.x_crs = ce_emcal_xcrs[i];
		  hit.y_crs = ce_emcal_ycrs[i];
		  hit.z_crs = ce_emcal_zcrs[i];
		  hit.Et_dep = ce_emcal_Etot_dep[i];
		  hit.E_digi = ce_emcal_ADC[i];
		  hit.time = ce_emcal_TDC[i];
		  hit.npe = ce_emcal_Npe[i];
	  
		  hhit.push_back(hit);
		}
	    }
	}

      Cluster cluster;
      cluster = ComputeCluster(hhit);

      Cl_seed_energy = cluster.C_seed_energy;
      Cl_seed_npe = cluster.C_seed_npe;
      Cl_energy = cluster.C_energy;

      Cl_seed_x = cluster.C_seed_x;
      Cl_seed_y = cluster.C_seed_y;
      Cl_seed_z = cluster.C_seed_z;
      Cl_x = cluster.C_x;
      Cl_y = cluster.C_y;
      Cl_radius = cluster.C_radius;
      Cl_theta = cluster.C_theta;
      Cl_phi = cluster.C_phi;
      Cl_size = cluster.C_size;
      Cl_Energy_tot_simul = cluster.C_Energy_tot_simul;
      Cl_size_simul = cluster.C_size_simul;

      g_px = gen_prt_dir_x[1] * gen_prt_tot_mom[1];
      g_py = gen_prt_dir_y[1] * gen_prt_tot_mom[1];
      g_pz = gen_prt_dir_z[1] * gen_prt_tot_mom[1];
      g_E = gen_prt_tot_e[1];
      //      cout << g_px << " " << g_py << " " << g_pz << " " << g_E << endl; 

      outTree->Fill();

      
    }// end for the read g4e_output

  fout->cd();
  outTree->Write();
  fout->Close();
  
  // ======================================
  // Draw the Results
  // ======================================
}
