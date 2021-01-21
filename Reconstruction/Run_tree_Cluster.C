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
    double Et_dep;
    double E_digi;
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
    double CRYS_ZPOS = 2110; //ECAL zpos from the center of EIC detector
    double Ethr = 10.; //MeV
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


    //  Loop all hits per event
    //  Cluster Seed;
    for (int i = 0; i < Size; i++) {
        if (hit.at(i).E_digi > Ethr && hit.at(i).E_digi > ClusSeed_Ene) {
            ClusSeed_Ene = hit.at(i).E_digi;
            ClusSeed_xcrs = hit.at(i).x_crs;
            ClusSeed_ycrs = hit.at(i).y_crs;
            ClusSeed_zcrs = hit.at(i).z_crs;
            ClusSeed_npe = hit.at(i).npe;
            //      cout << "SEED "<<ClusSeed_xcrs<< " "<<ClusSeed_ycrs<<endl;

        }
    }


    //ENERGY TOT simul starting fro Et_dep
    for (int i = 0; i < Size; i++) {
        if (hit.at(i).Et_dep > Ethr) {
            double Dx = hit.at(i).x_crs - ClusSeed_xcrs;
            double Dy = hit.at(i).y_crs - ClusSeed_ycrs;
            if (sqrt(Dx * Dx + Dy * Dy) <= 3 * Rmoliere) {
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
    }   // for energy tot

    // Cluster Center

    double w_tot = 0;
    double x, y;
    x = 0;
    y = 0;

    for (int i = 0; i < Size; i++) {
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

    double sigmax2 = Clus_xx - std::pow(Clus_x, 2.);
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

    cluster.C_seed_energy = ClusSeed_Ene;
    cluster.C_energy = Clus_Etot;
    cluster.C_seed_x = ClusSeed_xcrs;
    cluster.C_seed_y = ClusSeed_ycrs;
    cluster.C_seed_z = ClusSeed_zcrs;
    cluster.C_x = Clus_x;
    cluster.C_y = Clus_y;
    cluster.C_radius = Clus_Radius;
    cluster.C_theta = Clus_Theta;
    cluster.C_phi = Clus_phi;
    cluster.C_size = Clus_size;
    cluster.C_Energy_tot_simul = Clus_Energy_tot_simul;
    cluster.C_size_simul = Clus_size_simul;
    cluster.C_seed_npe = ClusSeed_npe;


    return cluster;


}


int main() {
  
    vector<int> *Npe = 0;
    vector<double> *ADC = 0;
    vector<double> *TDC = 0;
    vector<double> *Etot_dep = 0;
    vector<double> *xcrs = 0;
    vector<double> *ycrs = 0;
    vector<double> *zcrs = 0;

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
    int ene = 0;


    int Dene = 5; //in MeV


    for (int run = 0; run < 18; run++) {
        //  for(int run=3; run<4; run++){
        if (run < 3) ene = 200 * run + 100;
        if (run == 3) ene = 1;
        if (run > 3) ene += 1;

        //   if(run==6 || run== 8 || run== 10 || run== 11) continue;
        //std::string fileName = "/Users/Mariangela/work/simul_eic/rootfile/ecal_30cm_12x12sipm/out."+std::to_string(ene)+".ene.root";
	//        std::string fileName = "/Users/Mariangela/work/simul_eic/test_newVersion.root";
	std::string fileName = "/home/pu-kai/mnt/g4e/build/g4e_output_10k_events_crossing_angle.root";
        cout << "*********** " << fileName << " **************" << endl;
//    std::string fileName_out = "/Users/Mariangela/work/simul_eic/rootfile/cluster_30cm_12x12sipm_1Othr/outCluster."+std::to_string(ene)+".root";
        std::string fileName_out = "/home/pu-kai/Work/Reconstruction/data/outCluster."+std::to_string(ene)+".root";
        TFile *f1 = new TFile(fileName.c_str(), "Read");

	// Read the g4e output
	//
        TTree *events = (TTree *) f1->Get("events");
        events->SetBranchAddress("ce_emcal_Etot_dep", &Etot_dep);
        events->SetBranchAddress("ce_emcal_Npe", &Npe);
        events->SetBranchAddress("ce_emcal_ADC", &ADC);
        events->SetBranchAddress("ce_emcal_TDC", &TDC);
        events->SetBranchAddress("ce_emcal_xcrs", &xcrs);
        events->SetBranchAddress("ce_emcal_ycrs", &ycrs);
        events->SetBranchAddress("ce_emcal_zcrs", &zcrs);


	// Create the root file to save the result
	//
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
        int nbin_npe;
        if (run >= 3) ene = ene * 1000;
        double xmin = ene - ene * 0.5;
        double xmax = ene + ene * 0.4;
        if (run <= 3) {
            Dene = 1;
            nbin_npe = 230400 / 10;
        } else {
            Dene = 10;
            nbin_npe = 230400 / 100;
        }
        int nbin = (xmax - xmin) / Dene;
	
        TH1D *hEnergy_cluster = new TH1D("hEnergy_cluster", "hEnergy_cluster; Energy ;counts", nbin, xmin, xmax);
        TH1D *hEnergy_cluster_simul = new TH1D("hEnergy_cluster_simul", "hEnergy_cluster_simul; Energy ;counts", nbin,
                                               xmin, xmax);
        TH1D *hNpe_seed = new TH1D("hNpe_seed", "hNpe_seed; Npe; counts", nbin_npe, 0, 230400);

        if (run >= 3) ene = ene / 1000;



	
        //  cout <<"Number of events: "<< events->GetEntries()<<endl;
        for (int i = 0; i < events->GetEntries(); i++) {

            //  cout<< "********** Event "<<i<<" ****************"<<endl;

            events->GetEntry(i);
            //    cout<<"N di crs colpiti "<<Etot_dep->size()<<endl;
            int nsize = Etot_dep->size();
            vector<Hit> hhit;  // The Hit structure is defined at the begining
            for (int j = 0; j < Etot_dep->size(); j++) {
                Hit hit;
                hit.x_crs = xcrs->at(j);
                hit.y_crs = ycrs->at(j);
                hit.z_crs = zcrs->at(j);
                hit.Et_dep = Etot_dep->at(j);
                hit.E_digi = ADC->at(j);
                hit.time = TDC->at(j);
                hit.npe = Npe->at(j);
                hhit.push_back(hit);
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
            outTree->Fill();
            /*
         cout <<"Cluster Seed Ene xy "<<Cl_seed_energy<< " "<<Cl_seed_x<< " "<<Cl_seed_y<<endl;
         cout <<"Cluster Energy xy: "<<Cl_energy<<" "<<Cl_x<<" "<<Cl_y<<endl;
         cout <<"Cluster size: "<<Cl_size<<endl;
         cout <<"Cluster angle theta phi "<<Cl_theta<<" "<<Cl_phi<<endl;
         */
            hEnergy_cluster->Fill(Cl_energy);
            hEnergy_cluster_simul->Fill(Cl_Energy_tot_simul);
            hNpe_seed->Fill(Cl_seed_npe);

        } //for on event
        fout->cd();
        outTree->Write();
        hEnergy_cluster->Write();
        hEnergy_cluster_simul->Write();
        hNpe_seed->Write();
        fout->Close();

    }//end for run
} // end main()

