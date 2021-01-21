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

double GausM(double *x, double *par) {

    return par[0] * exp(-0.5 * TMath::Power(((x[0] - par[1]) / par[2]), 2)) + 0.11;
}


Double_t Eresolution_fit(double *x, double *par) {

    return sqrt(TMath::Power(par[0] / sqrt(x[0]), 2) + TMath::Power(par[1], 2));
    //  return sqrt(TMath::Power(par[0]/sqrt(x[0]),2));
    //   return (par[0]*(1/sqrt(x[0]))+par[1]);
}


int main() {
    double *Cl_seed_energy = 0;
    int *Cl_seed_npe = 0;
    double *Cl_energy = 0;
    double *Cl_seed_x = 0;
    double *Cl_seed_y = 0;
    double *Cl_seed_z = 0;
    double *Cl_x = 0;
    double *Cl_y = 0;
    double *Cl_radius = 0;
    double *Cl_theta = 0;
    double *Cl_phi = 0;
    int *Cl_size = 0;
    int Cl_size_simul = 0;
    double Cl_Energy_tot_simul = 0;
    int ene = 0;
    double Vene[18];
    double Rene_sim[18];
    double Rene[18];
    double npe[18];
    double npe_sigma[18];

    //    std::string path = "/Users/Mariangela/work/simul_eic/rootfile/cluster_30cm_12x12sipm_1Othr/";
    //  std::string path = "/Users/Mariangela/work/simul_eic/rootfile/";
    std::string path = "/home/pu-kai/Work/Reconstruction/data/";
    std::string fileName_out = path + "outEnergy_resolution_test.root";

    cout << fileName_out << endl;

    //   "/Users/Mariangela/work/simul_eic/rootfile/cluster_20cm_6x6sipm_Othr/outEnergy_resolution.root";

    double emin = -0.05;
    double emax = 15.05;
    int nbine = (emax - emin) / 0.1;
    TH1D *hEnergy_Resolution = new TH1D("hEnergy_Resolution", "hEnergy_Resolution; Energy ;#sigma_{E}/E", nbine, emin,
                                        emax);


    for (int run = 0; run < 18; run++) {         //run<18
        if (run < 3) ene = 200 * run + 100;
        if (run == 3) ene = 1;
        if (run > 3) ene += 1;

        //   if(run==6 || run==8 || run==10 || run==11) continue;


        std::string fileName = path + "outCluster." + std::to_string(ene) + ".root";

        //"/Users/Mariangela/work/simul_eic/rootfile/cluster_20cm_6x6sipm_Othr/outCluster."+std::to_string(ene)+".root";
        cout << "************ " << fileName << " ************" << endl;
        cout << run << " " << ene << endl;
        /*
        double xmin = ene - ene*0.4;
        double xmax = ene + ene*0.4;
        int nbin =(xmax-xmin)/1;
        
        TH1D* hEnergy_cluster2 = new TH1D("hEnergy_cluster2", "hEnergy_cluster2; Energy ;counts", nbin, xmin, xmax);
    */
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

        cout << "Number of events: " << outTree->GetEntries() << endl;
        /*
             for (int i=0; i<outTree->GetEntries(); i++){
                  outTree->GetEntry(i);
                 cout <<Cl_Energy_tot_simul<<endl;
            hEnergy_cluster2->Fill(Cl_Energy_tot_simul);
             }
       */
        TH1D *hCluster = (TH1D *) f1->Get("hEnergy_cluster");
        TH1D *hCluster_sim = (TH1D *) f1->Get("hEnergy_cluster_simul");
        TH1D *hSeed_npe = (TH1D *) f1->Get("hNpe_seed");


        TCanvas *a0 = new TCanvas("a0", "a0", 1000, 1000);
        hSeed_npe->Draw();

        TF1 *fit_npe = new TF1("fit_npe", "gaus", 0, 230400);
        hSeed_npe->Fit("fit_npe", "R");
        fit_npe->SetLineColor(2);
        fit_npe->Draw("same");

        double par_nep[3];
        fit_npe->GetParameters(par_nep);
        npe[run] = par_nep[1];

        cout << "********** npe " << npe[run] << " *************" << endl;
        //  double npe_sigma[18];



        TCanvas *b0 = new TCanvas("b0", "b0", 1000, 1000);
        hCluster->Draw();
        if (run >= 3) ene = ene * 1000;
        double xmax = ene + 0.3 * ene;
        double xmin = ene - 0.4 * ene;
        if (run >= 3) ene = ene / 1000;

        TF1 *fit1 = new TF1("fit1", "gaus", xmin, xmax);
        hCluster->Fit("fit1", "", "", xmin, xmax);
        fit1->SetLineColor(4);
        fit1->Draw("same");

        double par[3];
        fit1->GetParameters(par);

        float Eresolution = (par[2] / par[1]) * 100;


        float ene2;
        if (run < 3) {
            ene2 = ene * 1E-3;
        } else {
            ene2 = ene;
        }
        hEnergy_Resolution->Fill(ene2, Eresolution);
        cout << "Energy " << ene2 << " GeV - resolution(%) " << Eresolution << endl;

        Vene[run] = ene2;
        Rene[run] = Eresolution;

        TCanvas *b1 = new TCanvas("b1", "b1", 1000, 1000);
        hCluster_sim->Draw();

        TF1 *fit1S = new TF1("fit1S", "gaus", xmin, xmax);
        hCluster_sim->Fit("fit1S", "", "", xmin, xmax);
        fit1S->SetLineColor(4);
        fit1S->Draw("same");

        double parS[3];
        fit1S->GetParameters(parS);

        float EresolutionS = (parS[2] / parS[1]) * 100;

        Rene_sim[run] = EresolutionS;

        cout << "Energy simul " << ene2 << " GeV - resolution(%) " << EresolutionS << endl;
    }//end for run

    TCanvas *c0 = new TCanvas("c0", "c0", 1000, 1000);
    hEnergy_Resolution->SetLineColor(1);
    hEnergy_Resolution->SetMarkerStyle(20);
    hEnergy_Resolution->Draw("hist p");


    TF1 *fit2 = new TF1("fit2", Eresolution_fit, 0., 15., 2);
    //  hEnergy_Resolution->Fit("fit2","","",0, 15);

    //    fit2->FixParameter(0,2.2);
    // fit2->FixParameter(1,0.4);
    // fit2->SetParLimits(0,2.7,3.9);
    //fit2->SetParLimits(1,0.1,0.7);
    fit2->SetLineColor(4);
    hEnergy_Resolution->Fit("fit2");
    fit2->Draw("same");


    TGraphErrors *gr_resolution = new TGraphErrors(18, Vene, Rene, 0, 0);
    gr_resolution->SetName("gr_resolution");
    gr_resolution->SetTitle("gr_resolution");

    TFile *fout = new TFile(fileName_out.c_str(), "Recreate");
    fout->cd();
    gr_resolution->Write();
//    fout->Close();

    for (int run = 0; run < 18; run++) {
        cout << "energy " << Vene[run] << " resolution " << Rene[run] << endl;
        cout << "Energy simul " << Vene[run] << " GeV - resolution(%) " << Rene_sim[run] << endl;
    }

    TGraphErrors *gr_resolution_intrinsic = new TGraphErrors(18, Vene, Rene_sim, 0, 0);
    gr_resolution_intrinsic->SetName("gr_resolution_intrinsic");
    gr_resolution_intrinsic->SetTitle("gr_resolution_intrinsic");

    fout->cd();
    gr_resolution_intrinsic->Write();

    TCanvas *c2 = new TCanvas("c2", "c2", 1000, 1000);
    gr_resolution_intrinsic->SetMarkerStyle(20);
    gr_resolution_intrinsic->GetYaxis()->SetRangeUser(0, 12.5);
    gr_resolution_intrinsic->SetMarkerSize(2);
    gr_resolution_intrinsic->SetLineColor(2);
    gr_resolution_intrinsic->SetMarkerColor(2);
    gr_resolution_intrinsic->Draw("AP");


    TCanvas *c1 = new TCanvas("c1", "c1", 1000, 1000);

    gr_resolution->SetMarkerStyle(20);
    gr_resolution->GetYaxis()->SetRangeUser(0, 12.5);
    gr_resolution->SetMarkerSize(2);
    gr_resolution->SetLineColor(4);
    gr_resolution->SetMarkerColor(4);
    gr_resolution->Draw("AP");

    TF1 *fit3 = new TF1("fit3", Eresolution_fit, 0.1, 16., 2);

    //   fit3->SetParLimits(0,2.8,3.2);
    //  fit3->SetParLimits(1,0.1,0.5);
    gr_resolution->Fit("fit3", "", "", 0.1, 16);

    fit3->SetLineColor(2);
    fit3->Draw("same");

    gr_resolution_intrinsic->Draw("sameP");

    TGraphErrors *gr_npe = new TGraphErrors(18, Vene, npe, 0, 0);
    gr_npe->SetName("gr_npe");
    gr_npe->SetTitle("gr_npe");

    TCanvas *aa1 = new TCanvas("aa1", "aa1", 1000, 1000);

    gr_npe->SetMarkerStyle(20);
    //  gr_npe->GetYaxis()->SetRangeUser(0,12.5);
    gr_npe->SetMarkerSize(2);
    gr_npe->SetLineColor(4);
    gr_npe->SetMarkerColor(4);
    gr_npe->Draw("AP");
    gr_npe->Write();
    fout->Close();
    return 0;
} // end main()

