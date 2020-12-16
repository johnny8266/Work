#define Analysis_cxx
#include "Analysis.h"
#include <TH2.h>
#include <TH1.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include <iostream>
#include <fstream>

constexpr double mn=0.9396;
constexpr double mp=0.9383;
constexpr double mN=(mp+mn)/2.;
constexpr double e4he=-0.028295;
constexpr double m4he =2.*(mp+mn)+e4he;

void Analysis::Loop()
{

    TH1D * Q2c0 = new TH1D("Q2c0","Q2c0",50,2,10);
    TH1D * Q2c1 = new TH1D("Q2c1","Q2c1",50,2,10);
    TH1D * Q2c2 = new TH1D("Q2c2","Q2c2",50,2,10);
    TH1D * Q2c3 = new TH1D("Q2c3","Q2c3",50,2,10);

    TH1D * Xbc0 = new TH1D("Xbc0","Xbc0",50,0.,.5);
    TH1D * Xbc1 = new TH1D("Xbc1","Xbc1",50,0.,.5);
    TH1D * Xbc2 = new TH1D("Xbc2","Xbc2",50,0.,.5);
    TH1D * Xbc3 = new TH1D("Xbc3","Xbc3",50,0.,.5);

    TH1D * Phc0 = new TH1D("Phc0","Phc0",50,0,360);
    TH1D * Phc1 = new TH1D("Phc1","Phc1",50,0,360);
    TH1D * Phc2 = new TH1D("Phc2","Phc2",50,0,360);
    TH1D * Phc3 = new TH1D("Phc3","Phc3",50,0,360);

    TH1D * mtc0 = new TH1D("mtc0","mtc0",50,0,.5);
    TH1D * mtc1 = new TH1D("mtc1","mtc1",50,0,.5);
    TH1D * mtc2 = new TH1D("mtc2","mtc2",50,0,.5);
    TH1D * mtc3 = new TH1D("mtc3","mtc3",50,0,.5);

    if (fChain == 0) return;

    Long64_t nentries = fChain->GetEntriesFast();

    Long64_t nbytes = 0, nb = 0;
    for (Long64_t jentry=0; jentry<nentries;jentry++) {
        Long64_t ientry = LoadTree(jentry);
        if (ientry < 0) break;
        nb = fChain->GetEntry(jentry);   nbytes += nb;
        // if (Cut(ientry) < 0) continue;

        TLorentzVector *eB = new TLorentzVector(0,0,-ElBeam,ElBeam);
        TLorentzVector *hB = new TLorentzVector(0,0, std::sqrt(EhBeam*EhBeam-m4he*m4he),EhBeam);
        TLorentzVector *eS = new TLorentzVector(part_px[0],part_py[0],part_pz[0],part_e[0]);
        TLorentzVector *gS = new TLorentzVector(part_px[1],part_py[1],part_pz[1],part_e[1]);
        TLorentzVector *hS = new TLorentzVector(part_px[2],part_py[2],part_pz[2],part_e[2]);

        TLorentzVector *gV = new TLorentzVector(*eB-*eS);

        Q2c0->Fill(Q2); Xbc0->Fill(Xbj); Phc0->Fill(phih); mtc0->Fill(t);
        if (eS->Theta()>3.14159 -.0604) continue;
        Q2c1->Fill(Q2); Xbc1->Fill(Xbj); Phc1->Fill(phih); mtc1->Fill(t);
        if (gS->Theta()>3.14159 -.0222) continue;
        Q2c2->Fill(Q2); Xbc2->Fill(Xbj); Phc2->Fill(phih); mtc2->Fill(t);
        if (hS->Perp()<.4 || hS->Perp()>1.2) continue;
        Q2c3->Fill(Q2); Xbc3->Fill(Xbj); Phc3->Fill(phih); mtc3->Fill(t);
    }

    TCanvas *c0 = new TCanvas();
    Q2c2->Draw();
    TCanvas *c1 = new TCanvas();
    Xbc2->Draw();
    TCanvas *c2 = new TCanvas();
    Phc2->Draw();
    TCanvas *c3 = new TCanvas();
    mtc0->Draw(); mtc1->Draw("SAME"); mtc2->Draw("PFC SAME");
}

void Analysis::GenTxtFile(const char * name)
{
  ofstream myfile;
  myfile.open (name);
  myfile << "SIMPLE Event FILE\n";
  myfile << "============================================\n";
  myfile << "I, ievent, nParticles\n";
  myfile << "============================================\n";
  myfile << "I  K(I,1)  K(I,2)  K(I,3)  K(I,4)  K(I,5)  P(I,1)  P(I,2)  P(I,3)  P(I,4)  P(I,5)  V(I,1)  V(I,2)  V(I,3)\n";
  myfile << "============================================\n";

    if (fChain == 0) return;

    Long64_t nentries = fChain->GetEntriesFast();

    Long64_t nbytes = 0, nb = 0;
    for (Long64_t jentry=0; jentry<nentries;jentry++) {
        Long64_t ientry = LoadTree(jentry);
        if (ientry < 0) break;
        nb = fChain->GetEntry(jentry);   nbytes += nb;
        // if (Cut(ientry) < 0) continue;

        TLorentzVector *eB = new TLorentzVector(0,0,-ElBeam,ElBeam);
        TLorentzVector *hB = new TLorentzVector(0,0, std::sqrt(EhBeam*EhBeam-m4he*m4he),EhBeam);
        TLorentzVector *eS = new TLorentzVector(part_px[0],part_py[0],part_pz[0],part_e[0]);
        TLorentzVector *gS = new TLorentzVector(part_px[1],part_py[1],part_pz[1],part_e[1]);
        TLorentzVector *hS = new TLorentzVector(part_px[2],part_py[2],part_pz[2],part_e[2]);

        TLorentzVector *gV = new TLorentzVector(*eB-*eS);
  myfile << "0     " << jentry << "  3\n";
  myfile << "============================================\n";
  myfile << "1   " << "21   " << "11           " << "0   " << "3   " << "4   " << eB->Px() << "   " << eB->Py() << "   " << eB->Pz() << "   " << eB->E() << "   0.000511"   << "   0   " << "0   " << "0\n";
  myfile << "2   " << "21   " << "1000020040   " << "0   " << "6   " << "0   " << hB->Px() << "   " << hB->Py() << "   " << hB->Pz() << "   " << hB->E() << "   " << hB->M()<< "   0   " << "0   " << "0\n";
  myfile << "3   " << "21   " << "22           " << "1   " << "5   " << "0   " << gV->Px() << "   " << gV->Py() << "   " << gV->Pz() << "   " << gV->E() << "   " << gV->M()<< "   0   " << "0   " << "0\n";
  myfile << "4   " << "1    " << "11           " << "1   " << "0   " << "0   " << eS->Px() << "   " << eS->Py() << "   " << eS->Pz() << "   " << eS->E() << "   0.000511"   << "   0   " << "0   " << "0\n";
  myfile << "5   " << "1    " << "22           " << "3   " << "0   " << "0   " << gS->Px() << "   " << gS->Py() << "   " << gS->Pz() << "   " << gS->E() << "   0  "        << "   0   " << "0   " << "0\n";
  myfile << "6   " << "1    " << "1000020040   " << "2   " << "0   " << "0   " << hS->Px() << "   " << hS->Py() << "   " << hS->Pz() << "   " << hS->E() << "   " << hS->M()<< "   0   " << "0   " << "0\n";
  myfile << "=============== Event finished =============\n";

    if (jentry>=100) break;
    }


  myfile.close();
  return  ;
}
