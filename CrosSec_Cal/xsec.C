#include "TGenDVCS.h"
#include "TGVKelly.h"
#include "TSystem.h"
#include "TTree.h"
#include "TFile.h"
#include "TBranch.h"
#include "TH1F.h"
#include "TH2F.h"
#include <iostream>
#include <fstream>
using namespace std;

double xsec(void)
{
  //  ofstream myfile;
  //  myfile.open("weight.txt");
  
  //returns d4sigma/dQ2dxdtdphi in nb/GeV4
  gSystem->Load("libPhysics.so");
  gSystem->Load("./libTGVKelly.so");
  gSystem->Load("./libTGenGeo.so");
  gSystem->Load("./libTGenBase.so");
  gSystem->Load("./libTGenDVCS.so");

  //  TFile *hfile = new TFile("result.root", "update");
  TFile *hfile = new TFile("C1000_S200_B8_1M_2nd.root", "update");
  TTree *T = (TTree*)hfile->Get("T");
  Int_t N_events = (Int_t)T->GetEntries();
  Double_t Eb=2132.03, Q2, xb, t_var, phi, psf, phi_def, xsec, SLdt=10., NTOT;  //Eb is energy for fixed target.

  TBranch *add_br = T->Branch("xsec", &xsec, "xsec/D");
  T->SetBranchAddress("Q2", &Q2);
  T->SetBranchAddress("xb", &xb);
  T->SetBranchAddress("t_var", &t_var);
  T->SetBranchAddress("phi", &phi);
  T->SetBranchAddress("psf", &psf);
  //  T->SetBranchAddress("phi_def", &phi_def);

  TH1F *h_norm1 = new TH1F("h_norm1", "h_norm1", 90, 0., 45.);
  TH1F *h_norm2 = new TH1F("h_norm2", "h_norm2", 100, 0., 0.1);
  TH1F *h_norm3 = new TH1F("h_norm3", "h_norm3", 100, -2., 0.);
  TH1F *h_norm4 = new TH1F("h_norm4", "h_norm4", 63, 0., 6.3);
  TH2F *Q2_xsec = new TH2F("Q2_xsec", "Q2_xsec", 90, 0., 45., 200, 0., 100000.);
  TH2F *xb_xsec = new TH2F("xb_xsec", "xb_xsec", 100, 0., 0.1, 200, 0., 100000.);
  TH2F *t_var_xsec = new TH2F("t_var_xsec", "t_var_xsec", 100, -2., 0., 200, 0., 100000.);
  
  cout << "Number of events: " << N_events << endl << endl;
  cout << "Cross section  ||  Q2  ||  xb  ||  t  ||  phi" << endl;

  TGVKelly *tgv2 = new TGVKelly(Eb, kFALSE, kTRUE);
  TGenDVCS *gEv = new TGenDVCS(Eb, 0, 0, 0);
  Double_t ConvGeV2nbarn = 0.389379304e+6; // Unit conversion
  Double_t BHp, BHm, VCSp, VCSm, Ip, Im, SigmaTotPlus, SigmaTotMoins;

  
  for(int i = 0 ; i < N_events ; i++)
    {
      T->GetEntry(i);

      Eb = 2132.03;
      ConvGeV2nbarn = 0.389379304e+6; // Unit conversion
      
      Double_t* CFF = gEv->Interpol_CFF(Q2, xb, t_var);
      if(!CFF)
	{
	  cout << "Error !!" << endl;
	  return 0;
	}

      //      Double_t BHp, BHm, VCSp, VCSm, Ip, Im;
      //      Double_t SigmaTotPlus, SigmaTotMoins;
  
      BHp = tgv2->CrossSectionBH( Q2, xb, t_var, -phi, 1, 0, kTRUE );
      VCSp = tgv2->CrossSectionVCS( Q2, xb, t_var, -phi, 1, 0, CFF[0], CFF[1], CFF[2], CFF[3], CFF[4], CFF[5], CFF[6], CFF[7], kTRUE );
      Ip = tgv2->CrossSectionInterf( Q2, xb, t_var, -phi, 1, 0, -1, CFF[0], CFF[1], CFF[2], CFF[3], CFF[4], CFF[5], CFF[6], CFF[7], kTRUE );
      BHm = tgv2->CrossSectionBH( Q2, xb, t_var, -phi, -1, 0, kTRUE );
      VCSm = tgv2->CrossSectionVCS( Q2, xb, t_var, -phi, -1, 0, CFF[0], CFF[1], CFF[2], CFF[3], CFF[4], CFF[5], CFF[6], CFF[7], kTRUE );
      Im = tgv2->CrossSectionInterf( Q2, xb, t_var, -phi, -1, 0, -1, CFF[0], CFF[1], CFF[2], CFF[3], CFF[4], CFF[5], CFF[6], CFF[7], kTRUE );
      SigmaTotPlus = BHp + VCSp + Ip;
      SigmaTotMoins = BHm + VCSm + Im;
      //  if(opt==1) return TMath::Pi()*(BHp+BHm)* ConvGeV2nbarn;
      xsec = TMath::Pi() * ( SigmaTotPlus + SigmaTotMoins ) * ConvGeV2nbarn;

      NTOT = SLdt * xsec * psf * 1000000. / N_events;
      h_norm1->Fill(Q2, NTOT);
      h_norm2->Fill(xb, NTOT);
      h_norm3->Fill(t_var, NTOT);
      h_norm4->Fill(phi, NTOT);
      Q2_xsec->Fill(Q2, xsec);
      xb_xsec->Fill(xb, xsec);
      t_var_xsec->Fill(t_var, xsec);

      
      add_br->Fill();
      
      if(i % 10000 == 0)
	cout << i << "th value: " << xsec << ", " << Q2 << ", " << xb << ", " << t_var << ", " << phi << endl;
	//	cout << i << "th value: " << xsec << ", " << Q2 << ", " << xb << ", " << t_var << ", " << phi_def << endl;

      //    myfile << xsec << endl;
    }

  
  delete tgv2;  
  //  myfile.close();
  T->Write();
  hfile->Write();
  delete hfile;
  
  return xsec;  
}
