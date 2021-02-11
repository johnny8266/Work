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

  //  TFile *hfile = new TFile("./root_file/DVCS_4Pars.root", "update");
  TFile *hfile = new TFile("./root_file/DVCS_4Pars.root", "update");
  TTree *T = (TTree*)hfile->Get("T");
  Int_t N_events = (Int_t)T->GetEntries();
  Double_t Eb=2131.2132, Q2, xb, t_var, phi, psf, phi_def, xsec, xsec_inte;  //Eb is energy for fixed target.

  TBranch *add_br = T->Branch("xsec", &xsec, "xsec/D");
  //  T->Branch("xsec_inte", &xsec_inte, "xsec_inte/D");
  
  T->SetBranchAddress("Q2", &Q2);
  T->SetBranchAddress("xb", &xb);
  T->SetBranchAddress("t_var", &t_var);
  T->SetBranchAddress("phi", &phi);
  //  T->SetBranchAddress("psf", &psf);
  //  T->SetBranchAddress("phi_def", &phi_def);

  
  cout << "Number of events: " << N_events << endl << endl;
  cout << "              Cross section  ||  Q2  ||  xb  ||  t  ||  phi" << endl;


  TGenDVCS *gEv = new TGenDVCS(Eb, 0, 0, 0);
  Double_t ConvGeV2nbarn = 0.389379304e+6; // Unit conversion
  Double_t BHp, BHm, VCSp, VCSm, Ip, Im, SigmaTotPlus, SigmaTotMoins;

  
  for(int i = 0 ; i < N_events ; i++)
    {
      TGVKelly *tgv2 = new TGVKelly(Eb, kFALSE, kTRUE);
      ConvGeV2nbarn = 0.389379304e+6; // Unit conversion
      
      T->GetEntry(i);

      Double_t xb_0 = xb;
      xb = 0.1;
      
      Double_t* CFF = gEv->Interpol_CFF(Q2, xb, t_var);
      if(!CFF)
	{
	  cout << "Error !!" << endl;
	  return 0;
	}
  
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

      xsec = xsec * TMath::Sqrt(0.1 / xb_0);  // modify the Uniform cross section with the sqrt(1./xb)
      //      xsec = xsec * TMath::Sqrt(0.1 / xb);  // modify the Foam cross section with the sqrt(0.1/xb)
	    
      add_br->Fill();
      
      if(i % 2000 == 0)
	cout << i << "th value: " << xsec << ", " << Q2 << ", " << xb << ", " << t_var << ", " << phi << endl;

      delete tgv2;    
    }

  
  T->Write();
  hfile->Write();
  delete hfile;
  
  return xsec;  
}
