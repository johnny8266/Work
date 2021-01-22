/// \file
/// \ingroup tutorial_FOAM
/// \notebook -nodraw
/// Demonstrate the TFoam class.
///
///  To run this macro type from CINT command line
///
/// ~~~{.cpp}
///  root [0] gSystem->Load("libFoam.so")
///  root [1] .x foam_demo.C+
/// ~~~
///
/// \macro_code
///
/// \ author Stascek Jadach
/// \ modified by Wang, Pu-Kai 2020



#include "Riostream.h"
#include "TFile.h"
#include "TTree.h"
#include "TFoam.h"
#include "TH1F.h"
#include "TH1.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TFoamIntegrand.h"
#include "TRandom2.h"
#include "TRandom3.h"
#include "TFDISTR.h"
#include <time.h>
#include <vector>
#include <algorithm>
#include <fstream>
using namespace std;

Int_t main()
{
  //  gSystem->Load("libFoa.so");
  TFile RootFile("DVCS_4Pars.root","RECREATE");
  TTree *T = new TTree("T", "Fill simulated DVCS parameters");
  long   loop;
  Double_t MCresult, MCerror, MCwt;      
  Double_t Mass = 0.938271998, xb=0., xb_min=0., Q2_max=0., Q2=0., t_var=0., phi=0., psf=0., xsec_Integral=0., xsec_Integral_err=0.;

  T->Branch("Q2", &Q2, "Q2/D");
  T->Branch("xb", &xb, "xb/D");
  T->Branch("xb_min", &xb_min, "xb_min/D");
  T->Branch("t_var", &t_var, "t_var/D");
  T->Branch("phi", &phi, "phi/D");
  T->Branch("xsec_Integral", &xsec_Integral, "xsec_Integral/D");
  //  T->Branch("psf", &psf, "psf/D");

  /*
  TH1F *h1 = new TH1F("h1", "h1", 90, 0., 45.);
  TH1F *h2 = new TH1F("h2", "h2", 100, 0., 0.1);
  TH1F *h3 = new TH1F("h3", "h3", 100, -2., 0.);
  TH1F *h4 = new TH1F("h4", "h4", 63, 0., 6.3);
  */



  //-----------------------------------------
  long NevTot   =    100000;   // Total MC statistics
  Int_t  kDim   =         2;   // total dimension
  Int_t  nCells   =    5000;   // Number of Cells
  Int_t  nSampl   =     200;   // Number of MC events per cell in build-up
  Int_t  nBin     =       8;   // Number of bins in build-up
  Int_t  OptRej   =       1;   // Wted events for OptRej=0; wt=1 for OptRej=1 (default)
  Int_t  OptDrive =       2;   // (D=2) Option, type of Drive =0,1,2 for TrueVol,Sigma,WtMax
  Int_t  EvPerBin =      25;   // Maximum events (equiv.) per bin in buid-up
  Int_t  Chat     =       1;   // Chat level
  TRandom3 *PseRan   = new TRandom3();  // Create random number generator
  TFoam   *FoamX    = new TFoam("FoamX");   // Create Simulator
  TFoamIntegrand *rho = new TFDISTR();
  PseRan->SetSeed(0);
  //  long int tim = (long int)time(NULL);  cout << "time: " << tim << endl << endl;  PseRan->SetSeed(tim);
  //  long int tim = 1029384759;   PseRan->SetSeed(tim);
  Double_t *MCvect = new Double_t[kDim]; // vector generated in the MC run
      
  cout << "*****   Demonstration Program for Foam version " << FoamX->GetVersion() << "    *****" << endl;
  FoamX->SetkDim(        kDim);      // Mandatory!!!
  FoamX->SetnCells(      nCells);    // optional
  FoamX->SetnSampl(      nSampl);    // optional
  FoamX->SetnBin(        nBin);      // optional
  FoamX->SetOptRej(      OptRej);    // optional
  FoamX->SetOptDrive(    OptDrive);  // optional
  FoamX->SetEvPerBin(    EvPerBin);  // optional
  FoamX->SetChat(        Chat);      // optional
  FoamX->SetMaxWtRej(1.1);           // Maximum weight used to get w=1 MC events d=1.1	      
  //  FoamX->SetInhiDiv(1, 1);           // optional
  //  FoamX->SetInhiDiv(2, 1);           // optional
  //  FoamX->SetInhiDiv(3, 1);           // optional
  
  FoamX->SetRho(rho);
  FoamX->SetPseRan(PseRan);  
  FoamX->Initialize(); // Initialize simulator
  FoamX->Write("FoamX");     // Writing Foam on the disk, TESTING PERSISTENCY!!!
  
  long nCalls = FoamX->GetnCalls();
  cout << "====== Initialization done, entering MC loop" << endl;
  
  Double_t Eb = 2131.2132;
  Long_t n_effec=0;

  // Run the simulator
  //
  for(loop = 0 ; loop < NevTot ; loop++)
    {
      FoamX->MakeEvent();           // generate MC event
      FoamX->GetMCvect(MCvect);
      FoamX->GetIntegMC(xsec_Integral, xsec_Integral_err);
      //      FoamX->GetIntNorm(xsec_Integral, xsec_Integral_err);

      n_effec = FoamX->GetnEffev();

      Q2 = MCvect[0] * 98. + 2.;
      //      Q2 = TMath::Power(10., (0. + MCvect[0] * 2.));
      //      Q2 = 5.;
      
      xb_min = 2. * Eb * Q2 / (Mass * (4 * TMath::Power(Eb, 2)-Q2));
      xb = MCvect[1] * (0.1 - 0.0001) + 0.0001;	
      //      xb = TMath::Power(10., (-1. - MCvect[1] * 3.));
      //      xb = 0.05;

      //      t_var = -1. * MCvect[2];
      //      t_var = -1. * TMath::Power(10, (-4.*MCvect[2]));
      t_var = -0.5;
      
      //      phi = MCvect[3] * 2. * TMath::Pi();
      phi = 0.1;
      
      T->Fill();      
      /*
      h1->Fill(Q2);
      h2->Fill(xb);
      h3->Fill(t_var);
      h4->Fill(phi);
      */      
 
      if( ((loop) % 5000) == 0 )
	cout << n_effec << endl;
	//	cout << "loop = " << loop << ", " << Q2 << ", " << xb << ", " << t_var << ", " << phi << " || Simulation integral: " << xsec_Integral << " || xsec value: " << endl;
    }

  T->Write();

  /*  
  TCanvas* c1 = new TCanvas("c1", "c1", 800, 800);
  c1->Divide(2,2);
  c1->cd(1);
  h1->Draw();
  c1->cd(2);
  h2->Draw();
  c1->cd(3);
  h3->Draw();
  c1->cd(4); 
  h4->Draw();
  */
  
  Double_t eps = 0.0005;
  Double_t Effic, WtMax, AveWt, Sigma;
  Double_t IntNorm, Errel;
  FoamX->Finalize(   IntNorm, Errel);     // final printout
  FoamX->GetIntegMC( MCresult, MCerror );  // get MC intnegral
  FoamX->GetIntNorm( xsec_Integral, xsec_Integral_err );
  FoamX->GetWtParams(eps, AveWt, WtMax, Sigma); // get MC wt parameters
  Effic=0; if(WtMax>0) Effic=AveWt/WtMax;
  cout << "================================================================" << endl;
  cout << " MCresult= " << MCresult << " +- " << MCerror << " RelErr= "<< MCerror/MCresult << endl;
  cout << " Dispersion/<wt>= " << Sigma/AveWt << endl;
  cout << "      <wt>/WtMax= " << Effic <<",    for epsilon = "<< eps << endl;
  cout << " nCalls (initialization only) =   " << nCalls << endl;
  cout << "================================================================" << endl;

  // =======================================
  // creat a file to write in the data
  // =======================================
  fstream file;
  file.open("integral_vs_pars.txt", ios::app);
  file << nCells << " " << nSampl << " " << nBin << " " << MCresult << " " << MCerror << " " << xsec_Integral << " " << xsec_Integral_err << " " << Sigma/AveWt << " " << Effic << " " << eps << " " << nCalls << endl;
  file.close();
   
  delete [] MCvect;
  
  RootFile.Write();
  RootFile.Close();
  cout << "***** End of Demonstration Program  *****" << endl;
   
  return 0;
} 
