#include <iostream>
#include <time.h>
#include <fstream>
using namespace std;

Double_t Camel2(Int_t nDim, Double_t *Xarg);
Double_t sqr(Double_t x);

void foam_test()
{
  gSystem->Load("libFoam");

  //  TFile Rootfile("foam_test.root", "RECREATE");
  TH2D *hst_xy = new TH2D("hst_xy", "x-y plot", 50, 0, 1.0, 50, 0, 1.0);
  TH1F *h1 = new TH1F("h1", "h1", 100, 0., 1.);
  TH1F *h2 = new TH1F("h2", "h2", 100, 0., 1.);
  TH1F *h3 = new TH1F("h3", "h3", 100, 0., 1.);
  TH1F *h4 = new TH1F("h4", "h4", 100, 0., 1.);
  TH1F *h5 = new TH1F("h5", "h5", 100, 0., 1.);
  TH1F *h6 = new TH1F("h6", "h6", 100, 0., 1.);

  Int_t  nCells   =   10000;   // Number of Cells
  Int_t  nSampl   =     200;   // Number of MC events per cell in build-up
  Int_t  nBin     =       5;   // Number of bins in build-up
  
  Double_t *MCvect =new Double_t[2]; // 2-dim vector generated in the MC run
  TRandom3  *PseRan   = new TRandom3();  // Create random number generator
  long int tim = (long int)time(NULL);  cout << "time: " << tim << endl << endl;
  PseRan->SetSeed(tim);
  //  PseRan->SetSeed(4357);                // Set seed
  TFoam   *FoamX    = new TFoam("FoamX");   // Create Simulator
  FoamX->SetkDim(6);          // No. of dimensions, obligatory!
  FoamX->SetnCells(nCells);      // No. of cells, can be omitted, default=2000
  FoamX->SetnBin(nBin);           // Number of bins in build-up
  FoamX->SetnSampl(nSampl);         // Number of MC events per cell in build-up
  FoamX->SetRhoInt(Camel2);   // Set 2-dim distribution, included below
  FoamX->SetPseRan(PseRan);   // Set random number generator
  FoamX->Initialize();        // Initialize simulator, takes a few seconds...

  long nCalls = FoamX->GetnCalls();

  // From now on FoamX is ready to generate events according to Camel2(x,y)
  for(Long_t loop = 0 ; loop < 10000 ; loop++)
    {
      FoamX->MakeEvent();          // generate MC event

      FoamX->GetMCvect( MCvect);   // get generated vector (x,y)
      Double_t x=MCvect[0];
      Double_t y=MCvect[1];
      Double_t z=MCvect[0];
      Double_t w=MCvect[1];
      Double_t u=MCvect[0];
      Double_t v=MCvect[1];
      //      cout << "  (x,y) =  ( " << x << ", " << y << " )" << endl << endl;
      if( loop % 10000 == 0 )
	cout << loop << "th events done!" << endl;
      h1->Fill(x);
      h2->Fill(y);
      h3->Fill(z);
      h4->Fill(w);
      h5->Fill(u);
      h6->Fill(v);
      //      hst_xy->Fill(x,y);           // fill scattergram
    }

  Double_t eps = 0.0005;
  Double_t WtMax, AveWt, Sigma;
  Double_t IntNorm, Errel;
  Double_t mcResult, mcError, Effic;
  FoamX->Finalize(   IntNorm, Errel);     // final printout
  FoamX->GetIntegMC( mcResult, mcError);  // get MC integral, should be one
  FoamX->GetWtParams(eps, AveWt, WtMax, Sigma);  // get MC wt parameters
  Effic=0; if(WtMax>0) Effic=AveWt/WtMax;
  std::cout << " mcResult= " << mcResult << " +- " << mcError <<std::endl;
  // now hst_xy will be plotted visualizing generated distribution

  fstream file;
  file.open("integral_vs_pars.txt", ios::app);
  file << nCells << " " << nSampl << " " << nBin << " " << mcResult << " " << mcError << " " << Sigma/AveWt << " " << Effic << " " << eps << " " << nCalls << endl;
  file.close();
  
  TCanvas *cKanwa = new TCanvas("cKanwa","Canvas for plotting",1200,800);
  cKanwa->Divide(3,2);
  cKanwa->cd(1);
  h1->SetStats(0);
  h1->Draw();
  cKanwa->cd(2);
  h2->SetStats(0);
  h2->Draw();
  cKanwa->cd(3);
  h3->SetStats(0);
  h3->Draw();
  cKanwa->cd(4);
  h4->SetStats(0);
  h4->Draw();
  cKanwa->cd(5);
  h5->SetStats(0);
  h5->Draw();
  cKanwa->cd(6);
  h6->SetStats(0);
  h6->Draw();
  
  //  Rootfile.Write();
}//kanwa



//////////////////////////////////////////////////
//                                              //
// The input functions for the foam class       //
//                                              //
//////////////////////////////////////////////////

Double_t sqr(Double_t x){return x*x;};

Double_t Camel2(Int_t nDim, Double_t *Xarg)
{
  // 2-dimensional distribution for FOAM, normalized to one (within 1e-5)
  Double_t x = Xarg[0];
  Double_t y = Xarg[1];
  Double_t z = Xarg[2];
  Double_t w = Xarg[3];
  Double_t u = Xarg[4];
  Double_t v = Xarg[5];
  Double_t GamSq = sqr(0.100e0);
  Double_t Dist=exp( -( sqr(x-1./3)+sqr(y-1./6)+sqr(z-1./7)+sqr(w-1./2)+sqr(u-1./8)+sqr(v-1./4) ) / GamSq ) / GamSq / TMath::Pi()*2.;
  Dist        +=exp( -( sqr(x-2./3)+sqr(y-2./5)+sqr(z-2./3)+sqr(w-2./9)+sqr(u-2./8)+sqr(v-2./4) ) / GamSq) / GamSq / TMath::Pi();
  //  cout << "X: " << x << " Y: " << y << " Dist: " << Dist << " ";
  return 0.5*Dist;
}// Camel2




