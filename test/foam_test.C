#include <iostream>
using namespace std;

Double_t Camel2(Int_t nDim, Double_t *Xarg);
Double_t sqr(Double_t x);

void foam_test()
{
  gSystem->Load("libFoam");

  TFile Rootfile("foam_test.root", "RECREATE");
  TH2D *hst_xy = new TH2D("hst_xy", "x-y plot", 50, 0, 1.0, 50, 0, 1.0);
  TH1F *h1 = new TH1F("h1", "h1", 200, 0., 1.);
  TH1F *h2 = new TH1F("h2", "h2", 200, 0., 1.);

  Double_t *MCvect =new Double_t[2]; // 2-dim vector generated in the MC run
  TRandom3  *PseRan   = new TRandom3();  // Create random number generator
  PseRan->SetSeed(4357);                // Set seed
  TFoam   *FoamX    = new TFoam("FoamX");   // Create Simulator
  FoamX->SetkDim(2);          // No. of dimensions, obligatory!
  FoamX->SetnCells(500);      // No. of cells, can be omitted, default=2000
  FoamX->SetnBin(2);           // Number of bins in build-up
  FoamX->SetnSampl(50);         // Number of MC events per cell in build-up
  FoamX->SetRhoInt(Camel2);   // Set 2-dim distribution, included below
  FoamX->SetPseRan(PseRan);   // Set random number generator
  FoamX->Initialize();        // Initialize simulator, takes a few seconds...
  
  // From now on FoamX is ready to generate events according to Camel2(x,y)
  for(Long_t loop = 0 ; loop < 100000 ; loop++)
    {
      FoamX->MakeEvent();          // generate MC event

      FoamX->GetMCvect( MCvect);   // get generated vector (x,y)
      Double_t x=MCvect[0];
      Double_t y=MCvect[1];
      //      cout << "  (x,y) =  ( " << x << ", " << y << " )" << endl << endl;
      
      h1->Fill(x);
      h2->Fill(y);
      hst_xy->Fill(x,y);           // fill scattergram
    }// loop
  
  Double_t mcResult, mcError;
  FoamX->GetIntegMC( mcResult, mcError);  // get MC integral, should be one
  std::cout << " mcResult= " << mcResult << " +- " << mcError <<std::endl;
  // now hst_xy will be plotted visualizing generated distribution

  TCanvas *cKanwa = new TCanvas("cKanwa","Canvas for plotting",600,600);
  cKanwa->cd();
  hst_xy->SetStats(0);
  hst_xy->Draw("lego2");
  h1->SetStats(0);
  h1->Draw();
  h2->SetStats(0);
  h2->Draw();
  
  Rootfile.Write();
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
  Double_t GamSq = sqr(0.100e0);
  Double_t Dist=exp(-(sqr(x-1./3) +sqr(y-1./3))/GamSq)/GamSq/TMath::Pi()*2.;
  Dist        +=exp(-(sqr(x-2./3) +sqr(y-2./3))/GamSq)/GamSq/TMath::Pi();
  //  cout << "X: " << x << " Y: " << y << " Dist: " << Dist << " ";
  return 0.5*Dist;
}// Camel2
