#include "Riostream.h"
#include "TRandom.h"



void Save_root()
{
  Double_t x=0., y=0., z=0., w=0.;
  Int_t nlines = 0;
  auto f = TFile::Open("basic.root","RECREATE");
  TH1F *h1 = new TH1F("h1","x distribution",300,-1.,5.);
  TNtuple ntuple("ntuple","data from ascii file","x:y:z:w");
  TRandom *R = new TRandom();

  for(int i = 0 ; i < 10000 ; i++)
    {
      x = R->Gaus(3, 0.02);
      y = R->Gaus(2, 0.2);
      z = R->Gaus(1, 2);
      w = R->Poisson(0.3);
      ntuple.Fill(x,y,z,w);

      h1->Fill(x);
      h1->Fill(y);
      h1->Fill(z);
      nlines++;
    }

  
  printf(" found %d points\n",nlines);
  f->Write();

  h1->Draw();
}

