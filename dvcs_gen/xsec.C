#include "TGenDVCS.h"
#include "TGVKelly.h"
#include "TSystem.h"
#include "TGraph.h"
#include <iostream>

double xsec(void){
  //returns d4sigma/dQ2dxdtdphi in nb/GeV4
  gSystem->Load("libPhysics.so");
  gSystem->Load("./libTGVKelly.so");
  gSystem->Load("./libTGenGeo.so");
  gSystem->Load("./libTGenBase.so");
  gSystem->Load("./libTGenDVCS.so");
  
  Double_t Eb=5.75, Q2=1, xb=0.3, t=-0.444;

  TGenDVCS *gEv=new TGenDVCS(Eb,0,0,0);
  Double_t* CFF=gEv->Interpol_CFF(Q2,xb,t);
  if(!CFF) return 0;
  Double_t x[101],xsec[101];

  Double_t ConvGeV2nbarn = 0.389379304e+6; // Unit conversion
  Double_t BHp, BHm, VCSp, VCSm, Ip, Im;
  Double_t SigmaTotPlus, SigmaTotMoins;
  for(Int_t i=0;i<101;i++){// 100 points in phi
    TGVKelly *tgv2=new TGVKelly(Eb,kFALSE,kTRUE);
    Double_t phi=i*TMath::Pi()*2/100.;
    BHp = tgv2->CrossSectionBH( Q2, xb, t, -phi, 1, 0, kTRUE );
    VCSp = tgv2->CrossSectionVCS( Q2, xb, t, -phi, 1, 0, CFF[0], CFF[1], CFF[2], CFF[3], CFF[4], CFF[5], CFF[6], CFF[7], kTRUE );
    Ip = tgv2->CrossSectionInterf( Q2, xb, t, -phi, 1, 0, -1, CFF[0], CFF[1], CFF[2], CFF[3], CFF[4], CFF[5], CFF[6], CFF[7], kTRUE );
    BHm = tgv2->CrossSectionBH( Q2, xb, t, -phi, -1, 0, kTRUE );
    VCSm = tgv2->CrossSectionVCS( Q2, xb, t, -phi, -1, 0, CFF[0], CFF[1], CFF[2], CFF[3], CFF[4], CFF[5], CFF[6], CFF[7], kTRUE );
    Im = tgv2->CrossSectionInterf( Q2, xb, t, -phi, -1, 0, -1, CFF[0], CFF[1], CFF[2], CFF[3], CFF[4], CFF[5], CFF[6], CFF[7], kTRUE );
    SigmaTotPlus = BHp + VCSp + Ip;
    SigmaTotMoins = BHm + VCSm + Im;
    
    delete tgv2;
    //  if(opt==1) return TMath::Pi()*(BHp+BHm)* ConvGeV2nbarn;
    xsec[i]=TMath::Pi() * ( SigmaTotPlus + SigmaTotMoins ) * ConvGeV2nbarn;
    x[i]=phi*TMath::RadToDeg();
  }
  TGraph *g=new TGraph(101,x,xsec);
  g->Draw("A*");
}
