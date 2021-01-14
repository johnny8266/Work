#include "TFDISTR.h"
#include "TMath.h"
#include "TRandom.h"
#include <iostream>

using namespace std;

ClassImp(TFDISTR)

//_________________________________________________________

TFDISTR::TFDISTR()
{
  Eb = 2131.2132;
  ConvGeV2nbarn = 0.389379304e+6;
  //  tgv=new TGVKelly(Eb,kFALSE,kTRUE);
  gEv=new TGenDVCS(Eb,0,0,0);
}

//_________________________________________________________

Double_t TFDISTR::Density(int nDim, Double_t *Xarg)
{
  Eb = 2131.2132;
  tgv=new TGVKelly(Eb,kFALSE,kTRUE);
  
  // Integrand for mFOAM
  Double_t M = 0.938271998;

  Double_t Q2 = Xarg[0] * 90. + 2.;
  //  Double_t Q2 = 5.;

  Double_t xBMin = 2. * Eb * Q2 / (M * (4 * TMath::Power(Eb, 2)-Q2));  
  Double_t xBMax = Q2/(Q2-TMath::Power(M,2));
  //  Double_t xB = Xarg[1] * (0.1 - 0.001) + 0.001;
  //  if( (xB < xBMin) || (xB > xBMax) )
  //    return 0;
  Double_t xB = xBMin + 0.005;
    
  Double_t t = -Xarg[1];  
  //  Double_t t = -0.1;

  Double_t phi = Xarg[2] * 2. * TMath::Pi();
  //  Double_t phi = 1.4;



  Double_t ConvGeV2nbarn = 0.389379304e+6; // Unit conversion
  Double_t BHp, BHm, VCSp, VCSm, Ip, Im;
  Double_t SigmaTotPlus, SigmaTotMoins;
  Double_t* cffs = gEv->Interpol_CFF(Q2,xB,t);


  
  ////////////:
  BHp = tgv->CrossSectionBH( Q2, xB, t, -phi, 1, 0, kTRUE );
  VCSp = tgv->CrossSectionVCS( Q2, xB, t, -phi, 1, 0, cffs[0], cffs[1], cffs[2], cffs[3], cffs[4], cffs[5], cffs[6], cffs[7], kTRUE );
  Ip = tgv->CrossSectionInterf( Q2, xB, t, -phi, 1, 0, -1, cffs[0], cffs[1], cffs[2], cffs[3], cffs[4], cffs[5], cffs[6], cffs[7], kTRUE );
  BHm = tgv->CrossSectionBH( Q2, xB, t, -phi, -1, 0, kTRUE );
  VCSm = tgv->CrossSectionVCS( Q2, xB, t, -phi, -1, 0, cffs[0], cffs[1], cffs[2], cffs[3], cffs[4], cffs[5], cffs[6], cffs[7], kTRUE );
  Im = tgv->CrossSectionInterf( Q2, xB, t, -phi, -1, 0, -1, cffs[0], cffs[1], cffs[2], cffs[3], cffs[4], cffs[5], cffs[6], cffs[7], kTRUE );
  SigmaTotPlus = BHp + VCSp + Ip;
  SigmaTotMoins = BHm + VCSm + Im;
  DVCSxsec = TMath::Pi() * ( SigmaTotPlus + SigmaTotMoins ) * ConvGeV2nbarn;// Total DVCS cross section in nb/GeV4

  delete tgv;
  
  if(TMath::IsNaN(DVCSxsec))
    return 0.;
  return DVCSxsec;
}
