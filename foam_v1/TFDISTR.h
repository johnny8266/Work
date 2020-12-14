#include "TFoamIntegrand.h"
#include "TGVKelly.h"
#include "TGenDVCS.h"
#include <iostream>
using namespace std;

class TFDISTR: public TFoamIntegrand
{
  
private:
  Double_t Eb;//Electron beam energy
  double cffs[8];
  Double_t ConvGeV2nbarn; // Changement d'unit
  Double_t BHp, BHm, VCSp, VCSm, Ip, Im;
  Double_t SigmaTotPlus, SigmaTotMoins;
  Double_t DVCSxsec;
  //  Double_t Q2_foam, Xb_foam, t_foam, phi_foam;
  int* pFuncAddress;//!
  int* pArgsAddress;//!
  TGVKelly *tgv;
  TGenDVCS *gEv;

   
public:
  TFDISTR();
  Double_t Density(int nDim, Double_t *Xarg);
  Double_t Q2_foam, Xb_foam, t_foam, phi_foam, xsec_value;

  //  void Set_xsec(Double_t xs);
  Double_t Get_xsec();

  ClassDef(TFDISTR,1) //Class of testing functions for FOAM

};

