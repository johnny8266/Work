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
  /*
  inline void Set_Q2(Double_t q2){Q2_foam = q2;}
  inline void Set_Xb(Double_t xb){Xb_foam = xb;}
  inline void Set_t(Double_t tvar){t_foam = tvar;}
  inline void Set_phi(Double_t phivar){phi_foam = phivar;}
  */
   
public:
  TFDISTR();
  Double_t Density(int nDim, Double_t *Xarg);
  Double_t Q2_foam, Xb_foam, t_foam, phi_foam;
  /*
  inline Double_t Get_Q2() const {return Q2_foam;}
  inline Double_t Get_Xb() const {return Xb_foam;}
  inline Double_t Get_t() const {return t_foam;}
  inline Double_t Get_phi() const {return phi_foam;}
  */
  ClassDef(TFDISTR,1) //Class of testing functions for FOAM

};

