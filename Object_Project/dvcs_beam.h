#include "TLorentzVector.h"
#include "TVector3.h"
using namespace std;

class dvcs_beam
{
 public:
  dvcs_beam();
  
  void LeptonicProcess();
  void HadronicProcess();
  
  inline void SetElectron(TLorentzVector e){ e0 = e; }
  inline void SetProton(TLorentzVector p){ p0 = p; }
  inline void SetQ2(Double_t q2){ Q2 = q2; }
  inline void SetXb(Double_t xb){ Xb = xb; }
  inline void SetT(Double_t t){ t_var = t; }
  inline void Setphi(Double_t phi){ Phi = phi; }

  inline TLorentzVector GetProtonSC(){ return p1; }
  inline TLorentzVector GetElectronSC(){ return e1; }
  inline TLorentzVector GetPhotonSC(){ return photon; }
  inline TLorentzVector GetVPhotonSC(){ return Vphoton; }
  
 private:
  TLorentzVector e0, p0, e1, p1, photon, Vphoton; 
  TVector3 CM_frame_hadron, CM_frame_fix_beam, z_axis, rotate_axis;
  Double_t Q2, Xb, t_var, Phi, s_var, t0_min, t0_max;
  Double_t e1_S_angle_cos, e1_S_angle_sin, S_angle_cos_CMS, S_angle_sin_CMS;
  
};
