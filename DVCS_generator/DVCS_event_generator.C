#include "TLorentzVector.h"
#include "TMath.h" 
#include "TRandom.h"
using namespace std;

void DVCS_event_generator()
{
  // Declaration the variables and Create random function
  //
  TRandom *R = new TRandom();
  Double_t Q2=2.3, xb=0.36, Eb=5.75, M=0.938, s_var=0., t_var=-0.2, t0=0., u_var=0., u0=0.;
  TLorentzVector e0(0, 0., 5.75, 5.75), p0(0., 0., 0., 0.938);
  TLorentzVector CM_frame_4, Virtual_photon, e1, p1, photon, Cheat_vector(-0.947635, 0, 3.60565, 3.40559);
  TLorentzVector* saveTL;
  TVector3 CM_frame_3, z_axis(0, 0, 1.), rotate_axis;
  Double_t Virtual_photon_E=0., e1E=0., e1_S_angle_cos=0., e1_S_angle_sin=0., photon_mass=0.;
  Double_t S_angle_cos_CMS=0., S_angle_sin_CMS=0., E_CMS[4]={0.}, P_CMS[4]={0.}, M0_square[4]={-2.3, 0.879844, 0., 0.879844};
  Double_t test=0., phi=TMath::Pi()/18.*17.;
  Int_t Iteration=1;


  for(int i = 0 ; i < Iteration ; i++)
    {
      // Random the parameters 
      //

      
      // =================================
      // Calculate the leptonic reaction
      // =================================
      Virtual_photon_E = Q2 / (2. * M * xb);  // Energy of Virtual photon
      e1E = Eb - Virtual_photon_E;  // Energy of scattering electron
      e1_S_angle_cos = 1. - Q2 /(2. * Eb * e1E);  // Cos theta value of electron scattering
      e1_S_angle_sin = sqrt(1. - e1_S_angle_cos * e1_S_angle_cos);
      e1.SetE(e1E); e1.SetPz(e1E * e1_S_angle_cos); e1.SetPy(0.); e1.SetPx(e1E * e1_S_angle_sin);
      Virtual_photon = e0 - e1;
      cout << "Leptonic reaction: " << endl;  cout <<"e1: ";  e1.Print();  cout << "VP: ";  Virtual_photon.Print();
      cout << endl << endl;


      // =================================
      // Calculate the hadronic reaction
      // =================================
      s_var = (Virtual_photon + p0) * (Virtual_photon + p0);
      cout << "S Var: " << s_var << endl;
      CM_frame_4 = Virtual_photon + p0;
      CM_frame_3 = CM_frame_4.BoostVector();
      Virtual_photon.Boost(-CM_frame_3);  // Boost the Tlorentzvector to the CMS frame
      p0.Boost(-CM_frame_3);
      
      rotate_axis = (Virtual_photon.Vect()).Cross(z_axis);
      Double_t rotate_angle=(Virtual_photon.Vect()).Angle(z_axis);
      //      cout << "After rotate in CMS frame :" << endl;  Virtual_photon.Print();  p0.Print();
      Virtual_photon.Rotate(rotate_angle, rotate_axis);  // Align the particle with the beamline
      p0.Rotate(rotate_angle, rotate_axis);
      //      cout << "After rotate in CMS frame :" << endl;  Virtual_photon.Print();  p0.Print();

      // Calculate the Energy and momentum
      E_CMS[0] = (s_var + M0_square[0] - M0_square[1]) / (2. * sqrt(s_var));   E_CMS[1] = (s_var + M0_square[1] - M0_square[0]) / (2. * sqrt(s_var));
      E_CMS[2] = (s_var + M0_square[2] - M0_square[3]) / (2. * sqrt(s_var));   E_CMS[3] = (s_var + M0_square[3] - M0_square[2]) / (2. * sqrt(s_var));
      P_CMS[0] = sqrt(E_CMS[0] * E_CMS[0] - M0_square[0]);   P_CMS[1] = sqrt(E_CMS[1] * E_CMS[1] - M0_square[1]);
      P_CMS[2] = sqrt(E_CMS[2] * E_CMS[2] - M0_square[2]);   P_CMS[3] = sqrt(E_CMS[3] * E_CMS[3] - M0_square[3]); 
    
      t0 = (E_CMS[0] - E_CMS[2]) * (E_CMS[0] - E_CMS[2]) - (P_CMS[0] - P_CMS[2]) * (P_CMS[0] - P_CMS[2]);  
  
      S_angle_cos_CMS = 1. - ( (t0 - t_var) / (2. * P_CMS[0] * P_CMS[2]) );
      S_angle_sin_CMS = sqrt(1. - S_angle_cos_CMS * S_angle_cos_CMS);

      photon.SetE(E_CMS[2]);   photon.SetPz(E_CMS[2] * S_angle_cos_CMS);   photon.SetPy(0.);    photon.SetPx(E_CMS[2] * S_angle_sin_CMS);
      p1 = Virtual_photon + p0 - photon;

      photon.Rotate(-rotate_angle, rotate_axis);  // Rotate the particle back to origin direction
      p1.Rotate(-rotate_angle, rotate_axis);
      Virtual_photon.Rotate(-rotate_angle, rotate_axis);
            
      photon.Boost(CM_frame_3); // Boost CM frame back to Lab frame
      p1.Boost(CM_frame_3);
      Virtual_photon.Boost(CM_frame_3);
      //      Virtual_photon.Print();
      
      rotate_angle = (Virtual_photon.Vect()).Angle(z_axis);
      p1.RotateY(rotate_angle);  photon.RotateY(rotate_angle);
      //      cout << endl << "Before rotate phi" << endl;
      //      photon.Print();  p1.Print();
      p1.RotateZ(phi);  photon.RotateZ(phi);
      p1.RotateY(-rotate_angle);  photon.RotateY(-rotate_angle);
      cout << "After rotate phi" << endl;  cout << "photon: ";  photon.Print();  cout << "p1: ";  p1.Print();
      cout << endl;
      
      //      p1.RotateY(-rotate_angle);  photon.RotateY(-rotate_angle);
      //      photon.Print();  p1.Print();
    }
    

  /*
  // Creating the Tree, Draw and Save Data
  //
  auto f = TFile::Open("DVCS_evt_gen_result.root","RECREATE");
  TH1F *h1 = new TH1F("h1","x distribution",300,-1.,5.);
  TTree Save4Vec("Save4Vec", "Save the TLorentzVector");
  Save4Vec.SetBranchAddress("saveTL", &saveTL);
  Save4Vec.Write();
  //TNtuple ntuple("ntuple","data from ascii file","x:y:z:w");
  */
  
}
