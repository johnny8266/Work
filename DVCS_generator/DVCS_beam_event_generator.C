#include "TLorentzVector.h"
#include "TMath.h" 
#include "TRandom.h"
#include <ctime>
#include <cmath>
#include <iostream>
using namespace std;

void DVCS_beam_event_generator()
{
  // Declaration the variables and Create random function
  //
  
  TFile *hfile = new TFile("result.root", "RECREATE");
  TTree *T = new TTree("T", "Fill the DVCS data");
  TRandom *R = new TRandom();
  TLorentzVector CM_frame_HR_4, Virtual_photon, VP, e1, p1, photon;
  TLorentzVector* saveTL;
  TVector3 CM_frame_HR_3, CM_frame_fix_beam_3, z_axis(0, 0, 1.), rotate_axis, v1, v2;
  Double_t Q2, Q2_max, xb, xb_min, Eb, M, s_var, t_var, t0_min, t0_max, psf;
  Double_t Virtual_photon_E=0., e1E=0., e1_S_angle_cos=0., e1_S_angle_sin=0.;
  Double_t S_angle_cos_CMS=0., S_angle_sin_CMS=0., E_CMS[4]={0.}, P_CMS[4]={0.}, M0_square[4]={0.};
  Double_t phi, phi_def, e1_S_angle, p1_S_angle, photon_S_angle;
  Double_t e1_px, e1_py, e1_pz, e1_E, Vg_px, Vg_py, Vg_pz, Vg_E;
  Double_t p1_px, p1_py, p1_pz, p1_E, g_px, g_py, g_pz, g_E;
  Int_t Iteration=1000000;
  long int t = (long int)time(NULL);  R->SetSeed(t);  //Get current time & set the random seed

  
  T->Branch("Q2", &Q2, "Q2/D");
  T->Branch("xb", &xb, "xb/D");
  T->Branch("t_var", &t_var, "t_var/D");
  T->Branch("phi", &phi, "phi/D");
  T->Branch("phi_def", &phi_def, "phi_def/D");
  T->Branch("psf", &psf, "psf/D");
  T->Branch("e1_S_angle", &e1_S_angle, "e1_S_angle/D");
  T->Branch("p1_S_angle", &p1_S_angle, "p1_S_angle/D");
  T->Branch("photon_S_angle", &photon_S_angle, "photon_S_angle/D");
  T->Branch("e1_px", &e1_px, "e1_px/D");
  T->Branch("e1_py", &e1_py, "e1_py/D");
  T->Branch("e1_pz", &e1_pz, "e1_pz/D");
  T->Branch("e1_E", &e1_E, "e1_E/D");
  T->Branch("p1_px", &p1_px, "p1_px/D");
  T->Branch("p1_py", &p1_py, "p1_py/D");
  T->Branch("p1_pz", &p1_pz, "p1_pz/D");
  T->Branch("p1_E", &p1_E, "p1_E/D");
  T->Branch("Vg_px", &Vg_px, "Vg_px/D");
  T->Branch("Vg_py", &Vg_py, "Vg_py/D");
  T->Branch("Vg_pz", &Vg_pz, "Vg_pz/D");
  T->Branch("Vg_E", &Vg_E, "Vg_E/D");
  T->Branch("g_px", &g_px, "g_px/D");
  T->Branch("g_py", &g_py, "g_py/D");
  T->Branch("g_pz", &g_pz, "g_pz/D");
  T->Branch("g_E", &g_E, "g_E/D");
  

  
  // Run the Event Generator
  //
  for(int i = 0 ; i < Iteration ; i++)
    {
      if(i % 10000 == 0) cout << i << " events are generated ......" << endl;
      
      
      // =================================
      // Boost the beam-beam collider to fix target, then calculation
      // =================================
      TLorentzVector e0(0, 0., 10., 10.), p0(0., 0., -100., 100.0044);      
      CM_frame_fix_beam_3 = p0.BoostVector();
      p0.Boost(-CM_frame_fix_beam_3);
      e0.Boost(-CM_frame_fix_beam_3);
      //      p0.Print();  e0.Print();  cout << endl;
      Eb = e0.E();
      //      cout << Eb << endl;
      
      

      // =================================
      // Initialze all parameters
      // =================================
      Q2=0.; xb=0.; xb_min=0.01; M=0.938; s_var=0.; t_var=0.; t0_min=0.; t0_max=0.;
      //      Q2 = R->Uniform(5., 90.);
      //      xb_min = Q2 / (Eb * 2. * M); 
      //      xb = R->Uniform(xb_min, 0.1);
      xb = R->Uniform(xb_min, 0.1);
      Q2_max = xb_min * Eb * 2. * M;
      Q2 = R->Uniform(1., Q2_max);
      
      phi = R->Uniform(0., 2.*TMath::Pi());
      M0_square[0] = -Q2;  M0_square[1] = 0.879844;  M0_square[2] = 0;  M0_square[3] = 0.879844;  
      
      
      // =================================
      // Calculate the leptonic reaction
      // =================================
      Virtual_photon_E = Q2 / (2. * M * xb);  // Energy of Virtual photon
      //      cout << Eb << " " << Virtual_photon_E << endl;
      e1E = Eb - Virtual_photon_E; //cout << e1E << endl; // Energy of scattering electron
      e1_S_angle_cos = 1. - Q2 /(2. * Eb * e1E);  // Cos theta value of electron scattering
      e1_S_angle_sin = sqrt(1. - e1_S_angle_cos * e1_S_angle_cos);  
      e1_S_angle = TMath::ACos(e1_S_angle_cos);  
      //      cout << "Electron Scattering angle: " << e1_S_angle << endl; 
      e1.SetE(e1E); e1.SetPz(e1E * e1_S_angle_cos); e1.SetPy(0.); e1.SetPx(e1E * e1_S_angle_sin);
      Virtual_photon = e0 - e1;
      e1.Boost(CM_frame_fix_beam_3);
      e1_px = e1.Px(); e1_py = e1.Py(); e1_pz = e1.Pz(); e1_E = e1.E(); 

      /*
      cout << "Leptonic reaction: " << endl << "====================================" << endl;
      cout << "Scattering Electron: ";  e1.Print();  cout << "Virtual photon in fix target frame: ";  Virtual_photon.Print();
      cout << "====================================" << endl;
      */

      // =================================
      // Calculate the hadronic reaction
      // =================================
      s_var = (Virtual_photon + p0) * (Virtual_photon + p0);
      //      cout << "S Var: " << s_var << endl;
      CM_frame_HR_4 = Virtual_photon + p0;
      CM_frame_HR_3 = CM_frame_HR_4.BoostVector();
      Virtual_photon.Boost(-CM_frame_HR_3);  // Boost the Tlorentzvector to the CMS frame
      p0.Boost(-CM_frame_HR_3);
      //      cout << "Before rotate in CMS frame :" << endl << "Vp:";  Virtual_photon.Print();  cout << "p0:";  p0.Print();
      
      rotate_axis = (Virtual_photon.Vect()).Cross(z_axis);
      Double_t rotate_angle=(Virtual_photon.Vect()).Angle(z_axis);
      Virtual_photon.Rotate(rotate_angle, rotate_axis);  // Align the particle with the beamline
      p0.Rotate(rotate_angle, rotate_axis);
      //      cout << "After rotate in CMS frame :" << endl << "Vp:";  Virtual_photon.Print();  cout << "p0:";  p0.Print();

      // Calculate the Energy and momentum
      E_CMS[0] = (s_var + M0_square[0] - M0_square[1]) / (2. * sqrt(s_var));   E_CMS[1] = (s_var + M0_square[1] - M0_square[0]) / (2. * sqrt(s_var));
      E_CMS[2] = (s_var + M0_square[2] - M0_square[3]) / (2. * sqrt(s_var));   E_CMS[3] = (s_var + M0_square[3] - M0_square[2]) / (2. * sqrt(s_var));
      P_CMS[0] = sqrt(E_CMS[0] * E_CMS[0] - M0_square[0]);   P_CMS[1] = sqrt(E_CMS[1] * E_CMS[1] - M0_square[1]);
      P_CMS[2] = sqrt(E_CMS[2] * E_CMS[2] - M0_square[2]);   P_CMS[3] = sqrt(E_CMS[3] * E_CMS[3] - M0_square[3]); 

      //      for(int j = 0 ; j < 4 ; j++)
      //	cout << E_CMS[j] << " || " << P_CMS[j] << endl;
      
      t0_min = (E_CMS[0] - E_CMS[2]) * (E_CMS[0] - E_CMS[2]) - (P_CMS[0] - P_CMS[2]) * (P_CMS[0] - P_CMS[2]);
      t0_max = (E_CMS[0] - E_CMS[2]) * (E_CMS[0] - E_CMS[2]) - (P_CMS[0] + P_CMS[2]) * (P_CMS[0] + P_CMS[2]);
      //      t_var = R->Uniform(t0_max, t0_min);

      if( t0_min < -3. )  cout << "t0_min smaller than -3." << endl;
      t_var = R->Uniform(-3., t0_min);
      if( isnan(t_var) )  cout << "Error" << endl;

      //      psf = (0.1 - xb_min) * (Q2_max - 5.) * (t0_min - t0_max) * 2. * TMath::Pi();
      psf = (0.1 - xb_min) * (Q2_max - 5.) * (t0_min + 3.) * 2. * TMath::Pi();
      //      cout << "Q2: " << Q2 << " || xb: " << xb << " || t: " << t_var << endl;      
      //      cout << "t0 min: " << t0_min << " || t0 max: " << t0_max << " || t: " << t_var << endl << endl;

      S_angle_cos_CMS = 1. - ( (t0_min - t_var) / (2. * P_CMS[0] * P_CMS[2]) );
      S_angle_sin_CMS = sqrt(1. - S_angle_cos_CMS * S_angle_cos_CMS);

      //      cout << "Cos: " << S_angle_cos_CMS << " || Sin: " << S_angle_sin_CMS << endl << endl;
      
      photon.SetE(E_CMS[2]);   photon.SetPz(E_CMS[2] * S_angle_cos_CMS);   photon.SetPy(0.);    photon.SetPx(E_CMS[2] * S_angle_sin_CMS);
      p1 = Virtual_photon + p0 - photon;

      photon.Rotate(-rotate_angle, rotate_axis);  // Rotate the particle back to origin direction
      p1.Rotate(-rotate_angle, rotate_axis);
      Virtual_photon.Rotate(-rotate_angle, rotate_axis);
            
      photon.Boost(CM_frame_HR_3); // Boost CM frame back to Lab frame
      p1.Boost(CM_frame_HR_3);
      Virtual_photon.Boost(CM_frame_HR_3);

      //      rotate_angle = (Virtual_photon.Vect()).Angle(z_axis);  // Rotate in phi direction
      //      p1.Print();  photon.Print(); cout << endl;
      //      p1.RotateY(rotate_angle);  photon.RotateY(rotate_angle);
      //      p1.Print();  photon.Print(); cout << endl;
      p1.Rotate(phi, Virtual_photon.Vect());  photon.Rotate(phi, Virtual_photon.Vect());
      //      p1.RotateY(phi);  photon.RotateY(phi);
      //      p1.Print();  photon.Print();
      //      p1.RotateY(-rotate_angle);  photon.RotateY(-rotate_angle);
      
      photon.Boost(CM_frame_fix_beam_3); // Boost fix-target frame back to beam-beam frame
      p1.Boost(CM_frame_fix_beam_3);
      Virtual_photon.Boost(CM_frame_fix_beam_3);


      

      // Define the phi from scattering process and compare with random one
      //
      
      v1 = (Virtual_photon.Vect()).Cross(e1.Vect());
      v2 = (Virtual_photon.Vect()).Cross(photon.Vect());
      phi_def = v1.Angle(v2);
      if( (Virtual_photon.Vect()).Dot(v1.Cross(v2)) < 0. )
	{
	  phi_def = 2. * TMath::Pi() - phi_def;
	  //	  cout << "!!!" << "      ";
	}
      //      cout << "phi: " << phi << " || phi_def: " << phi_def << " || Difference: " << (phi - phi_def) << endl;

      
      // Calculate the interested physics quantity
      //
      p1_px = p1.Px(); p1_py = p1.Py(); p1_pz = p1.Pz(); p1_E = p1.E();
      Vg_px = Virtual_photon.Px(); Vg_py = Virtual_photon.Py(); Vg_pz = Virtual_photon.Pz(); Vg_E = Virtual_photon.E();
      g_px = photon.Px(); g_py = photon.Py(); g_pz = photon.Pz(); g_E = photon.E();
      
      p1_S_angle = TMath::ATan( (p1.Px() / p1.Pz()) );
      photon_S_angle = TMath::ATan( (photon.Px() / photon.Pz()) );

      T->Fill();
      
      //      cout << "Proton Scattering angle: " << p1_S_angle << endl << endl;

      /*
      cout << "Hadronic reaction: " << endl << "====================================" << endl;
      cout << "Virtual Photon: "; Virtual_photon.Print(); cout << "Scattering Photon: ";  photon.Print();  cout << "Scattering proton: ";  p1.Print();
      cout << "====================================" << endl << endl << endl;
      */
  
    }  
  cout << "Finish the event generated !" << endl;

  T->Write();

  delete hfile;

}
