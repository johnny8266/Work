#include "TLorentzVector.h"
#include "TMath.h" 
#include "TRandom.h"
#include <vector>
#include <ctime>
#include <cmath>
#include <iostream>
#include <algorithm>
using namespace std;

void DVCS_beam_event_gen_V2()
{
  TFile *rfile = new TFile("DVCS_4Pars.root");
  TTree *T = (TTree*)rfile->Get("T");
  Int_t Iteration = (Int_t)T->GetEntries();
  Double_t Q2, xb, Eb, M, s_var, t_var, t0_min, t0_max, phi, phi_def, xsec, psf;
  vector<Double_t> Q_2, X_B, T_Var, P_hi, X_sec;
  cout << Iteration << " set of pars..." << endl;

  T->SetBranchAddress("Q2", &Q2);
  T->SetBranchAddress("xb", &xb);
  T->SetBranchAddress("t_var", &t_var);
  T->SetBranchAddress("phi", &phi);
  T->SetBranchAddress("xsec", &xsec);
  
  for(int j = 0 ; j < Iteration ; j++)
    {
      T->GetEntry(j);
      
      Q_2.push_back(Q2);
      X_B.push_back(xb);
      T_Var.push_back(t_var);
      P_hi.push_back(phi);
      X_sec.push_back(xsec);
    }
  Q2 = 0.;  xb = 0.;  t_var = 0.;  phi = 0.;  xsec=0.;

  Double_t D_Q2, D_xb, D_t;

  D_Q2 = (*max_element(Q_2.begin(), Q_2.end())) - (*min_element(Q_2.begin(), Q_2.end()));
  D_xb = (*max_element(X_B.begin(), X_B.end())) - (*min_element(X_B.begin(), X_B.end()));
  D_t = (*max_element(T_Var.begin(), T_Var.end())) - (*min_element(T_Var.begin(), T_Var.end()));

  cout << "Q2 max - min: " << D_Q2 << " || XB max - min: " << D_xb << " || t max - min: " << D_t << endl;
  
  delete rfile;
  
  
  // Declaration the variables and Create random function
  //
  
  TFile *hfile = new TFile("result.root", "RECREATE");
  TTree *DVCS = new TTree("DVCS", "Fill the DVCS data");
  //  TRandom *R = new TRandom();
  TLorentzVector CM_frame_HR_4, Virtual_photon, VP, e1, p1, photon;
  TLorentzVector* saveTL;
  TVector3 CM_frame_HR_3, CM_frame_fix_beam_3, z_axis(0, 0, 1.), rotate_axis, v1, v2;
  Double_t Virtual_photon_E=0., e1E=0., e1_S_angle_cos=0., e1_S_angle_sin=0.;
  Double_t S_angle_cos_CMS=0., S_angle_sin_CMS=0., E_CMS[4]={0.}, P_CMS[4]={0.}, M0_square[4]={0.};
  Double_t e1_S_angle, p1_S_angle, photon_S_angle;
  Double_t e1_px, e1_py, e1_pz, e1_E, Vg_px, Vg_py, Vg_pz, Vg_E;
  Double_t p1_px, p1_py, p1_pz, p1_E, g_px, g_py, g_pz, g_E, g_p_amp;
  Int_t count=0;
  //  long int t = (long int)time(NULL);  R->SetSeed(t);  //Get current time & set the random seed

  
  DVCS->Branch("Q2", &Q2, "Q2/D");
  DVCS->Branch("xb", &xb, "xb/D");
  DVCS->Branch("t_var", &t_var, "t_var/D");
  DVCS->Branch("phi", &phi, "phi/D");
  DVCS->Branch("phi_def", &phi_def, "phi_def/D");
  DVCS->Branch("psf", &psf, "psf/D");
  DVCS->Branch("xsec", &xsec, "xsec/D");
  DVCS->Branch("e1_S_angle", &e1_S_angle, "e1_S_angle/D");
  DVCS->Branch("p1_S_angle", &p1_S_angle, "p1_S_angle/D");
  DVCS->Branch("photon_S_angle", &photon_S_angle, "photon_S_angle/D");
  DVCS->Branch("e1_px", &e1_px, "e1_px/D");
  DVCS->Branch("e1_py", &e1_py, "e1_py/D");
  DVCS->Branch("e1_pz", &e1_pz, "e1_pz/D");
  DVCS->Branch("e1_E", &e1_E, "e1_E/D");
  DVCS->Branch("p1_px", &p1_px, "p1_px/D");
  DVCS->Branch("p1_py", &p1_py, "p1_py/D");
  DVCS->Branch("p1_pz", &p1_pz, "p1_pz/D");
  DVCS->Branch("p1_E", &p1_E, "p1_E/D");
  DVCS->Branch("Vg_px", &Vg_px, "Vg_px/D");
  DVCS->Branch("Vg_py", &Vg_py, "Vg_py/D");
  DVCS->Branch("Vg_pz", &Vg_pz, "Vg_pz/D");
  DVCS->Branch("Vg_E", &Vg_E, "Vg_E/D");
  DVCS->Branch("g_px", &g_px, "g_px/D");
  DVCS->Branch("g_py", &g_py, "g_py/D");
  DVCS->Branch("g_pz", &g_pz, "g_pz/D");
  DVCS->Branch("g_E", &g_E, "g_E/D");
  

  
  // Run the Event Generator
  //
  //  for(int i = 0 ; i < Iteration ; i++)
  for(int i = 0 ; i < 1 ; i++)
    {
      if(i % 10000 == 0) cout << i << " events are generated ......" << endl;
      
      
      // =================================
      // Boost the beam-beam collider to fix target, then calculation
      // =================================
      TLorentzVector e0(0, 0., -10., 10.), p0( (100.*TMath::Sin(0.025)), 0., (100.*TMath::Cos(0.025)), 100.004402);
      //      TLorentzVector e0(0, 0., -10., 10.), p0(0., 0., 100., 100.004402);
      //      p0.Print();
      //      TLorentzVector e0(0, 0., 2132.03, 2132.03), p0(0., 0., 0., 0.938);
      CM_frame_fix_beam_3 = p0.BoostVector();
      p0.Boost(-CM_frame_fix_beam_3);
      e0.Boost(-CM_frame_fix_beam_3);
      //p0.Print();    e0.Print();   cout << endl;
      //      Eb = e0.E();
      cout << Eb << endl << endl;

      rotate_axis = (e0.Vect()).Cross(z_axis);
      Double_t rotate_angle=(e0.Vect()).Angle(z_axis);
      e0.Rotate(rotate_angle, rotate_axis);  // Align the particle with the beamline
      e0.Print();      Eb = e0.E();
      //      cout << Eb << endl;
      

      // =================================
      // Initialze all parameters
      // =================================
      Q2=0.; xb=0.; M=0.938271998; s_var=0.; t_var=0.; t0_min=0.; t0_max=0.;

      Q2 = Q_2[i];  xb = X_B[i];  t_var = T_Var[i];  phi = P_hi[i];  xsec = X_sec[i];
      
      M0_square[0] = -Q2;  M0_square[1] = M * M;  M0_square[2] = 0;  M0_square[3] = M * M;  
      
      
      // =================================
      // Calculate the leptonic reaction
      // =================================
      Virtual_photon_E = Q2 / (2. * M * xb);  // Energy of Virtual photon
      if ( Virtual_photon_E > Eb )
	{
	  cout << i << "th event exceed the range: " << Q2 << " " << xb << " " << t_var << " " << phi << " " << endl;
	  continue;
	}
      
	
      //      cout << Eb << " " << Virtual_photon_E << endl;
      e1E = Eb - Virtual_photon_E; //cout << e1E << endl; // Energy of scattering electron
      e1_S_angle_cos = 1. - Q2 /(2. * Eb * e1E);  // Cos theta value of electron scattering
      e1_S_angle_sin = sqrt(1. - e1_S_angle_cos * e1_S_angle_cos);  
      e1_S_angle = TMath::ACos(e1_S_angle_cos);  // !!! this angle is under fix target scenario, in collider kinematic the angle will be different
      //      cout << "Electron Scattering angle: " << e1_S_angle << endl; 
      e1.SetE(e1E);
      e1.SetPz(-e1E * e1_S_angle_cos);
      e1.SetPy(0.);
      e1.SetPx(e1E * e1_S_angle_sin);      
      Virtual_photon = e0 - e1;
      e1.Boost(CM_frame_fix_beam_3);
      e1_px = e1.Px(); e1_py = e1.Py(); e1_pz = e1.Pz(); e1_E = e1.E();

      //      cout << e1_S_angle_cos << "  " << e1_S_angle_sin << "  " << e1_S_angle << endl;

      
      cout << "Leptonic reaction: " << endl << "====================================" << endl;
      cout << "Scattering Electron: ";  e1.Print();  cout << "Virtual photon in fix target frame: ";  Virtual_photon.Print();
      cout << "====================================" << endl;
      

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
      rotate_angle = (Virtual_photon.Vect()).Angle(z_axis);
      Virtual_photon.Rotate(rotate_angle, rotate_axis);  // Align the particle with the beamline
      p0.Rotate(rotate_angle, rotate_axis);
      //      cout << "After rotate in CMS frame :" << endl << "Vp:";  Virtual_photon.Print();  cout << "p0:";  p0.Print();

      // Calculate the Energy and momentum
      E_CMS[0] = (s_var + M0_square[0] - M0_square[1]) / (2. * sqrt(s_var));   E_CMS[1] = (s_var + M0_square[1] - M0_square[0]) / (2. * sqrt(s_var));
      E_CMS[2] = (s_var + M0_square[2] - M0_square[3]) / (2. * sqrt(s_var));   E_CMS[3] = (s_var + M0_square[3] - M0_square[2]) / (2. * sqrt(s_var));
      P_CMS[0] = sqrt(E_CMS[0] * E_CMS[0] - M0_square[0]);   P_CMS[1] = sqrt(E_CMS[1] * E_CMS[1] - M0_square[1]);
      P_CMS[2] = sqrt(E_CMS[2] * E_CMS[2] - M0_square[2]);   P_CMS[3] = sqrt(E_CMS[3] * E_CMS[3] - M0_square[3]); 

      /*
      for(int j = 0 ; j < 4 ; j++)
	if( isnan(E_CMS[j]) ) 
	  cout << i << "th: E[" << E_CMS[j] << "] || P[" << P_CMS[j] << "]" << endl;
      */
      
      t0_min = (E_CMS[0] - E_CMS[2]) * (E_CMS[0] - E_CMS[2]) - (P_CMS[0] - P_CMS[2]) * (P_CMS[0] - P_CMS[2]);
      t0_max = (E_CMS[0] - E_CMS[2]) * (E_CMS[0] - E_CMS[2]) - (P_CMS[0] + P_CMS[2]) * (P_CMS[0] + P_CMS[2]);
      //      t_var = R->Uniform(t0_max, t0_min);

      psf = D_Q2 * D_xb * D_t * 2. * TMath::Pi();

      S_angle_cos_CMS = 1. - ( (t0_min - t_var) / (2. * P_CMS[0] * P_CMS[2]) );
      S_angle_sin_CMS = sqrt(1. - S_angle_cos_CMS * S_angle_cos_CMS);
      
      if( isnan( S_angle_sin_CMS ) ) 
	{
	  cout << i << "th: " << (t0_min - t_var) << " " << P_CMS[0] << " " << P_CMS[2] << " " << S_angle_cos_CMS << " " << S_angle_sin_CMS << endl;
	  count++;
	  continue;
	}

      //      cout << "Cos: " << S_angle_cos_CMS << " || Sin: " << S_angle_sin_CMS << endl << endl;
      
      photon.SetE(E_CMS[2]);   photon.SetPz(E_CMS[2] * S_angle_cos_CMS);   photon.SetPy(0.);    photon.SetPx(E_CMS[2] * S_angle_sin_CMS);
      p1 = Virtual_photon + p0 - photon;

      //      if( isnan( p1.Px() ) )
      //      	cout << i << " " << p1.Px() << " " << photon.Px() << endl;
      
      photon.Rotate(-rotate_angle, rotate_axis);  // Rotate the particle back to origin direction
      p1.Rotate(-rotate_angle, rotate_axis);
      Virtual_photon.Rotate(-rotate_angle, rotate_axis);
            
      photon.Boost(CM_frame_HR_3); // Boost CM frame back to Lab frame
      p1.Boost(CM_frame_HR_3);
      Virtual_photon.Boost(CM_frame_HR_3);

      
      rotate_angle = (Virtual_photon.Vect()).Angle(z_axis);  // Rotate in phi direction
      //      p1.Print();  photon.Print(); cout << endl;
      p1.RotateY(rotate_angle);  photon.RotateY(rotate_angle);
      //      p1.Print();  photon.Print(); cout << endl;
      //      p1.Rotate(phi, Virtual_photon.Vect());  photon.Rotate(phi, Virtual_photon.Vect());
      p1.RotateZ(phi);  photon.RotateZ(phi);
      //      p1.Print();  photon.Print();
      p1.RotateY(-rotate_angle);  photon.RotateY(-rotate_angle);


      
      photon.Boost(CM_frame_fix_beam_3); // Boost fix-target frame back to beam-beam frame
      p1.Boost(CM_frame_fix_beam_3);
      Virtual_photon.Boost(CM_frame_fix_beam_3);
      //      p1_S_angle = (p1.Vect()).Angle(z_axis);  cout << "Proton angle in fix target: " << (TMath::Pi() - p1_S_angle) << endl;

      

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


      
      p1_S_angle = TMath::Pi() - (p1.Vect()).Angle(z_axis); // cout << p1_S_angle << "  ";
      //      p1_S_angle = TMath::ACos( p1_pz / TMath::Sqrt(p1_px * p1_px + p1_py * p1_py + p1_pz * p1_pz) ); cout << TMath::Pi() - p1_S_angle << endl << endl;
      //      cout << "Proton angle in beam-beam: " << (TMath::Pi() - p1_S_angle) << endl << endl;
      photon_S_angle = (photon.Vect()).Angle(z_axis);  

      DVCS->Fill();
      
      //      cout << "Proton Scattering angle: " << p1_S_angle << endl << endl;

      /*
      cout << "Hadronic reaction: " << endl << "====================================" << endl;
      cout << "Virtual Photon: "; Virtual_photon.Print(); cout << "Scattering Photon: ";  photon.Print();  cout << "Scattering proton: ";  p1.Print();
      cout << "====================================" << endl << endl << endl;
      */
  
    }  
  cout << "Finish the event generated !" << endl;
  cout << "Event number of t_var smaller than the t_min: " << count << endl;
  
  DVCS->Write();

  delete hfile;

}
