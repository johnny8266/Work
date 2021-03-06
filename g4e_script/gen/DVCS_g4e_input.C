#include <iostream>
#include <fstream>
#include <ctime>
#include "TMath.h" 
#include "TRandom.h"

using namespace std;

void DVCS_g4e_input()
{
 
  TFile *hfile = new TFile("DVCS_foam_100k.root");
  TTree *DVCS = (TTree*)hfile->Get("DVCS");
  Double_t phi, phi_def, e1_S_angle, p1_S_angle, photon_S_angle, Q2, xb, t_var, xsec, psf;
  Double_t e1_px, e1_py, e1_pz, e1_E, Vg_px, Vg_py, Vg_pz, Vg_E;
  Double_t p1_px, p1_py, p1_pz, p1_E, g_px, g_py, g_pz, g_E;
  Double_t e_sc_theta;
  Int_t N_events = (Int_t)DVCS->GetEntries(), count=0;
  
  cout << "N_events: " << N_events << endl;
 

  DVCS->SetBranchAddress("Q2", &Q2);
  DVCS->SetBranchAddress("xb", &xb);
  DVCS->SetBranchAddress("t_var", &t_var);
  DVCS->SetBranchAddress("phi", &phi);
  DVCS->SetBranchAddress("psf", &psf);
  DVCS->SetBranchAddress("e1_S_angle", &e1_S_angle);
  DVCS->SetBranchAddress("p1_S_angle", &p1_S_angle);
  DVCS->SetBranchAddress("photon_S_angle", &photon_S_angle);
  DVCS->SetBranchAddress("e1_px", &e1_px);
  DVCS->SetBranchAddress("e1_py", &e1_py);
  DVCS->SetBranchAddress("e1_pz", &e1_pz);
  DVCS->SetBranchAddress("e1_E", &e1_E);
  DVCS->SetBranchAddress("p1_px", &p1_px);
  DVCS->SetBranchAddress("p1_py", &p1_py);
  DVCS->SetBranchAddress("p1_pz", &p1_pz);
  DVCS->SetBranchAddress("p1_E", &p1_E);
  DVCS->SetBranchAddress("Vg_px", &Vg_px);
  DVCS->SetBranchAddress("Vg_py", &Vg_py);
  DVCS->SetBranchAddress("Vg_pz", &Vg_pz);
  DVCS->SetBranchAddress("Vg_E", &Vg_E);
  DVCS->SetBranchAddress("g_px", &g_px);
  DVCS->SetBranchAddress("g_py", &g_py);
  DVCS->SetBranchAddress("g_pz", &g_pz);
  DVCS->SetBranchAddress("g_E", &g_E);
 

  ofstream myfile;
  myfile.open ("../../Data/g4e_input/dvcs_input.txt");

  cout << "Total events: " << N_events << endl;
  cout << "Half events: " << N_events << endl;
  
  for(int i = 0 ; i < N_events ; i++)
    {

      DVCS->GetEntry(i);
      
      myfile << "3 " << Q2 << " " << xb << " " << t_var << " " << phi_def << " " << psf << " " << e1_S_angle << " " << p1_S_angle << " " << photon_S_angle << " " << xsec << endl;
      myfile << "1 " << "-1 " << "1 " << "11 " << "1 " << "1 " << e1_px << " " << e1_py << " " << e1_pz << " " << e1_E << " " << "0.000511 " << "0.0 0.0 0.0" << endl;
      myfile << "2 " << "0 " << "1 " << "22 " << "1 " << "1 " << g_px << " " << g_py << " " << g_pz << " " << g_E << " " << "0.0 " << "0.0 0.0 0.0" << endl;
      myfile << "3 " << "1 " << "1 " << "2212 " << "1 " << "1 " << p1_px << " " << p1_py << " " << p1_pz << " " << p1_E << " " << "0.938272 " << "0.0 0.0 0.0" << endl;  
		  
      if( i % 500 == 0 ) cout << i << " th events finished......" << endl;
      
      if(i > 9999)
	break;
    }
  myfile.close();
 
  
  /*
  Int_t N_events = 0, loop = 0;

  auto gE_distribution = new TH1F("gE_distribution", "gE_distribution", 22., 0., 11.);

  auto pos_map = new TH2F("pos_map", "pos_map", 140, -1400., 1400., 140, -1400., 1400.);
  
  TRandom *R = new TRandom();
  R->SetSeed(0);
  ofstream myfile;
  myfile.open ("../Data/g4e_input/dvcs_input.txt");

  double gx = 0., gy = 0., gz = 0., gE = 0.;
  vector<double> g_x, g_y, g_z, g_E;

  
  while(1)
    {      
      int empty_flag = 0;
      g_x.clear();  g_y.clear();  g_z.clear();  g_E.clear();  
      
      for(int i = 0 ; i < 3 ; i++)
	{
	  gz = -1. * R->Uniform(0.1, 10.);
	  gx = -1. * gz * R->Uniform(-0.589, 0.589);
	  gy = -1. * gz * R->Uniform(-0.589, 0.589);
	  // gx = -1. * gz * R->Uniform(-0.335, 0.335);
	  // gy = -1. * gz * R->Uniform(-0.335, 0.335);
	  gE = TMath::Sqrt(gx * gx + gy * gy + gz * gz);

	  double ring_R = TMath::Sqrt( (gx / gz * 2240.) * (gx / gz * 2240.) + (gy / gz * 2240.) * (gy / gz * 2240.) );
	  
	  //	  if( ((((gx < (-0.067 * gz)) && (gx > 0.)) || ((gx > (0.067 * gz)) && (gx < 0.))) && (((gy < (-0.067 * gz)) && (gy > 0.)) || ((gy > (0.067 * gz)) && (gy < 0.)))) || (ring_R > 1320.) )
	  //	  if( (((gx < (-0.067 * gz)) && (gx > 0.)) || ((gx > (0.067 * gz)) && (gx < 0.))) && (((gy < (-0.067 * gz)) && (gy > 0.)) || ((gy > (0.067 * gz)) && (gy < 0.))) )
	  if( ((((gx < (-0.375 * gz)) && (gx > 0.)) || ((gx > (0.375 * gz)) && (gx < 0.))) && (((gy < (-0.375 * gz)) && (gy > 0.)) || ((gy > (0.375 * gz)) && (gy < 0.)))) || (ring_R > 1320.) )
	    {
	      //	      cout << ring_R << endl;
	      empty_flag = 1;
	      break;
	    }

	  g_x.push_back(gx);  g_y.push_back(gy);  g_z.push_back(gz);  g_E.push_back(gE);
	}
      
      if( empty_flag == 0 )
	{
	  myfile << "3 " << 0.1 << " " << 0.01 << " " << -0.1 << " " << 3.14 << " " << 11.1 << " " << 0.1 << " " << 0.1 << " " << 0.1 << " " << 0.1 << endl;
	  for(int j = 0 ; j < 3 ; j++)
	    {
	      gx = g_x[j];  gy = g_y[j];  gz = g_z[j];  gE = g_E[j];
	      myfile << j << " 0 " << "1 " << "22 " << "1 " << "1 " << gx << " " << gy << " " << gz << " " << gE << " " << "0.0 " << "0.0 0.0 0.0" << endl;
	      gE_distribution->Fill(gE);
	      pos_map->Fill((gx / gz * 2240.), (gy / gz * 2240.));
	    }
	  N_events++;
	}
	  		  
      if( N_events % 500 == 0 )
	cout << N_events << "th events finished......" << endl;
	//	cout << i << " th: [" << gx << " " << gy << " " << gz << " " << gE << "] events finished......" << endl;

      if( N_events > 20000)
	break;

      loop++;
    }

  cout << "while loop: " << loop << endl;
  
  auto c1 = new TCanvas("c1", "c1", 1000, 500);
  c1->Divide(2,1);
  c1->cd(1);
  gE_distribution->Draw();
  c1->cd(2);
  pos_map->Draw("colorz");
  */  



}
