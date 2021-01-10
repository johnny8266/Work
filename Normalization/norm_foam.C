#include <iostream>
#include <fstream>
using namespace std;

void norm_foam()
{
  TFile *hfile = new TFile("result_foam.root");
  TTree *DVCS = (TTree*)hfile->Get("DVCS");
  TCanvas *c1 = new TCanvas("c1", "c1", 800, 800);
  Double_t phi, phi_def, e1_S_angle, p1_S_angle, photon_S_angle, Q2, xb, t_var, xsec, psf;
  Double_t e1_px, e1_py, e1_pz, e1_E, Vg_px, Vg_py, Vg_pz, Vg_E;
  Double_t p1_px, p1_py, p1_pz, p1_E, g_px, g_py, g_pz, g_E;
  Double_t e_sc_theta;
  TVector3 z_axis(0., 0., 1.), electron;
  Int_t N_events = (Int_t)DVCS->GetEntries(), count=0;
  
  cout << "N_events: " << N_events << endl;
  //  cout << "count: " << count << endl;


  DVCS->SetBranchAddress("Q2", &Q2);
  DVCS->SetBranchAddress("xb", &xb);
  DVCS->SetBranchAddress("t_var", &t_var);
  DVCS->SetBranchAddress("phi", &phi);
  DVCS->SetBranchAddress("phi_def", &phi_def);
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
  DVCS->SetBranchAddress("xsec", &xsec);

  Double_t SLdt=10., Q2_weight=0., NTOT=0.;

  TH1F* Q2_norm = new TH1F("Q2_norm", "Q2_norm", 90, 0., 45.);
  ofstream myfile;
  myfile.open ("/home/pu-kai/mnt/g4e/DVCS/dvcs_input.txt");
  //  myfile.open ("/vol0/pwang-l/g4e/DVCS/dvcs_input.txt");
    
  for(int i = 0 ; i < N_events ; i++)
    {
      DVCS->GetEntry(i);
      Q2_weight = SLdt / N_events * xsec * psf;
      if( i % 1 == 0 )
      	{
	  //	  cout <<  Q2 << " " << xb << " " << t_var << " " << phi_def << " " << psf << endl;
	  //	  cout << xsec << "  " << psf << "  " << Q2_weight << endl;
	  myfile << "3 " << Q2 << " " << xb << " " << t_var << " " << phi_def << " " << psf << " " << e1_S_angle << " " << p1_S_angle << " " << photon_S_angle << " " << xsec << endl;
	  myfile << "1 " << "-1 " << "1 " << "11 " << "1 " << "1 " << e1_px << " " << e1_py << " " << e1_pz << " " << e1_E << " " << "0.000511 " << "0.0 0.0 0.0" << endl;
	  myfile << "2 " << "0 " << "1 " << "22 " << "1 " << "1 " << g_px << " " << g_py << " " << g_pz << " " << g_E << " " << "0.0 " << "0.0 0.0 0.0" << endl;
	  myfile << "3 " << "1 " << "1 " << "2212 " << "1 " << "1 " << p1_px << " " << p1_py << " " << p1_pz << " " << p1_E << " " << "0.938272 " << "0.0 0.0 0.0" << endl;  
		  
	}
      if( i % 500 == 0 ) cout << i << " th events finished......" << endl;
      NTOT = NTOT + Q2_weight;
    }
  myfile.close();
  cout << "NTOT: " << NTOT << endl;

  for(int i = 0 ; i < N_events ; i++)
    {
      DVCS->GetEntry(i);
      Q2_norm->Fill(Q2, (NTOT / N_events));
    }
  
  Q2_norm->SetTitle("Normalization of uniform distribution");
  Q2_norm->Draw();
}
