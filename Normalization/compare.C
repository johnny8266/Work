#include <iostream>
using namespace std;

void compare()
{
  //==============================================
  // Normalization of the uniform DVCS generator
  //==============================================
  
  TFile *hfile = new TFile("./root_file/uniform_20012021.root");
  TTree *T = (TTree*)hfile->Get("T");
  Double_t phi, phi_def, e1_S_angle, p1_S_angle, photon_S_angle, Q2, xb, t_var, xsec, psf;
  Double_t e1_px, e1_py, e1_pz, e1_E, Vg_px, Vg_py, Vg_pz, Vg_E;
  Double_t p1_px, p1_py, p1_pz, p1_E, g_px, g_py, g_pz, g_E;
  Double_t e_sc_theta;
  Int_t N_events = (Int_t)T->GetEntries(), count=0;
  
  TH1F* Q2_uni = new TH1F("Q2_uni", "Q2_uni", 90, 0., 45.);
  TH1F* xb_uni = new TH1F("xb_uni", "xb_uni", 110, 0., 0.11);
  TH1F* t_uni = new TH1F("t_uni", "t_uni", 64, -3.1, 0.1);
  TH1F* phi_uni = new TH1F("phi_uni", "phi_uni", 63, 0., 6.3);
  TH1F* xsec_uni = new TH1F("xsec_uni", "xsec_uni", 1000, 0., 1.);
  TH1F* psf_uni = new TH1F("psf_uni", "psf_uni", 20, 0., 100.);
  TH1F* Q2_norm_uni = new TH1F("Q2_norm_uni", "Q2_norm_uni", 90, 2., 20.);
  TH1F* Xb_norm_uni = new TH1F("Xb_norm_uni", "Xb_norm_uni", 99, 0.001, 0.1);
  TH1F* t_norm_uni = new TH1F("t_norm_uni", "t_norm_uni", 100, -2., 0.);
  TH1F* phi_norm_uni = new TH1F("phi_norm_uni", "phi_norm_uni", 63, 0., 6.3);
  
  TH1F* Q2_foam = new TH1F("Q2_foam", "Q2_foam", 90, 0., 45.);  
  TH1F* xb_foam = new TH1F("xb_foam", "xb_foam", 110, 0., 0.11);
  TH1F* t_foam = new TH1F("t_foam", "t_foam", 48, -1.1, 0.1);
  TH1F* phi_foam = new TH1F("phi_foam", "phi_foam", 63, 0., 6.3);
  TH1F* xsec_foam = new TH1F("xsec_foam", "xsec_foam", 1000, 0., 1.);
  TH1F* psf_foam = new TH1F("psf_foam", "psf_foam", 20, 0., 100.);
  TH1F* Q2_norm_foam = new TH1F("Q2_norm_foam", "Q2_norm_foam", 90, 2., 20.);
  TH1F* Xb_norm_foam = new TH1F("Xb_norm_foam", "Xb_norm_foam", 99, 0.001, 0.1);
  TH1F* t_norm_foam = new TH1F("t_norm_foam", "t_norm_foam", 100, -2., 0.);
  TH1F* phi_norm_foam = new TH1F("phi_norm_foam", "phi_norm_foam", 63, 0., 6.3);

  TH2F* foam_Q2_xsec = new TH2F("foam_Q2_xsec", "foam_Q2_xsec", 10, 0., 5., 500, 0., 5.);
  TH2F* uni_Q2_xsec = new TH2F("uni_Q2_xsec", "uni_Q2_xsec", 10, 0., 5., 500, 0., 5.);
  
  T->SetBranchAddress("Q2", &Q2);
  T->SetBranchAddress("xb", &xb);
  T->SetBranchAddress("t_var", &t_var);
  T->SetBranchAddress("phi", &phi);
  T->SetBranchAddress("phi_def", &phi_def);
  T->SetBranchAddress("psf", &psf);
  T->SetBranchAddress("e1_S_angle", &e1_S_angle);
  T->SetBranchAddress("p1_S_angle", &p1_S_angle);
  T->SetBranchAddress("photon_S_angle", &photon_S_angle);
  T->SetBranchAddress("e1_px", &e1_px);
  T->SetBranchAddress("e1_py", &e1_py);
  T->SetBranchAddress("e1_pz", &e1_pz);
  T->SetBranchAddress("e1_E", &e1_E);
  T->SetBranchAddress("p1_px", &p1_px);
  T->SetBranchAddress("p1_py", &p1_py);
  T->SetBranchAddress("p1_pz", &p1_pz);
  T->SetBranchAddress("p1_E", &p1_E);
  T->SetBranchAddress("Vg_px", &Vg_px);
  T->SetBranchAddress("Vg_py", &Vg_py);
  T->SetBranchAddress("Vg_pz", &Vg_pz);
  T->SetBranchAddress("Vg_E", &Vg_E);
  T->SetBranchAddress("g_px", &g_px);
  T->SetBranchAddress("g_py", &g_py);
  T->SetBranchAddress("g_pz", &g_pz);
  T->SetBranchAddress("g_E", &g_E);
  T->SetBranchAddress("xsec", &xsec);

  Double_t SLdt=10., Q2_weight=0.;

  cout << "Event numbers of uniform generator: " << N_events << endl;
  for(int i = 0 ; i < N_events ; i++)
    {
      T->GetEntry(i);
      Q2_weight = (SLdt * xsec * psf * 1000000.) / N_events;

      //      if ( (i % 50 == 0) && Q2 < 2. ) cout << xsec << endl;
      if(xsec > 0.1)
	count++;
      
      Q2_uni->Fill(Q2);
      xb_uni->Fill(xb);
      t_uni->Fill(t_var);
      phi_uni->Fill(phi);
      xsec_uni->Fill(xsec);
      psf_uni->Fill(psf);
      uni_Q2_xsec->Fill(Q2, xsec);
      
      Q2_norm_uni->Fill(Q2, Q2_weight);
      Xb_norm_uni->Fill(xb, Q2_weight);
      t_norm_uni->Fill(t_var, Q2_weight);
      phi_norm_uni->Fill(phi, Q2_weight);
      //      cout << Q2 << "  " << Q2_weight << endl;
    }

  cout << "Uniform part finish and xsec > 0.1: " << count << endl << endl;
  count = 0;
  
  
  //======================================
  // Normalization of the foam generator
  //======================================
  TFile *gfile = new TFile("./root_file/foam_imposed_19012021.root");
  TTree *DVCS = (TTree*)gfile->Get("DVCS");  
  N_events = (Int_t)DVCS->GetEntries();

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

  Double_t NTOT=0.;

  cout << "Event numbers of foam generator: " << N_events << endl;
  for(int i = 0 ; i < N_events ; i++)
    {
      DVCS->GetEntry(i);
      NTOT = (SLdt * xsec * psf * 1000000.) / N_events;
      //      if( i % 1000 == 0 ) cout << NTOT << endl;

      if(xsec > 0.1)
	count++;

      Q2_foam->Fill(Q2);
      xb_foam->Fill(xb);
      t_foam->Fill(t_var);
      phi_foam->Fill(phi);
      xsec_foam->Fill(xsec);
      psf_foam->Fill(psf);
      foam_Q2_xsec->Fill(Q2, xsec);
      
      Q2_norm_foam->Fill(Q2, NTOT);
      Xb_norm_foam->Fill(xb, NTOT);
      t_norm_foam->Fill(t_var, NTOT);
      phi_norm_foam->Fill(phi, NTOT);
    }

  cout << "Foam part finish and xsec > 0.1: " << count << endl << endl;
  count = 0;
  
 
  TCanvas *c = new TCanvas ("c","c", 1000, 1000);
  c->Divide(2,2);
  c->cd(1);
  Q2_norm_uni->GetYaxis()->SetMaxDigits(2);
  Q2_norm_uni->Draw();
  c->cd(2);
  Xb_norm_uni->GetYaxis()->SetMaxDigits(2);
  Xb_norm_uni->Draw();
  c->cd(3);
  t_norm_uni->GetYaxis()->SetMaxDigits(2);
  t_norm_uni->Draw();
  c->cd(4);
  phi_norm_uni->GetYaxis()->SetMaxDigits(2);
  phi_norm_uni->Draw();

  
  TCanvas *c1 = new TCanvas ("c1","c1", 1000, 1000);
  c1->Divide(2,2);
  c1->cd(1);
  Q2_norm_foam->GetYaxis()->SetMaxDigits(3);
  Q2_norm_foam->Draw();
  c1->cd(2);
  Xb_norm_foam->GetYaxis()->SetMaxDigits(3);
  Xb_norm_foam->Draw();
  c1->cd(3);
  t_norm_foam->GetYaxis()->SetMaxDigits(3);
  t_norm_foam->Draw();
  c1->cd(4);
  phi_norm_foam->GetYaxis()->SetMaxDigits(3);
  phi_norm_foam->Draw();

  /*
  TCanvas *c2 = new TCanvas ("c2","c2", 1000, 1000);
  c2->Divide(2,2);
  c2->cd(1);
  Q2_foam->Draw();
  c2->cd(2);
  xb_foam->Draw();
  c2->cd(3);
  t_foam->Draw();
  c2->cd(4);
  phi_foam->Draw();
  

  TCanvas *c3 = new TCanvas ("c3","c3", 1000, 1000);
  c3->Divide(2,2);
  c3->cd(1);
  Q2_uni->Draw();
  c3->cd(2);
  xb_uni->Draw();
  c3->cd(3);
  t_uni->Draw();
  c3->cd(4);
  phi_uni->Draw();
  */    

  TCanvas *c4 = new TCanvas ("c4","c4", 1000, 1000);
  c4->Divide(2,2);
  c4->cd(1);
  xsec_uni->Draw();
  c4->cd(2);
  psf_uni->Draw();
  c4->cd(3);
  xsec_foam->Draw();
  c4->cd(4);
  psf_foam->Draw();
  
  
  TCanvas *c5 = new TCanvas ("c5","c5", 1000, 500);
  c5->Divide(2,1);
  c5->cd(1);
  uni_Q2_xsec->Draw("colorz");
  c5->cd(2);
  foam_Q2_xsec->Draw("colorz");
  
  
}
