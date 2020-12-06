void Evt_gen_read()
{
  TFile *hfile = new TFile("result.root");
  TTree *T = (TTree*)hfile->Get("T");
  TCanvas *c1 = new TCanvas();
  Double_t phi, e1_S_angle, p1_S_angle, photon_S_angle, Q2, xb, t_var;
  Double_t e1_px, e1_py, e1_pz, e1_E, Vg_px, Vg_py, Vg_pz, Vg_E;
  Double_t p1_px, p1_py, p1_pz, p1_E, g_px, g_py, g_pz, g_E;
  Int_t N_events = (Int_t)T->GetEntries();


  T->SetBranchAddress("Q2", &Q2);
  T->SetBranchAddress("xb", &xb);
  T->SetBranchAddress("t_var", &t_var);
  /*  T->SetBranchAddress("e1_S_angle", &e1_S_angle);
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
  T->SetBranchAddress("g_E", &g_E);*/
  
  
  TH1F* Q2_distribution = new TH1F("Q2_distribution", "Q2_distribution", 51, 0., 51.);
  TH1F* xb_distribution = new TH1F("xb_distribution", "xb_distribution", 101, 0., 0.101);
  TH1F* t_distribution = new TH1F("t_distribution", "t_distribution", 200, -2000., 0.);
  TH1F* e_Scattering_angle = new TH1F("e_Scattering_angle", "e_Scattering_angle", 200, -0.01, 0.01);
  TH1F* p_Scattering_angle = new TH1F("p_Scattering_angle", "p_Scattering_angle", 630, -3.15, 3.15);
  TH1F* photon_Scattering_angle = new TH1F("photon_Scattering_angle", "photon_Scattering_angle", 630, -3.15, 3.15);
  TH2F* Q2_vs_xb = new TH2F("Q2_vs_xb", "Q2_vs_xb", 51, 0., 51., 101, 0., 0.101);
  TH2F* theta_p_vs_Q2 = new TH2F("theta_p_vs_Q2", "theta_p_vs_Q2", 320, -1.6, 1.6, 51, 0., 51.);
  TH2F* theta_p_vs_xb = new TH2F("theta_p_vs_xb", "theta_p_vs_xb", 320, -1.6, 1.6, 101, 0., 0.101);
  TH2F* theta_e_vs_Q2 = new TH2F("theta_e_vs_Q2", "theta_e_vs_Q2", 175, 0., 0.035, 51, 0., 51.);
  TH2F* theta_e_vs_xb = new TH2F("theta_e_vs_xb", "theta_e_vs_xb", 175, 0., 0.035, 101, 0., 0.101);
  TH2F* theta_photon_vs_Q2 = new TH2F("theta_photon_vs_Q2", "theta_photon_vs_Q2", 320, -1.6, 1.6, 51, 0., 51.);
  TH2F* theta_photon_vs_xb = new TH2F("theta_photon_vs_xb", "theta_photon_vs_xb", 320, -1.6, 1.6, 101, 0., 0.101);


  
  for(int i = 0 ; i < N_events ; i++)
    {
      T->GetEntry(i);

      cout << Q2 <<  "  " << xb << "  " << t_var << endl;
      /*
      p_Scattering_angle->Fill(p1_S_angle);
      e_Scattering_angle->Fill(e1_S_angle);
      photon_Scattering_angle->Fill(photon_S_angle);
      Q2_distribution->Fill(Q2);  xb_distribution->Fill(xb);  t_distribution->Fill(t_var);  Q2_vs_xb->Fill(Q2, xb);
      theta_e_vs_Q2->Fill(e1_S_angle, Q2);  theta_e_vs_xb->Fill(e1_S_angle, xb);
      theta_p_vs_Q2->Fill(p1_S_angle, Q2);  theta_p_vs_xb->Fill(p1_S_angle, xb);
      theta_photon_vs_Q2->Fill(photon_S_angle, Q2);  theta_photon_vs_xb->Fill(photon_S_angle, xb);
      */
    }


  /*
  TFile *histo_save = new TFile("./pics/histo.root", "RECREATE");
  
  Q2_vs_xb->SetStats(0);
  Q2_vs_xb->GetXaxis()->SetTitle("Q2");  Q2_vs_xb->GetYaxis()->SetTitle("Xb");
  Q2_vs_xb->GetXaxis()->CenterTitle(true);  Q2_vs_xb->GetYaxis()->CenterTitle(true);
  Q2_vs_xb->Draw("colorz");
  Q2_vs_xb->Write();
  c1->SaveAs("./pics/Q2_vs_xb.png");
  
  theta_p_vs_Q2->SetStats(0);
  theta_p_vs_Q2->GetXaxis()->SetTitle("Proton scattering angle [rad]");  theta_p_vs_Q2->GetYaxis()->SetTitle("Q2");
  theta_p_vs_Q2->GetXaxis()->CenterTitle();  theta_p_vs_Q2->GetYaxis()->CenterTitle();
  theta_p_vs_Q2->Draw("colorz");
  theta_p_vs_Q2->Write();
  c1->SaveAs("./pics/theta_p_vs_Q2.png");

  theta_p_vs_xb->SetStats(0);
  theta_p_vs_xb->GetXaxis()->SetTitle("Proton scattering angle [rad]");  theta_p_vs_xb->GetYaxis()->SetTitle("Xb");
  theta_p_vs_xb->GetXaxis()->CenterTitle();  theta_p_vs_xb->GetYaxis()->CenterTitle();
  theta_p_vs_xb->Draw("colorz");
  theta_p_vs_xb->Write();
  c1->SaveAs("./pics/theta_p_vs_xb.png");

  theta_e_vs_Q2->SetStats(0);
  theta_e_vs_Q2->GetXaxis()->SetTitle("Electron scattering angle [rad]");  theta_e_vs_Q2->GetYaxis()->SetTitle("Q2");
  theta_e_vs_Q2->GetXaxis()->CenterTitle();  theta_e_vs_Q2->GetYaxis()->CenterTitle();
  theta_e_vs_Q2->Draw("colorz");
  theta_e_vs_Q2->Write();
  c1->SaveAs("./pics/theta_e_vs_Q2.png");

  theta_e_vs_xb->SetStats(0);
  theta_e_vs_xb->GetXaxis()->SetTitle("Electron scattering angle [rad]");  theta_e_vs_xb->GetYaxis()->SetTitle("Xb");
  theta_e_vs_xb->GetXaxis()->CenterTitle();  theta_e_vs_xb->GetYaxis()->CenterTitle();
  theta_e_vs_xb->Draw("colorz");
  theta_e_vs_xb->Write();
  c1->SaveAs("./pics/theta_e_vs_xb.png");

  theta_photon_vs_Q2->SetStats(0);
  theta_photon_vs_Q2->GetXaxis()->SetTitle("Photon scattering angle [rad]");  theta_photon_vs_Q2->GetYaxis()->SetTitle("Q2");
  theta_photon_vs_Q2->GetXaxis()->CenterTitle();  theta_photon_vs_Q2->GetYaxis()->CenterTitle();
  theta_photon_vs_Q2->Draw("colorz");
  theta_photon_vs_Q2->Write();
  c1->SaveAs("./pics/theta_photon_vs_Q2.png");

  theta_photon_vs_xb->SetStats(0);
  theta_photon_vs_xb->GetXaxis()->SetTitle("Photon scattering angle [rad]");  theta_photon_vs_xb->GetYaxis()->SetTitle("Xb");
  theta_photon_vs_xb->GetXaxis()->CenterTitle();  theta_photon_vs_xb->GetYaxis()->CenterTitle();
  theta_photon_vs_xb->Draw("colorz");
  theta_photon_vs_xb->Write();
  c1->SaveAs("./pics/theta_photon_vs_xb.png");

  histo_save->Write();
  */
}
