void norm_uniform()
{
  TFile *hfile = new TFile("result_uniform_27_11_2020.root");
  TTree *T = (TTree*)hfile->Get("T");
  Double_t phi, phi_def, e1_S_angle, p1_S_angle, photon_S_angle, Q2, xb, t_var, xsec, psf;
  Double_t e1_px, e1_py, e1_pz, e1_E, Vg_px, Vg_py, Vg_pz, Vg_E;
  Double_t p1_px, p1_py, p1_pz, p1_E, g_px, g_py, g_pz, g_E;
  Double_t e_sc_theta;
  TVector3 z_axis(0., 0., 1.), electron;
  Int_t N_events = (Int_t)T->GetEntries(), count=0;
  
  //  cout << "N_events: " << N_events << endl;
  //  cout << "count: " << count << endl;


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

  TH1F* Q2_norm = new TH1F("Q2_norm", "Q2_norm", 90, 0., 45.);
  TH2F* Q2_xsec = new TH2F("Q2_xsec", "Q2_xsec", 700, 0., 7., 90, 0., 45.);
  TH2F* Q2_psf = new TH2F("Q2_psf", "Q2_psf", 800, 0., 80000., 90, 0., 45.);


  for(int i = 0 ; i < N_events ; i++)
    {
      T->GetEntry(i);

      Q2_xsec->Fill(xsec, Q2);
      Q2_psf->Fill(psf, Q2);
      Q2_weight = SLdt / N_events * xsec * psf;
      Q2_norm->Fill(Q2, Q2_weight);
    }

  TCanvas *c = new TCanvas ("c","c", 800, 800);
  //  c->Divide(1,2);
  
  Q2_norm->SetTitle("Normalization of uniform distribution");
  Q2_norm->Draw();

  /*
  Q2_xsec->SetTitle("Q2_vs_xsec");
  Q2_xsec->SetStats(0);
  Q2_xsec->GetXaxis()->SetTitle("xsec");    Q2_xsec->GetYaxis()->SetTitle("Q2");
  c->cd(1);
  Q2_xsec->Draw("colorz");
  Q2_psf->SetTitle("Q2_vs_psf");
  Q2_psf->SetStats(0);
  Q2_psf->GetXaxis()->SetTitle("psf");    Q2_psf->GetYaxis()->SetTitle("Q2");
  c->cd(2);
  Q2_psf->Draw("colorz");
  */

  
}
