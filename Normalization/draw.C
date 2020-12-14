void draw()
{
  TFile *hfile = new TFile("C1000_S200_B8_1M_2nd.root");
  TTree *T = (TTree*)hfile->Get("T");
  Double_t Q2, xb, t_var, phi, xsec, psf;
  Int_t N_events = (Int_t)T->GetEntries(), count=0;

  TH1F* Q2_foam = new TH1F("Q2_foam", "Q2_foam", 60, 0., 15.);  
  TH1F* xb_foam = new TH1F("xb_foam", "xb_foam", 80, 0., 0.04);
  TH1F* t_foam = new TH1F("t_foam", "t_foam", 48, -1.1, 0.1);
  TH1F* phi_foam = new TH1F("phi_foam", "phi_foam", 63, 0., 6.3);
  TH1F* xsec_foam = new TH1F("xsec_foam", "xsec_foam", 300, 0., 3000.);
  
  T->SetBranchAddress("Q2", &Q2);
  T->SetBranchAddress("xb", &xb);
  T->SetBranchAddress("t_var", &t_var);
  T->SetBranchAddress("phi", &phi);
  T->SetBranchAddress("xsec", &xsec);

  for(int i = 0 ; i < N_events ; i++)
    {
      T->GetEntry(i);

      Q2_foam->Fill(Q2);
      xb_foam->Fill(xb);
      t_foam->Fill(t_var);
      phi_foam->Fill(phi);
      xsec_foam->Fill(xsec);
    }

  TCanvas *c1 = new TCanvas ("c1","c1", 400, 400);
  xsec_foam->SetStats(0);
  xsec_foam->Draw();
  
  TCanvas *c2 = new TCanvas ("c2","c2", 800, 800);
  c2->Divide(2,2);
  c2->cd(1);
  Q2_foam->SetStats(0);
  Q2_foam->Draw();
  c2->cd(2);
  xb_foam->SetStats(0);
  xb_foam->Draw();
  c2->cd(3);
  t_foam->SetStats(0);
  t_foam->Draw();
  c2->cd(4);
  phi_foam->SetStats(0);
  phi_foam->Draw();
}
