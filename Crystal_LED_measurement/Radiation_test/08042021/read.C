double GausM(double *x, double *par)
{
  return par[0] * exp(-0.5 * TMath::Power(((x[0] - par[1]) / par[2]), 2)) + par[3];
}

void read()
{
  
  string File_num[4] = {"background_only_crystal", "Co60", "Cs137", "Na22"};
  TH1D *h1[4], *h2[4], *background;
  
  for(int i = 0 ; i < 4 ; i++)
    {
      //      string root_file_name = "./080671_Crystal/" + File_num[i] + "/plot.root";
      string root_file_name = "./160669_Crystal/" + File_num[i] + "/plot.root";
      const char *rfn = root_file_name.c_str();
      
      TFile* gfile = TFile::Open(rfn);
      TDirectory* dir = gFile->GetDirectory("Energy");
      
      dir->GetObject("_R_EnergyCH0@DT5730_1204", h1[i]);
    }

  TFile* gfile = TFile::Open("./background_without_crystal/background.root");
  TDirectory* dir = gFile->GetDirectory("Energy");
      
  dir->GetObject("_R_EnergyCH0@DT5730_1204", background);
      

  TCanvas *c1 = new TCanvas("c1", "c1", 900, 900);
  c1->Divide(2,2);

  auto *pe_peek_Cs = new TF1("pe_peek_Cs", GausM, 0., 500., 4);
  pe_peek_Cs->SetParLimits(0, 100., 800.);
  pe_peek_Cs->SetParLimits(1, 240., 260.);
  pe_peek_Cs->SetParLimits(2, 10., 50.);
  pe_peek_Cs->SetParLimits(3, -200, 1000.);

  auto *pe_peek_general = new TF1("pe_peek_general", GausM, 0., 500., 4);
  pe_peek_general->SetParLimits(0, 100., 25000.);
  pe_peek_general->SetParLimits(1, 330., 450.);
  pe_peek_general->SetParLimits(2, 10., 100.);
  pe_peek_general->SetParLimits(3, -400, 2000.);

  
  for(int j = 0 ; j < 4 ; j++)
    {
      TLegend *leg1;
      leg1 = new TLegend(0.55, 0.6, 0.85, 0.8);
      leg1->SetBorderSize(0);

      c1->cd(j+1);

      h2[j] = new TH1D(Form("h2_%d", j), Form("h2_%d", j), 4095, 0, 4095);

      if(j == 0)
	*h2[j] = *h1[j] - *background;
      else if(j > 0)
	*h2[j] = *h1[j] - *h1[0];
      
      h2[j]->RebinX(5);
      const char *title = File_num[j].c_str();
      h2[j]->SetTitle(title);

      double h_mean = h2[j]->GetMean(), h_std = h2[j]->GetStdDev(), fit_mean = 0.;
      if(h_mean < 0.)
	h_mean = -h_mean;
      
      if(j==2)
	{
	  h2[j]->Fit("pe_peek_Cs", "R", "", (h_mean - 1.*h_std), (h_mean + 1.*h_std));
	  fit_mean = pe_peek_Cs->GetParameter(1);
	  cout << " ??? " << endl;
	}
      else
	{
	  h2[j]->Fit("pe_peek_general", "R", "", (h_mean - 1.*h_std), (h_mean + 1.*h_std));
	  fit_mean = pe_peek_general->GetParameter(1);
	}
      
      h2[j]->SetStats(0);
      h2[j]->GetXaxis()->SetRangeUser(0, 1000);
      h2[j]->Draw();

      leg1->AddEntry(Form("h2_%d", j), Form("Mean: %.1f", fit_mean), "");
      leg1->Draw("same");
      
    }

  
  
}
