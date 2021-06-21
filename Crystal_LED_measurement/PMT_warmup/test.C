double GausM(double *x, double *par)
{
  return par[0] * exp(-0.5 * TMath::Power(((x[0] - par[1]) / par[2]), 2)) + par[3];
}



void test()
{
  TH1D* h1[6];
  TH1F *h2 = new TH1F("h2", "h2", 35, 0., 35.);
  string File_num[6] = {"5", "10", "15", "20", "25", "30"};
  double x[6], y[6];
  
  TCanvas *c1 = new TCanvas("c1", "c1", 1200, 800);
  c1->Divide(3,2);


  auto *pe_peek_Co = new TF1("pe_peek_Co", GausM, 800., 1200., 4);
  pe_peek_Co->SetParLimits(0, 2000., 4000.);
  pe_peek_Co->SetParLimits(1, 850., 1100.);
  pe_peek_Co->SetParLimits(2, 100., 200.);
  pe_peek_Co->SetParLimits(3, -500, 500.); 

  
  for(int i = 0 ; i < 6 ; i++)
    {
      c1->cd(i+1);

      string root_file_name = "./" + File_num[i] + "mins.root";
      const char *rfn = root_file_name.c_str();
      
      TFile* gfile = TFile::Open(rfn);
      TDirectory* dir = gFile->GetDirectory("Energy");

      dir->GetObject("_R_EnergyCH0@DT5730_1204", h1[i]);
  
      h1[i]->RebinX(9);
      h1[i]->Fit("pe_peek_Co", "R", "", 800., 1200.);

      
      string plot_name = "PMT warmed up " + File_num[i] + " mins";
      const char *p_name = plot_name.c_str();
      
      h1[i]->SetTitle(p_name);
      h1[i]->SetStats(0);
      h1[i]->GetXaxis()->SetRangeUser(400, 1600);
      //      h1[i]->GetXaxis()->SetRangeUser(800, 1100);
      h1[i]->Draw();


      double mean_value = pe_peek_Co->GetParameter(1);
      double t = 5. * (i+1);
      x[i] = t;
      y[i] = mean_value;
      //      h2->Fill(t, mean_value);
      
    }

  TCanvas *c2 = new TCanvas("c2", "c2", 600, 600);
  TGraph* gr = new TGraph(6,x,y);
  gr->SetMarkerStyle(2);
  gr->SetMarkerSize(3);
  gr->SetTitle("PMT warm up");
  gr->GetXaxis()->SetTitle("Warming time[m]");
  gr->GetYaxis()->SetTitle("Collected charge");
  gr->Draw("ap");
  
}
