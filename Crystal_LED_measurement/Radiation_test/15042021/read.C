#include "TMath.h"


double Co60_dis(double *x, double *par)
{
  double peak = par[0] * exp(-0.5 * TMath::Power(((x[0] - par[1]) / par[2]), 2));
  double background = par[3] * TMath::Exp(-(x[0] - par[4]) / par[5]) + par[6];
  return peak + background;
}

double Cs137_dis(double *x, double *par)
{
  double peak = par[0] * exp(-0.5 * TMath::Power(((x[0] - par[1]) / par[2]), 2));
  double background = par[3] * TMath::Exp(-(x[0] - par[4]) / par[5]) + par[6];
  return peak + background;
}

double Na22_dis(double *x, double *par)
{
  double peak_1 = par[0] * exp(-0.5 * TMath::Power(((x[0] - par[1]) / par[2]), 2));
  double peak_2 = par[3] * exp(-0.5 * TMath::Power(((x[0] - par[4]) / par[5]), 2));
  double background = par[6];
  return peak_1 + peak_2 + background;
}




void read()
{
  
  string File_num[4] = {"background_only_crystal", "Co60", "Cs137", "Na22"};
  TH1D *h1[4], *h2[4], *background;
  
  for(int i = 0 ; i < 4 ; i++)
    {
      //      string root_file_name = "./080671_Crystal/" + File_num[i] + "/plot.root";
      string root_file_name = "./160669_Crystal/200sec_" + File_num[i] + ".root";
      const char *rfn = root_file_name.c_str();
      
      TFile* gfile = TFile::Open(rfn);
      TDirectory* dir = gFile->GetDirectory("Energy");
      
      dir->GetObject("_R_EnergyCH1@DT5730_1204", h1[i]);
    }
  

  TCanvas *c1 = new TCanvas("c1", "c1", 900, 900);
  c1->Divide(2,2);


  auto *self_fit = new TF1("self_fit", Co60_dis, 500., 4000., 7);
  self_fit->SetParLimits(0, 200., 800.);
  self_fit->SetParLimits(1, 1500., 2000.);
  self_fit->SetParLimits(2, 100., 500.);

  self_fit->SetParLimits(3, 100., 8000.);
  self_fit->SetParLimits(4, 100., 1000.);
  self_fit->SetParLimits(5, 0.1, 1000.);
  self_fit->SetParLimits(6, -100., 1000.);

  
  auto *Co60_fit = new TF1("Co60_fit", Co60_dis, 500., 4000., 7);
  Co60_fit->SetParLimits(0, 4000., 8000.);
  Co60_fit->SetParLimits(1, 1250., 1600.);
  Co60_fit->SetParLimits(2, 100., 500.);

  Co60_fit->SetParLimits(3, 100., 8000.);
  Co60_fit->SetParLimits(4, 100., 1000.);
  Co60_fit->SetParLimits(5, 0.1, 1000.);
  Co60_fit->SetParLimits(6, -100., 1000.);
  
  
  auto *Cs137_fit = new TF1("Cs137_fit", Cs137_dis, 200., 3000., 7);
  Cs137_fit->SetParLimits(0, 1000., 12000.);
  Cs137_fit->SetParLimits(1, 600., 720.);
  Cs137_fit->SetParLimits(2, 150., 400.);

  Cs137_fit->SetParLimits(3, 100., 5000.);
  Cs137_fit->SetParLimits(4, 100., 1000.);
  Cs137_fit->SetParLimits(5, 0.1, 5000.);
  Cs137_fit->SetParLimits(6, -100., 1000.);
    

  auto *Na22_fit = new TF1("Na22_fit", Na22_dis, 200., 3000., 7);
  Na22_fit->SetParLimits(0, 5000., 12000.);
  Na22_fit->SetParLimits(1, 400., 600.);
  Na22_fit->SetParLimits(2, 100., 500.);
  
  Na22_fit->SetParLimits(3, 800., 2000.);
  Na22_fit->SetParLimits(4, 1100., 1600.);
  Na22_fit->SetParLimits(5, 300., 1000.);

  Na22_fit->SetParLimits(6, -500., 5000.);
  
  
  
  for(int j = 0 ; j < 4 ; j++)
    {
      TLegend *leg1;
      leg1 = new TLegend(0.55, 0.55, 0.9, 0.8);
      leg1->SetBorderSize(0);

      c1->cd(j+1);

      h2[j] = new TH1D(Form("h2_%d", j), Form("h2_%d", j), 4095, 0, 4095);

      if(j == 0)
      	*h2[j] = *h1[j];
      if(j > 0)
	*h2[j] = *h1[j] - *h1[0];
      
      h2[j]->RebinX(3);
      const char *title = File_num[j].c_str();
      h2[j]->SetTitle(title);

      double h_mean = h2[j]->GetMean(), h_std = h2[j]->GetStdDev(), fit_mean = 0., fit_width = 0.;
      if(h_mean < 0.)
	h_mean = -h_mean;


      if(j==0)
	{
	  h2[j]->Fit("self_fit", "R", "", 500., 4000.);
	  fit_mean = self_fit->GetParameter(1);
	  fit_width = self_fit->GetParameter(2);
	  leg1->AddEntry(Form("h2_%d", j), Form("First peak: %.1f", fit_mean), "");
	  leg1->AddEntry(Form("h2_%d", j), Form("Width: %.1f", fit_width), "");
	}
      else if(j==1)
	{
	  h2[j]->Fit("Co60_fit", "R", "", 500., 4000.);
	  fit_mean = Co60_fit->GetParameter(1);
	  fit_width = Co60_fit->GetParameter(2);
	  leg1->AddEntry(Form("h2_%d", j), Form("First peak: %.1f", fit_mean), "");
	  leg1->AddEntry(Form("h2_%d", j), Form("Width: %.1f", fit_width), "");
	}
      else if(j==2)
	{
	  h2[j]->Fit("Cs137_fit", "R", "", 200., 3000.);
	  fit_mean = Cs137_fit->GetParameter(1);
	  fit_width = Cs137_fit->GetParameter(2);
	  leg1->AddEntry(Form("h2_%d", j), Form("First peak: %.1f", fit_mean), "");
	  leg1->AddEntry(Form("h2_%d", j), Form("Width: %.1f", fit_width), "");
	}
      else if(j==3)
	{
	  //	  h2[j]->Fit("Na22_fit", "R", "", (h_mean - 1.*h_std), (h_mean + 1.*h_std));
	  h2[j]->Fit("Na22_fit", "R", "", 200., 3000.);
	  fit_mean = Na22_fit->GetParameter(1);
	  fit_width = Na22_fit->GetParameter(2);
	  leg1->AddEntry(Form("h2_%d", j), Form("First peak: %.1f", fit_mean), "");
	  leg1->AddEntry(Form("h2_%d", j), Form("Width: %.1f", fit_width), "");
	  fit_mean = Na22_fit->GetParameter(4);
	  fit_width = Na22_fit->GetParameter(5);
	  leg1->AddEntry(Form("h2_%d", j), Form("Second peak: %.1f", fit_mean), "");
	  leg1->AddEntry(Form("h2_%d", j), Form("Width: %.1f", fit_width), "");
	}
      
      
      h2[j]->SetStats(0);
      h2[j]->GetXaxis()->SetRangeUser(0, 4000);
      h2[j]->Draw();
      leg1->Draw("same");
      
    }

  
  
}
