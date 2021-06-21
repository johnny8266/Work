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




void fitting()
{
  
  string File_num[4] = {"background_only_crystal", "Co60", "Cs137", "Na22"};
  string Crystals_num[10] = {"010667", "030651", "050658", "050666", "010647", "040668", "994724", "994763", "080671", "160669"};
  string channel_num[10] = {"6", "6", "6", "6", "6", "6", "6", "6", "1", "1"};
  double T_420[10] = {0.6769, 0.6668, 0.6961, 0.6839, 0.6813, 0.7117, 0.6913, 0.6918, 0.7057, 0.687};
  TH1D *h1[40], *h2[40];
  vector<double> v1;


  TCanvas *c_4[4];
  for(int c = 0 ; c < 4 ; c++)
    {
      c_4[c] = new TCanvas(Form("c_4_%d", (c+1)), Form("c_4_%d", (c+1)), 1200, 900);
      c_4[c]->Divide(4,3);
    }


  // Set the fitting parameter
  //
  auto *self_fit = new TF1("self_fit", Co60_dis, 500., 4000., 7);
  self_fit->SetParLimits(0, 200., 800.);
  self_fit->SetParLimits(1, 1300., 2000.);
  self_fit->SetParLimits(2, 100., 500.);

  self_fit->SetParLimits(3, 100., 8000.);
  self_fit->SetParLimits(4, 100., 1000.);
  self_fit->SetParLimits(5, 0.1, 1000.);
  self_fit->SetParLimits(6, -100., 1000.);

  
  auto *Co60_fit = new TF1("Co60_fit", Co60_dis, 500., 4000., 7);
  Co60_fit->SetParLimits(0, 4000., 8000.);
  Co60_fit->SetParLimits(1, 1150., 1600.);
  Co60_fit->SetParLimits(2, 100., 500.);

  Co60_fit->SetParLimits(3, 100., 8000.);
  Co60_fit->SetParLimits(4, 100., 1000.);
  Co60_fit->SetParLimits(5, 0.1, 1000.);
  Co60_fit->SetParLimits(6, -100., 1000.);
  
  
  auto *Cs137_fit = new TF1("Cs137_fit", Cs137_dis, 200., 3000., 7);
  Cs137_fit->SetParLimits(0, 1000., 12000.);
  Cs137_fit->SetParLimits(1, 500., 720.);
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
  Na22_fit->SetParLimits(4, 900., 1600.);
  Na22_fit->SetParLimits(5, 300., 1000.);

  Na22_fit->SetParLimits(6, -500., 5000.);
  


  // Read the plot and fit
  //
  int count = 0;
  double pe_1 = 53.5, Npe_a = 0., Npe_b = 0., Np = 0., Opf, gain_PMT = 167;
  
  for(int Cn = 0 ; Cn < 10 ; Cn++)
    {
      for(int i = 0 ; i < 4 ; i++)
	{
	  c_4[i]->cd(Cn+1);
	  
	  string root_file_name = "./" + Crystals_num[Cn] + "/200sec_" + File_num[i] + ".root";
	  const char *rfn = root_file_name.c_str();
      
	  TFile* gfile = TFile::Open(rfn);
	  TDirectory* dir = gFile->GetDirectory("Energy");

	  string plot_table = "_R_EnergyCH" + channel_num[Cn] + "@DT5730_1204";
	  const char *p_t = plot_table.c_str();
	  dir->GetObject(p_t, h1[count]);
	  h1[count]->Draw();

	  TLegend *leg1;
	  leg1 = new TLegend(0.55, 0.6, 0.85, 0.8);
	  leg1->SetBorderSize(0);
	  
	  h2[count] = new TH1D(Form("h2_%d", count), Form("h2_%d", count), 4095, 0, 4095);

	  int background_n = 4 * Cn;
	  
	  if( (count % 4) == 0 )
	    *h2[count] = *h1[count];
	  if( (count % 4) > 0 )
	    *h2[count] = *h1[count] - *h1[background_n];
      
	  h2[count]->RebinX(3);

  
	  //	  if(i == 0)
	    
	  const char *title;
	  if(i == 0)
	    title = ("Bg  " + Crystals_num[Cn]).c_str();
	  else if(i == 1)
	    title = ("Co60  " + Crystals_num[Cn]).c_str();
	  else if(i == 2)
	    title = ("Cs137  " + Crystals_num[Cn]).c_str();
	  else if(i == 3)
	    title = ("Na22  " + Crystals_num[Cn]).c_str();

	  h2[count]->SetTitle(title);

	  double h_mean = h2[count]->GetMean(), h_std = h2[count]->GetStdDev(), fit_mean = 0., fit_width = 0.;
	  if(h_mean < 0.)
	    h_mean = -h_mean;

	  
	  if( (count % 4) == 0 )
	    {
	      h2[count]->Fit("self_fit", "R", "", 500., 4000.);
	      fit_mean = self_fit->GetParameter(1);
	      fit_width = self_fit->GetParameter(2);
	      leg1->AddEntry(Form("h2_%d", count), Form("Bg peak: %.1f", fit_mean), "");
	      //	      leg1->AddEntry(Form("h2_%d", count), Form("Width: %.1f", fit_width), "");
	    }
	  else if( (count % 4) == 1 )
	    {
	      h2[count]->Fit("Co60_fit", "R", "", 500., 4000.);
	      fit_mean = Co60_fit->GetParameter(1);
	      fit_width = Co60_fit->GetParameter(2);
	      leg1->AddEntry(Form("h2_%d", count), Form("Co60: %.1f", fit_mean), "");
	      //	      leg1->AddEntry(Form("h2_%d", count), Form("Width: %.1f", fit_width), "");
	    }
	  else if( (count % 4) == 2 )
	    {
	      h2[count]->Fit("Cs137_fit", "R", "", 200., 3000.);
	      fit_mean = Cs137_fit->GetParameter(1);
	      fit_width = Cs137_fit->GetParameter(2);

	      Npe_a = fit_mean / pe_1;
	      Npe_b = (fit_mean * 2.5 * 2. / 1.609) / gain_PMT / T_420[Cn] / 0.9 / 0.6617;
	      Np = Npe_b / 0.25;
	      v1.push_back(Npe_a);  v1.push_back(Npe_b);  v1.push_back(Np);
	      leg1->AddEntry(Form("h2_%d", count), Form("Cs137: %.1f", fit_mean), "");
	      //	      leg1->AddEntry(Form("h2_%d", count), Form("#pe / MeV : %.1f", Npe_b), "");
	      //	      leg1->AddEntry(Form("h2_%d", count), Form("N pe (a): %.1f", Npe_a), "");
	      //	      leg1->AddEntry(Form("h2_%d", count), Form("Width: %.1f", fit_width), "");
	    }
	  else if( (count % 4) == 3 )
	    {
	      //	  h2[count]->Fit("Na22_fit", "R", "", (h_mean - 1.*h_std), (h_mean + 1.*h_std));
	      h2[count]->Fit("Na22_fit", "R", "", 200., 3000.);

	      fit_mean = Na22_fit->GetParameter(1);
	      fit_width = Na22_fit->GetParameter(2);
	      Npe_a = fit_mean / pe_1;
	      Npe_b = (fit_mean * 2.5 * 2. / 1.609) / gain_PMT / T_420[Cn] / 0.9 / 0.511;
	      Np = Npe_b / 0.25;
	      v1.push_back(Npe_a);  v1.push_back(Npe_b);  v1.push_back(Np);
	      leg1->AddEntry(Form("h2_%d", count), Form("Na22_1: %.1f", fit_mean), "");
	      //	      leg1->AddEntry(Form("h2_%d", count), Form("#pe / MeV : %.1f", Npe_b), "");		   
	      //	      leg1->AddEntry(Form("h2_%d", count), Form("Width: %.1f", fit_width), "");

	      fit_mean = Na22_fit->GetParameter(4);
	      fit_width = Na22_fit->GetParameter(5);
	      Npe_a = fit_mean / pe_1;
	      Npe_b = (fit_mean * 2.5 * 2. / 1.609) / gain_PMT / T_420[Cn] / 0.9 / 1.275;
	      Np = Npe_b / 0.25;
	      Opf = T_420[Cn] * 0.9;
	      v1.push_back(Npe_a);  v1.push_back(Npe_b);  v1.push_back(Np);  v1.push_back(Opf);
	      leg1->AddEntry(Form("h2_%d", count), Form("Na22_2: %.1f", fit_mean), "");
	      //	      leg1->AddEntry(Form("h2_%d", count), Form("#pe / MeV : %.1f", Npe_b), "");
	      //	      leg1->AddEntry(Form("h2_%d", count), Form("Width: %.1f", fit_width), "");
	    }

	  
	  h2[count]->SetStats(0);
	  h2[count]->GetXaxis()->SetRangeUser(0, 4000);
	  h2[count]->Draw();
	  leg1->Draw("same");
	  
	  count++;
	  
	}
    }

  int print_l = 0;

  //  cout << "Npe_a  "
  cout << endl << endl << endl;
  for(int i = 0 ; i < 10 ; i++)
    {
      //      cout << Crystals_num[i] << ": ";
      for(int j = 0 ; j < 10 ; j++)
	{
	  cout << setprecision(4) << v1[print_l] << "  ";
	  print_l++;
	}
      cout << endl;
    }

  TCanvas *c5 = new TCanvas("c5", "c5", 900, 900);
  TGraph *gr1[3];
  TMultiGraph *mg = new TMultiGraph();
  TLegend *legend, *legend_mark[3];
  legend = new TLegend(0.2, 0.2, 0.6, 0.3);
  legend->SetBorderSize(0);

  
  for(int c = 0 ; c < 3 ; c++)
    {
      
      double x[10], y[10];
      print_l = 0;
      int Count = 0;
      
      for(int i = 0 ; i < 10 ; i++)
	{
	  print_l = 1 + c * 3 + 10 * i;
	  x[Count] = i+1;
	  y[Count] = v1[print_l];
	  Count++;
	}
      gr1[c] = new TGraph(10, x, y);
      gr1[c]->SetName(Form("gr1_%d", c));
      gr1[c]->GetXaxis()->SetTitle("No. Crystal");
      gr1[c]->GetYaxis()->SetTitle("N photons/MeV");
      gr1[c]->GetXaxis()->SetRangeUser(0., 12.);	  
      gr1[c]->GetHistogram()->SetMaximum(40.);
      gr1[c]->GetHistogram()->SetMinimum(0.);
      gr1[c]->SetMarkerStyle(c+2);
      gr1[c]->SetMarkerSize(3);
      gr1[c]->SetMarkerColor(c+2);
      //      legend->AddEntry(Form("gr1_%d", c), Form("gr1_%d", c), "P");
      
      mg->Add( (gr1[c]), "AP");

    }
  mg->GetXaxis()->SetTitle("No. Crystal");
  mg->GetYaxis()->SetTitle("N photons/MeV");
  mg->GetXaxis()->SetRangeUser(0., 12.);	  
  mg->GetHistogram()->SetMaximum(40.);
  mg->GetHistogram()->SetMinimum(0.);
  mg->Draw("a");

  legend->AddEntry("gr1_0", "Cs137 [661 keV]", "P");
  legend->AddEntry("gr1_1", "Na22 [511 keV]", "P");
  legend->AddEntry("gr1_2", "Na22 [1275 keV]", "P");
  legend->Draw();
  
  /*
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
  */
  
  
}
