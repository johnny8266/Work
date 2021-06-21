#include <iostream>
#include <string>
#include <stdlib.h>


double GausM(double *x, double *par)
{
  return par[0] * exp(-0.5 * TMath::Power(((x[0] - par[1]) / par[2]), 2)) + par[3];
}


double Single_pe(double *x, double *par)
{
  double pedestal = par[0] * exp(-0.5 * TMath::Power(((x[0] - par[1]) / par[2]), 2));
  double pe_1 = par[3] * exp(-0.5 * TMath::Power(((x[0] - par[4]) / par[5]), 2));
  //  double background = par[6] / x[0];
  return pedestal + pe_1;
}

double Multi_pe_2(double *x, double *par)
{
  double pedestal = par[0] * exp(-0.5 * TMath::Power(((x[0] - par[1]) / par[2]), 2));
  double pe_1 = par[3] * exp(-0.5 * TMath::Power(((x[0] - par[4]) / par[5]), 2));
  double pe_2 = par[6] * exp(-0.5 * TMath::Power(((x[0] - par[7]) / par[8]), 2));
  //  double background = par[9] / x[0];
  return pedestal + pe_1 + pe_2;
}

double Multi_pe_3(double *x, double *par)
{
  double pedestal = par[0] * exp(-0.5 * TMath::Power(((x[0] - par[1]) / par[2]), 2));
  double pe_1 = par[3] * exp(-0.5 * TMath::Power(((x[0] - par[4]) / par[5]), 2));
  double pe_2 = par[6] * exp(-0.5 * TMath::Power(((x[0] - par[7]) / par[8]), 2));
  double pe_3 = par[9] * exp(-0.5 * TMath::Power(((x[0] - par[10]) / par[11]), 2));
  //  double background = par[9] / x[0];
  return pedestal + pe_1 + pe_2 + pe_3;
}

double Multi_pe_4(double *x, double *par)
{
  double pedestal = par[0] * exp(-0.5 * TMath::Power(((x[0] - par[1]) / par[2]), 2));
  double pe_1 = par[3] * exp(-0.5 * TMath::Power(((x[0] - par[4]) / par[5]), 2));
  double pe_2 = par[6] * exp(-0.5 * TMath::Power(((x[0] - par[7]) / par[8]), 2));
  double pe_3 = par[9] * exp(-0.5 * TMath::Power(((x[0] - par[10]) / par[11]), 2));
  double pe_4 = par[12] * exp(-0.5 * TMath::Power(((x[0] - par[13]) / par[14]), 2));
  double pe_5 = par[15] * exp(-0.5 * TMath::Power(((x[0] - par[16]) / par[17]), 2));
  //  double background = par[9] / x[0];
  return pedestal + pe_1 + pe_2 + pe_3 + pe_4 + pe_5;
}


void pe()
{
  TH1D* h1[5];
  TH1F *h2 = new TH1F("h2", "h2", 35, 0., 35.);
  string File_num[5] = {"only_pedestal", "2.4Vpp_700mV_offset", "2.4Vpp_760mV_offset", "2.4Vpp_820mV_offset", "2.4Vpp_905mV_offset"};
  double x[6], y[6];
  int i = 0, p_i = 0;
  
  
  TCanvas *c1 = new TCanvas("c1", "c1", 900, 900);
  c1->Divide(2,2);
  TCanvas *c2 = new TCanvas("c2", "c2", 900, 900);
  c2->Divide(1,1);

  TLegend *legend[4];
  
  
  auto *pe_peek_Co = new TF1("pe_peek_Co", GausM, 800., 1200., 4);
  pe_peek_Co->SetParLimits(0, 2500., 5000.);
  pe_peek_Co->SetParLimits(3, -500, 500.); 


  auto *pedestal_1_pe = new TF1("pedestal_1_pe", Single_pe, 40., 200., 6);
  pedestal_1_pe->SetParLimits(0, 10000., 15000.);
  pedestal_1_pe->SetParLimits(1, 50., 80.);
  pedestal_1_pe->SetParLimits(2, 2., 20.);

  pedestal_1_pe->SetParLimits(3, 10., 1000.);
  pedestal_1_pe->SetParLimits(4, 90., 130.);
  pedestal_1_pe->SetParLimits(5, 2., 30.);
  //  pedestal_1_pe->SetParLimits(6, -10., 10.);


  auto *pedestal_multi_pe_2 = new TF1("pedestal_multi_pe_2", Multi_pe_3, 40., 300., 9);
  pedestal_multi_pe_2->SetLineColor(3); 
  pedestal_multi_pe_2->SetParLimits(0, 2000., 15000.);
  pedestal_multi_pe_2->SetParLimits(1, 50., 80.);
  pedestal_multi_pe_2->SetParLimits(2, 2., 20.);
  
  pedestal_multi_pe_2->SetParLimits(3, 10., 1000.);
  pedestal_multi_pe_2->SetParLimits(4, 90., 120.);
  pedestal_multi_pe_2->SetParLimits(5, 2., 30.);

  pedestal_multi_pe_2->SetParLimits(6, 10., 1000.);
  pedestal_multi_pe_2->SetParLimits(7, 120., 170.);
  pedestal_multi_pe_2->SetParLimits(8, 2., 100.);

  
  auto *pedestal_multi_pe_3 = new TF1("pedestal_multi_pe_3", Multi_pe_3, 40., 500., 12);
  pedestal_multi_pe_3->SetLineColor(3);
  pedestal_multi_pe_3->SetParLimits(0, 2000., 15000.);
  pedestal_multi_pe_3->SetParLimits(1, 55., 70.);
  pedestal_multi_pe_3->SetParLimits(2, 2., 20.);
  
  pedestal_multi_pe_3->SetParLimits(3, 10., 1000.);
  pedestal_multi_pe_3->SetParLimits(4, 90., 120.);
  pedestal_multi_pe_3->SetParLimits(5, 2., 30.);

  pedestal_multi_pe_3->SetParLimits(6, 10., 1000.);
  pedestal_multi_pe_3->SetParLimits(7, 120., 185.);
  pedestal_multi_pe_3->SetParLimits(8, 2., 100.);

  pedestal_multi_pe_3->SetParLimits(9, 10., 1000.);
  pedestal_multi_pe_3->SetParLimits(10, 180., 300.);
  pedestal_multi_pe_3->SetParLimits(11, 10., 150.);


  auto *pedestal_multi_pe_4 = new TF1("pedestal_multi_pe_4", Multi_pe_4, 40., 500., 18);
  pedestal_multi_pe_4->SetLineColor(3);
  pedestal_multi_pe_4->SetParLimits(0, 300., 1200.);
  pedestal_multi_pe_4->SetParLimits(1, 55., 70.);
  pedestal_multi_pe_4->SetParLimits(2, 2., 20.);
  
  pedestal_multi_pe_4->SetParLimits(3, 300., 1000.);
  pedestal_multi_pe_4->SetParLimits(4, 90., 120.);
  pedestal_multi_pe_4->SetParLimits(5, 2., 30.);

  pedestal_multi_pe_4->SetParLimits(6, 300., 1000.);
  pedestal_multi_pe_4->SetParLimits(7, 120., 200.);
  pedestal_multi_pe_4->SetParLimits(8, 2., 100.);

  pedestal_multi_pe_4->SetParLimits(9, 100., 1000.);
  pedestal_multi_pe_4->SetParLimits(10, 180., 250.);
  pedestal_multi_pe_4->SetParLimits(11, 10., 150.);

  pedestal_multi_pe_4->SetParLimits(12, 10., 1000.);
  pedestal_multi_pe_4->SetParLimits(13, 220., 300.);
  pedestal_multi_pe_4->SetParLimits(14, 10., 150.);

  pedestal_multi_pe_4->SetParLimits(15, 10., 1000.);
  pedestal_multi_pe_4->SetParLimits(16, 260., 500.);
  pedestal_multi_pe_4->SetParLimits(17, 10., 200.);
  

  
  
  for(int j = 0 ; j < 5 ; j++)
    {

      string root_file_name = "./" + File_num[j] + ".root";
      const char *rfn = root_file_name.c_str();
      
      TFile* gfile = TFile::Open(rfn);
      TDirectory* dir = gFile->GetDirectory("Energy");

      if(j==4)
	dir->GetObject("_R_EnergyCH5@DT5730_1204", h1[j]);  // This is only validate for root file 2.4Vpp_905mV_offset
      else 
	dir->GetObject("_R_EnergyCH7@DT5730_1204", h1[j]);
      
      string plot_name = File_num[j];
      const char *p_name = plot_name.c_str();
      /*
      if(j > 0)
	{
	  c2->cd(j);
	  h1[j]->RebinX(3);
	  h1[j]->GetXaxis()->SetRangeUser(0, 500);
	  h1[j]->SetStats(0);
	  h1[j]->SetTitle(p_name);
	  h1[j]->Draw();
	}
      */  
      
      if(j > 0)
	{
	  if(j==4)
	    c2->cd(1);
	  else
	    c1->cd(j);
	  //	  h1[j]->SetTitle(p_name);
	  //	  h1[j]->RebinX(3);

	  legend[j-1] = new TLegend(0.6, 0.5, 0.88, 0.8);
	  legend[j-1]->SetBorderSize(0);
	  
	  if( j==1 )
	    {
	      h1[j]->Fit("pedestal_1_pe", "R", "", 40., 200.);
	      h1[j]->GetXaxis()->SetRangeUser(0, 500);
	      h1[j]->SetStats(0);
	      h1[j]->Draw();

	      for(int x = 1 ; x < 5 ; x = x + 3)
		{
		  double mean_peak = pedestal_1_pe->GetParameter(x);
		  legend[j-1]->AddEntry("Mean: ", Form("%.1f", mean_peak), "");
		}
	      legend[j-1]->Draw("same");
	    }
	  else if( j==2 )
	    {
	      h1[j]->Fit("pedestal_multi_pe_2", "R", "", 40., 300.);
	      h1[j]->GetXaxis()->SetRangeUser(0, 500);
	      h1[j]->SetStats(0);
	      h1[j]->Draw();

	      for(int x = 1 ; x < 8 ; x = x + 3)
		{
		  double mean_peak = pedestal_multi_pe_2->GetParameter(x);
		  legend[j-1]->AddEntry("Mean: ", Form("%.1f", mean_peak), "");
		}
	      legend[j-1]->Draw("same");
	      
	      TF1 *pe_gaus[2];
	      int i0 = 3;
	      for(int a = 0 ; a < 2 ; a++)
		{
		  pe_gaus[a] = new TF1("pe_gaus","gaus(0)", 40., 300.);
		  for(int b = 0 ; b < 3 ; b++)
		    {
		      pe_gaus[a]->SetParameter(b, (pedestal_multi_pe_2->GetParameter(i0)) );
		      i0++;
		    }
		  pe_gaus[a]->Draw("same");
		}
	    }
	  else if( j==3 )
	    {
	      h1[j]->Fit("pedestal_multi_pe_3", "R", "", 40., 500.);
	      h1[j]->GetXaxis()->SetRangeUser(0, 500);
	      h1[j]->SetStats(0);
	      h1[j]->Draw();

	      for(int x = 1 ; x < 11 ; x = x + 3)
		{
		  double mean_peak = pedestal_multi_pe_3->GetParameter(x);
		  legend[j-1]->AddEntry("Mean: ", Form("%.1f", mean_peak), "");
		}
	      legend[j-1]->Draw("same");
	      
	      TF1 *pe_gaus[3];
	      int i0 = 3;
	      for(int a = 0 ; a < 3 ; a++)
		{
		  pe_gaus[a] = new TF1("pe_gaus","gaus(0)", 40., 500.);
		  for(int b = 0 ; b < 3 ; b++)
		    {
		      pe_gaus[a]->SetParameter(b, (pedestal_multi_pe_3->GetParameter(i0)) );
		      i0++;
		    }
		  pe_gaus[a]->Draw("same");
		}
	    }
	  else if( j==4 )
	    {
	      h1[j]->Scale(1./20.);
	      h1[j]->Fit("pedestal_multi_pe_4", "R", "", 40., 500.);
	      h1[j]->GetXaxis()->SetRangeUser(0, 700);
	      h1[j]->SetStats(0);
	      h1[j]->Draw();

	      for(int x = 1 ; x < 14 ; x = x + 3)
		{
		  double mean_peak = pedestal_multi_pe_4->GetParameter(x);
		  legend[j-1]->AddEntry("Mean: ", Form("%.1f", mean_peak), "");
		}
	      legend[j-1]->Draw("same");
	      
	      TF1 *pe_gaus[5];
	      int i0 = 3;
	      for(int a = 0 ; a < 5 ; a++)
		{
		  pe_gaus[a] = new TF1("pe_gaus","gaus(0)", 40., 500.);
		  for(int b = 0 ; b < 3 ; b++)
		    {
		      pe_gaus[a]->SetParameter(b, (pedestal_multi_pe_4->GetParameter(i0)) );
		      i0++;
		    }
		  pe_gaus[a]->Draw("same");
		}
	    }   
	}
      

      /*
	double mean_value = pe_peek_Co->GetParameter(1);
	double t = std::stod(File_num[i]);
	x[i] = t;
	y[i] = mean_value;
	//      h2->Fill(t, mean_value);
	*/
    }

  /*
  TCanvas *c2 = new TCanvas("c2", "c2", 600, 600);
  TGraph* gr = new TGraph(6,x,y);
  gr->SetMarkerStyle(2);
  gr->SetMarkerSize(3);
  gr->SetTitle("PMT warm up");
  gr->GetXaxis()->SetTitle("Warming time[m]");
  gr->GetYaxis()->SetTitle("Collected charge");
  gr->Draw("ap");
  */
}
