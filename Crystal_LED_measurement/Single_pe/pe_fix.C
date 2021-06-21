#include <iostream>
#include <string>
#include <stdlib.h>


// The fitting functions in this code
// are Gaussian, but their position are not free.
// They interval between each Gaussian are fix!

double GausM(double *x, double *par)
{
  return par[0] * exp(-0.5 * TMath::Power(((x[0] - par[1]) / par[2]), 2)) + par[3];
}


double Single_pe(double *x, double *par)
{
  double pedestal = par[0] * exp(-0.5 * TMath::Power(((x[0] - par[1]) / par[2]), 2));
  double pe_1 = par[3] * exp(-0.5 * TMath::Power(((x[0] - par[4]) / par[5]), 2));
  return pedestal + pe_1;
}

double Multi_pe_2(double *x, double *par)
{
  // double pedestal = par[0] * exp(-0.5 * TMath::Power(((x[0] - par[1]) / par[2]), 2));
  // double pe_1 = par[3] * exp(-0.5 * TMath::Power(((x[0] - par[4]) / par[5]), 2));
  // double pe_2 = par[6] * exp(-0.5 * TMath::Power(((x[0] - par[7]) / par[8]), 2));

  double pedestal = par[0] * exp(-0.5 * TMath::Power(((x[0] - par[1]) / par[2]), 2));
  double pe_1 = par[3] * exp(-0.5 * TMath::Power(((x[0] - (par[1] + par[4]) ) / par[5]), 2));
  double pe_2 = par[6] * exp(-0.5 * TMath::Power(((x[0] - (par[1] + 2. * par[4]) ) / par[7]), 2));
  return pedestal + pe_1 + pe_2;
}

double Multi_pe_3(double *x, double *par)
{
  // double pedestal = par[0] * exp(-0.5 * TMath::Power(((x[0] - par[1]) / par[2]), 2));
  // double pe_1 = par[3] * exp(-0.5 * TMath::Power(((x[0] - par[4]) / par[5]), 2));
  // double pe_2 = par[6] * exp(-0.5 * TMath::Power(((x[0] - par[7]) / par[8]), 2));
  // double pe_3 = par[9] * exp(-0.5 * TMath::Power(((x[0] - par[10]) / par[11]), 2));

  double pedestal = par[0] * exp(-0.5 * TMath::Power(((x[0] - par[1]) / par[2]), 2));
  double pe_1 = par[3] * exp(-0.5 * TMath::Power(((x[0] - (par[1] + par[4]) ) / par[5]), 2));
  double pe_2 = par[6] * exp(-0.5 * TMath::Power(((x[0] - (par[1] + 2. * par[4]) ) / par[7]), 2));
  double pe_3 = par[8] * exp(-0.5 * TMath::Power(((x[0] - (par[1] + 3. * par[4]) ) / par[9]), 2));
  return pedestal + pe_1 + pe_2 + pe_3;
}

double Multi_pe_4(double *x, double *par)
{
  // double pedestal = par[0] * exp(-0.5 * TMath::Power(((x[0] - par[1]) / par[2]), 2));
  // double pe_1 = par[3] * exp(-0.5 * TMath::Power(((x[0] - par[4]) / par[5]), 2));
  // double pe_2 = par[6] * exp(-0.5 * TMath::Power(((x[0] - par[7]) / par[8]), 2));
  // double pe_3 = par[9] * exp(-0.5 * TMath::Power(((x[0] - par[10]) / par[11]), 2));
  // double pe_4 = par[12] * exp(-0.5 * TMath::Power(((x[0] - par[13]) / par[14]), 2));
  // double pe_5 = par[15] * exp(-0.5 * TMath::Power(((x[0] - par[16]) / par[17]), 2));

  double pedestal = par[0] * exp(-0.5 * TMath::Power(((x[0] - par[1]) / par[2]), 2));
  double pe_1 = par[3] * exp(-0.5 * TMath::Power(((x[0] - (par[1] + par[4]) ) / par[5]), 2));
  double pe_2 = par[6] * exp(-0.5 * TMath::Power(((x[0] - (par[1] + 2. * par[4]) ) / par[7]), 2));
  double pe_3 = par[8] * exp(-0.5 * TMath::Power(((x[0] - (par[1] + 3. * par[4]) ) / par[9]), 2));
  double pe_4 = par[10] * exp(-0.5 * TMath::Power(((x[0] - (par[1] + 4. * par[4]) ) / par[11]), 2));
  double pe_5 = par[12] * exp(-0.5 * TMath::Power(((x[0] - (par[1] + 5. * par[4]) ) / par[13]), 2));
  //  double pe_6 = par[14] * exp(-0.5 * TMath::Power(((x[0] - (par[1] + 6. * par[4]) ) / par[15]), 2));
  //  double background = par[9] / x[0];
  return pedestal + pe_1 + pe_2 + pe_3 + pe_4 + pe_5;
}


//.......ooooo0000000ooooo..............ooooo0000000ooooo..............ooooo0000000ooooo..............ooooo0000000ooooo.......


void pe_fix()
{
  const int n_files = 3;
  TH1D* h1[n_files];
  TH1F *h2 = new TH1F("h2", "h2", 35, 0., 35.);
  string File_num[n_files] = {"10_degree", "30_degree", "50_degree"};
  double x[6], y[6];
  int i = 0, p_i = 0;
    
  TCanvas *c1 = new TCanvas("c1", "c1", 1800, 600);
  c1->Divide(3,1);

  TLegend *legend[4];
  
  
  auto *pedestal_multi_pe_2 = new TF1("pedestal_multi_pe_2", Multi_pe_3, 40., 300., 8);
  pedestal_multi_pe_2->SetLineColor(3); 
  pedestal_multi_pe_2->SetParLimits(0, 2000., 15000.);
  pedestal_multi_pe_2->SetParLimits(1, 50., 80.);
  pedestal_multi_pe_2->SetParLimits(2, 2., 20.);
  
  pedestal_multi_pe_2->SetParLimits(3, 10., 1000.);
  pedestal_multi_pe_2->SetParLimits(4, 40., 60.);
  pedestal_multi_pe_2->SetParLimits(5, 2., 30.);

  pedestal_multi_pe_2->SetParLimits(6, 10., 1000.);
  pedestal_multi_pe_2->SetParLimits(7, 10., 100.);

  
  auto *pedestal_multi_pe_3 = new TF1("pedestal_multi_pe_3", Multi_pe_3, 40., 500., 10);
  pedestal_multi_pe_3->SetLineColor(3);
  pedestal_multi_pe_3->SetParLimits(0, 600., 1000.);
  pedestal_multi_pe_3->SetParLimits(1, 55., 70.);
  pedestal_multi_pe_3->SetParLimits(2, 2., 20.);
  
  pedestal_multi_pe_3->SetParLimits(3, 10., 1000.);
  pedestal_multi_pe_3->SetParLimits(4, 40., 60.);
  pedestal_multi_pe_3->SetParLimits(5, 2., 50.);

  pedestal_multi_pe_3->SetParLimits(6, 10., 1000.);
  pedestal_multi_pe_3->SetParLimits(7, 10., 100.);

  pedestal_multi_pe_3->SetParLimits(8, 10., 1000.);
  pedestal_multi_pe_3->SetParLimits(9, 10., 100.);


  auto *pedestal_multi_pe_4 = new TF1("pedestal_multi_pe_4", Multi_pe_4, 40., 600., 14);
  pedestal_multi_pe_4->SetLineColor(3);
  pedestal_multi_pe_4->SetParLimits(0, 5000., 10000.);
  pedestal_multi_pe_4->SetParLimits(1, 250., 270.);
  pedestal_multi_pe_4->SetParLimits(2, 2., 20.);
  
  pedestal_multi_pe_4->SetParLimits(3, 1000., 10000.);
  pedestal_multi_pe_4->SetParLimits(4, 40., 60.);
  pedestal_multi_pe_4->SetParLimits(5, 10., 100.);

  pedestal_multi_pe_4->SetParLimits(6, 500., 7000.);
  pedestal_multi_pe_4->SetParLimits(7, 10., 200.);

  pedestal_multi_pe_4->SetParLimits(8, 300., 7000.);
  pedestal_multi_pe_4->SetParLimits(9, 10., 200.);

  pedestal_multi_pe_4->SetParLimits(10, 100., 5000.);
  pedestal_multi_pe_4->SetParLimits(11, 10., 200.);

  pedestal_multi_pe_4->SetParLimits(12, 100., 3000.);
  pedestal_multi_pe_4->SetParLimits(13, 10., 200.);

  

  
  
  for(int j = 0 ; j < n_files ; j++)
    {
      c1->cd(j+1);
      
      string root_file_name = "./temperature_dependance/" + File_num[j] + ".root";
      const char *rfn = root_file_name.c_str();
      
      TFile* gfile = TFile::Open(rfn);
      TDirectory* dir = gFile->GetDirectory("Energy");
      
      dir->GetObject("EnergyCH0@DT5730_1204", h1[j]);  // This is only validate for root file
      string plot_name = File_num[j];
      const char *p_name = plot_name.c_str();
      

      h1[j]->RebinX(3);
      h1[j]->Fit("pedestal_multi_pe_4", "R", "", 220., 600.);
      //      gstyle->SetOptFit(1);
      h1[j]->GetXaxis()->SetRangeUser(200, 700);
      h1[j]->SetStats(0);
      h1[j]->Draw();

      cout << endl << endl << plot_name << " interval between each pe peaks: " << pedestal_multi_pe_4->GetParameter(4);
      cout << endl << pedestal_multi_pe_4->GetChisquare() << endl;
      cout << pedestal_multi_pe_4->GetNDF();
      cout << endl << endl;
      
      TF1 *pe_gaus[5];
      int i0 = 3;

      /*
      for(int a = 0 ; a < 5 ; a++)
	{
	  pe_gaus[a] = new TF1("pe_gaus","gaus(0)", 40., 600.);
	  for(int b = 0 ; b < 3 ; b++)
	    {
	      if(b==0)
		{
		  if(a==0)
		    pe_gaus[a]->SetParameter(b, (pedestal_multi_pe_4->GetParameter(3)) );
		  else
		    {
		      int n_amp = 4 + a * 2;
		      pe_gaus[a]->SetParameter(b, (pedestal_multi_pe_4->GetParameter(n_amp)) );
		    }
		}
	      else if(b==1)
		{
		  double peinterval = pedestal_multi_pe_4->GetParameter(4);
		  pe_gaus[a]->SetParameter(b, (pedestal_multi_pe_4->GetParameter(1) + (a+1)*peinterval) );			  
		}
	      else if(b==2)
		{
		  int n_width = 5 + a * 2;
		  pe_gaus[a]->SetParameter(b, (pedestal_multi_pe_4->GetParameter(n_width)) );
		}
	    }
	  
	  pe_gaus[a]->Draw("same");
	}
      */
      
    }

      /*
      if(j >= 0)
	{
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

	      for(int x = 0 ; x < 3 ; x++)
		{
		  double mean_peak = pedestal_multi_pe_2->GetParameter(1) + (pedestal_multi_pe_2->GetParameter(4)) * x;
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
		      if(b==0)
			{
			  if(a==0)
			    pe_gaus[a]->SetParameter(b, (pedestal_multi_pe_2->GetParameter(3)) );
			  else
			    {
			      int n_amp = 4 + a * 2;
			      pe_gaus[a]->SetParameter(b, (pedestal_multi_pe_2->GetParameter(n_amp)) );
			    }
			}
		      else if(b==1)
			{
			  double peinterval = pedestal_multi_pe_2->GetParameter(4);
			  pe_gaus[a]->SetParameter(b, (pedestal_multi_pe_2->GetParameter(1) + (a+1)*peinterval) );			  
			}
		      else if(b==2)
			{
			  int n_width = 5 + a * 2;
			  pe_gaus[a]->SetParameter(b, (pedestal_multi_pe_2->GetParameter(n_width)) );
			}
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

	      for(int x = 0 ; x < 4 ; x++)
		{
		  double mean_peak = pedestal_multi_pe_3->GetParameter(1) + (pedestal_multi_pe_3->GetParameter(4)) * x;
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
		      if(b==0)
			{
			  if(a==0)
			    pe_gaus[a]->SetParameter(b, (pedestal_multi_pe_3->GetParameter(3)) );
			  else
			    {
			      int n_amp = 4 + a * 2;
			      pe_gaus[a]->SetParameter(b, (pedestal_multi_pe_3->GetParameter(n_amp)) );
			    }
			}
		      else if(b==1)
			{
			  double peinterval = pedestal_multi_pe_3->GetParameter(4);
			  pe_gaus[a]->SetParameter(b, (pedestal_multi_pe_3->GetParameter(1) + (a+1)*peinterval) );			  
			}
		      else if(b==2)
			{
			  int n_width = 5 + a * 2;
			  pe_gaus[a]->SetParameter(b, (pedestal_multi_pe_3->GetParameter(n_width)) );
			}
		    }
		  pe_gaus[a]->Draw("same");
		}
	    }
	  else if( j==4 )
	    {
	      h1[j]->Scale(1./20.);
	      h1[j]->Fit("pedestal_multi_pe_4", "R", "", 40., 600.);
	      h1[j]->GetXaxis()->SetRangeUser(0, 700);
	      h1[j]->SetStats(0);
	      h1[j]->Draw();

	      for(int x = 0 ; x < 6 ; x++)
		{
		  double mean_peak = pedestal_multi_pe_4->GetParameter(1) + (pedestal_multi_pe_4->GetParameter(4)) * x;
		  legend[j-1]->AddEntry("Mean: ", Form("%.1f", mean_peak), "");
		}
	      legend[j-1]->Draw("same");
	      
	      TF1 *pe_gaus[5];
	      int i0 = 3;
	      
	      for(int a = 0 ; a < 5 ; a++)
		{
		  pe_gaus[a] = new TF1("pe_gaus","gaus(0)", 40., 600.);
		  for(int b = 0 ; b < 3 ; b++)
		    {
		      if(b==0)
			{
			  if(a==0)
			    pe_gaus[a]->SetParameter(b, (pedestal_multi_pe_4->GetParameter(3)) );
			  else
			    {
			      int n_amp = 4 + a * 2;
			      pe_gaus[a]->SetParameter(b, (pedestal_multi_pe_4->GetParameter(n_amp)) );
			    }
			}
		      else if(b==1)
			{
			  double peinterval = pedestal_multi_pe_4->GetParameter(4);
			  pe_gaus[a]->SetParameter(b, (pedestal_multi_pe_4->GetParameter(1) + (a+1)*peinterval) );			  
			}
		      else if(b==2)
			{
			  int n_width = 5 + a * 2;
			  pe_gaus[a]->SetParameter(b, (pedestal_multi_pe_4->GetParameter(n_width)) );
			}

		      
		      //		      i0++;
		    }
		  pe_gaus[a]->Draw("same");
		}
	      
	    }   
	}
      */

     




  
  /*  
  TCanvas *c3 = new TCanvas("c3", "c3", 900, 900);
  TH1D* h1_clear;
  TFile* lfile = TFile::Open("./pe_distribution_clear/Really_clear_pe_5minutes_per_intensity_930_1018mV_8mv_interval.root");
  TDirectory* dir = lfile->GetDirectory("Energy");
  
  dir->GetObject("_R_EnergyCH0@DT5730_1204", h1_clear);  // This is only validate for root file 2.4Vpp_905mV_offset
  TLine *add_line_1 = new TLine(65, 0., 65, 81000.);
  TLine *add_line_2 = new TLine(122, 0., 122, 81000.);
  TLine *add_line_3 = new TLine(179, 0., 179, 81000.);

  h1_clear->SetStats(0);
  h1_clear->GetYaxis()->SetTitle("");
  h1_clear->GetYaxis()->SetRange(0., 80000.);
  h1_clear->RebinX(3);
  h1_clear->GetXaxis()->SetRangeUser(0, 500);
  h1_clear->Draw();
  add_line_1->SetLineWidth(2);
  add_line_1->SetLineColor(2);
  add_line_1->Draw("same");
  add_line_2->SetLineWidth(2);
  add_line_2->SetLineColor(2);
  add_line_2->Draw("same");
  add_line_3->SetLineWidth(2);
  add_line_3->SetLineColor(2);
  add_line_3->Draw("same");
  */
}
