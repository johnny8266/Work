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


void pe_temp()
{
  const int n_files = 3;
  TH1D* h1[n_files];
  string File_num[n_files] = {"10_degree", "30_degree", "50_degree"};
  int i = 0;
    
  TCanvas *c1 = new TCanvas("c1", "c1", 1800, 600);
  c1->Divide(3,1);

  

  /*
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
  */


  auto *pedestal_multi_pe_4 = new TF1("pedestal_multi_pe_4", Multi_pe_4, 40., 600., 14);
  pedestal_multi_pe_4->SetLineColor(3);
  pedestal_multi_pe_4->SetParLimits(0, 5000., 10000.);
  pedestal_multi_pe_4->SetParLimits(1, 250., 270.);
  pedestal_multi_pe_4->SetParLimits(2, 2., 20.);
  
  pedestal_multi_pe_4->SetParLimits(3, 1000., 10000.);
  pedestal_multi_pe_4->SetParLimits(4, 40., 60.);  // interval between each pe peak
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
      h1[j]->GetXaxis()->SetRangeUser(200, 700);
      h1[j]->SetStats(0);
      h1[j]->SetTitle(p_name);
      h1[j]->Draw();

      TLegend *legend_temp;
      legend_temp = new TLegend(0.55, 0.6, 0.85, 0.85);
      legend_temp->SetBorderSize(0);
      legend_temp->AddEntry((TObject*)0, Form(" #chi^{2}:  %.3f", (pedestal_multi_pe_4->GetChisquare() / pedestal_multi_pe_4->GetNDF()) ), "");
      legend_temp->AddEntry((TObject*)0, Form(" pe interval:  %.3f", pedestal_multi_pe_4->GetParameter(4) ), "");
      legend_temp->Draw("same");
      
      cout << endl << endl << plot_name << " interval between each pe peaks: " << pedestal_multi_pe_4->GetParameter(4) << endl;


      
      TF1 *pe_gaus[5];
      int i0 = 3;

      // Draw the each pe peak except pedestal
      for(int a = 0 ; a < 5 ; a++)
	{
	  pe_gaus[a] = new TF1("pe_gaus","gaus(0)", 220., 600.);
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
            
    }


}
