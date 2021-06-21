#include <string>
using namespace std;

double GausM(double *x, double *par)
{
  return par[0] * exp(-0.5 * TMath::Power(((x[0] - par[1]) / par[2]), 2)) + par[3];
}



void test()
{
  TH1D* h1[56];
  TH1F *h2 = new TH1F("h2", "h2", 35, 0., 35.);
  
  //  string File_num[6] = {"5", "10", "15", "20", "25", "30"};
  
  string File_num[2] = {"before", "after"};
  int start_voltage[4] = {600, 600, 590, 610};
  
  double x[7], y[7];
  int h_count = 0, LED_count = 0;

  TLegend *legend[11];
  TGraph *LED_linearity[8];

  TCanvas *c1 = new TCanvas("c1", "c1", 1400, 700);
  c1->Divide(4,2);
  
  for(int i = 0 ; i < 4 ; i++)
    {
      for(int j = 0 ; j < 2 ; j++)
	{
	  for(int k = 0 ; k < 7 ; k++)
	    {
	      // legend[h_count] = new TLegend(0.6, 0.65, 0.8, 0.8);
	      // legend[h_count]->SetBorderSize(0);
		
	      string th = to_string(i), S_V = to_string( (start_voltage[i] + k * 5) );
	      	  
	      string root_file_name = "./LED_0" + th + "/" + File_num[j] + "/" + S_V + "_60sec.root";
	      const char *rfn = root_file_name.c_str();
      
	      TFile* gfile = TFile::Open(rfn);
	      TDirectory* dir = gFile->GetDirectory("Energy");

	      dir->GetObject("_R_EnergyCH0@DT5730_1204", h1[h_count]);
  
	      h1[h_count]->RebinX(5);


	      double h_entry = h1[h_count]->GetEntries(), h_mean = h1[h_count]->GetMean();
	      auto *pe_peek_Co = new TF1("pe_peek_Co", GausM, (h_mean-150.), (h_mean+150.), 4);
 
	      pe_peek_Co->SetParLimits(0, 50., h_entry);
	      pe_peek_Co->SetParLimits(1, (h_mean-250.), (h_mean+250.));
	      pe_peek_Co->SetParLimits(2, 50., 1000.);
	      pe_peek_Co->SetParLimits(3, -2000, 1000.); 
	      h1[h_count]->Fit("pe_peek_Co", "R", "", (h_mean-150.), (h_mean+150.));

	      if( (i == 0) && (j == 0) )
		{
		  c1->cd(k+1);
		  h1[h_count]->Draw();
		  //		  pe_peek_Co->Draw("same");		  
		}

	      
	      double mean_value = pe_peek_Co->GetParameter(1);
	      double sigma = pe_peek_Co->GetParameter(2);
	      y[k] = sigma;
	      x[k] = (start_voltage[i] + k * 5);
	      /*
	      h1[h_count]->SetStats(0);
	      h1[h_count]->GetXaxis()->SetRangeUser(30, 1000);
	      h1[h_count]->SetTitle(Form("LED_%d  %dth plug", i, (j+1)));
	      h1[h_count]->Draw();
	      	      
	      double mean_value = pe_peek_Co->GetParameter(1);	    
	      legend[h_count]->AddEntry("", Form("peak: %.1f", mean_value), "");
	      legend[h_count]->Draw("same");
	      */
	      
	      h_count++;	      
	    }
	  
	  LED_linearity[LED_count] = new TGraph(7, x, y);
	  LED_linearity[LED_count]->SetMarkerStyle(21);
	  LED_linearity[LED_count]->SetMarkerSize(1);
	  LED_count++;
	}
    }


  /*
  TCanvas *c2 = new TCanvas("c2", "c2", 1400, 700);
  c2->Divide(4,2);

  for( int i = 0 ; i < 8 ; i++ )
    {
      c2->cd(i+1);
      LED_linearity[i]->Draw("AP");
      
    }
  */
  
}
