#include <string>
using namespace std;

double GausM(double *x, double *par)
{
  return par[0] * exp(-0.5 * TMath::Power(((x[0] - par[1]) / par[2]), 2)) + par[3];
}



void test()
{
  TH1D* h1[8];
  TH1F *h2 = new TH1F("h2", "h2", 35, 0., 35.);
  //  string File_num[6] = {"5", "10", "15", "20", "25", "30"};
  
  string File_num[2] = {"plug_0th_after_warmup.root", "plug_1th_after_warmup.root"};
  
  double x[6], y[6];
  int count = 0;
  
  TCanvas *c1 = new TCanvas("c1", "c1", 1400, 700);
  c1->Divide(4,2);

  TLegend *legend[11];
  
  for(int i = 0 ; i < 4 ; i++)
    {
      for(int j = 0 ; j < 2 ; j++)
	{
	  c1->cd(i + 1 + j * 4);
	  legend[count] = new TLegend(0.6, 0.65, 0.8, 0.8);
	  legend[count]->SetBorderSize(0);
	  string th = to_string(i);
	  
	  string root_file_name = "./LED_0" + th + "/" + File_num[j];
	  const char *rfn = root_file_name.c_str();
      
	  TFile* gfile = TFile::Open(rfn);
	  TDirectory* dir = gFile->GetDirectory("Energy");

	  dir->GetObject("_R_EnergyCH0@DT5730_1204", h1[count]);
  
	  h1[count]->RebinX(5);

	  double h_entry = h1[count]->GetEntries(), h_mean = h1[count]->GetMean();
	  auto *pe_peek_Co = new TF1("pe_peek_Co", GausM, (h_mean-150.), (h_mean+150.), 4);
 
	  pe_peek_Co->SetParLimits(0, 500., 5000.);
	  pe_peek_Co->SetParLimits(1, (h_mean-50.), (h_mean+50.));
	  pe_peek_Co->SetParLimits(2, 50., 300.);
	  pe_peek_Co->SetParLimits(3, -2000, 1000.); 
	  h1[count]->Fit("pe_peek_Co", "R", "", (h_mean-150.), (h_mean+150.));

	  h1[count]->SetStats(0);
	  h1[count]->GetXaxis()->SetRangeUser(30, 1000);
	  h1[count]->SetTitle(Form("LED_%d  %dth plug", i, (j+1)));
	  h1[count]->Draw();

	  double mean_value = pe_peek_Co->GetParameter(1);	    
	  legend[count]->AddEntry("", Form("peak: %.1f", mean_value), "");
	  legend[count]->Draw("same");

	  count++;
	}
    }

      /*
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
      */
      

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
