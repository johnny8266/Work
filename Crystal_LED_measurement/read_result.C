#include "TLorentzVector.h"
#include "TMath.h" 
#include "TRandom.h"
#include <vector>
#include <ctime>
#include <cmath>
#include <iostream>
#include <algorithm>
using namespace std;


double GausM(double *x, double *par)
{
  return par[0] * exp(-0.5 * TMath::Power(((x[0] - par[1]) / par[2]), 2)) + par[3];
}



void read_result()
{

  //  TH1I *ped[3], *Co_ped[3], *Co_pure[3], *Cs_ped[3], *Cs_pure[3];

  const int N_source = 7; 
  //  string File_num[N_source] = {"compassF_run-nosource.root", "compassF_run-Co60.root", "compassF_run-Cs137.root"};
  string File_num[N_source] = {"compassF_run_Am241_5mins.root", "compassF_run_Co60_5mins.root", "compassF_run_Fe55_5mins.root", "compassF_run_Ba133_5mins.root", "compassF_run_Cs137_5mins.root", "compassF_run_Na22_5mins.root", "compassF_run_pedestal_5mins.root"};

  TH1I *ped = new TH1I("ped", "ped", 200, 0, 200);
  TH1I *Co_ped = new TH1I("Co_ped", "Co_ped", 200, 0, 200);
  TH1I *Co_pure = new TH1I("Co_pure", "Co_pure", 200, 0, 200);
  TH1I *Cs_ped = new TH1I("Cs_ped", "Cs_ped", 200, 0, 200);
  TH1I *Cs_pure = new TH1I("Cs_pure", "Cs_pure", 200, 0, 200);
  TH1I *Am_ped = new TH1I("Am_ped", "Am_ped", 200, 0, 200);
  TH1I *Am_pure = new TH1I("Am_pure", "Am_pure", 200, 0, 200);
  TH1I *Fe_ped = new TH1I("Fe_ped", "Fe_ped", 200, 0, 200);
  TH1I *Fe_pure = new TH1I("Fe_pure", "Fe_pure", 200, 0, 200);
  TH1I *Ba_ped = new TH1I("Ba_ped", "Ba_ped", 200, 0, 200);
  TH1I *Ba_pure = new TH1I("Ba_pure", "Ba_pure", 200, 0, 200);
  TH1I *Na_ped = new TH1I("Na_ped", "Na_ped", 200, 0, 200);
  TH1I *Na_pure = new TH1I("Na_pure", "Na_pure", 200, 0, 200);

  
  /*
  for(int i = 0 ; i < 3 ; i++)
    {
      ped[i] = new TH1I(Form("ped_%d", i), Form("ped_%d", i), 500, 0, 500);
      Co_ped[i] = new TH1I(Form("Co_ped_%d", i), Form("Co_ped_%d", i), 500, 0, 500);
      Co_pure[i] = new TH1I(Form("Co_pure_%d", i), Form("Co_pure_%d", i), 500, 0, 500);
      Cs_ped[i] = new TH1I(Form("Cs_ped_%d", i), Form("Cs_ped_%d", i), 500, 0, 500);
      Cs_pure[i] = new TH1I(Form("Cs_pure_%d", i), Form("Cs_pure_%d", i), 500, 0, 500);
    }
  */

  for(int j = 0 ; j < N_source ; j++)
    {

      //      string root_file_name = "./22-03-2021/second/" + File_num[j];
      string root_file_name = "./23-03-2021/source/" + File_num[j];
      const char *rfn = root_file_name.c_str();

      TFile *rfile = TFile::Open(rfn);
      TTree *Data_F = (TTree*)rfile->Get("Data_F");

      TTreeReader fReader("Data_F", rfile);

      TTreeReaderArray<UShort_t>  Board = {fReader, "Board"};
      TTreeReaderArray<UShort_t>  Channel = {fReader, "Channel"};
      TTreeReaderArray<ULong64_t> Timestamp = {fReader, "Timestamp"};
      TTreeReaderArray<UShort_t>  Energy = {fReader, "Energy"};
      TTreeReaderArray<UShort_t>  EnergyShort = {fReader, "EnergyShort"};
      TTreeReaderArray<UInt_t>    Flags = {fReader, "Flags"};

      int events_numer = 0;
  
      while (fReader.Next())
	{
	  //	  if(++events_numer > 5)
	  //	    break;
	  ++events_numer;
	  

	  if(j == 0)
	    Am_ped->Fill(Energy[0]);

	  else if(j == 1)
      	    Co_ped->Fill(Energy[0]);

	  else if(j == 2)
	    Fe_ped->Fill(Energy[0]);

	  else if(j == 3)
	    Ba_ped->Fill(Energy[0]);

	  else if(j == 4)
      	    Cs_ped->Fill(Energy[0]);

	  else if(j == 5)
	    Na_ped->Fill(Energy[0]);

	  else if(j == 6)
	    ped->Fill(Energy[0]);
	  //      cout << Board[0] << " " << Channel[0] << " " << Timestamp[0] << " " << Energy[0] << " " << EnergyShort[0] << " " << Flags[0] << endl;
      
	}

      cout << events_numer << endl;
      
      //      ped->Draw();
    }

  *Co_pure = *Co_ped - *ped;
  *Cs_pure = *Cs_ped - *ped;
  *Am_pure = *Am_ped - *ped;
  *Fe_pure = *Fe_ped - *ped;
  *Na_pure = *Na_ped - *ped;
  *Ba_pure = *Ba_ped - *ped;
  

  TLegend *legend_1;

  auto c3 = new TCanvas("c3", "c3", 1350, 900);
  c3->Divide(3,2);
  c3->cd(1);
  Co_pure->Draw();
  c3->cd(2);
  Cs_pure->Draw();
  c3->cd(3);
  Am_pure->Draw();
  c3->cd(4);
  Fe_pure->Draw();
  c3->cd(5);
  Na_pure->Draw();
  c3->cd(6);
  Ba_pure->Draw();

  /*
  auto c1 = new TCanvas("c1", "c1", 900, 900);
  c1->Divide(2,2);
  c1->cd(1);
  Co_ped->SetStats(0);
  Co_ped->Draw();
  c1->cd(2);
  Cs_ped->SetStats(0);
  Cs_ped->Draw();

  
  auto *pe_peek_Co = new TF1("pe_peek_Co", GausM, 40., 95., 4);
  pe_peek_Co->SetParLimits(0, 5., 200000.);
  pe_peek_Co->SetParLimits(1, 50., 80.);
  pe_peek_Co->SetParLimits(2, 15., 25.);
  pe_peek_Co->SetParLimits(3, 1., 1000.);
  c1->cd(3);
  Co_pure->Fit("pe_peek_Co", "R", "", 40., 95.);
  Co_pure->SetStats(0);
  Co_pure->GetXaxis()->SetTitle("ADC");
  Co_pure->SetTitle("Pure Co signal");
  Co_pure->Draw();
  legend_1 = new TLegend(0.5, 0.6, 0.8, 0.8);
  legend_1->SetBorderSize(0);
  legend_1->AddEntry((TObject*)0, Form("ADC peak: %.3f", (pe_peek_Co->GetParameter(1)) ), "");
  legend_1->AddEntry((TObject*)0, Form("Sigma: %.3f", (pe_peek_Co->GetParameter(2)) ), "");
  legend_1->AddEntry((TObject*)0, Form("Chi / NDF: %.3f", (pe_peek_Co->GetChisquare() / pe_peek_Co->GetNDF()) ), "");
  legend_1->Draw("same");


  auto *pe_peek_Cs = new TF1("pe_peek_Cs", GausM, 30., 65., 4);
  pe_peek_Cs->SetParLimits(0, 5., 200000.);
  pe_peek_Cs->SetParLimits(1, 30., 55.);
  pe_peek_Cs->SetParLimits(2, 10., 25.);
  pe_peek_Cs->SetParLimits(3, 1., 10000.);
  c1->cd(4);
  Cs_pure->Fit("pe_peek_Cs", "R", "", 30., 65.);
  Cs_pure->SetTitle("Pure Cs signal");
  Cs_pure->SetStats(0);
  Cs_pure->GetXaxis()->SetTitle("ADC");
  Cs_pure->Draw();
  legend_1 = new TLegend(0.5, 0.6, 0.8, 0.8);
  legend_1->SetBorderSize(0);
  legend_1->AddEntry((TObject*)0, Form("ADC peak: %.3f", (pe_peek_Cs->GetParameter(1)) ), "");
  legend_1->AddEntry((TObject*)0, Form("Sigma: %.3f", (pe_peek_Cs->GetParameter(2)) ), "");
  legend_1->AddEntry((TObject*)0, Form("Chi / NDF: %.3f", (pe_peek_Cs->GetChisquare() / pe_peek_Co->GetNDF()) ), "");
  legend_1->Draw("same");


  auto c2 = new TCanvas("c2", "c2", 450, 450);
  auto *pedstal = new TF1("pedstal", GausM, 20., 35., 4);
  pedstal->SetParLimits(0, 1000., 20000.);
  pedstal->SetParLimits(1, 20., 35.);
  pedstal->SetParLimits(2, 3., 10.);
  pedstal->SetParLimits(3, 1., 10000.);
  ped->Fit("pedstal", "R", "", 20., 35.);
  ped->SetStats(0);
  ped->GetXaxis()->SetTitle("ADC");
  ped->Draw();
  legend_1 = new TLegend(0.5, 0.6, 0.8, 0.8);
  legend_1->SetBorderSize(0);
  legend_1->AddEntry((TObject*)0, Form("ADC peak: %.3f", (pedstal->GetParameter(1)) ), "");
  legend_1->AddEntry((TObject*)0, Form("Sigma: %.3f", (pedstal->GetParameter(2)) ), "");
  legend_1->AddEntry((TObject*)0, Form("Chi / NDF: %.3f", (pedstal->GetChisquare() / pe_peek_Co->GetNDF()) ), "");
  legend_1->Draw("same");
  */
  
}
