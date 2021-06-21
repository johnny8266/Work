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



void read_waveform()
{

  const int N_source = 1; 
  //  string File_num[4] = {"LED_signal_pedestal_to_2pe.root", "20sec_crystal_040668_no_source.root", "3sec_crystal_040668_Na22.root", "3sec_crystal_040668_Cs137.root"};
  
  string File_num[1] = {"0.5Vpp_TR_300lsb_Co60_1secs.root"};
  //  string File_num[1] = {"2Vpp_TR_75lsb_Co60_1secs.root"};

  auto *c1 = new TCanvas("c1", "c1", 1350, 900);
  c1->Divide(3,2);
  
  for(int j = 0 ; j < N_source ; j++)
    {
      string root_file_name  = "./" + File_num[j];
      const char *rfn = root_file_name.c_str();
      TFile *rfile = TFile::Open(rfn);
      TTree *Data_R = (TTree*)rfile->Get("Data_R");
      //      double fArray[];
      //      vector<double> fArray;
      auto *chain = new TChain ("Data_R");
      chain->Add(rfn);
      
      UShort_t Channel;
      UShort_t Energy;

      TArrayS mWaveform;
      TArrayS *mWaveformPnt = &mWaveform;
      
      // Data_R->Print();
      // Data_R->SetBranchAddress("Channel", &Channel);
      // Data_R->SetBranchAddress("Energy", &Energy);
      chain->SetBranchAddress("Samples", &mWaveformPnt);
      chain->SetBranchAddress("Energy", &Energy);
      //      Samples->Browse();
      //      Samples->FindBranch("fArray");
      //      Data_R->SetBranchAddress("fArray", &fArray);

      int events_numer = (int)Data_R->GetEntries();
      cout << events_numer << endl;
      int pe_flag = 0, ped_flag = 0;
      TH1I *pe_h[6], *ped_h[3];
      
      for(int i = 0; i < events_numer; i = i + 20)
	{
	  chain->GetEvent(i);
	  //	  cout << Energy << endl;

	  if( (Energy > 1000) && (Energy < 1200) && (pe_flag < 6) )
	    {
	      c1->cd(1+pe_flag);
	      int a = mWaveformPnt->GetSize();
	      pe_h[pe_flag] = new TH1I(Form("Cs_137_%d", pe_flag), Form("Cs_137_%d", pe_flag), 512, 0, 512);
	      
	      for(int j = 0 ; j < a ; j++)
		{
		  int ADC = mWaveformPnt->GetAt(j);
		  pe_h[pe_flag]->SetBinContent(j, ADC);
		}
	      pe_h[pe_flag]->SetStats(0);
	      pe_h[pe_flag]->SetMaximum(13500);
	      pe_h[pe_flag]->SetMinimum(11000);
	      pe_h[pe_flag]->GetXaxis()->SetRangeUser(0, 200);
	      pe_h[pe_flag]->GetXaxis()->SetTitle("N Sample");
	      pe_h[pe_flag]->GetYaxis()->SetTitle("ADC");
	      pe_h[pe_flag]->Draw();
	      pe_flag++;
	      cout << " ??? " << endl;
	    }
	  /*
	  if( (Energy > 600) && (Energy < 700) && (ped_flag < 3) )
	    {
	      c1->cd(4+ped_flag);
	      int a = mWaveformPnt->GetSize();
	      ped_h[ped_flag] = new TH1I(Form("No_source_sg_%d", ped_flag), Form("No_source_sg_%d", ped_flag), 512, 0, 512);	      
	      
	      for(int j = 0 ; j < a ; j++)
		{
		  int ADC = mWaveformPnt->GetAt(j);
		  ped_h[ped_flag]->SetBinContent(j, ADC);
		}
	      ped_h[ped_flag]->SetStats(0);
	      ped_h[ped_flag]->SetMaximum(14000);
	      ped_h[ped_flag]->SetMinimum(11000);
	      ped_h[ped_flag]->GetXaxis()->SetRangeUser(0, 200);
	      ped_h[ped_flag]->GetXaxis()->SetTitle("N Sample");
	      ped_h[ped_flag]->GetYaxis()->SetTitle("ADC");
	      ped_h[ped_flag]->Draw();
	      ped_flag++;	      
	    }
	  */
	  //	  if( (pe_flag == 3) && (ped_flag == 3) )
	  if( pe_flag == 6 )
	    break;
	   
	}
      
    }// End read file


  
  
}
