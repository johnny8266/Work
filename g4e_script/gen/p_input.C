#include <iostream>
#include <fstream>
#include <ctime>
#include "TMath.h" 
#include "TRandom.h"

using namespace std;

void p_input()
{

  Int_t N_events = 0, loop = 0;

  auto gE_distribution = new TH1F("gE_distribution", "gE_distribution", 20, 100., 102.);

  auto pos_map = new TH2F("pos_map", "pos_map", 20, 2.3, 2.7, 30, -0.15, 0.15);
  
  TRandom *R = new TRandom();
  R->SetSeed(0);
  ofstream myfile;
  myfile.open ("../../Data/g4e_input/p_input.txt");

  double gx = 0., gy = 0., gz = 0., gE = 0.;
  vector<double> g_x, g_y, g_z, g_E;

  
  while(1)
    {      
      int empty_flag = 0;
      g_x.clear();  g_y.clear();  g_z.clear();  g_E.clear();  

      myfile << "3 " << 0.1 << " " << 0.01 << " " << -0.1 << " " << 3.14 << " " << 11.1 << " " << 0.1 << " " << 0.1 << " " << 0.1 << " " << 0.1 << endl;
      
      for(int i = 0 ; i < 3 ; i++)
	{
	  //	  gz = 1. * R->Uniform(90., 100.);
	  double theta = (0.025 + 1. * R->Uniform(-0.001, 0.001)), phi = 1. * R->Uniform(-0.001, 0.001);
	  gz = 100. * TMath::Cos(0.025);
	  gx = 100. * TMath::Sin(theta) * TMath::Cos(phi);
	  gy = 100. * TMath::Sin(phi);

	  gE = TMath::Sqrt(gx * gx + gy * gy + gz * gz + 0.938271998 * 0.938271998);
	  
	  myfile << i << " 0 " << "1 " << "2212 " << "1 " << "1 " << gx << " " << gy << " " << gz << " " << gE << " " << "0.938271998 " << "0.0 0.0 0.0" << endl;
	  gE_distribution->Fill(gE);
	  pos_map->Fill(gx, gy);
	  //	  N_events++;
	}
	  		  
      //      if( N_events % 500 == 0 )
      //	cout << N_events << "th events finished......" << endl;
	//	cout << i << " th: [" << gx << " " << gy << " " << gz << " " << gE << "] events finished......" << endl;

      if( loop > 20000)
	break;

      loop++;
    }

  cout << "while loop: " << loop << endl;
  
  auto c1 = new TCanvas("c1", "c1", 1000, 500);
  c1->Divide(2,1);
  c1->cd(1);
  gE_distribution->Draw();
  c1->cd(2);
  pos_map->Draw("colorz");

  
}
