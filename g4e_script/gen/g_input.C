#include <iostream>
#include <fstream>
#include <ctime>
#include "TMath.h" 
#include "TRandom.h"

using namespace std;

void g_input()
{

  Int_t N_events = 0, loop = 0;

  auto gE_distribution = new TH1F("gE_distribution", "gE_distribution", 22., 0., 11.);

  auto pos_map = new TH2F("pos_map", "pos_map", 140, -1400., 1400., 140, -1400., 1400.);
  
  TRandom *R = new TRandom();
  R->SetSeed(0);
  ofstream myfile;
  myfile.open ("../../Data/g4e_input/g_input.txt");

  double gx = 0., gy = 0., gz = 0., gE = 0.;
  vector<double> g_x, g_y, g_z, g_E;

  
  while(1)
    {      
      int empty_flag = 0;
      g_x.clear();  g_y.clear();  g_z.clear();  g_E.clear();  
      
      for(int i = 0 ; i < 3 ; i++)
	{
	  gz = -1. * R->Uniform(0.1, 10.);
	  gx = -1. * gz * R->Uniform(-0.589, 0.589);
	  gy = -1. * gz * R->Uniform(-0.589, 0.589);
	  // gx = -1. * gz * R->Uniform(-0.335, 0.335);
	  // gy = -1. * gz * R->Uniform(-0.335, 0.335);
	  gE = TMath::Sqrt(gx * gx + gy * gy + gz * gz);

	  double ring_R = TMath::Sqrt( (gx / gz * 2240.) * (gx / gz * 2240.) + (gy / gz * 2240.) * (gy / gz * 2240.) );
	  
	  //	  if( ((((gx < (-0.067 * gz)) && (gx > 0.)) || ((gx > (0.067 * gz)) && (gx < 0.))) && (((gy < (-0.067 * gz)) && (gy > 0.)) || ((gy > (0.067 * gz)) && (gy < 0.)))) || (ring_R > 1320.) )
	  //	  if( (((gx < (-0.067 * gz)) && (gx > 0.)) || ((gx > (0.067 * gz)) && (gx < 0.))) && (((gy < (-0.067 * gz)) && (gy > 0.)) || ((gy > (0.067 * gz)) && (gy < 0.))) )
	  if( ((((gx < (-0.4 * gz)) && (gx > 0.)) || ((gx > (0.4 * gz)) && (gx < 0.))) && (((gy < (-0.4 * gz)) && (gy > 0.)) || ((gy > (0.4 * gz)) && (gy < 0.)))) || (ring_R > 1100.) ) // Glass region range
	    {
	      //	      cout << ring_R << endl;
	      empty_flag = 1;
	      break;
	    }

	  g_x.push_back(gx);  g_y.push_back(gy);  g_z.push_back(gz);  g_E.push_back(gE);
	}
      
      if( empty_flag == 0 )
	{
	  myfile << "3 " << 0.1 << " " << 0.01 << " " << -0.1 << " " << 3.14 << " " << 11.1 << " " << 0.1 << " " << 0.1 << " " << 0.1 << " " << 0.1 << endl;
	  for(int j = 0 ; j < 3 ; j++)
	    {
	      gx = g_x[j];  gy = g_y[j];  gz = g_z[j];  gE = g_E[j];
	      myfile << j << " 0 " << "1 " << "22 " << "1 " << "1 " << gx << " " << gy << " " << gz << " " << gE << " " << "0.0 " << "0.0 0.0 0.0" << endl;
	      gE_distribution->Fill(gE);
	      pos_map->Fill((gx / gz * 2240.), (gy / gz * 2240.));
	    }
	  N_events++;
	}
	  		  
      if( N_events % 500 == 0 )
	cout << N_events << "th events finished......" << endl;
	//	cout << i << " th: [" << gx << " " << gy << " " << gz << " " << gE << "] events finished......" << endl;

      if( N_events > 20000)
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
