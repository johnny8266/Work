#include <fstream>
#include <iostream>


void read()
{

  const int channel = 30, rsheet = 5, ratio = 5; 
  //  const int N_plot = rsheet * ratio;
  TH2I *si_response[rsheet][ratio];
  double R_value[rsheet] = {100., 1000., 2000., 5000., 10000.};
  double ratio_value[ratio] = {0., 0.25, 0.5, 0.75, 1.};

  for(int i = 0 ; i < rsheet ; i++)
    for(int j = 0 ; j < ratio ; j++)
      si_response[i][j] = new TH2I(Form("si_response_%d_%d", i, j), Form("si_response_%d_%d", i, j), 6, 0, 6, 5, 0, 5);

  
  Double_t V[channel][rsheet][ratio] = {};
  

  
  ifstream f("out.txt");
  //  ifstream f("test.txt");
  Double_t dum1,dum2,dum3;
  for(int ichannel = 0 ; ichannel < channel ; ichannel++)
  //  for(int ichannel = 0 ;  ichannel < 1 ; ichannel++ )
    {
      for(int irsheet = 0 ; irsheet < rsheet ; irsheet++)
      //      for(int irsheet = 0 ; irsheet < 1 ; irsheet++ )
	{
	  for(int iratio = 0 ; iratio < ratio ; iratio++)
	    {
	      //	      f >> dum1;
	      f >> V[ichannel][irsheet][iratio];
	      //	      cout << dum1 << " "; 
	    }
	}
    }


  /*
  for(int ichannel = 0 ; ichannel < channel ; ichannel++)
    {
      for(int irsheet = 0 ; irsheet < rsheet ; irsheet++)
      	{
	  for(int iratio = 0 ; iratio < ratio ; iratio++)
	    {
	      int x = ichannel % 6, y = ichannel / 6;
	      si_response[irsheet][iratio]->Fill(x, y, (V[ichannel][irsheet][iratio]));
	    }
	}
    }
  */



  for(int irsheet = 0 ; irsheet < rsheet ; irsheet++)
    {
      for(int iratio = 0 ; iratio < ratio ; iratio++)
	{
	  double W = 0.;
	  double base_v = 0.;
	  
	  for(int ichannel = 0 ; ichannel < channel ; ichannel++)
	    {
	      int x = ichannel % 6, y = ichannel / 6;

	      if( ichannel == 0 )
		base_v = V[ichannel][irsheet][iratio];

	      // if( (ichannel == 14) || (ichannel == 15) )
	      // 	continue;
	      // else
		{
		  si_response[irsheet][iratio]->Fill(x, y, (V[ichannel][irsheet][iratio]));		  
		  W = W + ((V[ichannel][irsheet][iratio]) - base_v);
		}


	      //	      W = W + (V[ichannel][irsheet][iratio]) * (V[ichannel][irsheet][iratio]) / R_value[irsheet];

	    }
	  cout << R_value[irsheet] << "   " << ratio_value[iratio] << "   " << W << endl;
	  //	  cout << W << endl;	 
	}
      cout << endl << endl;
    }






  


  
  int count = 1;
  TCanvas *c1 = new TCanvas("c1", "c1", 900, 900);
  c1->Divide(5,5);
  for(int irsheet = 0 ; irsheet < rsheet ; irsheet++)
    {      
      for(int iratio = 0 ; iratio < ratio ; iratio++)
	{
	  c1->cd(count);
	  si_response[irsheet][iratio]->SetStats(0);
	  si_response[irsheet][iratio]->Draw("colorz");
	  count++;
	}
    }
  



  /*
  for(int ichannel = 0 ;  ichannel < channel ; ichannel++ )
    for(int irsheet = 0 ; irsheet < rsheet ; irsheet++ )
      for(int iratio = 0 ; iratio < ratio ; iratio++ )
	cout << V[ichannel][irsheet][iratio] << endl;
  */

}
