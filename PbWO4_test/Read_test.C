#include <vector> 
#include <sstream> 
#include <string> 
#include <iostream> 
#include <iomanip> 
#include <fstream>
#include <sstream>
using namespace std; 

void Read_test()
{ 
  TFile RootFile("plot.root","RECREATE");
  TTree *Tree = new TTree("Tree", "Fill simulated DVCS parameters");
  string line, whichbook; //for storing words 
  vector<string> words; //unspecified size vector 
  vector<double> wave_length, transmittance;
  stringstream line_string;
  int count=0;
  double L, T;
  TGraph *gr1[11];

  //    cout << "Welcome to the book analysis program. Please input the filename of the book you would like to analyze: "; 
  //    cin >> whichbook; 
  //    cout << endl; 

  ifstream bookread; 

  //    bookread.open(whichbook.c_str());
  for(int i = 1 ; i < 12 ; i++)
    {
      bookread.open(Form("./PbWO4_test_09122020/Sample%d.Sample.Raw.asc", i));

      if(bookread.is_open())
	{ 
	  while(!(bookread.eof()))
	    { 
	      line = ""; 
	      getline(bookread, line); 

	      string lineToAdd = "";
	      if( (line.size() == 21) || (line.size() == 20) )
		{
		  count++;

		  line_string << line;
		  line_string >> L >> T;
		  //		cout << L << " " << T << endl;
		  wave_length.push_back(L);
		  transmittance.push_back(T);
		}
	    }
	  double x[count], y[count];
	  for(int j = 0 ; j < count ; j++)
	    {
	      x[j] = wave_length[j];
	      y[j] = transmittance[j];
	    }
	  int a = i - 1;
	  gr1[a] = new TGraph (count, x, y);
	  //	  gr1[i]->Draw("same");
	}
      count = 0;
      wave_length.clear();
      transmittance.clear();
      bookread.close();
    }

  TCanvas *c1, *c2, *c3;
  
  for(int k = 0 ; k < 11 ; k++)
    {
      if( k == 0 )
	continue;
      //      gr1[k]->Draw();
      else if( k == 1 )
	{
	  c1 = new TCanvas("c1", "c1", 600, 600);
	  gr1[k]->SetLineColor(k);
	  gr1[k]->SetLineWidth(2);
	  gr1[k]->Draw();
	}
      else if( k < 5 )
	{
	  gr1[k]->SetLineColor(k);
	  gr1[k]->SetLineWidth(2);
	  gr1[k]->Draw("same");
	}
      else if( k == 5 )
	{
	  c2 = new TCanvas("c2", "c2", 600, 600);
	  gr1[k]->SetLineColor(k);
	  gr1[k]->SetLineWidth(2);
	  gr1[k]->Draw();
	}
      else if( k < 9 )
	{
	  gr1[k]->SetLineColor(k);
	  gr1[k]->SetLineWidth(2);
	  gr1[k]->Draw("same");
	}
      else if( k == 9 )
	{
	  c3 = new TCanvas("c3", "c3", 600, 600);
	  gr1[k]->SetLineColor(k);
	  gr1[k]->SetLineWidth(2);
	  gr1[k]->Draw();
	}
      else if( k < 11 )
	{
	  gr1[k]->SetLineColor(k);
	  gr1[k]->SetLineWidth(2);
	  gr1[k]->Draw("same");
	}
      gr1[k]->Write();
    }

  RootFile.Write();

}
