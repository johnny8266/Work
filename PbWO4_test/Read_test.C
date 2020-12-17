#include <vector> 
#include <sstream> 
#include <string>
#include <cstring>
#include <iostream> 
#include <iomanip> 
#include <fstream>
#include <sstream>
using namespace std; 

void Read_test()
{ 
  string line, Crtstal_ID, root_file_name, data_file_path, file_0_str; //for storing words 

  cout << "Please input the Crystal ID you would like to analyze: "; 
  cin >> Crtstal_ID; 
  cout << endl << "Please input the start file number: ";
  cin >> file_0_str;
  cout << endl << endl;

  root_file_name = "Crystal_" + Crtstal_ID + ".root";

  const char *rfn = root_file_name.c_str();
  
  TFile RootFile(rfn,"RECREATE");
  TTree *Tree = new TTree("Tree", "PbWO4 test result");

  vector<string> words; //unspecified size vector 
  vector<double> wave_length, transmittance;
  stringstream line_string;
  int count=0, file_0;
  double L, T;
  TGraph *gr1[8];

  
  ifstream bookread; 
  file_0 = atoi(file_0_str.c_str());

  //Read the raw files
  //
  for(int i = file_0 ; i < (file_0 + 8) ; i++)
    {
      data_file_path = "./PbWO4_test_15122020_New_crystal/Crystal_" + Crtstal_ID + "/Sample"+ i +".Sample.Raw.asc";
      const char *dfp = data_file_path.c_str();
      bookread.open(dfp);

      if(bookread.is_open())
	{ 
	  while(!(bookread.eof()))
	    { 
	      line = ""; 
	      getline(bookread, line); 

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
	  int a = i - file_0;
	  gr1[a] = new TGraph (count, x, y);
	  //	  gr1[i]->Draw("same");
	}
      count = 0;
      wave_length.clear();
      transmittance.clear();
      bookread.close();
    }


  // Draw the plots
  //
  TCanvas *c1, *c2;  
  for(int k = 0 ; k < 8 ; k++)
    {
      int linecolor = 1;
      if( k < 4 )
	{
	  c1 = new TCanvas("c1", "c1", 600, 600);
	  linecolor+=k;
	  gr1[k]->SetLineColor(linecolor);
	  gr1[k]->SetLineWidth(2);
	  gr1[k]->Draw();
	}
      else 
	{
	  c2 = new TCanvas("c2", "c2", 600, 600);
	  linecolor+=k;
	  gr1[k]->SetLineColor(linecolor);
	  gr1[k]->SetLineWidth(2);
	  gr1[k]->Draw();
	}
      gr1[k]->Write();
    }

  RootFile.Write();

}
