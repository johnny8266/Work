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
  string line, root_file_name, data_file_path, file_0_str; //for storing words
  string Crtstal_ID_str[10] = {"010647", "010667", "030651", "040668", "050658", "050666", "080671", "160669", "994724", "994763"};
  /*
  cout << "Please input the Crystal ID you would like to analyze: "; 
  cin >> Crtstal_ID_str; 
  cout << endl << "Please input the start file number: ";
  cin >> file_0_str;
  cout << endl << endl;
  */

  //  root_file_name = "Crystal_" + Crtstal_ID_str + ".root";
  //  const char *rfn = root_file_name.c_str();
  //  TFile RootFile(rfn,"RECREATE");
  //  TTree *Tree = new TTree("Tree", "PbWO4 test result");

  vector<string> words; //unspecified size vector 
  vector<double> wave_length, transmittance;
  stringstream line_string;
  int count=0, file_0[10] = {68, 77, 33, 115, 41, 97, 88, 52, 106, 60}, pic_count=0;
  double L, T;
  TGraph *gr1[90];

  
  ifstream bookread; 
  //  file_0 = atoi(file_0_str.c_str());

  //Read the raw files
  //
  for(int j = 0 ; j < 10 ; j++)
    {
      int file_start = file_0[j];
      string ID = Crtstal_ID_str[j];
      for(int i = file_start ; i < (file_start + 9) ; i++)
	{
	  data_file_path = "./PbWO4_test_15122020_New_crystal/Crystal_" + ID + "/Sample"+ i +".Sample.Raw.asc";
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
	      gr1[pic_count] = new TGraph (count, x, y);
	      pic_count++;
	    }
	  count = 0;
	  wave_length.clear();
	  transmittance.clear();
	  bookread.close();
	}
    }
  cout << "pic_count: " << pic_count << endl;

  
  // Draw the plots
  //

  TCanvas *c1 = new TCanvas("c1", "c1", 1200, 900);
  c1->Divide(4,3);
  for(int k = 0 ; k < 10 ; k++)
    {
      c1->cd((k+1));
      int linecolor = 1;
      for(int l = 0 ; l < 9 ; l++)
	{
	  int m = k * 9 + l;
	  gr1[m]->SetTitle("");
	  gr1[m]->GetXaxis()->SetTitle("#lambda [nm]");
	  gr1[m]->GetYaxis()->SetTitle("T [%]");
	  if( m % 9 == 0 )
	    continue;
	  else if( m % 9 == 1 )
	    {
	      linecolor+=l;
	      gr1[m]->SetLineColor(linecolor);
	      gr1[m]->SetLineWidth(2);
	      gr1[m]->Draw();
	    }
	  else if( m % 9 == 8 )
	    {
	      double x0, x1, y0, y1;
	      x0 = gr1[m]->GetXaxis()->GetXmin();
	      x1 = gr1[m]->GetXaxis()->GetXmax();
	      y0 = gr1[m]->GetYaxis()->GetXmin();
	      y1 = gr1[m]->GetYaxis()->GetXmax();
	      TLine *line1 = new TLine(x0, 70., 620, 70.);
	      TLine *line2 = new TLine(x0, 60., 420, 60.);
	      TLine *line3 = new TLine(x0, 35., 360, 35.);
	      TLine *line4 = new TLine(620, y0, 620, 70.);
	      TLine *line5 = new TLine(420, y0, 420, 60.);
	      TLine *line6 = new TLine(360, y0, 360, 35.);
	      linecolor+=l;
	      gr1[m]->SetLineColor(linecolor);
	      gr1[m]->SetLineWidth(2);
	      gr1[m]->Draw("same");
	      line4->SetLineStyle(2);
	      line5->SetLineStyle(2);
	      line6->SetLineStyle(2);
	      line1->Draw("same");
	      line2->Draw("same");
	      line3->Draw("same");
	      line4->Draw("same");
	      line5->Draw("same");
	      line6->Draw("same");
	    }
	  else
	    {
	      linecolor+=l;
	      gr1[m]->SetLineColor(linecolor);
	      gr1[m]->SetLineWidth(2);
	      gr1[m]->Draw("same");
	    }
	  
	  //      gr1[k]->Write();
	}
    }

  //  RootFile.Write();

}
