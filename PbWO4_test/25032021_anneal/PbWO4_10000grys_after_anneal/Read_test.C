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
  string line, root_file_name, data_file_path, file_0_str, crystal; //for storing words
  string Crtstal_ID_str[10] = {"010647", "010667", "030651", "040668", "050658", "050666", "080671", "160669", "994724", "994763"};

  vector<string> words; //unspecified size vector 
  vector<double> wave_length, transmittance;
  stringstream line_string;
  int count=0, file_0[10] = {69, 78, 34, 116, 42, 98, 89, 53, 107, 61}, pic_count=0;
  double L, T;
  TGraph *gr1[20], *gr_mark[10];
  TLegend *legend[11], *legend_mark[10];

  
  ifstream bookread; 
  //  file_0 = atoi(file_0_str.c_str());

  //Read the raw files
  //
  for(int j = 0 ; j < 10 ; j++)
    {
      int file_start = file_0[j];
      string ID = Crtstal_ID_str[j];
      for(int i = file_start ; i < (file_start + 2) ; i++)
	{
	  data_file_path = "./" + ID + "/Sample"+ i +".Sample.Raw.asc";
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
  int linecolor;
  double mark_x[3] = {360, 420, 620}, mark_y[3] = {35, 60, 70};
  TCanvas *c1 = new TCanvas("c1", "c1", 1200, 900);
  const char *CID;
  c1->Divide(4,3);
  for(int k = 0 ; k < 10 ; k++)
    {
      c1->cd((k+1));
      
      for(int l = 0 ; l < 2 ; l++)
	{
	  int m = k * 2 + l;
	  crystal = "PbWO4 Crystal " + Crtstal_ID_str[k];
	  CID = (crystal).c_str();
	  gr1[m]->SetName(Form("gr_%d", m));
	  gr1[m]->SetTitle(CID);
	  gr1[m]->GetXaxis()->SetTitle("#lambda [nm]");
	  gr1[m]->GetYaxis()->SetTitle("T [%]");
	  gr1[m]->GetXaxis()->SetRangeUser(260., 800.);
	  gr1[m]->GetHistogram()->SetMaximum(80.);

	  if( m % 2 == 0 )
	    {
	      gr1[m]->SetLineColor(1);
	      gr1[m]->SetLineWidth(1);
	      gr1[m]->Draw();
	    }
	  else if( m % 2 == 1 )
	    {
	      gr1[m]->SetLineColor(2);
	      gr1[m]->SetLineWidth(1);
	      gr1[m]->Draw("same");
	    }
	}
    }


  /*
  
  c1->cd(11);
  linecolor = 1;
  legend[10] = new TLegend(0.5, 0.3, 0.8, 0.65);
  legend[10]->SetBorderSize(0);
  for(int i = 0 ; i < 90 ; i++)
    {
      if( i == 0 )
	{
	  gr1[i]->GetXaxis()->SetRangeUser(260., 800.);
	  gr1[i]->GetHistogram()->SetMaximum(80.);
	  gr1[i]->SetTitle("First setup of 10 Crystals");
	  gr1[i]->SetLineColor(linecolor);
	  gr1[i]->SetLineWidth(1);
	  gr1[i]->Draw();
	  linecolor++;
	  int k = i / 9;
	  crystal = "Crystal " + Crtstal_ID_str[k];
	  CID = (crystal).c_str();
	  legend[10]->AddEntry(Form("gr_%d",i), CID, "l");
	}
      else if( (i % 9 == 0) && (i > 1) )
	{
	  gr1[i]->SetLineColor(linecolor);
	  gr1[i]->SetLineWidth(1);
	  gr1[i]->Draw("same");
	  linecolor++;
	  if(linecolor == 10) linecolor++;

	  int k = i / 9;
	  crystal = "Crystal " + Crtstal_ID_str[k];
	  CID = (crystal).c_str();
	  legend[10]->AddEntry(Form("gr_%d",i), CID, "l");
	}      
    }
  gr_mark[9]->Draw("same p");
  legend[10]->Draw("same");
  legend_mark[9]->Draw("same");
  */


  //  RootFile.Write();

}
