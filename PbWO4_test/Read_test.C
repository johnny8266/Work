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
  TGraph *gr1[90], *gr_mark[10];
  TLegend *legend[10], *legend_mark[10];

  
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
  
  double mark_x[3] = {360, 420, 620}, mark_y[3] = {35, 60, 70};
  TCanvas *c1 = new TCanvas("c1", "c1", 1400, 1050);
  c1->Divide(4,3);
  for(int k = 0 ; k < 10 ; k++)
    {
      c1->cd((k+1));
      legend[k] = new TLegend(0.5, 0.3, 0.8, 0.6);
      legend[k]->SetBorderSize(0);
      legend_mark[k] = new TLegend(0.5, 0.2, 0.8, 0.25);
      legend_mark[k]->SetBorderSize(0);
      int linecolor = 1;
      for(int l = 0 ; l < 9 ; l++)
	{
	  int m = k * 9 + l;
	  crystal = "PbWO4 Crystal " + Crtstal_ID_str[k];
	  const char *CID = (crystal).c_str();
	  gr1[m]->SetName(Form("gr_%d", m));
	  gr1[m]->SetTitle(CID);
	  gr1[m]->GetXaxis()->SetTitle("#lambda [nm]");
	  gr1[m]->GetYaxis()->SetTitle("T [%]");
	  gr1[m]->GetXaxis()->SetRangeUser(260., 800.);
	  gr1[m]->GetHistogram()->SetMaximum(80.);
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
	      /*
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
	      */
	      gr_mark[k] = new TGraph (3, mark_x, mark_y);
	      gr_mark[k]->SetName(Form("mark_%d", k));
	      gr_mark[k]->SetMarkerStyle(2);
	      gr_mark[k]->SetMarkerSize(5);

	      linecolor+=l;
	      gr1[m]->SetLineColor(linecolor);
	      gr1[m]->SetLineWidth(1);
	      gr1[m]->Draw("same");
	      gr_mark[k]->Draw("same p");
	      /*
	      line4->SetLineStyle(2);
	      line5->SetLineStyle(2);
	      line6->SetLineStyle(2);
	      line1->Draw("same");
	      line2->Draw("same");
	      line3->Draw("same");
	      line4->Draw("same");
	      line5->Draw("same");
	      line6->Draw("same");
	      */
	    }
	  else
	    {
	      linecolor+=l;
	      gr1[m]->SetLineColor(linecolor);
	      gr1[m]->SetLineWidth(2);
	      gr1[m]->Draw("same");
	    }
	  if(m % 9 == 1)
	    legend[k]->AddEntry(Form("gr_%d",m), "rotate 0 deg", "l");
	  if(m % 9 == 2)
	    legend[k]->AddEntry(Form("gr_%d",m), "rotate 90 deg", "l");
	  if(m % 9 == 3)
	    legend[k]->AddEntry(Form("gr_%d",m), "rotate 180 deg", "l");
	  if(m % 9 == 4)
	    legend[k]->AddEntry(Form("gr_%d",m), "rotate 270 deg", "l");
	  if(m % 9 == 5)
	    legend[k]->AddEntry(Form("gr_%d",m), "rotate 0 deg", "l");
	  if(m % 9 == 6)
	    legend[k]->AddEntry(Form("gr_%d",m), "rotate 90 deg", "l");
	  if(m % 9 == 7)
	    legend[k]->AddEntry(Form("gr_%d",m), "rotate 180 deg", "l");
	  if(m % 9 == 8)
	    legend[k]->AddEntry(Form("gr_%d",m), "rotate 270 deg", "l");
	}
      legend_mark[k]->AddEntry(Form("mark_%d", k), "     Criteria", "p");
      legend[k]->Draw("same");
      legend_mark[k]->Draw("same");
      
    }

  //  RootFile.Write();

}
