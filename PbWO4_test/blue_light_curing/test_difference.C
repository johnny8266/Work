#include <vector> 
#include <sstream> 
#include <string>
#include <cstring>
#include <iostream> 
#include <iomanip> 
#include <fstream>
#include <sstream>
using namespace std; 

void test_difference()
{
  const int crystal_nums = 2;
  string line, root_file_name, data_file_path, file_0_str, crystal; //for storing words
  //  string Crtstal_ID_str[crystal_nums] = {"080671"};
  string Crtstal_ID_str[crystal_nums] = {"050658", "994724"};
  stringstream line_string;
  double L, T;
  int count=0, pic_count=0, flag=0;
  int file_0[crystal_nums] = {257};
  //  int file_0[crystal_nums] = {269};
  const char *CID;
  vector<double> Tb, Ta, Tlambda;

  
  TH1F *h1[crystal_nums], *h2[crystal_nums], *h3[crystal_nums],
       *h4[crystal_nums], *h5[crystal_nums], *h6[crystal_nums];
  TLine *add_line_flat = new TLine(260, 1.1, 800, 1.1);
  TLine *add_line_vert = new TLine(340, -3., 340, 3.);
  TCanvas *c1 = new TCanvas("c1", "c1", 800, 800);
  //  c1->Divide(2,1);

  ifstream bookread; 
    
  for(int i = 0 ; i < crystal_nums ; i++)
    {
      //      c1->cd(i+1);
      pic_count = 0;
      
      int file_start = file_0[i];
      string ID = Crtstal_ID_str[i];

      h1[count] = new TH1F(Form("h1_%d", i),"test_before", 270, 260, 800);
      h2[count] = new TH1F(Form("h2_%d", i),"test_after_1", 270, 260, 800);
      h3[count] = new TH1F(Form("h3_%d", i),"test_after_2", 270, 260, 800);
      h4[count] = new TH1F(Form("h4_%d", i),"difference of [test1, before]", 270, 260, 800);
      h5[count] = new TH1F(Form("h5_%d", i),"difference of [test2, before]", 270, 260, 800);
      h6[count] = new TH1F(Form("h6_%d", i),"dk", 270, 260, 800);

      crystal = "PbWO4 Crystal " + Crtstal_ID_str[i];
	
      CID = (crystal).c_str();
      
      h1[count]->SetStats(0);
      h1[count]->GetXaxis()->SetTitle("#lambda [nm]");
      h1[count]->GetYaxis()->SetTitle("T difference [%]");
      h1[count]->SetTitle(CID);
      
      h2[count]->SetStats(0);
      h2[count]->GetXaxis()->SetTitle("#lambda [nm]");
      h2[count]->GetYaxis()->SetTitle("T difference [%]");
      h2[count]->SetTitle(CID);
      
      h3[count]->SetStats(0);
      h3[count]->GetXaxis()->SetTitle("#lambda [nm]");
      h3[count]->GetYaxis()->SetTitle("T difference [%]");
      h3[count]->SetTitle(CID);

      h4[count]->SetStats(0);
      h4[count]->GetXaxis()->SetTitle("#lambda [nm]");
      h4[count]->GetYaxis()->SetTitle("T difference [%]");
      h4[count]->SetTitle(CID);

      h5[count]->SetStats(0);
      h5[count]->GetXaxis()->SetTitle("#lambda [nm]");
      h5[count]->GetYaxis()->SetTitle("T difference [%]");
      h5[count]->SetTitle(CID);

      h6[count]->SetStats(0);
      h6[count]->GetXaxis()->SetTitle("#lambda [nm]");
      h6[count]->GetYaxis()->SetTitle("dk [1/m]");
      h6[count]->SetTitle(CID);
      h6[count]->SetMaximum(1.5);
      h6[count]->SetMinimum(0.);
      
            
      for(int j = file_start ; j < (file_start + 3) ; j++)
	{
	  
	  data_file_path = "./" + ID + "/Sample"+ j +".Sample.Raw.asc";
	  cout << data_file_path << endl;
	  
	  const char *dfp = data_file_path.c_str();
	  bookread.open(dfp);
	  if(bookread.is_open())
	    { 
	      while(!(bookread.eof()))
		{ 
		  line = ""; 
		  getline(bookread, line); 
		  
		  //		  if( (line.size() == 21) || (line.size() == 20) || (line.size() == 13) || (line.size() == 14) )
		  if( (line.size() == 21) || (line.size() == 20) )
		    {

		      line_string << line;
		      line_string >> L >> T;  // L:wavelength, T:transmittance
		      cout << L << "  " << T << endl;

		      if( (pic_count % 3) == 0 )
			{
			  int n_bin = h1[count]->FindBin(L);
			  h1[count]->SetBinContent(n_bin, T);     // measurements before irradiation
			  Tb.push_back(T);
			  Tlambda.push_back(L);
			}
		      if( (pic_count % 3) == 1 )
			{
			  int n_bin = h2[count]->FindBin(L);
			  h2[count]->SetBinContent(n_bin, T);     // 1st. measurement after irradiation
			  Ta.push_back(T);
			}
		      if( (pic_count % 3) == 2 )
			{
			  int n_bin = h3[count]->FindBin(L);
			  h3[count]->SetBinContent(n_bin, T);     // 2nd. measurement after irradiation
			}
		      
		    } // select the word
		} // read file
	    } // open n file

	  bookread.close();
	  pic_count++;
	} // loop the test results I want to draw

      /*      
      for(int tc = 0 ; tc < Tb.size() ; tc++)
	{
	  L = Tlambda[tc];
	  int n_bin = h6[count]->FindBin(L);
	  if( (Tb[tc] > 0.) && (Ta[tc] > 0.) && L > 340. )
	    {
	      double dk = TMath::Log(Tb[tc] / Ta[tc]) / 0.2;;
	      if( (dk < 0.) && Tb[tc] > 1. ) // examine the abnormal of the transmittance
		cout << Tb[tc] << " " << Ta[tc] << " " << L << endl;
	      if( (dk > 1.1) && (flag == 0) )
		flag = 1;
	      h6[count]->SetBinContent(n_bin, dk);
	    }
	}
      */
      //      h6[count]->Draw();
      //      add_line_flat->SetLineStyle(2);
      //      add_line_flat->SetLineColor(2);
      //      add_line_flat->Draw("same");
      //      add_line_vert->Draw("same");
      
      
      //      *h4[count] = *h2[count] - *h1[count];
      *h4[count] = *h3[count] - *h1[count];
      //      h4[count]->SetLineColor(1);
      //      h5[count]->SetLineColor(2);

      h1[count]->SetMarkerColor(1);
      h2[count]->SetMarkerColor(2);
      //      h3[count]->SetMarkerColor(3);

      h1[count]->SetMarkerStyle(6);
      h2[count]->SetMarkerStyle(6);
      h3[count]->SetMarkerStyle(6);
      h4[count]->SetMarkerStyle(6);
      
      h1[count]->SetMarkerSize(5);
      h2[count]->SetMarkerSize(5);
      h3[count]->SetMarkerSize(5);
      h4[count]->SetMarkerSize(5);
      
      //      h1[count]->Draw("p");
      //      h2[count]->Draw("p same");
      h4[count]->Draw("p");
      

      cout << "count: " << count << " " << flag << " finish ..." << endl << endl;
      count++;  flag = 0;  Ta.clear();  Tb.clear();  Tlambda.clear();
    } 

}
