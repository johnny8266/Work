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
  const int crystal_nums = 9;
  string line, root_file_name, data_file_path, file_0_str, crystal; //for storing words
  //  string Crtstal_ID_str[10] = {"010647", "010667", "030651", "040668", "050658", "050666", "080671", "160669", "994724", "994763"};
  string Crtstal_ID_str[crystal_nums] = {"040668", "050658", "994724", "050666", "010647", "030651", "080671", "994763", "010667"};
  stringstream line_string;
  double L, T;
  int count=0, pic_count=0, flag=0;
  //  int file_0[10] = {67, 76, 32, 114, 40, 96, 87, 51, 105, 59};
  int file_0[crystal_nums] = {114, 40, 105, 96, 67, 32, 87, 59, 76};
  const char *CID;
  vector<double> Tb, Ta, Tlambda;

  
  TH1F *h1[crystal_nums], *h2[crystal_nums], *h3[crystal_nums],
       *h4[crystal_nums], *h5[crystal_nums], *h6[crystal_nums];
  TCanvas *c1 = new TCanvas("c1", "c1", 900, 1200);
  TLine *add_line_flat = new TLine(260, 1.1, 800, 1.1);
  TLine *add_line_vert = new TLine(340, -3., 340, 3.);
  c1->Divide(3,4);

  ifstream bookread; 
    
  for(int i = 0 ; i < crystal_nums ; i++)
    {
      if( i == 0 )
	c1->cd(9);
      
      if( i == 1 )
	c1->cd(6);
      
      if( i == 2 )
	c1->cd(11);
      
      if( i == 3 )
	c1->cd(8);
      
      if( i == 4 )
	c1->cd(5);
      
      if( i == 5 )
	c1->cd(10);
      
      if( i == 6 )
	c1->cd(7);
      
      if( i == 7 )
	c1->cd(4);
      
      if( i == 8 )
	c1->cd(2);
      
      pic_count = 0;
      
      int file_start = file_0[i];
      string ID = Crtstal_ID_str[i];

      crystal = "PbWO4 Crystal " + Crtstal_ID_str[i];
      CID = (crystal).c_str();

      h1[count] = new TH1F(Form("h1_%d", i),"test_before", 270, 260, 800);
      h2[count] = new TH1F(Form("h2_%d", i),"test_after_1", 270, 260, 800);
      h3[count] = new TH1F(Form("h3_%d", i),"test_after_2", 270, 260, 800);
      h4[count] = new TH1F(Form("h4_%d", i),"difference of [test1, before]", 270, 260, 800);
      h5[count] = new TH1F(Form("h5_%d", i),"difference of [test2, before]", 270, 260, 800);
      h6[count] = new TH1F(Form("h6_%d", i),"dk", 270, 260, 800);

      //      h1[count]->SetStats(0);
      h1[count]->GetXaxis()->SetTitle("#lambda [nm]");
      h1[count]->GetYaxis()->SetTitle("T difference [%]");
      h1[count]->SetTitle(CID);
      
      //      h2[count]->SetStats(0);
      h2[count]->GetXaxis()->SetTitle("#lambda [nm]");
      h2[count]->GetYaxis()->SetTitle("T difference [%]");
      h2[count]->SetTitle(CID);
      
      //      h3[count]->SetStats(0);
      h3[count]->GetXaxis()->SetTitle("#lambda [nm]");
      h3[count]->GetYaxis()->SetTitle("T difference [%]");
      h3[count]->SetTitle(CID);

      //      h4[count]->SetStats(0);
      h4[count]->GetXaxis()->SetTitle("#lambda [nm]");
      h4[count]->GetYaxis()->SetTitle("T difference [%]");
      h4[count]->SetTitle(CID);

      h5[count]->SetStats(0);
      h5[count]->GetXaxis()->SetTitle("#lambda [nm]");
      h5[count]->GetYaxis()->SetTitle("T difference [%]");
      h5[count]->SetTitle(CID);

      h6[count]->SetStats(0);
      h6[count]->GetXaxis()->SetTitle("#lambda [nm]");
      h6[count]->GetYaxis()->SetTitle("dk");
      h6[count]->SetTitle(CID);
      h6[count]->SetMaximum(3.);
      h6[count]->SetMinimum(-3.);
      
            
      for(int j = file_start ; j < (file_start + 3) ; j++)
	{
	  //	  if(pic_count == 2)
	  //	    continue;
	  
	  data_file_path = "./irradiation_test_27012021/" + ID + "/Sample"+ j +".Sample.Raw.asc";
	  //	  cout << data_file_path << endl;
	  
	  const char *dfp = data_file_path.c_str();
	  bookread.open(dfp);
	  if(bookread.is_open())
	    { 
	      while(!(bookread.eof()))
		{ 
		  line = ""; 
		  getline(bookread, line); 

		  //		  cout << "???";
		  
		  if( (line.size() == 21) || (line.size() == 20) )
		    {

		      line_string << line;
		      line_string >> L >> T;  // L:wavelength, T:transmittance
		      //		      cout << L << "  " << T << endl;

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

      
      for(int tc = 0 ; tc < Tb.size() ; tc++)
	{
	  L = Tlambda[tc];
	  int n_bin = h6[count]->FindBin(L);
	  if( (Tb[tc] > 0.) && (Ta[tc] > 0.) )
	    {
	      double dk = TMath::Log(Tb[tc] / Ta[tc]) / 0.2;;
	      if( (dk < 0.) && Tb[tc] > 1. ) // examine the abnormal of the transmittance
		cout << Tb[tc] << " " << Ta[tc] << " " << L << endl;
	      if( (dk > 1.1) && (flag == 0) )
		flag = 1;
	      h6[count]->SetBinContent(n_bin, dk);
	    }
	}

      
      //      *h4[count] = *h2[count] - *h1[count];
      //      *h5[count] = *h3[count] - *h1[count];
      //      *h4[count]  = TMath::Log(*h1[count] / *h2[count]) / 20.;
      
      //      h4[count]->SetLineColor(1);
      //      h5[count]->SetLineColor(2);
      
      h6[count]->Draw();
      //      if(flag == 1)
      add_line_flat->SetLineStyle(2);
      add_line_flat->SetLineColor(2);
      add_line_flat->Draw("same");
      add_line_vert->Draw("same");

      cout << "count: " << count << " " << flag << " finish ..." << endl << endl;
      count++;  flag = 0;  Ta.clear();  Tb.clear();  Tlambda.clear();
    } 

}
