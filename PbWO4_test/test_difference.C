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
  string line, root_file_name, data_file_path, file_0_str, crystal; //for storing words
  string Crtstal_ID_str[10] = {"010647", "010667", "030651", "040668", "050658", "050666", "080671", "160669", "994724", "994763"};
  stringstream line_string;
  double L, T;
  int count=0, pic_count=0;
  int file_0[10] = {67, 76, 32, 114, 40, 96, 87, 51, 105, 59};
  const char *CID;
    

  
  TH1F *h1[10], *h2[10], *h3[10], *h3_sud[10];
  TCanvas *c1 = new TCanvas("c1", "c1", 1500, 1000);
  c1->Divide(4,3);

  ifstream bookread; 
    
  for(int i = 0 ; i < 10 ; i++)
    {
      c1->cd((i + 1));
      pic_count = 0;
      
      int file_start = file_0[i];
      string ID = Crtstal_ID_str[i];

      crystal = "PbWO4 Crystal " + Crtstal_ID_str[i];
      CID = (crystal).c_str();

      h1[count] = new TH1F(Form("h1_%d", i),"test1", 270, 260, 800);
      h2[count] = new TH1F(Form("h2_%d", i),"test2", 270, 260, 800);
      h3[count] = new TH1F(Form("h3_%d", i),"difference of [test1, test2]", 270, 260, 800);
      
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
            
      for(int j = file_start ; j < (file_start + 3) ; j++)
	{
	  //	  if(pic_count == 2)
	  //	    continue;
	  
	  data_file_path = "./PbWO4_test_15122020_New_crystal/Crystal_" + ID + "/Sample"+ j +".Sample.Raw.asc";
	  //	  cout << data_file_path << endl;
	  
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

		      line_string << line;
		      line_string >> L >> T;  // L:wavelength, T:transmittance
		      //		      cout << L << "  " << T << endl;

		      if( (pic_count % 3) == 0 )
			{
			  int n_bin = h1[count]->FindBin(L);
			  h1[count]->SetBinContent(n_bin, T);
			}
		      if( (pic_count % 3) == 2 )
			{
			  int n_bin = h2[count]->FindBin(L);
			  h2[count]->SetBinContent(n_bin, T);
			}
		    }
		}
	    }
	  //	  h1[count]->Draw();
	  //	  h2[count]->Draw("same");

	  bookread.close();
	  pic_count++;
	} // loop the test results I want to draw 
      *h3[count] = *h2[count] - *h1[count];
      //      h3[count]->Draw("same");
      h3[count]->Draw();
      cout << "count: " << count << endl;
      count++;

    } // loop 10 crystal


}
