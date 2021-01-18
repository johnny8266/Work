

void test_difference()
{
  string line, root_file_name, data_file_path, file_0_str, crystal; //for storing words
  string Crtstal_ID_str[10] = {"010647", "010667", "030651", "040668", "050658", "050666", "080671", "160669", "994724", "994763"};
  int count=0, file_0[10] = {68, 77, 33, 115, 41, 97, 88, 52, 106, 60}, pic_count=0;

  
  TH1F *h1[10], *h2[10], *h3_sud[10];
  TCanvas *c1 = new TCanvas("c1", "c1", 1200, 800);
  c1->Divide(4,3);

  
  for(int i = 0 ; i < 10 ; i++)
    {
      c1->cd((i + 1));
      
      int file_start = file_0[i];
      string ID = Crtstal_ID_str[i];
      
      h3_sud[i] = new TH1F(Form("h3_sud_%d", i), "test2", 100, -2, 4);
      h3_sud[i]->SetStats(0);
      h3_sud[i]->SetTitle("");
      
      for(int j = file_start ; j < (file_start + 2) ; j++)
	{
	  data_file_path = "./PbWO4_test_15122020_New_crystal/Crystal_" + ID + "/Sample"+ j +".Sample.Raw.asc";
	  h1[count] = new TH1F(Form("h1_%d", i),"test1", 200, 300, 4);
	  h2[count] = new TH1F(Form("h2_%d", i),"test2", 200, -4, 4);
	  
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
		      line_string >> L >> T;

		      wave_length.push_back(L);
		      transmittance.push_back(T);
		    }
		}
	    }

	  h1[j]->Draw();
	  count++;
	  
	  TPad *subpad = new TPad("subpad", "", 0., 0., 1., 0.2);
	  subpad->Draw();
	  subpad->cd();

	}
      
      h3_sud[i]->Draw();
      count++;
    }


}
