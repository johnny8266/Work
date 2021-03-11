#include <vector> 
#include <sstream> 
#include <string> 
#include <iostream> 
#include <iomanip> 
#include <fstream>
#include <sstream>
#include <stdlib.h>
using namespace std; 

void Read_test()
{ 
  string line, whichbook; //for storing words 
  vector<string> words; //unspecified size vector 
  vector<double> wave_length, transmittance;
  stringstream line_string;
  int count=0;
  double L, T;

  TH1F *h1 = new TH1F("h1", "h1", 100, 0., 2000.);
  //    cout << "Welcome to the book analysis program. Please input the filename of the book you would like to analyze: "; 
  //    cin >> whichbook; 
  //    cout << endl; 

  ifstream bookread; 

  //    bookread.open(whichbook.c_str());
  bookread.open("record.txt");

  if(bookread.is_open())
   { 
      while(!(bookread.eof()))
	{ 
	  line = ""; 
	  getline(bookread, line); 

	  count++;

	  //	  line_string << line;
	  //	  line_string >> L;

	  L = atof(line.c_str());;
	  cout << L << endl;
	  h1->Fill(L);
	  //		cout << L << " " << T << endl;
	  
	}
    }

  cout << count << endl;
  h1->Draw();
}
