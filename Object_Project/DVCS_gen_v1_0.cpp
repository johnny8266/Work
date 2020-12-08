#include "calculation.h"
#include "dvcs_beam.h"
#include <iostream>
using namespace std;

int main(int argc, char* argv[])
{
  calculation c1;
  dvcs_beam test;
  string a="";
  int e_E=0., p_E=0., add=0.;
  cout << "Input the electron beam energy[GeV]: " << endl;
  cin >> e_E;
  cout << "Input the proton beam energy[GeV]: " << endl;
  cin >> p_E;

  test.SetQ2(e_E);

  test.LeptonicProcess();
  test.HadronicProcess();
  /*
  add = e_E + p_E;
  c1.set_add(add);
  cout << c1.get_add() << endl;
  
  for(int i = 0 ; i < argc ; i++)
    cout << argv[i] << "  ";
  cout << endl;
  */
  
  return 0;
}
