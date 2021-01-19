#include <iostream>
using namespace std;

void test()
{
  TH1F *h = new TH1F("h", "h", 270, 260, 800);
  cout << h->FindBin(278) << endl;;
}
