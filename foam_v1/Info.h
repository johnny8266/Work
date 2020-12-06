//#include <TFDISTR.h>
#include <iostream>

using namespace std;

class Info : public TFDISTR
{
 public:

  inline Double_t Get_Q2() const {return Q2_foam;}
  inline Double_t Get_Xb() const {return Xb_foam;}
  inline Double_t Get_t() const {return t_foam;}
  inline Double_t Get_phi() const {return phi_foam;}
};
