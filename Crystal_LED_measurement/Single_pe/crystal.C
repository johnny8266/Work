#include <TMath.h>
#include <TF1NormSum.h>
#include <TF1.h>
#include <cmath>

// see math/mathcore/src/PdfFuncMathCore.cxx in ROOT 6.x
double crystalball_function(double x, double alpha, double n, double sigma, double mean) {
  // evaluate the crystal ball function
  if (sigma < 0.)
    return 0.;
  double z = (x - mean)/sigma; 
  if (alpha < 0)
    z = -z; 
  double abs_alpha = std::abs(alpha);
  // double C = n/abs_alpha * 1./(n-1.) * std::exp(-alpha*alpha/2.);
  // double D = std::sqrt(M_PI/2.)*(1.+ROOT::Math::erf(abs_alpha/std::sqrt(2.)));
  // double N = 1./(sigma*(C+D));
  if ( z > abs_alpha )
    {
      //double A = std::pow(n/abs_alpha,n) * std::exp(-0.5*abs_alpha*abs_alpha);
      double nDivAlpha = n/abs_alpha;
      double AA =  std::exp(-0.5*abs_alpha*abs_alpha);
      double B = nDivAlpha -abs_alpha;
      double arg = nDivAlpha/(B-z);
      
      return AA * std::pow(arg,n);
    }
  else
    {
      return std::exp(- 0.5 * z * z);
    }
}

double crystalball_function(const double *x, const double *p)
{
  // if ((!x) || (!p)) return 0.; // just a precaution
  // [Constant] * ROOT::Math::crystalball_function(x, [Alpha], [N], [Sigma], [Mean])
  double pedestal = p[0] * exp(-0.5 * TMath::Power(((x[0] - p[1]) / p[2]), 2));
  double pe_1 = (p[3] * crystalball_function(x[0], p[6], p[7], p[5], p[4]));
  double pe_2 = (p[8] * crystalball_function(x[0], p[11], p[12], p[10], p[9]));
  double pe_3 = (p[13] * crystalball_function(x[0], p[16], p[17], p[15], p[14]));
  return pedestal + pe_1 + pe_2 + pe_3;
}

double Multi_pe_3(double *x, double *par)
{
  double pedestal = par[0] * exp(-0.5 * TMath::Power(((x[0] - par[1]) / par[2]), 2));
  double pe_1 = par[3] * exp(-0.5 * TMath::Power(((x[0] - par[4]) / par[5]), 2));
  double pe_2 = par[6] * exp(-0.5 * TMath::Power(((x[0] - par[7]) / par[8]), 2));
  double pe_3 = par[9] * exp(-0.5 * TMath::Power(((x[0] - par[10]) / par[11]), 2));
  //  double background = par[9] / x[0];
  return pedestal + pe_1 + pe_2 + pe_3;
}




void crystal()
{

  TCanvas *c1 = new TCanvas("c1", "c1", 800, 800);
  //  c1->Divide(2, 1);

  TF1 *pede = new TF1("pede", "gaus(x)", 40., 80.);
  pede->SetParLimits(0, 10000., 13000.);
  pede->SetParLimits(1, 50., 70.);
  pede->SetParLimits(2, 1., 20.);


  auto *pedestal_multi_pe_3 = new TF1("pedestal_multi_pe_3", Multi_pe_3, 40., 500., 12);
  pedestal_multi_pe_3->SetLineColor(3);
  pedestal_multi_pe_3->SetParLimits(0, 2000., 15000.);
  pedestal_multi_pe_3->SetParLimits(1, 55., 70.);
  pedestal_multi_pe_3->SetParLimits(2, 2., 20.);
  
  pedestal_multi_pe_3->SetParLimits(3, 10., 1000.);
  pedestal_multi_pe_3->SetParLimits(4, 90., 120.);
  pedestal_multi_pe_3->SetParLimits(5, 2., 30.);

  pedestal_multi_pe_3->SetParLimits(6, 10., 1000.);
  pedestal_multi_pe_3->SetParLimits(7, 120., 185.);
  pedestal_multi_pe_3->SetParLimits(8, 2., 100.);

  pedestal_multi_pe_3->SetParLimits(9, 10., 1000.);
  pedestal_multi_pe_3->SetParLimits(10, 180., 300.);
  pedestal_multi_pe_3->SetParLimits(11, 10., 150.);
  

  
  TF1 *f_cry = new TF1("f_cry", crystalball_function, 40., 400., 18);
  f_cry->SetLineColor(2);
  f_cry->SetParLimits(0, 2000., 2500.);
  f_cry->SetParLimits(1, 50., 70.);
  f_cry->SetParLimits(2, 1., 20.);

  f_cry->SetParLimits(3, 500., 1500.);
  f_cry->SetParLimits(4, 100., 150.);
  f_cry->SetParLimits(5, 10., 100.);
  f_cry->SetParLimits(6, 10., 100.);
  f_cry->SetParLimits(7, 1., 100.);
  
  f_cry->SetParLimits(8, 100., 800.);
  f_cry->SetParLimits(9, 150., 300.);
  f_cry->SetParLimits(10, 10., 100.);
  f_cry->SetParLimits(11, 1., 30.);
  f_cry->SetParLimits(12, 1., 30.);

  f_cry->SetParLimits(13, 100., 800.);
  f_cry->SetParLimits(14, 200., 300.);
  f_cry->SetParLimits(15, 10., 100.);
  f_cry->SetParLimits(16, 1., 30.);
  f_cry->SetParLimits(17, 1., 30.);


  

  /*
  TF1 *f_cry = new TF1("f_cry", "crystalball(x)", 80., 200.);
  f_cry->SetParLimits(0, 400., 800.);
  f_cry->SetParLimits(1, 100., 150.);
  f_cry->SetParLimits(2, 1., 50.);
  f_cry->SetParLimits(3, 0.1, 10.);
  f_cry->SetParLimits(4, 0.1, 30.);
  */
  
  //  TF1 *f_sum = new TF1("f_sum", "pede(0)+f_cry(3)", 40, 200);
  Double_t par[9];
  
  TH1D* h1;
      
  TFile* gfile = TFile::Open("./2.4Vpp_820mV_offset.root");
  TDirectory* dir = gFile->GetDirectory("Energy");

  dir->GetObject("_R_EnergyCH7@DT5730_1204", h1);
  h1->RebinX(3);

  /*
  h1->Fit(pede, "R");
  pede->GetParameters(&par[0]);
  h1->Fit(f_cry, "R+");
  pede->GetParameters(&par[3]);
  
  f_sum->SetParameters(par);
  h1->Fit("f_sum", "R+");
  */

  /*
  h1->Fit("f_cry", "R", "", 40, 400);
  h1->GetXaxis()->SetRangeUser(0, 500);
  h1->SetStats(0);
  h1->Draw();
  */
  
  //  h1->Fit("f_cry", "R", "", 40, 400);
  h1->Fit("pedestal_multi_pe_3", "R", "", 40, 400);
  h1->GetXaxis()->SetRangeUser(0, 500);
  h1->SetStats(0);
  h1->Draw();
  
  // TF1 *f_cry = new TF1("f_cry", "crystalball(x)+crystalball(x)", 3., 10.);
  // f_cry->SetParNames("f_cry");
  // f_cry->SetParameters(1, 6., 0.3, 2, 1.5);

  //    f_cry_1->Draw();
  //  f_cry->Draw();

  // TF1 *f_sum = new TF1("f_sum", "f_cry_1+f_cry", -2., 10.);
  // f_sum->Draw();
  
  /*
      TF1 *fcos = new TF1 ("fcos", "[0]*cos(2*x)", 0., 10.);
      fcos->SetParNames( "cos");
      fcos->SetParameter( 0, 0.5);
 
      TF1 *fsin = new TF1 ("fsin", "[0]*sin(x)", 0., 10.);
      fsin->SetParNames( "sin");
      fsin->SetParameter( 0, 2.);
 
      TF1 *fsincos = new TF1 ("fsc", "fcos+fsin", 0., 10.);
           fsincos->Draw();
      //      TF1 *fs2 = new TF1 ("fs2", "fsc+fsc", 0., 10.);
      //      fs2->Draw();
      */
}
