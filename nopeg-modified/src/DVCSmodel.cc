#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <iostream>
#include "TMath.h"
#include "DVCSmodel.h"

using namespace std;

ClassImp(TFDISTR)

int DVCSmodel::counter = 0;
time_t DVCSmodel::TimeC = 0;
int DVCSmodel::ModelNb = 0;
int DVCSmodel::DimNb = 4;
double DVCSmodel::Eb = 2132.03;
double DVCSmodel::M = 0.938271998;
double DVCSmodel::xB = 0.;
double DVCSmodel::Q2 = 0.;
double DVCSmodel::t = 0.;
double DVCSmodel::phi = 0.;
double DVCSmodel::EBene = 0.;
double DVCSmodel::HBene = 0.;
double DVCSmodel::HBmom = 0.;
double DVCSmodel::EBpol = 0.;
double DVCSmodel::s = 0.;
double DVCSmodel::ym{0.};
double DVCSmodel::yM{1.};
double DVCSmodel::qm{1.};
double DVCSmodel::qM{10.};
double DVCSmodel::Wm{15.};
double DVCSmodel::TM{3.14};
double DVCSmodel::tmi{0.};
double DVCSmodel::tR{1.};
double DVCSmodel::xM{1.5}; // Not set in Config !
double DVCSmodel::mass=0.;
double *DVCSmodel::MCV=NULL;
double *DVCSmodel::ExV=NULL;

DVCSmodel::DVCSmodel(int MNb, double b1, double b2, double a1, double a2, double a3, 
        double a4, double a5, double a6, double a7, double a8,double a9){
    // TODO add safety variable for static functions
    ModelNb = MNb;
    EBene = b1;
    HBmom = b2;
  
    if(MNb==1 ||MNb ==2 || MNb ==3){
      mass =m4he;
    }
    else{
      mass = mp;
    }
    
    HBene = std::sqrt(b2*b2+mass*mass);
    EBpol = a9;
    s = mass*mass + 2*EBene*(HBene+HBmom);


    SetKinLimits(a1, a2, a3, a4, a5, a6, a7, a8);

    ExV = new double[DimNb];
    MCV = new double[DimNb];

    for(int i=0; i<DimNb; i++){
        ExV[i]=0;
        MCV[i]=0;
    }
}

DVCSmodel::~DVCSmodel(){
    delete []ExV;
    delete []MCV;
}

void DVCSmodel::SetKinLimits(double a1, double a2, double a3, double a4, double a5, double a6, 
        double a7, double a8)        {

    ym=a1; yM=a2; qm=a3; qM=a4; Wm=a5; TM=a6; tmi=a7; tR=a8;

    // Reset proper y limits if in conflict with other cuts
    double epsilon = 0.01;
    double Q2a = 4*EBene*EBene*(1-yM)*std::pow(std::tan(TM/2.),2);
    if(Q2a<qm){
        yM = 1 - qm/4/EBene/EBene/std::pow(std::tan(TM/2.),2)-epsilon;
        std::cout << "ERROR y cut are too wide in DVCSmodel::SetKinLimits" << std::endl;
        std::cout << "y max is reset to " << yM << std::endl;
        if(ym>yM) exit(EXIT_FAILURE);
    }
    double Q2b = mass*mass+2.*EBene*(HBene+HBmom)*ym-Wm;
    if(Q2b<qm){
        ym = (qm - mass*mass + Wm)/(2.*EBene*(HBene+HBmom)) +epsilon;
        std::cout << "ERROR y cut are too wide in DVCSmodel::SetKinLimits" << std::endl;
        std::cout << "y min is reset to " << ym << std::endl;
        if(ym>yM) exit(EXIT_FAILURE);
    }
    return;
}


double DVCSmodel::Model(int ND,double* V){
    if (ND!=DimNb) {
        std::cout << "ERROR Nb of dimension mismatch in DVCSmodel::Model" << std::endl;
        return 0;
    }
    for (int i=0; i<DimNb; i++){
        if (V[i]<0. || V[i]>1. ) {
            std::cout << "ERROR Variable out of bound in DVCSmodel::Model" << std::endl;
            return 0;
        }
        MCV[i]=V[i];
    }
    DVCSmodel::Transform();

    if (counter==0) TimeC=time(NULL);
    counter++;
    time_t TimeP = time(NULL)-TimeC;
    if (counter%10000 == 100) cout << "Time elapsed is "
        <<  int(TimeP/3600) << " h "
            <<  int(TimeP/60)%60 << " min "  
            <<  TimeP%60 << " sec, for " << counter << " calls." << std::endl;

    switch(ModelNb){
        case 1: return DVCSmodel::Model1(); 
        case 2: return DVCSmodel::Model2();
        case 3: return DVCSmodel::Model3();
	case 4: return DVCSmodel::Model4();
    }
    std::cout << "ERROR code should never get here ! DVCSmodel::Model" << std::endl;
    return 0;
}



// Determine Q2 max
double DVCSmodel::GetqM(){
    double epsilon = 0.01;
    double Q2a = 4*EBene*EBene*(1-ExV[0])*std::pow(std::tan(TM/2.),2);
    double Q2b = mass*mass+(s-mass*mass)*ExV[0]-Wm;
    double Q2c = xM*mp/mass*(s-mass*mass)*ExV[0];

    // Select the stringent Q2 max
    double Q2M = qM;
    if (Q2M>Q2a) Q2M=Q2a;
    if (Q2M>Q2b) Q2M=Q2b;
    if (Q2M>Q2c) Q2M=Q2c;
    Q2M -= epsilon;

    if(Q2M<qm){
        std::cout << "ERROR y cut are too wide, problem in DVCSmodel::GetqM" << std::endl;
        std::cout << ExV[0] << " " << Q2a << " " << Q2b << " " << Q2c  << " " << Q2M  << std::endl;
        return 0.;
    }

    return Q2M;
}

// Determine xb
double DVCSmodel::Getxa(){
    double xb=  mass/mp*(ExV[1]/((s-mass*mass)*ExV[0]));
    if (xb >xM){
        std::cout <<"WARNING: strange kinematics! Fix xb value in DVCSmodel::Getxa"<<std::endl;
        std::cout << xb <<" "<< ExV[1] <<" "<< ExV[0] <<" "<<std::endl;
        xb = xM;
    }

    double xa = mp/mass*xb ;
    return xa;
}

// Determine t min
double DVCSmodel::Gettm(){
    double xa = Getxa();
    double csi = xa/(2.-xa);
    double tmin = 4.*mass*mass*csi*csi/(1.-csi*csi);
    tmin = std::max(tmin,tmi);
    return tmin;
}


//Transform MC kinematic to Experimental
void DVCSmodel::Transform(double* event){
    if(ModelNb==1 || ModelNb==2 || ModelNb==3 || ModelNb==4){
        // Rescale var 0 to y
        event[0]=(event[0]*(yM-ym))+ym;
        ExV[0]=event[0];

        // Rescale var 1 to Q2
        double Q2M = GetqM();
        event[1]=(event[1]*(Q2M-qm))+qm;
        ExV[1]=event[1];

        // Rescale var 2 to t
        double tmin = Gettm();
        event[2]=(event[2]*tR)+tmin;
        ExV[2]=event[2];

        // Rescale var 3 to phi
        event[3]=event[3]*360.0;
        ExV[3]=event[3];
    }
    else
      {
	xB = event[0]*(0.03-0.005)+0.005;
	double Q2max = 2. * M * Eb * xB;
	if( Q2max > 15. ) Q2max = 13.;
	Q2 = event[1] * Q2max + 2.;
	t = -event[2];
	phi = event[3] * 2. * TMath::Pi();
      }
    return ;
}

//Transform MC kinematic to Experimental
void DVCSmodel::Transform(){
    if(ModelNb==1 || ModelNb==2 ||ModelNb ==3 || ModelNb==4){
        // Rescale var 0 to y
        ExV[0]=(MCV[0]*(yM-ym))+ym;

        // Rescale var 1 to Q2
        double Q2M = GetqM();
        ExV[1]=(MCV[1]*(Q2M-qm))+qm;

        // Rescale var 2 to t
        double tmin = Gettm();
        ExV[2]=(MCV[2]*tR)+tmin;

        // Rescale var 3 to phi
        ExV[3]=MCV[3]*360. ;
    }
    return;
}

// Model basyed on Mohammad basic parametrization of CLAS data
// Warning! This has been modified and will not provide results 
// coherent with the original work.
double DVCSmodel::Model1(){

    // Parameters from table 4.1 in Moh.thesis
    const double alp  = 2.5;
    const double q0   = 1.0;
    const double d    = 0.4;
    const double xc   = 0.2;
    const double c    = 0.2;
    const double b    = 11.0;
    const double beta = 12.0;

    double Dist = 0;

    // 4-dimensional distribution for Foam
    double q  = ExV[1];
    double xb = q/(2*mN*EBene*ExV[0]);
    double t  = ExV[2];
    double f  = ExV[3]/57.2957795;

    Dist =1/(1.0+pow((xb-xc)/c,2))*pow(q0/q,alp)*1/(pow((1+b*t),beta))*(1-d*(1-TMath::Cos(f)));  

    return Dist;
} // dist Moh

// Model Fucini et al. with the real part of CFF
double DVCSmodel::Model2(){

    double ris=0;

    double ebeam = EBene;
    double ehadr = HBene;
    double lambda = EBpol;

    double q = ExV[1];
    double y = ExV[0];
    double t = -ExV[2];
    double f = ExV[3];

    double xa = Getxa();
    double xb = mass/mp*xa;

    cross_4he_(&xb,&t,&q,&ebeam,&ehadr,&f,&lambda,&ris);
    //std::cout<<"RISU"<<"  "<<xb<<"  "<< t<<"  "<<q<<"  "<<ebeam<<"  "<<f<<"  "<< lambda<<"  "<<ris<<std::endl;
    return ris;
} 

// Model Fucini et al. without the real part of CFF
double DVCSmodel::Model3(){

    double risre0=0;

    double ebeam = EBene;
    double ehadr = HBene;
    double lambda = EBpol;

    double q = ExV[1];
    double y = ExV[0];
    double t = -ExV[2];
    double f = ExV[3];

    double xa = Getxa();
    double xb = mass/mp*xa;


    cross_4he_re0_(&xb,&t,&q,&ebeam,&ehadr,&f,&lambda,&risre0);
    //std::cout<<"RISU"<<"  "<<xb<<"  "<< t<<"  "<<q<<"  "<<ebeam<<"  "<<f<<"  "<< lambda<<"  "<<risre0<<std::endl;

    return risre0;
}

double DVCSmodel::Model4(){

    double crossfree=0.;

    double ebeam = EBene;
    double ehadr = HBene;
    double lambda = EBpol;

    double q = ExV[1];
    double y = ExV[0];
    double t = -ExV[2];
    double f = ExV[3];

    double xa = Getxa();
    double xb = mass/mp*xa;


    xsection_freeproton_(&xb,&t,&q,&ebeam,&ehadr,&f,&lambda,&crossfree);
    //     std::cout<<"RISU"<<"  "<<xb<<"  "<< t<<"  "<<q<<"  "<<ebeam<<"  "<<f<<"  "<< lambda<<"  "<<crossfree<<std::endl;

    return crossfree;
}

double DVCSmodel::Model5()
{
  double ConvGeV2nbarn = 0.389379304e+6; // Unit conversion
  double BHp, BHm, VCSp, VCSm, Ip, Im;
  double SigmaTotPlus, SigmaTotMoins, DVCSxsec;
  double* cffs = gEv->Interpol_CFF(Q2,xB,t);

  BHp = tgv->CrossSectionBH( Q2, xB, t, -phi, 1, 0, kTRUE );
  VCSp = tgv->CrossSectionVCS( Q2, xB, t, -phi, 1, 0, cffs[0], cffs[1], cffs[2], cffs[3], cffs[4], cffs[5], cffs[6], cffs[7], kTRUE );
  Ip = tgv->CrossSectionInterf( Q2, xB, t, -phi, 1, 0, -1, cffs[0], cffs[1], cffs[2], cffs[3], cffs[4], cffs[5], cffs[6], cffs[7], kTRUE );
  BHm = tgv->CrossSectionBH( Q2, xB, t, -phi, -1, 0, kTRUE );
  VCSm = tgv->CrossSectionVCS( Q2, xB, t, -phi, -1, 0, cffs[0], cffs[1], cffs[2], cffs[3], cffs[4], cffs[5], cffs[6], cffs[7], kTRUE );
  Im = tgv->CrossSectionInterf( Q2, xB, t, -phi, -1, 0, -1, cffs[0], cffs[1], cffs[2], cffs[3], cffs[4], cffs[5], cffs[6], cffs[7], kTRUE );
  SigmaTotPlus = BHp + VCSp + Ip;
  SigmaTotMoins = BHm + VCSm + Im;
  DVCSxsec = TMath::Pi() * ( SigmaTotPlus + SigmaTotMoins ) * ConvGeV2nbarn;// Total DVCS cross section in nb/GeV4

  return DVCSxsec;
}
