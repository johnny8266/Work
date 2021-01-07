#ifndef DVCSmodel_h
#define DVCSmodel_h
#include <iostream>
#include <time.h>
#include "TGVKelly.h"
#include "TGenDVCS.h"
#include "TFoam.h"
#include "TRandom3.h"

extern "C" { void cross_4he_(double*,double*,double*,double*,double*,double*,double*,double*);}

extern "C" { void cross_4he_re0_(double*,double*,double*,double*,double*,double*,double*,double*);}

extern "C" { void xsection_freeproton_(double*,double*,double*,double*,double*,double*,double*,double*);}

constexpr double mn=0.9396;
constexpr double mp=0.9383;
constexpr double mN=(mp+mn)/2.;
constexpr double e4he=-0.028295;
constexpr double m4he =2.*(mp+mn)+e4he;

class DVCSmodel {
    private:
        static int counter;
        static time_t TimeC;
        static int ModelNb;
        static int DimNb;
        static double Eb, M, xB, Q2, t, phi;
        static double EBene;
        static double HBene;
        static double HBmom;
	static double EBpol;
	static double s;
        static double mass;

        static double *MCV;
        static double *ExV;

        static double ym, yM, qm, qM, Wm, TM, tmi, tR, xM;

        // Set the kinematic limits in order: y min, y Max, Q2 min, Q2 Max, 
        // W2 min, theta Max, t min and t Max
        void SetKinLimits(double, double, double, double, double, double, double, double);

        static double Gettm();
        static double GetqM();
	static double Getxa();

        static double Model1();
        static double Model2();
	static double Model3();
	static double Model4();
        static double Model5();
        //...

    public:
        DVCSmodel(int, double, double, double, double, double, double, double, double, double, double, double);
        ~DVCSmodel();

        // Gives the number of dimensions of the model
        static int GetDim(){return DimNb;}
        static double Model(int,double*);
        static void Transform();
        static void Transform(double*);

        static int GetCount(){return counter;}
        static int GetModNb(){return ModelNb;}
        static int GetDimNb(){return DimNb;}
        static double GetEBene(){return EBene;}
        static double GetHBene(){return HBene;}
        static double GetHBmom(){return HBmom;}
	static double GetEBpol(){return EBpol;}
	static double GetS(){return s;}
	static double Getmass(){return mass;}
};

#endif  // DVCSmodel_h
