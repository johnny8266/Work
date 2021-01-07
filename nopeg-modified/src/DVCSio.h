#ifndef DVCSio_h
#define DVCSio_h

#include <iostream>
#include <TFile.h>
#include <TTree.h>

class DVCSio {
    protected:
        TFile *file;
        TTree *tree;

        void branching();
        void InitVars();

        // Variables to be recorded in root
        int ievent=0, PIDlBeam, PIDhBeam;
        float ElBeam, EhBeam  ;
        float Q2,W,Nu,XBj,y,t,phih;
        int Nb,PID[10];
        float E[10],px[10],py[10],pz[10];

    public:

        DVCSio(char*);
        ~DVCSio();

        void Fill(double,double,double,double);

};

#endif //DVCSio_h
