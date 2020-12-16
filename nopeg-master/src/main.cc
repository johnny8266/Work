#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include "DVCSio.h"
#include "DVCSmodel.h"
#include "TFoam.h"
#include "TRandom3.h"

using namespace std;

int main(int argc, char * argv[])
{
    cout << "The Orsay Perugia Event Generator (TOPEG)" << endl;
    cout << "R. Dupre, S. Fucini, S. Scopetta" << endl << endl;
    string OptFileName;

    // Get the name of the configuration file
    if ( argv[1] == NULL ) OptFileName += "Config.txt";
    else OptFileName = (string) argv[1];

    // Read the configuration file
    ifstream OptFile;
    string line;
    if (argc > 1) OptFile.open(OptFileName);
    else OptFile.open("Config.txt");

    if (OptFile.is_open())
    {
        string opt0;
        double opt1,opt2;
        int opt3,opt4,opt5,opt6,opt7,opt8;
        double opt9,opt10,opt11,opt12,opt13,opt14,opt15,opt16;
        int opt17;
        double opt18,opt19;

        for (int i=0; i<20; i++  )
        {
            getline (OptFile,line);
            int found = line.find_first_of(" %");
            line.erase(found, string::npos);
            switch(i){
                case  0: opt0 = line;                  // Output file name 
                case  1: stringstream(line) >>  opt1;  // Electron beam energy
                case  2: stringstream(line) >>  opt2;  // Nuclear beam energy
                case  3: stringstream(line) >>  opt3;  // Z of nucleus
                case  4: stringstream(line) >>  opt4;  // A of nucleus
                case  5: stringstream(line) >>  opt5;  // Nb of event to generate
                case  6: stringstream(line) >>  opt6;  // Nb of TFoam cells
                case  7: stringstream(line) >>  opt7;  // Nb of samples per cell
                case  8: stringstream(line) >>  opt8;  // Seed
                case  9: stringstream(line) >>  opt9;  // y min
                case 10: stringstream(line) >> opt10;  // y Max
                case 11: stringstream(line) >> opt11;  // Q2 min
                case 12: stringstream(line) >> opt12;  // Q2 Max
                case 13: stringstream(line) >> opt13;  // W2 min
                case 14: stringstream(line) >> opt14;  // theta Max electron (degrees)
                case 15: stringstream(line) >> opt15;  // t min
                case 16: stringstream(line) >> opt16;  // t Range
                case 17: stringstream(line) >> opt17;  // Model nb (1: Param 2: Fucini et al)
                case 18: stringstream(line) >> opt18;  // Electron beam helicity
                case 19: stringstream(line) >> opt19;  // non used double
            }
        }
        OptFile.close();

        // TODO Add a safety function to check parameters
        cout << "These are the input parameters" << endl;
        cout << " Output file name                  " <<  opt0 << endl; 
        cout << " Electron beam energy              " <<  opt1 << endl; 
        cout << " Nuclear beam energy               " <<  opt2 << endl; 
        cout << " Z of nucleus                      " <<  opt3 << endl; 
        cout << " A of nucleus                      " <<  opt4 << endl; 
        cout << " Nb of event to generate           " <<  opt5 << endl; 
        cout << " Nb of TFoam cells                 " <<  opt6 << endl; 
        cout << " Nb of samples per cell            " <<  opt7 << endl; 
        cout << " Seed                              " <<  opt8 << endl; 
        cout << " y min                             " <<  opt9 << endl; 
        cout << " y Max                             " << opt10 << endl; 
        cout << " Q2 min                            " << opt11 << endl; 
        cout << " Q2 Max                            " << opt12 << endl; 
        cout << " W2 min                            " << opt13 << endl; 
        cout << " theta Max electron (degrees)      " << opt14 << endl; 
        cout << " t min                             " << opt15 << endl; 
        cout << " t Range                           " << opt16 << endl; 
        cout << " Model nb (1: Param 2: Full 3:Re=0)" << opt17 << endl; 
        cout << " Electron beam helicity            " << opt18 << endl; 
        cout << " non used double                   " << opt19 << endl; 
        cout << endl;

        // Get the name of the root file for the cells 
        int cutl1 = OptFileName.find_first_of("/\n");
        if (cutl1>0) OptFileName.erase(0, cutl1);
        string RooFileName = "grids";
        RooFileName += OptFileName;
        int found = RooFileName.find_first_of(". %\n");
        if (found>0) RooFileName.erase(found, string::npos);
        RooFileName += ".root";

        // Setup the root output class
        char *cstr = new char[opt0.length() + 1];
        strcpy(cstr, opt0.c_str());
        DVCSio * output = new DVCSio(cstr);
        delete [] cstr;

        // 4-dim vector generated in the MC run
        double *MCvect =new double[4];

        // Create random number generator
        TRandom *PseRan = new TRandom3();
        UInt_t s=opt8;
        PseRan->SetSeed(0);

        // Setup model
        DVCSmodel *M = new DVCSmodel(opt17,opt1,opt2,opt9,opt10,opt11,opt12,opt13,opt14,opt15,opt16,opt18);

        // Setup Foam
        TFoam *FoamX    = new TFoam("TOPEG");   // Create Simulator
        FoamX->SetPseRan(PseRan);  // Set random number generator

        TFile *GridFile = new TFile(RooFileName.c_str(),"READ");
        if ( GridFile->IsOpen() ) {
            // TODO for safety it could be interesting to store the opt* in the file to check them
            // before running
            printf("Grid file opened successfully.\n");
            FoamX = (TFoam*)GridFile->Get("FoamData");
            FoamX->SetRhoInt(DVCSmodel::Model);  // Set 4-dim distribution, included below
            PseRan->SetSeed(0);                  // Reset seed 
            FoamX->ResetPseRan(PseRan);
        }
        else {
            printf("Grid file not opened, starting initialize... \n");
            FoamX->SetkDim(DVCSmodel::GetDim());       // No. of dimensions, obligatory!
            FoamX->SetnCells(opt6);  // No. of cells, default=2000
            FoamX->SetnSampl(opt7);   // No. of MC events in the cell MC exploration d=200
            FoamX->SetnBin(8);          // No. of bins in edge-histogram in cell exploration d=8
            FoamX->SetOptRej(1);   // Wted events for OptRej=0; wt=1 for OptRej=1 (default)
            FoamX->SetOptDrive(2);   // Maximum weight reduction, =1 for variance reduction d=2
            FoamX->SetEvPerBin(25);   // Maximum number of the effective wt=1 events/bin
            FoamX->SetMaxWtRej(1.1);   // Maximum weight used to get w=1 MC events d=1.1
            FoamX->SetRhoInt(DVCSmodel::Model);  // Set 4-dim distribution, included below
            FoamX->Initialize();

            // Write the initialize info
            delete GridFile;
            GridFile = new TFile(RooFileName.c_str(),"NEW");
            if ( GridFile->IsOpen() ) printf("Grid file is created.\n");
            FoamX->Write("FoamData");
        }
        FoamX->CheckAll(1);
        cout << "Model was called " << M->GetCount() << " times" << endl;

        cout << "Event generation starts" << endl;
        for(long loop=0; loop<opt5; loop++){
            FoamX->MakeEvent();           // generate MC event
            FoamX->GetMCvect(MCvect);     // get generated vector
            M->Transform(MCvect);         // Calculate experimental values from MCV
            output->Fill(MCvect[0], MCvect[1], MCvect[2],  MCvect[3]);
            if (loop%(opt5/20)==0) cout << "Events generated " << loop << endl;
        }
        cout << "Event generation ends" << endl;

        cout << "Model was called " << M->GetCount() << " times" << endl;

        double IntNorm, Errel, MCresult, MCerror;
        FoamX->Finalize(IntNorm, Errel);
	FoamX->GetIntegMC( MCresult, MCerror);
        double XSecFactor = (opt10 - opt9) * (opt12 - opt11) * opt16 * 360 * 6.2832;
        cout << "Cross section : " << IntNorm * XSecFactor << " nb" << endl;
        cout << "Precision : " << Errel * XSecFactor << " nb" << endl;
        cout << "Events produced correspond to " << opt5 / IntNorm / XSecFactor << " nb^-1" << endl;
	cout << "MC integration: " << MCresult << ", error: " << Errel << endl;


        delete GridFile;
        delete []MCvect,PseRan;
        delete output;
        delete M;
        delete FoamX;

    }
    else cout << "Unable to open file"; 

    return 0;
}


