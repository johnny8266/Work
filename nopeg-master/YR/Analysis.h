//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Sep  1 12:10:26 2020 by ROOT version 6.22/00
// from TTree TOPEG/Tree with simulated events from TOPEG
// found on file: outputEIChigha.root
//////////////////////////////////////////////////////////

#ifndef Analysis_h
#define Analysis_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class Analysis {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           ievent;
   Float_t         ElBeam;
   Int_t           PIDlBeam;
   Float_t         EhBeam;
   Int_t           PIDhBeam;
   Float_t         Q2;
   Float_t         W;
   Float_t         Gamnu;
   Float_t         Xbj;
   Float_t         y;
   Float_t         t;
   Float_t         phih;
   Int_t           Nb_part;
   Int_t           part_id[3];   //[Nb_part]
   Float_t         part_px[3];   //[Nb_part]
   Float_t         part_py[3];   //[Nb_part]
   Float_t         part_pz[3];   //[Nb_part]
   Float_t         part_e[3];   //[Nb_part]

   // List of branches
   TBranch        *b_ievent;   //!
   TBranch        *b_ElBeam;   //!
   TBranch        *b_PIDlBeam;   //!
   TBranch        *b_EhBeam;   //!
   TBranch        *b_PIDhBeam;   //!
   TBranch        *b_Q2;   //!
   TBranch        *b_W;   //!
   TBranch        *b_Gamnu;   //!
   TBranch        *b_Xbj;   //!
   TBranch        *b_y;   //!
   TBranch        *b_t;   //!
   TBranch        *b_phih;   //!
   TBranch        *b_Nb_part;   //!
   TBranch        *b_part_id;   //!
   TBranch        *b_part_px;   //!
   TBranch        *b_part_py;   //!
   TBranch        *b_part_pz;   //!
   TBranch        *b_part_e;   //!

   Analysis(TTree *tree=0);
   virtual ~Analysis();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual void     GenTxtFile(const char *);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef Analysis_cxx
Analysis::Analysis(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("outputEIChigha.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("outputEIChigha.root");
      }
      f->GetObject("TOPEG",tree);

   }
   Init(tree);
}

Analysis::~Analysis()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t Analysis::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t Analysis::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void Analysis::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("ievent", &ievent, &b_ievent);
   fChain->SetBranchAddress("ElBeam", &ElBeam, &b_ElBeam);
   fChain->SetBranchAddress("PIDlBeam", &PIDlBeam, &b_PIDlBeam);
   fChain->SetBranchAddress("EhBeam", &EhBeam, &b_EhBeam);
   fChain->SetBranchAddress("PIDhBeam", &PIDhBeam, &b_PIDhBeam);
   fChain->SetBranchAddress("Q2", &Q2, &b_Q2);
   fChain->SetBranchAddress("W", &W, &b_W);
   fChain->SetBranchAddress("Gamnu", &Gamnu, &b_Gamnu);
   fChain->SetBranchAddress("Xbj", &Xbj, &b_Xbj);
   fChain->SetBranchAddress("y", &y, &b_y);
   fChain->SetBranchAddress("t", &t, &b_t);
   fChain->SetBranchAddress("phih", &phih, &b_phih);
   fChain->SetBranchAddress("Nb_part", &Nb_part, &b_Nb_part);
   fChain->SetBranchAddress("part_id", part_id, &b_part_id);
   fChain->SetBranchAddress("part_px", part_px, &b_part_px);
   fChain->SetBranchAddress("part_py", part_py, &b_part_py);
   fChain->SetBranchAddress("part_pz", part_pz, &b_part_pz);
   fChain->SetBranchAddress("part_e", part_e, &b_part_e);
   Notify();
}

Bool_t Analysis::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void Analysis::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t Analysis::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef Analysis_cxx
