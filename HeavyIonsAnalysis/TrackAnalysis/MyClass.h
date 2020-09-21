//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Jul 17 21:43:29 2020 by ROOT version 6.12/07
// from TTree trackTree/v1
// found on file: outputMiniAOD.root
//////////////////////////////////////////////////////////

#ifndef MyClass_h
#define MyClass_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
#include "vector"
#include "vector"
#include "vector"
#include "vector"
#include "vector"
#include "vector"

class MyClass {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           nRun;
   Int_t           nEv;
   Int_t           nLumi;
   vector<float>   *xVtx;
   vector<int>     *nTracksVtx;
   vector<float>   *ptSumVtx;
   vector<float>   *trkPt;
   vector<float>   *trkEta;
   vector<float>   *trkPhi;
   vector<char>    *trkCharge;
   vector<int>     *trkPDFId;
   vector<char>    *trkNHits;
   vector<bool>    *highPurity;
   vector<int>     *trkAssociatedVtxIndx;
   vector<int>     *jetNumDaughters;
   vector<float>   *jetEta;
   vector<float>   *jetPt;
   vector<float>   *jetPhi;
   vector<float>   *jetTheta;
   vector<float>   *jetMass;
   vector<int>     *muonMultiplicity;
   vector<int>     *chargedMultiplicity;
   Int_t           jetN;
   vector<float>   *dau_pt_sum;
   vector<vector<int> > *dau_chg;
   vector<vector<int> > *dau_pid;
   vector<vector<unsigned int> > *dau_vref;
   vector<vector<float> > *dau_pt;
   vector<vector<float> > *dau_eta;
   vector<vector<float> > *dau_phi;
   vector<vector<float> > *dau_theta;
   vector<vector<float> > *dau_vz;
   vector<vector<float> > *dau_vy;
   vector<vector<float> > *dau_vx;
   vector<vector<float> > *dau_vrefz;
   vector<vector<float> > *dau_vrefy;
   vector<vector<float> > *dau_vrefx;
   vector<vector<float> > *dau_vp_difZ;
   vector<vector<float> > *dau_vp_difY;
   vector<vector<float> > *dau_vp_difX;

   // List of branches
   //TBranch        *b_nRun;   //!
   TBranch        *b_nRun;   //!
   TBranch        *b_nLumi;   //!
   TBranch        *b_xVtx;   //!
   TBranch        *b_nTracksVtx;   //!
   TBranch        *b_ptSumVtx;   //!
   TBranch        *b_trkPt;   //!
   TBranch        *b_trkEta;   //!
   TBranch        *b_trkPhi;   //!
   TBranch        *b_trkCharge;   //!
   TBranch        *b_trkPDFId;   //!
   TBranch        *b_trkNHits;   //!
   TBranch        *b_highPurity;   //!
   TBranch        *b_trkAssociatedVtxIndx;   //!
   TBranch        *b_jetNumDaughters;   //!
   TBranch        *b_jetEta;   //!
   TBranch        *b_jetPt;   //!
   TBranch        *b_jetPhi;   //!
   TBranch        *b_jetTheta;   //!
   TBranch        *b_jetMass;   //!
   TBranch        *b_muonMultiplicity;   //!
   TBranch        *b_chargedMultiplicity;   //!
   TBranch        *b_jetN;   //!
   TBranch        *b_dau_pt_sum;   //!
   TBranch        *b_dau_chg;   //!
   TBranch        *b_dau_pid;   //!
   TBranch        *b_dau_vref;   //!
   TBranch        *b_dau_pt;   //!
   TBranch        *b_dau_eta;   //!
   TBranch        *b_dau_phi;   //!
   TBranch        *b_dau_theta;   //!
   TBranch        *b_dau_vz;   //!
   TBranch        *b_dau_vy;   //!
   TBranch        *b_dau_vx;   //!
   TBranch        *b_dau_vrefz;   //!
   TBranch        *b_dau_vrefy;   //!
   TBranch        *b_dau_vrefx;   //!
   TBranch        *b_dau_vp_difZ;   //!
   TBranch        *b_dau_vp_difY;   //!
   TBranch        *b_dau_vp_difX;   //!

   MyClass(TTree *tree=0);
   virtual ~MyClass();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef MyClass_cxx
MyClass::MyClass(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("total_output_ak8.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("total_output_ak8.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("total_output_ak8.root:/analyzer");
      dir->GetObject("trackTree",tree);

   }
   Init(tree);
}

MyClass::~MyClass()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t MyClass::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t MyClass::LoadTree(Long64_t entry)
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

void MyClass::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   xVtx = 0;
   nTracksVtx = 0;
   ptSumVtx = 0;
   trkPt = 0;
   trkEta = 0;
   trkPhi = 0;
   trkCharge = 0;
   trkPDFId = 0;
   trkNHits = 0;
   highPurity = 0;
   trkAssociatedVtxIndx = 0;
   jetNumDaughters = 0;
   jetEta = 0;
   jetPt = 0;
   jetPhi = 0;
   jetTheta = 0;
   jetMass = 0;
   muonMultiplicity = 0;
   chargedMultiplicity = 0;
   dau_pt_sum = 0;
   dau_chg = 0;
   dau_pid = 0;
   dau_vref = 0;
   dau_pt = 0;
   dau_eta = 0;
   dau_phi = 0;
   dau_theta = 0;
   dau_vz = 0;
   dau_vy = 0;
   dau_vx = 0;
   dau_vrefz = 0;
   dau_vrefy = 0;
   dau_vrefx = 0;
   dau_vp_difZ = 0;
   dau_vp_difY = 0;
   dau_vp_difX = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("nRun", &nRun, &b_nRun);
   fChain->SetBranchAddress("nEv", &nEv, &b_nRun);
   fChain->SetBranchAddress("nLumi", &nLumi, &b_nLumi);
   fChain->SetBranchAddress("xVtx", &xVtx, &b_xVtx);
   fChain->SetBranchAddress("nTracksVtx", &nTracksVtx, &b_nTracksVtx);
   fChain->SetBranchAddress("ptSumVtx", &ptSumVtx, &b_ptSumVtx);
   fChain->SetBranchAddress("trkPt", &trkPt, &b_trkPt);
   fChain->SetBranchAddress("trkEta", &trkEta, &b_trkEta);
   fChain->SetBranchAddress("trkPhi", &trkPhi, &b_trkPhi);
   fChain->SetBranchAddress("trkCharge", &trkCharge, &b_trkCharge);
   fChain->SetBranchAddress("trkPDFId", &trkPDFId, &b_trkPDFId);
   fChain->SetBranchAddress("trkNHits", &trkNHits, &b_trkNHits);
   fChain->SetBranchAddress("highPurity", &highPurity, &b_highPurity);
   fChain->SetBranchAddress("trkAssociatedVtxIndx", &trkAssociatedVtxIndx, &b_trkAssociatedVtxIndx);
   fChain->SetBranchAddress("jetNumDaughters", &jetNumDaughters, &b_jetNumDaughters);
   fChain->SetBranchAddress("jetEta", &jetEta, &b_jetEta);
   fChain->SetBranchAddress("jetPt", &jetPt, &b_jetPt);
   fChain->SetBranchAddress("jetPhi", &jetPhi, &b_jetPhi);
   fChain->SetBranchAddress("jetTheta", &jetTheta, &b_jetTheta);
   fChain->SetBranchAddress("jetMass", &jetMass, &b_jetMass);
   fChain->SetBranchAddress("muonMultiplicity", &muonMultiplicity, &b_muonMultiplicity);
   fChain->SetBranchAddress("chargedMultiplicity", &chargedMultiplicity, &b_chargedMultiplicity);
   fChain->SetBranchAddress("jetN", &jetN, &b_jetN);
   fChain->SetBranchAddress("dau_pt_sum", &dau_pt_sum, &b_dau_pt_sum);
   fChain->SetBranchAddress("dau_chg", &dau_chg, &b_dau_chg);
   fChain->SetBranchAddress("dau_pid", &dau_pid, &b_dau_pid);
   fChain->SetBranchAddress("dau_vref", &dau_vref, &b_dau_vref);
   fChain->SetBranchAddress("dau_pt", &dau_pt, &b_dau_pt);
   fChain->SetBranchAddress("dau_eta", &dau_eta, &b_dau_eta);
   fChain->SetBranchAddress("dau_phi", &dau_phi, &b_dau_phi);
   fChain->SetBranchAddress("dau_theta", &dau_theta, &b_dau_theta);
   fChain->SetBranchAddress("dau_vz", &dau_vz, &b_dau_vz);
   fChain->SetBranchAddress("dau_vy", &dau_vy, &b_dau_vy);
   fChain->SetBranchAddress("dau_vx", &dau_vx, &b_dau_vx);
   fChain->SetBranchAddress("dau_vrefz", &dau_vrefz, &b_dau_vrefz);
   fChain->SetBranchAddress("dau_vrefy", &dau_vrefy, &b_dau_vrefy);
   fChain->SetBranchAddress("dau_vrefx", &dau_vrefx, &b_dau_vrefx);
   fChain->SetBranchAddress("dau_vp_difZ", &dau_vp_difZ, &b_dau_vp_difZ);
   fChain->SetBranchAddress("dau_vp_difY", &dau_vp_difY, &b_dau_vp_difY);
   fChain->SetBranchAddress("dau_vp_difX", &dau_vp_difX, &b_dau_vp_difX);
   Notify();
}

Bool_t MyClass::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void MyClass::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t MyClass::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef MyClass_cxx
