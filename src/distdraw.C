//blh
#include "TH2D.h"
#include "TH1D.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TProfile.h"
#include "TGraph.h"
#include <vector>
#include "math.h"
#include "TRandom3.h"
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TCutG.h"
#include <numeric>
//#include "coordinateTools.h"

using TMath::ATan;
using TMath::Exp;

void distdraw()
{
    TFile *f1 = new TFile("SumS678_full_eta.root");
    TFile *f2 = new TFile("SumS678_full_eta.root");
    TFile *f3 = new TFile("SumS678_full_eta.root");
    TH1::SetDefaultSumw2(kFALSE);
    TH2::SetDefaultSumw2(kFALSE);

    const int trackbin 	= 6;
    const int ptbin 	= 15;
    const int   trackbinbounds[trackbin]    = {0, 85,86,87,88,89,90};
    const float ptbinbounds_lo[ptbin]       = {0,   4,  5,  6, 10};
    const float ptbinbounds_hi[ptbin]       = {200, 30, 30, 30, 30};

    int ptwant 	= 9;//manually change through 13... or create a for loop that does this but then you will get a 78 canvases.

    int Xscale 	= 1900;
    int Yscale 	= 1000;

    TCanvas* c1 = new TCanvas("c1", "c1", Xscale, Yscale);
    c1->Divide(2,1);
    //Defining Canvas Arrays
    //TCanvas* cA[trackbin];
    //for(int i = 0; i < trackbin; i++){
    //    cA[i] = new TCanvas(Form("cA_%d",i),Form("cA_%d",i), Xscale, Yscale);
    //    cA[i] ->Divide(2,1);
    //}

    int colors[trackbin] = {1,2,3,4,6,7};


    //Main Loops
    //for(int wtrk = 0; wtrk < 2; wtrk++){
    for(int wtrk = 0; wtrk < trackbin; wtrk++){

        TH2D* h1;
        h1 = (TH2D*)f->Get(Form("hEPDrawPluss_%d_%d",wtrk+1, ptwant))->Clone();

        TH2D* h2;
        h2 = (TH2D*)f->Get(Form("hEPDrawMinus_%d_%d",wtrk+1, ptwant))->Clone();
        
        int maxbinX = h1->GetNbinsX();
        int maxbinY = h1->GetNbinsY();

        
        TH1D *h1DPlusY = (TH1D*) h1->ProjectionY("",1,maxbinY)->Clone();
        TH1D *h1DMinusY = (TH1D*) h2->ProjectionY("",1,maxbinY)->Clone();
        //h1DPlusY->Rebin(6);
        //h1DMinusY->Rebin(6);

        h1DPlusY->Divide(h1DMinusY);

        TH1D *h1DPlusX = (TH1D*) h1->ProjectionX("",1,maxbinX)->Clone();
        TH1D *h1DMinusX = (TH1D*) h2->ProjectionX("",1, maxbinX)->Clone();
        //h1DPlusX->Rebin(6);
        //h1DMinusX->Rebin(6);

        h1DPlusX->Divide(h1DMinusX);


        h1DPlusX->SetMarkerStyle(20);
        h1DPlusX->SetMarkerColor(colors[wtrk]); 

        h1DPlusY->SetMarkerStyle(21);
        h1DPlusY->SetMarkerColor(colors[wtrk]); 

        c1->cd(1);
        h1DPlusX->Draw("SAME");

        c1->cd(2);
        h1DPlusY->Draw("SAME");

    }
}
