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

void Fitter(TH2D* h1, TH2D* h2, int YPlo, int YPhi,int color, int boolsame, int name)
{

        TH1::SetDefaultSumw2(kTRUE);
        TH2::SetDefaultSumw2(kTRUE);

        TH1D *histfit1 = (TH1D*) h1->ProjectionY("",YPlo,YPhi)->Clone();
        TH1D *histfit2 = (TH1D*) h2->ProjectionY("",YPlo,YPhi)->Clone();
        histfit1->Divide(histfit2);
        histfit1->Scale(h2->GetBinContent(h2->FindBin(0,0)));
        histfit1->SetMarkerStyle(20);
        std::string function = "[0]/(TMath::Pi()*2)*(1+2*([1]*TMath::Cos(x)+[2]*TMath::Cos(2*x)+[3]*TMath::Cos(3*x)+[4]*TMath::Cos(4*x)+[5]*TMath::Cos(5*x)))";
        TF1 func1(Form("deltaPhi_%d", name), function.c_str(), -0.5*TMath::Pi(), 1.5*TMath::Pi());
        func1.SetParameter(0, histfit1->GetMaximum());
        func1.SetParameter(1, 0.1);
        func1.SetParameter(2, 0.1);
        func1.SetParameter(3, 0.1);
        func1.SetParameter(4, 0.1);
        func1.SetParameter(5, 0.1);
        histfit1->Fit(&func1, "q 0");
        histfit1->Fit(&func1, "q 0");
        histfit1->Fit(&func1, "m q 0");
        histfit1->Fit(&func1, "m q 0");
        histfit1->Fit(&func1, "m q E 0");
        //histfit1->Fit(&func1, "m S E q");
        histfit1->Fit(&func1, "m E q 0");
        histfit1->SetStats(kFALSE);
        //histfit1->GetFunction("deltaPhi1")->SetLineColor(color);
        histfit1->SetMarkerColor(color);
        histfit1->SetMarkerStyle(33);
        histfit1->SetMarkerSize(1);
        if(boolsame == 1){
        histfit1->DrawCopy("PE");
        }
        if(boolsame == 2){
        histfit1->DrawCopy("PE SAME");
        }
        func1.SetLineColor(color);
        func1.DrawCopy("SAME");
        delete histfit1;
        delete histfit2;
}



void CentrCompareRG()
{
    TFile *f1 = new TFile("P470_Primary_Reco.root");
    TFile *f2 = new TFile("P470_Primary_Gen.root");

    TH1::SetDefaultSumw2(kTRUE);
    TH2::SetDefaultSumw2(kTRUE);

    const int trackbin 	= 3;
    const int trackbinC = 3;
    const int ptbin 	= 4;
    //const int   trackbinbounds[trackbin]    = {85,86,87,88,89,90,91,92};
    const int   trackbinbounds[trackbin]    = {85,87,89};
    const float ptbinbounds_lo[ptbin]       = { 4,  5,  6, 10};
    const float ptbinbounds_hi[ptbin]       = {30, 30, 30, 30};

    int ptwant 	= 2;//manually change through 13... or create a for loop that does this but then you will get a 78 canvases.

    int pt_lo 	= ptbinbounds_lo[ptwant-1];
    int pt_hi 	= ptbinbounds_hi[ptwant-1];

    int Xscale 	= 1600;
    int Yscale 	= 800;

    //Defining Canvas Arrays
    //TCanvas* c1 = new TCanvas("c1", "c1", Xscale, Yscale);
    //TCanvas* c2 = new TCanvas("c2", "c2", Xscale, Yscale);

    TCanvas* cA[trackbinC];
    for(int i = 0; i < trackbinC; i++){
        cA[i] = new TCanvas(Form("cA_%d",i),Form("Central Bin Comparison Between Reco and Gen_%d",i), Xscale, Yscale);
    }


    //Defining Jet count for normalization
    TH1D* hJ1_f1;
    hJ1_f1 = (TH1D*)f1->Get("hJet_Pass")->Clone();

    TH1D* hJ1_f2;
    hJ1_f2 = (TH1D*)f2->Get("hJet_Pass")->Clone();
    //Main Loops

    float BW1 = 0.3;

    for(int wtrk = 0; wtrk < trackbinC; wtrk++){
        cA[wtrk]->Divide(3,1);

        TH2D* hSC1;
        hSC1 = (TH2D*)f1->Get(Form("hSigC_%d_%d_to_%d",trackbinbounds[wtrk], pt_lo, pt_hi))->Clone();
        hSC1->Scale(1/(hJ1_f1->GetBinContent(wtrk+1)));
        hSC1->Scale(1./(BW1));

        TH2D* hBC1;
        hBC1 = (TH2D*)f1->Get(Form("hBckC_%d_%d_to_%d",trackbinbounds[wtrk], pt_lo, pt_hi))->Clone();

        TH2D* hSC2;
        hSC2 = (TH2D*)f2->Get(Form("hSigC_%d_%d_to_%d",trackbinbounds[wtrk], pt_lo, pt_hi))->Clone();
        hSC2->Scale(1/(hJ1_f2->GetBinContent(wtrk+1)));
        hSC2->Scale(1./(BW1));

        TH2D* hBC2;
        hBC2 = (TH2D*)f2->Get(Form("hBckC_%d_%d_to_%d",trackbinbounds[wtrk], pt_lo, pt_hi))->Clone();


        cA[wtrk]->cd(1);
        Fitter(hSC1,hBC1,27,35,3, 1, 1);       
        Fitter(hSC2,hBC2,27,35,2, 2, 2);       

        hSC1->Scale(hBC1->GetBinContent(hBC1->FindBin(0,0)));
        hSC2->Scale(hBC2->GetBinContent(hBC2->FindBin(0,0)));

        cA[wtrk]->cd(2);
        hSC1->Divide(hBC1);
        hSC1->Draw("FB BB SURF1");

        cA[wtrk]->cd(3);
        hSC2->Divide(hBC2);
        hSC2->Draw("FB BB SURF1");

   }
}

