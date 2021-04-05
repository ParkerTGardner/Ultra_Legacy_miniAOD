#include "TH2D.h"
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

using TMath::ATan;
using TMath::Exp;

void WORK()
{
    TFile *f = new TFile("all.root");

    TH1::SetDefaultSumw2(kFALSE);
    TH2::SetDefaultSumw2(kFALSE);

    int scale = 1000;

    TCanvas* cS1          = new TCanvas("cS1","cS1" , scale, scale);
    TCanvas* cS2          = new TCanvas("cS2","cS2" , scale, scale);
    TCanvas* cS3          = new TCanvas("cS3","cS3" , scale, scale);

    TCanvas* c1          = new TCanvas("c1","c1" , scale, scale);

    TH2D* hs1;
    TH2D* hb1;

    hs1 = (TH2D*)f->Get("hSig_4_3");
    hb1 = (TH2D*)f->Get("hBck_4_3");

    TH2D* hs2;
    TH2D* hb2;

    hs2 = (TH2D*)f->Get("hSig_4_4");
    hb2 = (TH2D*)f->Get("hBck_4_4");

    TH2D* hs3;
    TH2D* hb3;

    hs3 = (TH2D*)f->Get("hSig_5_3");
    hb3 = (TH2D*)f->Get("hBck_5_3");

    TH2D* hs4;
    TH2D* hb4;

    hs4 = (TH2D*)f->Get("hSig_5_4");
    hb4 = (TH2D*)f->Get("hBck_5_4");

    TH2D* hs5;
    TH2D* hb5;

    hs5 = (TH2D*)f->Get("hSig_6_3");
    hb5 = (TH2D*)f->Get("hBck_6_3");

    TH2D* hs6;
    TH2D* hb6;

    hs6 = (TH2D*)f->Get("hSig_6_4");
    hb6 = (TH2D*)f->Get("hBck_6_4");


    hs1->Add(hs2);
    hs1->Add(hs3);
    hs1->Add(hs4);
    hs1->Add(hs5);
    hs1->Add(hs6);

    hb1->Add(hb2);
    hb1->Add(hb3);
    hb1->Add(hb4);
    hb1->Add(hb5);
    hb1->Add(hb6);

    cout << "sig entries is: " << hs1->GetEntries() << "  and bck entries is: " << hb1->GetEntries() << endl;

    TH1D *histfit1 = (TH1D*) hs1->ProjectionY("",40,54)->Clone();
    TH1D *histfit2 = (TH1D*) hb1->ProjectionY("",40,54)->Clone();

    histfit1->Divide(histfit2);

    std::string function = "[0]/(TMath::Pi()*2)*(1+2*([1]*TMath::Cos(x)+[2]*TMath::Cos(2*x)+[3]*TMath::Cos(3*x)+[4]*TMath::Cos(4*x)+[5]*TMath::Cos(5*x)   ))";

    TF1 func1("deltaPhi1", function.c_str(), -0.5*3.14157, 1.5*3.14157);
    func1.SetParameter(0, histfit1->GetMaximum());
    func1.SetParameter(1, 0.1);
    func1.SetParameter(2, 0.1);
    func1.SetParameter(3, 0.1);
    func1.SetParameter(4, 0.1);
    func1.SetParameter(5, 0.1);

    histfit1->SetMarkerStyle(20);

    histfit1->Fit(&func1, "q");
    histfit1->Fit(&func1, "q");
    histfit1->Fit(&func1, "m q");
    histfit1->Fit(&func1, "m q");
    histfit1->Fit(&func1, "m q E");
    auto fitResult1 = histfit1->Fit(&func1, "m S E q");
    //fitResult1->Print();            
    //histo->DrawCopy("PE");

    TF1 *fv1 = new TF1("fv1","([0]*TMath::Cos(x))"  ,-0.5*3.14157, 1.5*3.14157);
    TF1 *fv2 = new TF1("fv2","([0]*TMath::Cos(2*x))",-0.5*3.14157, 1.5*3.14157);
    TF1 *fv3 = new TF1("fv3","([0]*TMath::Cos(3*x))",-0.5*3.14157, 1.5*3.14157);
    fv1->SetParameter(0,func1.GetParameter(1));
    fv2->SetParameter(0,func1.GetParameter(2));
    fv3->SetParameter(0,func1.GetParameter(3));

    fv1->SetLineColor(1);
    fv2->SetLineColor(2);
    fv3->SetLineColor(3);

    cS1->cd();
    fv1->DrawCopy();
    fv2->DrawCopy("SAME");
    fv3->DrawCopy("SAME");

    cS2->cd();
    histfit1->DrawCopy("PE");

    cS3->cd();
    hs1->Divide(hb1);
    hs1->DrawCopy("COLZ");

    c1->cd();
    hs1->DrawCopy("SURF1 fb");
}


















/*


   TCanvas* c1          = new TCanvas("c1","c1" , 1000, 1000);
   c1->Divide(ptbin,trackbin);
   TCanvas* c2          = new TCanvas("c2","c2" , 1000, 1000);
   c2->Divide(ptbin,trackbin);

//TCanvas* c3          = new TCanvas("c3","c3" , 1000, 1000);
//c3->Divide(ptbin,trackbin);
//TCanvas* c4          = new TCanvas("c4","c4" , 1000, 1000);
//c4->Divide(ptbin,trackbin);


TCanvas* c1Y         = new TCanvas("c1Y","c1Y" , 1000, 1000);
c1Y->Divide(ptbin,trackbin);
//TCanvas* c1Y_sub     = new TCanvas("c1Y_sub","c1Y_sub" , 1000, 1000);
//c1Y_sub->Divide(ptbin,trackbin);
//TCanvas* c1X         = new TCanvas("c1X","c1X" , 1000, 1000);
//c1X->Divide(5,3);
//TCanvas* c1X_sub     = new TCanvas("c1X_sub","c1X_sub" , 1000, 1000);
//c1X_sub->Divide(5,3);

for(int wtrk =1; wtrk <trackbin+1; wtrk++){
for(int wppt =1; wppt <ptbin+1; wppt++){
TH2D* h1 = (TH2D*)f->Get(Form("hSig_%d_%d",wtrk,wppt));
TH2D* h2 = (TH2D*)f->Get(Form("hBck_%d_%d",wtrk,wppt));
//TH2D* h3 = (TH2D*)f->Get(Form("hSig_sub_%d_%d",wtrk,wppt));
//TH2D* h4 = (TH2D*)f->Get(Form("hBck_sub_%d_%d",wtrk,wppt));
//h1->SetDirectory(nullptr);
//h2->SetDirectory(nullptr);
//h3->SetDirectory(nullptr);
//h4->SetDirectory(nullptr);

h1->SetStats(kFALSE);
h2->SetStats(kFALSE);
//h3->SetStats(kFALSE);
//h4->SetStats(kFALSE);

int where = wppt + (wtrk-1)*(ptbin);

//h1->Divide(h2);
//c1->cd(where);
//h1->ProjectionY()->DrawCopy();


c1->cd(where);
h1->DrawCopy("COLZ");

c2->cd(where);
h2->DrawCopy("COLZ");

//h1->Divide(h2);
//c3->cd(where);
//h1->DrawCopy("COLZ");

//h3->Divide(h4);
//c4->cd(where);
//h3->DrawCopy("COLZ");
*/
/*
   c1->cd(where);
   h1->DrawCopy("COLZ");

   c2->cd(where);
   h2->DrawCopy("COLZ");
   */
/*

   TH1D* hh1 = h1->ProjectionY("",1,39);
   TH1D* hh2 = h2->ProjectionY("",1,39);
   hh1->Divide(hh2);
   hh1->SetName(Form("hh1_%d_%d",wtrk,wppt));
   hh1->SetStats(kFALSE);


//h1->Divide(h2);
//h1->SetName(Form("h1_%d_%d",wtrk,wppt));
//h1->SetStats(kFALSE);

//h3->Divide(h4);
//h3->SetName(Form("h3_%d_%d",wtrk,wppt));
//h3->SetStats(kFALSE);

//c1->cd(where);
//h1->DrawCopy("COLZ");

c1Y->cd(where);
//h1->ProjectionY()->DrawCopy();
hh1->DrawCopy();            


//c1Y_sub->cd(where);
//h1->ProjectionY("",27,39)->DrawCopy();

//c1X->cd(where);
//h1->ProjectionX()->DrawCopy();

//c1X_sub->cd(where);
//h1->ProjectionX("",1,15)->DrawCopy();

}
}
}
*/
