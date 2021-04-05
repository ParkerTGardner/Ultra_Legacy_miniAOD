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


void OO()
{
    TFile *f = new TFile("ALL_trk4_pt3.root");

    TH1::SetDefaultSumw2(kFALSE);
    TH2::SetDefaultSumw2(kFALSE);


    int Xscale = 1900;
    int Yscale = 1000;

    TCanvas* c3          = new TCanvas("c3","c3", Xscale, Yscale);
    //TCanvas* c4          = new TCanvas("c4","c4", Xscale, Yscale);
    //c4->Divide(ptbin,trackbin);
    //TCanvas* c5          = new TCanvas("c5","c5", Xscale, Yscale);
    //c5->Divide(ptbin,trackbin);
    //TCanvas* c6          = new TCanvas("c6","c6", Xscale, Yscale);
    //c6->Divide(ptbin,trackbin);

    //c6->SetLeftMargin(0.2);
    //c6->SetTheta(60.839);
    //c6->SetPhi(38.0172);

    TH1D* hJ;
    hJ = (TH1D*)f->Get("hPass_Trig_Jet")->Clone();

    for(int i = 0; i < trackbin+2; i++){
    int jets = hJ->GetBinContent(i);
    cout << "for i = " << i << " the fill is " << jets << endl;
    }
//int wppt =1;
//int wtrk =1;

            TH2D* h1;
            h1 = (TH2D*)f->Get(Form("hSig_%d_%d",wtrk,wppt))->Clone();
            h1->SetStats(kFALSE);
            h1->Scale(1/(hJ->GetBinContent(wtrk)));
            h1->Scale(1./(           h1->GetXaxis()->GetBinWidth(3)  ));

            TH2D* h2;
            h2 = (TH2D*)f->Get(Form("hBck_%d_%d",wtrk,wppt))->Clone();
            h2->SetStats(kFALSE);
            h2->Scale(1/(hJ->GetBinContent(wtrk)));
            h2->Scale(1./(           h2->GetXaxis()->GetBinWidth(3)  ));


            cout << "wtrk is " << wtrk << "  and wppt is " << wppt << "  and sig entries is " << h1->GetEntries() << "  and bck entries is " << h2->GetEntries() << endl;
            int YPmax = h1->GetNbinsX();
            int YPlow = floor(0.75*YPmax);
            cout << YPmax << " " << YPlow << endl;

            TH1D *histfit1 = (TH1D*) h1->ProjectionY("",YPlow,YPmax)->Clone();
            TH1D *histfit2 = (TH1D*) h2->ProjectionY("",YPlow,YPmax)->Clone();


            c3->cd(where);

            histfit1->Divide(histfit2);
            histfit1->Scale(h2->GetMaximum());
            histfit1->SetMarkerStyle(20);

            std::string function = "[0]/(TMath::Pi()*2)*(1+2*([1]*TMath::Cos(x)+[2]*TMath::Cos(2*x)+[3]*TMath::Cos(3*x)+[4]*TMath::Cos(4*x)+[5]*TMath::Cos(5*x)   ))";

            TF1 func1("deltaPhi1", function.c_str(), -0.5*3.14159, 1.5*3.14159);
            func1.SetParameter(0, histfit1->GetMaximum());
            func1.SetParameter(1, 0.1);
            func1.SetParameter(2, 0.1);
            func1.SetParameter(3, 0.1);
            func1.SetParameter(4, 0.1);
            func1.SetParameter(5, 0.1);
            func1.SetLineColor(1);


            histfit1->Fit(&func1, "q");
            histfit1->Fit(&func1, "q");
            histfit1->Fit(&func1, "m q");
            histfit1->Fit(&func1, "m q");
            histfit1->Fit(&func1, "m q E");
            auto fitResult1 = histfit1->Fit(&func1, "m S E q");

            histfit1->SetStats(kFALSE);
            histfit1->DrawCopy("PE");

            //c4->cd(where);

            TF1 *fv1 = new TF1("fv1","[0]/(TMath::Pi()*2)*(1+   2*([1]*TMath::Cos(x)))"  ,-0.5*3.14157, 1.5*3.14157);
            TF1 *fv2 = new TF1("fv2","[0]/(TMath::Pi()*2)*(1+   2*([1]*TMath::Cos(2*x)))",-0.5*3.14157, 1.5*3.14157);
            TF1 *fv3 = new TF1("fv3","[0]/(TMath::Pi()*2)*(1+   2*([1]*TMath::Cos(3*x)))",-0.5*3.14157, 1.5*3.14157);

            fv1->SetParameter(0,func1.GetParameter(0));
            fv2->SetParameter(0,func1.GetParameter(0));
            fv3->SetParameter(0,func1.GetParameter(0));

            fv1->SetParameter(1,func1.GetParameter(1));
            fv2->SetParameter(1,func1.GetParameter(2));
            fv3->SetParameter(1,func1.GetParameter(3));

            fv1->SetLineStyle(7);
            fv2->SetLineStyle(7);
            fv3->SetLineStyle(7);
            fv1->SetLineColor(2);
            fv2->SetLineColor(3);
            fv3->SetLineColor(4);

            fv1->DrawCopy("SAME");
            fv2->DrawCopy("SAME");
            fv3->DrawCopy("SAME");

            //c5->cd(where);
/*
            h1->Divide(h2);
            h1->Scale(h2->GetMaximum());
            h1->SetStats(kFALSE);
            h1->DrawCopy("COLZ");

            int   minBin = h1->GetMinimumBin();
            float min    = h1->GetBinContent(minBin);

            int   maxBin = h1->GetMaximumBin();
            float max    = h1->GetBinContent(maxBin);

            if(wppt > 1){
            h1->GetXaxis()->SetRangeUser(-3.5,3.5);
            h1->GetZaxis()->SetRangeUser(min,(min + (max/6)));
            }
            if(wppt == 1){
            h1->GetZaxis()->SetRangeUser(min,(min + (max/4)));
            }

            c6->cd(where);
            c6->SetLeftMargin(0.2);
            c6->SetTheta(60.839);
            c6->SetPhi(38.0172);
            h1->Draw("SURF1 fb");
            c6->SetLeftMargin(0.2);
            c6->SetTheta(60.839);
            c6->SetPhi(38.0172);
*/
            //delete fv1;
            //delete fv2;
            //delete fv3;
            //delete histfit1;
            //delete histfit2;
            //delete h2;
            //delete h1;

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
