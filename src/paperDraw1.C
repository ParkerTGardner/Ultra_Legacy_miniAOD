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

void Fitter(TH2D* h1, TH2D* h2, int YPlo, int YPhi,int color, int boolsame)
{

    TH1::SetDefaultSumw2(kTRUE);
    TH2::SetDefaultSumw2(kTRUE);

    TH1D *histfit1 = (TH1D*) h1->ProjectionY("",YPlo,YPhi)->Clone();
    TH1D *histfit2 = (TH1D*) h2->ProjectionY("",YPlo,YPhi)->Clone();
    histfit1->Divide(histfit2);
    histfit1->Scale(h2->GetMaximum());
    histfit1->SetMarkerStyle(20);
    std::string function = "[0]/(TMath::Pi()*2)*(1+2*([1]*TMath::Cos(x)+[2]*TMath::Cos(2*x)+[3]*TMath::Cos(3*x)+[4]*TMath::Cos(4*x)+[5]*TMath::Cos(5*x)))";
    TF1 func1("deltaPhi1", function.c_str(), -0.5*TMath::Pi(), 1.5*TMath::Pi());
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
}



void paperDraw1()
{
    TFile *f1 = new TFile("All_P_HM.root");
    //TFile *f1 = new TFile("paper/draft/highM/P_highM.root");
//P_highM_correct.C
    TH1::SetDefaultSumw2(kTRUE);
    TH2::SetDefaultSumw2(kTRUE);

    const int trackbin 	= 1;
    const int trackbinC = 1;
    const int ptbin 	= 1;
    const int   trackbinbounds[trackbin]         = {80};
    const int   trackbinboundsUpper[trackbin]    = {1000};
    const float ptbinbounds_lo[ptbin]       = { 3};
    const float ptbinbounds_hi[ptbin]       = {30};

    int ptwant 	= 1;//manually change through 13... or create a for loop that does this but then you will get a 78 canvases.

    int pt_lo 	= ptbinbounds_lo[ptwant-1];
    int pt_hi 	= ptbinbounds_hi[ptwant-1];

    int Xscale 	= 800;
    int Yscale 	= 800;

    //Defining Canvas Arrays
    //TCanvas* c1 = new TCanvas("c1", "c1", Xscale, Yscale);
    //TCanvas* c2 = new TCanvas("c2", "c2", Xscale, Yscale);

    TCanvas* c1 = new TCanvas("c1","c1", Xscale, Yscale);
    TCanvas* c2 = new TCanvas("c2","c2", Xscale, Yscale);
    TCanvas* c3 = new TCanvas("c3_fake","c3_fake", Xscale, Yscale);
    TCanvas* c4 = new TCanvas("c4_fake","c4_fake", Xscale, Yscale);


    //Defining Jet count for normalization
    TH1D* hJ1_f1;
    hJ1_f1 = (TH1D*)f1->Get("hJet_Pass")->Clone();

    //Main Loops

    float BW1 = 0.3;
    float V1  = .05;
    float VV  = V1*V1;

    for(int wtrk = 0; wtrk < 1; wtrk++){

        TH2D* h1;
        //h1 = (TH2D*)f1->Get(Form("hSigS_%d_%d_to_%d",trackbinbounds[wtrk], pt_lo, pt_hi))->Clone();
        h1 = (TH2D*)f1->Get(Form("hSigS_%d_to_%d_and_%d_to_%d",trackbinbounds[wtrk],trackbinboundsUpper[wtrk], pt_lo, pt_hi))->Clone();
        h1->Scale(1/(hJ1_f1->GetBinContent(5)));
        h1->Scale(1./(BW1));

        TH2D* h2;
        h2 = (TH2D*)f1->Get(Form("hBckS_%d_to_%d_and_%d_to_%d",trackbinbounds[wtrk],trackbinboundsUpper[wtrk], pt_lo, pt_hi))->Clone();
        //h2 = (TH2D*)f1->Get(Form("hBckS_%d_%d_to_%d",trackbinbounds[wtrk], pt_lo, pt_hi))->Clone();

        c1->cd();
        Fitter(h1,h2,28,34,2, 1);       

        c2->cd();
        h1->Divide(h2);
        h1->SetStats(kFALSE);
        c2->SetLeftMargin(0.2);
        c2->SetTheta(60.839);
        c2->SetPhi(38.0172); 
        h1->Scale(h2->GetMaximum());
        h1->Draw("FB BB SURF1");
        c2->SetLeftMargin(0.2);
        c2->SetTheta(60.839);
        c2->SetPhi(38.0172); 


        c3->cd();
        TH2D* h3;
        //h3 = (TH2D*)f1->Get(Form("hSigS_%d_%d_to_%d",trackbinbounds[wtrk], pt_lo, pt_hi))->Clone();
        h3 = (TH2D*)f1->Get(Form("hSigS_%d_to_%d_and_%d_to_%d",trackbinbounds[wtrk],trackbinboundsUpper[wtrk], pt_lo, pt_hi))->Clone();
        h3->Scale(1/(hJ1_f1->GetBinContent(5)));
        h3->Scale(1./(BW1));

        TH2D* h4;
        h4 = (TH2D*)f1->Get(Form("hBckS_%d_to_%d_and_%d_to_%d",trackbinbounds[wtrk],trackbinboundsUpper[wtrk], pt_lo, pt_hi))->Clone();
        //h4 = (TH2D*)f1->Get(Form("hBckS_%d_%d_to_%d",trackbinbounds[wtrk], pt_lo, pt_hi))->Clone();

        TH1D *histfit3 = (TH1D*) h3->ProjectionY("",28,34)->Clone();
        TH1D *histfit4 = (TH1D*) h4->ProjectionY("",28,34)->Clone();
        histfit3->Divide(histfit4);
        histfit3->Scale(h4->GetMaximum());

        std::string function2 = "(TMath::Cos(2*x))";
        TF1 efunc("efunc",function2.c_str(),-0.5*3.14159,1.5*3.14159);
        float   projMax = (histfit3->GetMaximum())/(2*3.14159);
        int     loopMax = 400000;

        for (int i = 0; i<loopMax; i++) {
            float YY = efunc.GetRandom();

            if( (YY > -3.14159/4 && YY<3.14159/4) || (YY>3*3.14159/4 && YY<5*3.14159/4)){

                //histfit3->Fill(YY, (projMax*(VV*VV)/(loopMax/41)) );
                histfit3->Fill(YY, ((projMax*2*VV)/(loopMax/41)) );

            } else {

                //histfit3->Fill(YY, -(projMax*(VV*VV)/(loopMax/41)) );
                histfit3->Fill(YY, -((projMax*2*VV)/(loopMax/41)) );
            }
        }

        histfit3->SetMarkerStyle(20);

        std::string function = "[0]/(TMath::Pi()*2)*(1+2*([1]*TMath::Cos(x)+[2]*TMath::Cos(2*x)+[3]*TMath::Cos(3*x)+[4]*TMath::Cos(4*x)+[5]*TMath::Cos(5*x)   ))";

        TF1 func3("deltaPhi3", function.c_str(), -0.5*3.14159, 1.5*3.14159);
        func3.SetParameter(0, histfit3->GetMaximum());
        func3.SetParameter(1, 0.1);
        func3.SetParameter(2, 0.1);
        func3.SetParameter(3, 0.1);
        func3.SetParameter(4, 0.1);
        func3.SetParameter(5, 0.1);

        histfit3->Fit(&func3, "q");
        histfit3->Fit(&func3, "q");
        histfit3->Fit(&func3, "m q");
        histfit3->Fit(&func3, "m q");
        histfit3->Fit(&func3, "m q E");
        auto fitResult3 = histfit3->Fit(&func3, "m S E q");

        histfit3->SetStats(kFALSE);
        histfit3->DrawCopy("PE");

        h3->Divide(h4);
        h3->Scale(h4->GetMaximum());
        h3->SetStats(kFALSE);
        c4->cd();

        for (int j = 0; j<42; j++){
            for (int i = 0; i<loopMax; i++) {
                float YY = efunc.GetRandom();
                if( (YY > -3.14159/4 && YY<3.14159/4) || (YY>3*3.14159/4 && YY<5*3.14159/4)){
                    h3->Fill(-6.15+(j*0.299), YY, ((projMax*2*VV)/(loopMax/41)) );
                } else {
                    h3->Fill(-6.15+(j*0.299), YY, -((projMax*2*VV)/(loopMax/41)) );
                }
            }
        }
        c4->SetLeftMargin(0.2);
        c4->SetTheta(60.839);
        c4->SetPhi(38.0172);
        h3->Draw("SURF1 FB BB");
        c4->SetLeftMargin(0.2);
        c4->SetTheta(60.839);
        c4->SetPhi(38.0172);

    }
}

