//blh
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TLegend.h"
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
#include "TLatex.h"
//#include "coordinateTools.h"

using TMath::ATan;
using TMath::Exp;


double VOut(TH2D* h1, TH2D* h2, int YPlo, int YPhi,int vn)
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
    histfit1->Fit(&func1, "m E q 0");
    return func1.GetParameter(vn);
}
//GetParError
double VEOut(TH2D* h1, TH2D* h2, int YPlo, int YPhi, int vn)
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
    histfit1->Fit(&func1, "m E q 0");
    return func1.GetParError(vn);
}

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



void paperDraw2()
{
    TFile *f1 = new TFile("All_P_LM.root");
    TFile *f2 = new TFile("All_P_HM.root");
    TH1::SetDefaultSumw2(kTRUE);
    TH2::SetDefaultSumw2(kTRUE);

    const int Ltrackbin 	= 6;
    const int Htrackbin 	= 4;

    const int ptbin 	= 1;

    const int   LMtrackbinbounds[Ltrackbin]         = { 0,10,20,30,40,50};
    const int   LMtrackbinboundsUpper[Ltrackbin]    = {10,20,30,40,50,60};

    const int   HMtrackbinbounds[Htrackbin]         = {60,70,78,86};
    const int   HMtrackbinboundsUpper[Htrackbin]    = {70,78,86,1000};

    const float ptbinbounds_lo[ptbin]       = {3};
    const float ptbinbounds_hi[ptbin]       = {30};

    int ptwant 	= 1;

    int pt_lo 	= ptbinbounds_lo[ptwant-1];
    int pt_hi 	= ptbinbounds_hi[ptwant-1];

    int Xscale 	= 800;
    int Yscale 	= 800;

    TCanvas* c1 = new TCanvas("c1","c1", Xscale, Yscale);
    //TCanvas* c2 = new TCanvas("c2","c2", Xscale, Yscale);
    //TCanvas* c3 = new TCanvas("c3","c3", Xscale, Yscale);


    //Defining Jet count for normalization
    TH1D* hJ1_f1;
    hJ1_f1 = (TH1D*)f1->Get("hJet_Pass")->Clone();

    TH1D* hJ1_f2;
    hJ1_f2 = (TH1D*)f2->Get("hJet_Pass")->Clone();
    //Main Loops

    float BW1 = 0.3;

    double Lx_1[Ltrackbin] = {5,15,25,35,45,55};
    double Ly_1[Ltrackbin];
    double Lxe_1[Ltrackbin] = {5,5,5,5,5,5};
    double Lye_1[Ltrackbin];

    double Lx_2[Ltrackbin] = {5,15,25,35,45,55};
    double Ly_2[Ltrackbin];
    double Lxe_2[Ltrackbin] = {5,5,5,5,5,5};
    double Lye_2[Ltrackbin];

    double Lx_3[Ltrackbin] = {5,15,25,35,45,55};
    double Ly_3[Ltrackbin];
    double Lxe_3[Ltrackbin] = {5,5,5,5,5,5};
    double Lye_3[Ltrackbin];

    double Hx_1[Ltrackbin] = {65,75,85,95};
    double Hy_1[Ltrackbin];
    double Hxe_1[Ltrackbin] = {5,4,4,4};
    double Hye_1[Ltrackbin];

    double Hx_2[Ltrackbin] = {65,75,85,95};
    double Hy_2[Ltrackbin];
    double Hxe_2[Ltrackbin] = {5,4,4,4};
    double Hye_2[Ltrackbin];

    double Hx_3[Ltrackbin] = {65,75,85,95};
    double Hy_3[Ltrackbin];
    double Hxe_3[Ltrackbin] = {5,4,4,4};
    double Hye_3[Ltrackbin];

    for(int wtrk = 0; wtrk < Ltrackbin; wtrk++){

        TH2D* h1;
        h1 = (TH2D*)f1->Get(Form("hSigS_%d_%d_to_%d",LMtrackbinbounds[wtrk], pt_lo, pt_hi))->Clone();
        h1->Scale(1/(hJ1_f1->GetBinContent(wtrk+1)));
        h1->Scale(1./(BW1));

        TH2D* h2;
        h2 = (TH2D*)f1->Get(Form("hBckS_%d_%d_to_%d",LMtrackbinbounds[wtrk], pt_lo, pt_hi))->Clone();

        Ly_1[wtrk]  = VOut(h1,h2,28,34, 1);       
        Lye_1[wtrk] = VEOut(h1,h2,28,34, 1);       

        Ly_2[wtrk]  = VOut(h1,h2,28,34, 2);       
        Lye_2[wtrk] = VEOut(h1,h2,28,34, 2);       

        Ly_3[wtrk]  = VOut(h1,h2,28,34, 3);       
        Lye_3[wtrk] = VEOut(h1,h2,28,34, 3);       

    }

    for(int wtrk = 0; wtrk < Htrackbin; wtrk++){

        TH2D* h1;
        h1 = (TH2D*)f2->Get(Form("hSigS_%d_to_%d_and_%d_to_%d",HMtrackbinbounds[wtrk], HMtrackbinboundsUpper[wtrk], pt_lo, pt_hi))->Clone();
        h1->Scale(1/(hJ1_f2->GetBinContent(wtrk+1)));
        h1->Scale(1./(BW1));

        TH2D* h2;
        h2 = (TH2D*)f2->Get(Form("hBckS_%d_to_%d_and_%d_to_%d",HMtrackbinbounds[wtrk], HMtrackbinboundsUpper[wtrk], pt_lo, pt_hi))->Clone();

        Hy_1[wtrk]  = VOut(h1,h2,28,34, 1);       
        Hye_1[wtrk] = VEOut(h1,h2,28,34, 1);       

        Hy_2[wtrk]  = VOut(h1,h2,28,34, 2);       
        Hye_2[wtrk] = VEOut(h1,h2,28,34, 2);       

        Hy_3[wtrk]  = VOut(h1,h2,28,34, 3);       
        Hye_3[wtrk] = VEOut(h1,h2,28,34, 3);       

    }
    c1->cd();
    auto grL1 = new TGraphErrors(Ltrackbin,Lx_1,Ly_1,Lxe_1,Lye_1);
    grL1->SetTitle("v_{1}{2} ");
    grL1->SetMarkerColor(2);
    grL1->SetMarkerStyle(24);
    grL1->Draw("P");
    grL1->GetXaxis()->SetRangeUser(0,100);
    //grL1->Draw("ALP");
    //c1->Update();


    //c2->cd();
    auto grL2 = new TGraphErrors(Ltrackbin,Lx_2,Ly_2,Lxe_2,Lye_2);
    grL2->SetTitle("v_{2}{2} ");
    grL2->SetMarkerColor(4);
    grL2->SetMarkerStyle(24);
    //grL2->Draw("P");
    //grL2->GetXaxis()->SetRangeUser(0,100);
    //grL2->Draw("ALP");
    //c2->Update();

    //c3->cd();
    auto grL3 = new TGraphErrors(Ltrackbin,Lx_3,Ly_3,Lxe_3,Lye_3);
    grL3->SetTitle("v_{3}{2} ");
    grL3->SetMarkerColor(3);
    grL3->SetMarkerStyle(24);
    //grL3->Draw("P");
    //grL3->GetXaxis()->SetRangeUser(0,100);
    //grL3->Draw("ALP");
    //c3->Update();

    //c1->cd();
    auto grH1 = new TGraphErrors(Htrackbin,Hx_1,Hy_1,Hxe_1,Hye_1);
    grH1->SetTitle("v_{1}{2} (filtered)");
    grH1->SetMarkerColor(2);
    grH1->SetMarkerStyle(20);
    //grH1->Draw("P SAME");

    //c2->cd();
    auto grH2 = new TGraphErrors(Htrackbin,Hx_2,Hy_2,Hxe_2,Hye_2);
    grH2->SetTitle("v_{2}{2} (filtered)");
    grH2->SetMarkerColor(4);
    grH2->SetMarkerStyle(20);
    //grH2->Draw("P SAME");

    //c3->cd();
    auto grH3 = new TGraphErrors(Htrackbin,Hx_3,Hy_3,Hxe_3,Hye_3);
    grH3->SetTitle("v_{3}{2} (filtered)");
    grH3->SetMarkerColor(3);
    grH3->SetMarkerStyle(20);
    //grH3->Draw("P SAME");

    auto mg1  = new TMultiGraph();
    //auto mg2  = new TMultiGraph();
    //auto mg3  = new TMultiGraph();


    grL1->SetLineWidth(1);
    grL2->SetLineWidth(1);
    grL3->SetLineWidth(1);

    grH1->SetLineWidth(1);
    grH2->SetLineWidth(1);
    grH3->SetLineWidth(1);

    grL1->SetLineColor(2);
    grL2->SetLineColor(4);
    grL3->SetLineColor(3);

    grH1->SetLineColor(2);
    grH2->SetLineColor(4);
    grH3->SetLineColor(3);

    mg1->Add(grL1,"PC");
    mg1->Add(grH1,"Pc");
    mg1->Add(grL2,"PC");
    mg1->Add(grH2,"PC");
    mg1->Add(grL3,"PC");
    mg1->Add(grH3,"PC");
    c1->cd();
    mg1->Draw("AP");
    c1->BuildLegend(0.64, 0.16, 0.9, 0.34);

//auto legend = new TLegend(0.1,0.7,0.48,0.9);
//   legend->SetHeader("The Legend Title","C"); // option "C" allows to center the header
//   legend->AddEntry("grL1","Graph with error bars","lep");
//   legend->Draw();





   //c2->cd();
    //mg2->Draw("AP");
    //c3->cd();
    //mg3->Draw("AP");


}

