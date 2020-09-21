#include "TH2D.h"
#include "TMultiGraph.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TProfile.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include <vector>
#include "math.h"
#include "TRandom3.h"
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TCutG.h"
#include <numeric>
#include "coordinateTools.h"
#include "Advanced_Header_Ndivide.h"

using TMath::ATan;
using TMath::Exp;

void Advanced_Drawer()
{
    TFile *f = new TFile("narroow_ak8_full.root");

    TH1::SetDefaultSumw2(kTRUE);
    TH2::SetDefaultSumw2(kTRUE);

    auto mg1 = new TMultiGraph("mg1","mg1");
    auto mg2 = new TMultiGraph("mg2","mg2");
    auto mg3 = new TMultiGraph("mg3","mg3");

    for(int i = 0; i < trackbin; i++){
        xntrk[i]  =(90/trackbin)*(i+1)-(90/trackbin)/2;
        //xntrke[i] =(80/trackbin)/2;
        xntrke[i] = 0;
    }

    for(int i = 0; i < dEta_N; i++){

        double v1[trackbin];//  = {0};
        double v1e[trackbin];// = {0};
        double v2[trackbin];//  = {0};
        double v2e[trackbin];// = {0};
        double v3[trackbin];//  = {0};
        double v3e[trackbin];// = {0};

        int i_dEta_bin_cut = (int) dEta_bin_cut[i];
        int i_color_eta        = dEta_bin_color[i];

        for(int wppt =1; wppt <ptbin+1; wppt++){
            int i_style_pt = dPt_bin_style[wppt-1];

            for(int wtrk =1; wtrk <trackbin+1; wtrk++){
                int where = wppt + (wtrk-1)*(ptbin);

                TH2D* h1 = (TH2D*)f->Get(Form("hSig_%d_%d",wtrk,wppt));
                TH2D* h2 = (TH2D*)f->Get(Form("hBck_%d_%d",wtrk,wppt));
                
                if (divide_pt == false){
                TH2D* h3 = (TH2D*)f->Get(Form("hSig_%d_2",wtrk));
                TH2D* h4 = (TH2D*)f->Get(Form("hBck_%d_2",wtrk));
                TH2D* h5 = (TH2D*)f->Get(Form("hSig_%d_3",wtrk));
                TH2D* h6 = (TH2D*)f->Get(Form("hBck_%d_3",wtrk));
                h1->Add(h3);
                h1->Add(h5);               
                h2->Add(h4);
                h2->Add(h6);
                }

                h1->Divide(h2);
               
                TH1D *histfit1 = (TH1D*) h1->ProjectionY("",i_dEta_bin_cut,80)->Clone();

                std::string function = "[0]/(TMath::Pi()*2)*(1+2*([1]*TMath::Cos(x)+[2]*TMath::Cos(2*x)+[3]*TMath::Cos(3*x)+[4]*TMath::Cos(4*x)+[5]*TMath::Cos(5*x)   ))";
                TF1 func1("deltaPhi1", function.c_str(), -0.5*3.14157, 1.5*3.14157);
                func1.SetParameter(0, histfit1->GetMaximum());
                func1.SetParameter(1, 0.1);
                func1.SetParameter(2, 0.1);
                func1.SetParameter(3, 0.1);
                func1.SetParameter(4, 0.1);
                func1.SetParameter(5, 0.1);

                histfit1->Fit(&func1, "q");
                histfit1->Fit(&func1, "q");
                histfit1->Fit(&func1, "m q");
                histfit1->Fit(&func1, "m q");
                histfit1->Fit(&func1, "m q E");
                auto fitResult1 = histfit1->Fit(&func1, "m S E q");
 
                TF1 *fv1 = new TF1("fv1","([0]*TMath::Cos(x))"  ,-0.5*3.14157, 1.5*3.14157);
                TF1 *fv2 = new TF1("fv2","([0]*TMath::Cos(2*x))",-0.5*3.14157, 1.5*3.14157);
                TF1 *fv3 = new TF1("fv3","([0]*TMath::Cos(3*x))",-0.5*3.14157, 1.5*3.14157);
                fv1->SetParameter(0,func1.GetParameter(1));
                fv2->SetParameter(0,func1.GetParameter(2));
                fv3->SetParameter(0,func1.GetParameter(3));
                
                v1[wtrk-1]  = {func1.GetParameter(1)};
                v1e[wtrk-1] = {func1.GetParError(1)};
                v2[wtrk-1]  = {func1.GetParameter(2)};
                v2e[wtrk-1] = {func1.GetParError(2)};
                v3[wtrk-1]  = {func1.GetParameter(3)};
                v3e[wtrk-1] = {func1.GetParError(3)};

                delete fv1;
                delete fv2;
                delete fv3;
                delete histfit1;
                delete h2;
                delete h1;
            }
        
            float dEtaUnit = (((float)i_dEta_bin_cut/(float)80)*4);

            auto gre1 = new TGraphErrors(trackbin, xntrk, v1, xntrke, v1e);
            auto gre2 = new TGraphErrors(trackbin, xntrk, v2, xntrke, v2e);
            auto gre3 = new TGraphErrors(trackbin, xntrk, v3, xntrke, v3e);

            gre1->SetMarkerColor(i_color_eta);
            gre1->SetMarkerStyle(i_style_pt);
            gre2->SetMarkerColor(i_color_eta);
            gre2->SetMarkerStyle(i_style_pt);
            gre3->SetMarkerColor(i_color_eta);
            gre3->SetMarkerStyle(i_style_pt);

            gre1->SetMarkerSize(MS);
            gre2->SetMarkerSize(MS);
            gre3->SetMarkerSize(MS);

            if(wppt < ptbin){
            gre1->SetTitle(Form("dEta > %.1f, Pt: %.1f to %.0f", dEtaUnit, pt_edges[wppt-1], pt_edges[wppt] ));
            gre2->SetTitle(Form("dEta > %.1f, Pt: %.1f to %.0f", dEtaUnit, pt_edges[wppt-1], pt_edges[wppt] ));
            gre3->SetTitle(Form("dEta > %.1f, Pt: %.1f to %.0f", dEtaUnit, pt_edges[wppt-1], pt_edges[wppt] ));
            }
            if(wppt == ptbin){
            gre1->SetTitle(Form("dEta > %.1f, Pt: %.1f to inf", dEtaUnit, pt_edges[wppt-1]));
            gre2->SetTitle(Form("dEta > %.1f, Pt: %.1f to inf", dEtaUnit, pt_edges[wppt-1]));
            gre3->SetTitle(Form("dEta > %.1f, Pt: %.1f to inf", dEtaUnit, pt_edges[wppt-1]));
            }

            mg1->Add(gre1);
            mg2->Add(gre2);
            mg3->Add(gre3);
            
            //delete gre1;
            //delete gre2;
            //delete gre3;
        }
    }
 
    auto c2 = new TCanvas("c2","c2",800, 800);
    c2->cd();
    mg1->Draw("APL");
    mg1->SetTitle("AK8: V1 coefficient in pt (marker style) and dEta region (color).");
    mg1->GetXaxis()->SetTitle("N Chgd Daughters");
    mg1->GetYaxis()->SetTitle("V1");
    mg1->GetXaxis()->CenterTitle(true);
    mg1->GetYaxis()->SetTitleOffset(0.6);
    //mg1->GetYaxis()->CenterTitle(true);
    c2->BuildLegend(.6,.1,.9,.4);
    
    auto c3 = new TCanvas("c3","c3",800, 800);
    c3->cd();
    mg2->Draw("APL");
    mg2->SetTitle("AK8: V2 coefficient in pt (marker style) and dEta region (color).");
    mg2->GetXaxis()->SetTitle("N Chgd Daughters");
    mg2->GetYaxis()->SetTitle("V2");
    mg2->GetXaxis()->CenterTitle(true);
    //mg2->GetYaxis()->CenterTitle(true);
    c3->BuildLegend(.6,.6,.9,.9);
    
    auto c4 = new TCanvas("c4","c4",800, 800);
    c4->cd();
    mg3->Draw("APL");
    mg3->SetTitle("AK8: V3 coefficient in pt (marker style) and dEta region (color).");
    mg3->GetXaxis()->SetTitle("N Chgd Daughters");
    mg3->GetYaxis()->SetTitle("V3");
    mg3->GetYaxis()->SetTitleOffset(0.6);
    mg3->GetXaxis()->CenterTitle(true);
    //mg3->GetYaxis()->CenterTitle(true);
    c4->BuildLegend(.6,.1,.9,.4);

// *******************************************************************************
//      *********************************************************************
//           **********************************************************
//                ************************************************




}



