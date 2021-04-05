//blh
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
//#include "coordinateTools.h"

using TMath::ATan;
using TMath::Exp;


void Drawer()
{
    //TFile *f = new TFile("primary_full_save.root");
    TFile *f = new TFile("Total_Sum.root");

  /*  
    TCanvas* cA = new TCanvas("cA","cA" , 1000, 1000);
    TCanvas* cB = new TCanvas("cB","cB" , 1000, 1000);
    TCanvas* cC = new TCanvas("cC","cC" , 1000, 1000);


    TH2D* h1 = (TH2D*)f->Get(Form("hSig_%d_%d",1,1));
    TH2D* h2 = (TH2D*)f->Get(Form("hBck_%d_%d",1,1));
    
    cA->cd();
    h1->DrawCopy("COLZ");
    cB->cd();
    h2->Draw("COLZ");

    h1->Divide(h2);
    cC->cd();
    h1->Draw("COLZ");
  */

    TH1::SetDefaultSumw2(kFALSE);
    TH2::SetDefaultSumw2(kFALSE);

	int trackbin = 7;
	int ptbin = 3;
    int scale = 1000;

    TCanvas* c1          = new TCanvas("c1","c1" , scale, scale);
    c1->Divide(ptbin,trackbin);
    
    TCanvas* c2          = new TCanvas("c2","c2" , scale, scale);
    c2->Divide(ptbin,trackbin);

    TCanvas* c3a          = new TCanvas("c3a","c3a", scale, scale);
    c3a->Divide(ptbin,trackbin);
/*
    TCanvas* c3b          = new TCanvas("c3b","c3b", scale, scale);
    c3b->Divide(ptbin,trackbin);

    TCanvas* c3c          = new TCanvas("c3c","c3c", scale, scale);
    c3c->Divide(ptbin,trackbin);
*/


    //TH2D* hBckrnd[trackbin][ptbin];
/*

    const int dEta_N = 5;
    double dEta_bin_cut[dEta_N]  = {0};

    double v2[dEta_N][trackbin]  = {0};
    double v2e[dEta_N][trackbin] = {0};
    double v3[dEta_N][trackbin]  = {0};
    double v3e[dEta_N][trackbin] = {0};

    for(int i = 0; i < dEta_N; i++){
    dEtabin_cut[i]=10*(i+1);

    }
*/


    for(int wtrk =1; wtrk <trackbin+1; wtrk++){
        for(int wppt =1; wppt <ptbin+1; wppt++){
            int where = wppt + (wtrk-1)*(ptbin);
            c1->cd();

            TH2D* h1 = (TH2D*)f->Get(Form("hSig_%d_%d",wtrk,wppt));
            TH2D* h2 = (TH2D*)f->Get(Form("hBck_%d_%d",wtrk,wppt));

            h1->Divide(h2);
           
            TH1D *histfit1 = (TH1D*) h1->ProjectionY("",40,80)->Clone();

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

            fv1->DrawCopy();
            fv2->DrawCopy("SAME");
            fv3->DrawCopy("SAME");

            c2->cd(where);
            histfit1->DrawCopy("PE");

            c1->cd(where);
            h1->DrawCopy("COLZ");


/*
            c1full->cd(where);
            TF1 func3("deltaPhi3", function.c_str(), -0.5*3.14157, 1.5*3.14157);
            func3.SetParameter(0, histfit3->GetMaximum());
            func3.SetParameter(1, 0.1);
            func3.SetParameter(2, 0.1);
            func3.SetParameter(3, 0.1);

            histfit3->SetMarkerStyle(20);

            histfit3->Fit(&func3, "q");
            histfit3->Fit(&func3, "q");
            histfit3->Fit(&func3, "m q");
            histfit3->Fit(&func3, "m q E");
            auto fitResult3 = histfit3->Fit(&func3, "m S E q");

            histfit3->DrawCopy("PE");
            c1->cd(where);
            h3->DrawCopy("COLZ");
*/
            //cout<< "track bin " << wtrk << "   and pt bin " << wppt <<endl;
            //cout << "v2 is " << func1.GetParameter(2) << endl;
            //cout << "v3 is " << func1.GetParameter(3) << endl;


            delete fv1;
            delete fv2;
            delete fv3;
            delete histfit1;
            delete h2;
            delete h1;

        }
    }
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
