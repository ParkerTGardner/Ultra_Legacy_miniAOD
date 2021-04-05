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


void Drawer()
{
    TFile *f = new TFile("new_wei.root");

    TH1::SetDefaultSumw2(kFALSE);
    TH2::SetDefaultSumw2(kFALSE);

    int trackbin = 1;
    int ptbin = 1;

    int Xscale = 1900;
    int Yscale = 1000;

    //TCanvas* c1a          = new TCanvas("c1a","c1a", Xscale, Yscale);
    //c1a->Divide(ptbin,trackbin);
    //TCanvas* c1b          = new TCanvas("c1b","c1b", Xscale, Yscale);
    //c1b->Divide(ptbin,trackbin);
    TCanvas* c2          = new TCanvas("c2","c2", Xscale, Yscale);
    c2->Divide(ptbin,trackbin);
    TCanvas* c3          = new TCanvas("c3","c3", Xscale, Yscale);
    c3->Divide(ptbin,trackbin);
    TCanvas* c4          = new TCanvas("c4","c4", Xscale, Yscale);
    c4->Divide(ptbin,trackbin);
    //TCanvas* c4b          = new TCanvas("c4b","c4b", Xscale, Yscale);
    //c4b->Divide(ptbin,trackbin);
    //TCanvas* c5          = new TCanvas("c5","c5", Xscale, Yscale);
    //c5->Divide(ptbin,trackbin);
    TCanvas* c6          = new TCanvas("c6","c6", Xscale, Yscale);
    c6->Divide(ptbin,trackbin);

    c6->SetLeftMargin(0.2);
    c6->SetTheta(60.839);
    c6->SetPhi(38.0172);

    //TCanvas* c0a          = new TCanvas("c0a","c0a", Xscale, Yscale);
    //TCanvas* c0b          = new TCanvas("c0b","c0b", Xscale, Yscale);
    //TCanvas* c0c          = new TCanvas("c0c","c0c", Xscale, Yscale);
    //TCanvas* c0d          = new TCanvas("c0d","c0d", Xscale, Yscale);
    /*
       c0a->cd();
       TH1D* h0a;
       h0a = (TH1D*)f->Get("hN_ChgDAU");
    //h0a->SetStats(kFALSE);
    h0a->Draw("COLZ");

    c0b->cd();
    TH2D* h0b;
    h0b = (TH2D*)f->Get("hN_ChgDAUvsJPT");
    h0b->SetStats(kFALSE);
    h0b->Draw("COLZ");

    c0c->cd();
    TH1D* h0c;
    h0c = (TH1D*)f->Get("hEvent_Pass");
    //h0c->SetStats(kFALSE);
    h0c->Draw();

    c0d->cd();
    TH1D* h0d;
    h0d = (TH1D*)f->Get("hPass_Trig_Jet");
    //h0d->SetStats(kFALSE);
    h0d->Draw();
    */
    TH1D* hJ;
    hJ = (TH1D*)f->Get("hPass_Trig_Jet")->Clone();

    float VV = 0.25;

    for(int wtrk =1; wtrk < trackbin+1; wtrk++){
        /*
           c1a->cd(wtrk);
           TH2D* ha;
           ha = (TH2D*)f->Get(Form("hdNdEvsPT_%d",wtrk));
           ha->SetStats(kFALSE);
           ha->Draw("COLZ");

           c1b->cd(wtrk);
           TH1D* hb;
           hb = (TH1D*)f->Get(Form("hdNdE_%d",wtrk));
           hb->SetStats(kFALSE);
           hb->Draw();
           */
        for(int wppt =1; wppt <ptbin+1; wppt++){
            int where = wppt + (wtrk-1)*(ptbin);

            c2->cd(where);
            TH2D* h1;
            h1 = (TH2D*)f->Get(Form("hSig_%d_%d",wtrk,wppt))->Clone();
            h1->SetStats(kFALSE);
            h1->Scale(1/(hJ->GetEntries()));
            h1->Scale(1./(           h1->GetXaxis()->GetBinWidth(3)  ));
            h1->DrawCopy("SURF1 fb");

            c3->cd(where);
            TH2D* h2;
            h2 = (TH2D*)f->Get(Form("hBck_%d_%d",wtrk,wppt))->Clone();
            h2->SetStats(kFALSE);
            h2->Scale(1/(hJ->GetEntries()));
            h2->Scale(1./(           h2->GetXaxis()->GetBinWidth(3)  ));
            h2->DrawCopy("SURF1 fb");


            cout << "wtrk is " << wtrk << "  and wppt is " << wppt << "  and sig entries is " << h1->GetEntries() << "  and bck entries is " << h2->GetEntries() << endl;
            int YPmax = h1->GetNbinsX();
            int YPlow = floor(0.75*YPmax);
            cout << YPmax << " " << YPlow << endl;

            TH1D *histfit1 = (TH1D*) h1->ProjectionY("",YPlow,YPmax)->Clone();
            TH1D *histfit2 = (TH1D*) h2->ProjectionY("",YPlow,YPmax)->Clone();

            //c1a->cd(where);
            //histfit1->DrawCopy("PE");
            //c1b->cd(where);
            //histfit2->DrawCopy("PE");

            c4->cd(where);
            histfit1->Divide(histfit2);

            histfit1->Scale(h2->GetMaximum());

            std::string function2 = "(TMath::Cos(2*x))";
            TF1 efunc("efunc",function2.c_str(),-0.5*3.14159,1.5*3.14159);
            float   projMax = (histfit1->GetMaximum())/(2*3.14159);
            int     loopMax = 400000;

                for (int i = 0; i<loopMax; i++) {
                    float YY = efunc.GetRandom();

                    if( (YY > -3.14159/4 && YY<3.14159/4) || (YY>3*3.14159/4 && YY<5*3.14159/4)){

                        histfit1->Fill(YY, (projMax*(VV*VV)/(loopMax/40)) );

                    } else {

                        histfit1->Fill(YY, -(projMax*(VV*VV)/(loopMax/40)) );
                    }
                }




            histfit1->SetMarkerStyle(20);

            std::string function = "[0]/(TMath::Pi()*2)*(1+2*([1]*TMath::Cos(x)+[2]*TMath::Cos(2*x)+[3]*TMath::Cos(3*x)+[4]*TMath::Cos(4*x)+[5]*TMath::Cos(5*x)   ))";

            TF1 func1("deltaPhi1", function.c_str(), -0.5*3.14159, 1.5*3.14159);
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
            /*
               std::string functionFAKE = "[0]/(TMath::Pi()*2)*(1+2*([1]*TMath::Cos(x)+[2]*TMath::Cos(2*x)+[3]*TMath::Cos(3*x)+[4]*TMath::Cos(4*x)+[5]*TMath::Cos(5*x)   ))";

               TF1 func2F("deltaPhi1", functionFAKE.c_str(), -0.5*3.14157, 1.5*3.14157);
               func2F.SetParameter(0, 1.00*(func1.GetParameter(0)));
               func2F.SetParameter(1, 1.00*(func1.GetParameter(1)));
               func2F.SetParameter(2, 1.10*(func1.GetParameter(2)));
               func2F.SetParameter(3, 1.00*(func1.GetParameter(3)));
               func2F.SetParameter(4, func1.GetParameter(4));
               func2F.SetParameter(5, func1.GetParameter(5));

               func2F.SetLineColor(3);

            //cout << "f1  " << func2F.GetParameter(1) << "  and f2  " << func2F.GetParameter(2) << endl;
            //TF1 *fv1 = new TF1("fv1","([0]*TMath::Cos(x))"  ,-0.5*3.14157, 1.5*3.14157);
            //TF1 *fv2 = new TF1("fv2","([0]*TMath::Cos(2*x))",-0.5*3.14157, 1.5*3.14157);
            //TF1 *fv3 = new TF1("fv3","([0]*TMath::Cos(3*x))",-0.5*3.14157, 1.5*3.14157);

            //fv1->SetLineColor(1);
            //fv2->SetLineColor(2);
            //fv3->SetLineColor(3);

            //fv1->DrawCopy();
            //fv2->DrawCopy("SAME");
            //fv3->DrawCopy("SAME");
            */
            histfit1->SetStats(kFALSE);
            histfit1->DrawCopy("PE");

            //c4b->cd(where);
            //histfit1->SetStats(kFALSE);
            //histfit1->DrawCopy("PE");
            //func2F.DrawCopy("SAME PE");

            h1->Divide(h2);
            h1->Scale(h2->GetMaximum());
            h1->SetStats(kFALSE);
            //c5->cd(where);
            //h1->DrawCopy("COLZ");

            c6->cd();

            for (int j = 0; j<81; j++){ 
                for (int i = 0; i<loopMax; i++) {

                    float YY = efunc.GetRandom();


                    if( (YY > -3.14159/4 && YY<3.14159/4) || (YY>3*3.14159/4 && YY<5*3.14159/4)){

                        h1->Fill(-6+(j*0.149), YY, (projMax*(VV*VV)/(loopMax/40)) );

                    } else {

                        h1->Fill(-6+(j*0.149), YY, -(projMax*(VV*VV)/(loopMax/40)) );
                        //continue;
                    }
                }
            }







            c6->SetLeftMargin(0.2);
            c6->SetTheta(60.839);
            c6->SetPhi(38.0172);

            h1->Draw("SURF1 fb");

            c6->SetLeftMargin(0.2);
            c6->SetTheta(60.839);
            c6->SetPhi(38.0172);

            //h1->GetXaxis()->SetRangeUser(-3,3);
            //int voo = h1->FindBin(0,0,0);
            //int qqq = floor((h1->GetBinContent(voo))/2);
            //int spy = h1->GetMinimumBin();
            //int dji = h1->GetBinContent(spy);
            //if( qqq < dji) qqq = 2*qqq;
            //if( dji < 0) dji = 0;
            //h1->GetZaxis()->SetRangeUser(dji,qqq);


            //cout << "voo is " << voo << endl;
            //cout << "qqq is " << qqq << endl;
            //cout << "spy is " << spy << endl;
            //cout << "dji is " << dji << endl;

            //delete fv1;
            //delete fv2;
            //delete fv3;
            //delete histfit1;
            //delete histfit2;
            //delete h2;
            //delete h1;

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
