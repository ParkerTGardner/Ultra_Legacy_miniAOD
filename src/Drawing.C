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

void Drawing()
{
    //TFile *f = new TFile("SumP678_full_eta.root");
    //TFile *f = new TFile("safeSumP678.root");
    TFile *f = new TFile("MC_SumP600.root");
    //TFile *f = new TFile("MC_SumPHerwig.root");

    TH1::SetDefaultSumw2(kTRUE);
    TH2::SetDefaultSumw2(kTRUE);

    const int trackbin 	= 6;
    const int trackbinC = 6;
    const int ptbin 	= 15;
    //const int   trackbinbounds[trackbin]    = {85,86,87,88,89,90,91,92};
    const int   trackbinbounds[trackbin]    = {85,86,87,88,89,90};
    const float ptbinbounds_lo[ptbin]       = {0,   3,  4,  5,  6, 10,  0,  0,  0,  3,  3,  5,  5, 10, 10};
    const float ptbinbounds_hi[ptbin]       = {30, 30, 30, 30, 30, 30, 10, 20, 40, 20, 40, 20, 40, 20, 40};

    int ptwant 	= 4;//manually change through 13... or create a for loop that does this but then you will get a 78 canvases.

    int pt_lo 	= ptbinbounds_lo[ptwant-1];
    int pt_hi 	= ptbinbounds_hi[ptwant-1];

    float etalo     = .9;
    float etahi     = 3.5;

    int Xscale 	= 1900;
    int Yscale 	= 1000;

    //Defining Canvas Arrays
    TCanvas* c1 = new TCanvas("c1", "c1", Xscale, Yscale);
    TCanvas* c2 = new TCanvas("c2", "c2", Xscale, Yscale);

    TCanvas* cA[trackbinC];
    for(int i = 0; i < trackbinC; i++){
        cA[i] = new TCanvas(Form("cA_%d",i),Form("cA_%d",i), Xscale, Yscale);
    }
    //TCanvas* cB[trackbin];
    //for(int i = 0; i < trackbin; i++){
    //    cB[i] = new TCanvas(Form("cB_%d",i),Form("cB_%d",i), Xscale, Yscale);
    //}

    //Defining Jet count for normalization
    TH1D* hJ;
    hJ = (TH1D*)f->Get("hJet_Pass")->Clone();
    for(int i = 0; i < trackbin+2; i++){
        int jets = hJ->GetBinContent(i);
        cout << "for i = " << i << " the fill is " << jets << endl;
    }

    //Main Loops
    //for(int wtrk = 0; wtrk < 2; wtrk++){
    for(int wtrk = 0; wtrk < trackbinC; wtrk++){

        cA[wtrk]->Divide(2,1);
        //cB[wtrk]->Divide(2,1);

        TH2D* h1;
        //hSig_90_3_to_40
        h1 = (TH2D*)f->Get(Form("hSig_%d_%d_to_%d",trackbinbounds[wtrk], pt_lo, pt_hi))->Clone();
        float BW1 = h1->GetXaxis()->GetBinWidth(3);
        h1->Scale(1/(hJ->GetBinContent(wtrk+1)));
        h1->Scale(1./(BW1));

        TH2D* h2;
        h2 = (TH2D*)f->Get(Form("hBck_%d_%d_to_%d",trackbinbounds[wtrk], pt_lo, pt_hi))->Clone();
        float BW2 = h2->GetXaxis()->GetBinWidth(3);
        h2->Scale(1/(hJ->GetBinContent(wtrk+1)));
        h2->Scale(1./(BW2));

        c1->cd();
        h1->DrawCopy("SURF1 FB BB");
        c2->cd();
        h2->DrawCopy("SURF1 FB BB");

        //h1->GetXaxis->SetRangeUser(-4,4)
        //h2->GetXaxis->SetRangeUser(-4,4)


        cA[wtrk]->cd(1);

        int YPmax = h1->GetNbinsX();
        int YPlo = floor(0.5*YPmax+floor(2.0/BW1));
        int YPhi = floor(0.5*YPmax+floor(4.3/BW1));

        //int YPlo = 22;
        //int YPhi = 22;

        /*
        cout << " *** WTRK " << wtrk << " *** " << endl;
        cout << "IN BIN: YPmax is: " 	<< YPmax
            << " , YPlo is: " 
            << YPlo 	
            << " , YPhi is: " 
            << YPhi 	
            << endl;
        cout << "IN ETA: YPmax is: " 	<< (YPmax*BW1)/2 
            << " , YPlo is: " 
            << YPlo*BW1 - ((YPmax*BW1)/2) 
            << " , YPhi is: " 
            << YPhi*BW1 - ((YPmax*BW1)/2)   
            << endl;
        cout << "Guess?: YPmax is: " 	<< 6 
            << "? , YPlo is: " 
            << 2 
            << "? , YPhi is: " 
            << 4 
            << "?"   
            << endl;
        */
        cout << "@@@@@@@@@@@" << endl;
        cout << "h2->GetMaximum() : ... " << h2->GetMaximum() << endl;
        cout << "@@@@@@@@@@@" << endl;

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
        histfit1->Fit(&func1, "q");
        histfit1->Fit(&func1, "q");
        histfit1->Fit(&func1, "m q");
        histfit1->Fit(&func1, "m q");
        histfit1->Fit(&func1, "m q E");
        auto fitResult1 = histfit1->Fit(&func1, "m S E q");

        cout << "v2 with coeff = " 	<< func1.GetParameter(0)/(TMath::Pi()*2)*(1+2*(func1.GetParameter(2))) 		<< endl;
        cout << "v2 Ffit param = " 	<< func1.GetParameter(2) 		<< endl;
        cout << "0/(2pi) param = " 	<< func1.GetParameter(0)/(TMath::Pi()*2) 	<< endl;

        cA[wtrk]->cd(1);
        histfit1->SetStats(kFALSE);
        histfit1->SetTitle(Form("dPhi Profile for dEta %.2f to %.2f", YPlo*BW1-6.15,YPhi*BW1-6.15 ));
        histfit1->GetYaxis()->SetTitle("dN(chg) / dP per jet");
        histfit1->GetYaxis()->SetTitleOffset(1.3);
        histfit1->GetXaxis()->SetTitle("d Phi");
        histfit1->DrawCopy("PE");

        cA[wtrk]->cd(2);
        cA[wtrk]->SetTheta(60.839);
        cA[wtrk]->SetPhi(38.0172);
        //figure out the rotation
        h1->Divide(h2);
        h1->GetXaxis()->SetRangeUser(-3.0,3.0);
        h1->SetStats(kFALSE);
        h1->SetTitle(Form("Charged Mult > %d, pt from %.2f to %.2f", trackbinbounds[wtrk], (float)pt_lo/10, (float)pt_hi/10));

        //h1->GetZaxis()->SetLabelSize();
        h1->GetZaxis()->SetTitle("dN(chg) / dEdP per jet");
        h1->GetZaxis()->SetTitleOffset(1.45);
        //h1->GetZaxis()->SetTitleSize();

        //h1->GetXaxis()->SetLabelSize();
        h1->GetXaxis()->SetTitle("d Eta");
        h1->GetXaxis()->SetTitleOffset(2);
        //h1->GetXaxis()->SetTitleSize();

        //h1->GetYaxis()->SetLabelSize();
        h1->GetYaxis()->SetTitle("d Phi");
        h1->GetYaxis()->SetTitleOffset(2);
        //h1->GetYaxis()->SetTitleSize();
        int binx2 = h1->GetXaxis()->FindBin(.1);
        int biny2 = h1->GetYaxis()->FindBin(3.1415);

        int binx1 = h1->GetXaxis()->FindBin(-2.5);
        int biny1 = h1->GetYaxis()->FindBin(1);

        float r2 = (h1->GetBinContent(binx2,biny2));
        float r1 = (h1->GetBinContent(binx1,biny1)) - .4;

        h1->GetZaxis()->SetRangeUser(r1, r1 + (.65*(r2-r1)) );
        cA[wtrk]->SetTheta(60.839);
        cA[wtrk]->SetPhi(38.0172);
        h1->DrawCopy("SURF1 FB BB");
        cA[wtrk]->SetTheta(60.839);
        cA[wtrk]->SetPhi(38.0172);
/*
        //vary the delta eta range! function here outputs a plot of the v2 varried for various eta lower limit 
        ////each of the 6 wtrk plots will have a correspodning canvas showing the "transition" plateau, both coeff and param 
        ////12 canvases per run, 24 tpads
        cB[wtrk]->cd(1);

        int varyYPlo1  = 1 + ceil(0.5*YPmax+floor(etalo/BW1));
        int varyYPlo2  = ceil(0.5*YPmax+floor(etahi/BW1));

        cout << varyYPlo1 <<endl;
        cout << varyYPlo2 <<endl;
        cout << YPmax <<endl;
        cout << BW1 <<endl;
        int varylength = varyYPlo2 - varyYPlo1;
        cout << "varylength " << varylength <<endl;
        float v2param[8] = {0};
        float v2coeff[8] = {0};

        for(int i=0; i < varylength; i++){

            TH2D* h1n;
            //hSig_90_3_to_40
            h1n = (TH2D*)f->Get(Form("hSig_%d_%d_to_%d",trackbinbounds[wtrk], pt_lo, pt_hi))->Clone();
            h1n->Scale(1/(hJ->GetBinContent(wtrk+1)));
            h1n->Scale(1./(BW1));

            TH2D* h2n;
            h2n = (TH2D*)f->Get(Form("hBck_%d_%d_to_%d",trackbinbounds[wtrk], pt_lo, pt_hi))->Clone();
            h2n->Scale(1/(hJ->GetBinContent(wtrk+1)));
            h2n->Scale(1./(BW2));

            TH1D *varyfit1 = (TH1D*) h1n->ProjectionY("",varyYPlo1+i,YPhi)->Clone();
            TH1D *varyfit2 = (TH1D*) h2n->ProjectionY("",varyYPlo1+i,YPhi)->Clone();

            std::string varyfunction = "[0]/(TMath::Pi()*2)*(1+2*([1]*TMath::Cos(x)+[2]*TMath::Cos(2*x)+[3]*TMath::Cos(3*x)+[4]*TMath::Cos(4*x)+[5]*TMath::Cos(5*x)))";

            TF1 varyfunc1("varydeltaPhi1", varyfunction.c_str(), -0.5*TMath::Pi(), 1.5*TMath::Pi());
            varyfunc1.SetParameter(0, varyfit1->GetMaximum());
            varyfunc1.SetParameter(1, 0.1);
            varyfunc1.SetParameter(2, 0.1);
            varyfunc1.SetParameter(3, 0.1);
            varyfunc1.SetParameter(4, 0.1);
            varyfunc1.SetParameter(5, 0.1);

            varyfit1->Fit(&varyfunc1, "q");
            varyfit1->Fit(&varyfunc1, "q");
            varyfit1->Fit(&varyfunc1, "m q");
            varyfit1->Fit(&varyfunc1, "m q");
            varyfit1->Fit(&varyfunc1, "m q E");
            //varyfit1->Fit(&varyfunc1, "m S E q");
            auto varyfitResult1 = varyfit1->Fit(&varyfunc1, "m S E q");

            v2param[i] = varyfunc1.GetParameter(2);
            v2coeff[i] = 1; //pow((varyfunc1.GetParameter(2))*(2*TMath::Pi()/(varyfunc1.GetParameter(0))),0.5);

            delete h1n;
            delete h2n;

        }
        //new histogram
        TH1D* h3   = new TH1D("h3","h3", varylength,varyYPlo1*BW1-6,varyYPlo2*BW1-6);
        TH1D* h4   = new TH1D("h4","h4", varylength,0,varylength);
        for(int i=0; i < varylength; i++){
            h3->Fill((varyYPlo1*BW1)-6+(i*BW1)+(0.5*BW1), v2param[i]);
            h4->Fill(i, v2coeff[i]);
        }
        cB[wtrk]->cd(1);
        h3->SetMarkerStyle(20);
        h3->SetMarkerColor(2); 
        h3->SetTitle(Form("v2 fit param for lower range dE > (%.2f to %.2f) ", varyYPlo1*BW1-6, varyYPlo2*BW1-6)); 
        h3->GetYaxis()->SetTitle("v2 fit param");
        h3->GetYaxis()->SetTitleOffset(1.4);
        h3->GetXaxis()->SetTitle(Form(" Chg Mult>%d, pt from %.2f to %.2f ", trackbinbounds[wtrk], (float)pt_lo/10, (float)pt_hi/10));
        h3->SetStats(kFALSE);
        h3->DrawCopy("HIST P");
        cB[wtrk]->cd(2);
        h4->DrawCopy("HIST P");
        delete h3;
        delete h4;
*/

    }
}




/*


   for(int i = LLow; i < HHi; i++){
   int YPlow = i;
   cout << YPmax << " " << YPlow << endl;

   TH1D *histfit1 = (TH1D*) h1->ProjectionY("",YPlow,YPmax)->Clone();
   TH1D *histfit2 = (TH1D*) h2->ProjectionY("",YPlow,YPmax)->Clone();

//c3->cd(ii);

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

histfit1->Fit(&func1, "q");
histfit1->Fit(&func1, "q");
histfit1->Fit(&func1, "m q");
histfit1->Fit(&func1, "m q");
histfit1->Fit(&func1, "m q E");
auto fitResult1 = histfit1->Fit(&func1, "m S E q");


cout << func1.GetParameter(2) << endl;
cout << func1.GetParameter(0)/(TMath::Pi()) << endl;

c[wtrk]->cd(1);
histfit1->SetStats(kFALSE);
histfit1->GetYaxis->SetLabel();
histfit1->GetXaxis->SetLabel();
histfit1->Draw("PE");

h1->Divide(h2);

h1->SetStats(kFALSE);
h1->GetZaxis->SetLabel();
h1->GetYaxis->SetLabel();
h1->GetXaxis->SetLabel();
c[wtrk]->cd(2);
h1->Draw("SURF1");
}
*/
