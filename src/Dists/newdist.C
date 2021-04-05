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


void newdist()
{
    //TFile *fF1 = new TFile("full876_p_feb.root");
    //TFile *fF2 = new TFile("Full678DistS.root");

    TFile *fP1 = new TFile("pythia470_f_p_feb.root");
    //TFile *fP2 = new TFile("P470DistS.root");

    TFile *fH1 = new TFile("herwig_P_feb.root");
    //TFile *fH2 = new TFile("HerwigDistS.root");


    TH1::SetDefaultSumw2(kTRUE);
    TH2::SetDefaultSumw2(kTRUE);

    const int trackbin 	= 7;
    const int trackbinC = 2;
    const int ptbin 	= 45;
    //const int   trackbinbounds[trackbin]    = {85,86,87,88,89,90,91,92};
    const int   trackbinbounds[trackbin]    = {0, 85,86,87,88,89,90};
    const float ptbinbounds_lo[ptbin]       = { 4,  5,  6, 10};
    const float ptbinbounds_hi[ptbin]       = {200, 30, 30, 30, 30};

    int ptwant 	= 2;//manually change through 13... or create a for loop that does this but then you will get a 78 canvases.

    int pt_lo 	= ptbinbounds_lo[ptwant-1];
    int pt_hi 	= ptbinbounds_hi[ptwant-1];


    int Xscale 	= 1400;
    int Yscale 	= 800;

    //Defining Canvas Arrays
    //TCanvas* c1 = new TCanvas("c1", "c1", Xscale, Yscale);
    //TCanvas* c2 = new TCanvas("c2", "c2", Xscale, Yscale);

    TCanvas* cA[trackbinC];
    for(int i = 0; i < trackbinC; i++){
        cA[i] = new TCanvas(Form("cA_%d",i),Form("cA_%d",i), Xscale, Yscale);
    }

    TCanvas* cB[trackbinC];
    for(int i = 0; i < trackbinC; i++){
        cB[i] = new TCanvas(Form("cB_%d",i),Form("cB_%d",i), Xscale, Yscale);
    }

    TCanvas* cC[trackbinC];
    for(int i = 0; i < trackbinC; i++){
        cC[i] = new TCanvas(Form("cC_%d",i),Form("cC_%d",i), Xscale, Yscale);
    }

    TCanvas* cD[trackbinC];
    for(int i = 0; i < trackbinC; i++){
        cD[i] = new TCanvas(Form("cD_%d",i),Form("cD_%d",i), Xscale, Yscale);
    }


    for(int wtrk = 0; wtrk < trackbinC; wtrk++){

        cA[wtrk]->Divide(2,1);
        cB[wtrk]->Divide(2,1);
        cC[wtrk]->Divide(2,1);
        cD[wtrk]->Divide(2,1);
        


/////////
/*
        TH1D* hLeadJTanyF;
        hLeadJTanyF = (TH1D*)fF1->Get("hLeadJTany")->Clone();
        TH2D* hLeadJTany2dF;
        hLeadJTany2dF = (TH2D*)fF1->Get("hLeadJTany2d")->Clone();

        TH1D* hLeadJTchgF;
        hLeadJTchgF = (TH1D*)fF1->Get("hLeadJTchg")->Clone();
        TH2D* hLeadJTchg2dF;
        hLeadJTchg2dF = (TH2D*)fF1->Get("hLeadJTchg2d")->Clone();

        cA[wtrk]->cd(1);
        hLeadJTanyF->Draw();

        cB[wtrk]->cd(1);
        hLeadJTchgF->Draw();

        cC[wtrk]->cd(1);
        hLeadJTany2dF->Draw("COLZ");

        cD[wtrk]->cd(1);
        hLeadJTchg2dF->Draw("COLZ");
*/
//////////

        TH1D* hLeadJTanyH;
        hLeadJTanyH = (TH1D*)fH1->Get("hLeadJTany")->Clone();
        TH2D* hLeadJTany2dH;
        hLeadJTany2dH = (TH2D*)fH1->Get("hLeadJTany2d")->Clone();

        TH1D* hLeadJTchgH;
        hLeadJTchgH = (TH1D*)fH1->Get("hLeadJTchg")->Clone();
        TH2D* hLeadJTchg2dH;
        hLeadJTchg2dH = (TH2D*)fH1->Get("hLeadJTchg2d")->Clone();

        cA[wtrk]->cd(1);
        hLeadJTanyH->Draw();

        cB[wtrk]->cd(1);
        hLeadJTchgH->Draw();

        cC[wtrk]->cd(1);
        hLeadJTany2dH->Draw("COLZ");

        cD[wtrk]->cd(1);
        hLeadJTchg2dH->Draw("COLZ");
//////////

        TH1D* hLeadJTanyP;
        hLeadJTanyP = (TH1D*)fP1->Get("hLeadJTany")->Clone();
        TH2D* hLeadJTany2dP;
        hLeadJTany2dP = (TH2D*)fP1->Get("hLeadJTany2d")->Clone();

        TH1D* hLeadJTchgP;
        hLeadJTchgP = (TH1D*)fP1->Get("hLeadJTchg")->Clone();
        TH2D* hLeadJTchg2dP;
        hLeadJTchg2dP = (TH2D*)fP1->Get("hLeadJTchg2d")->Clone();

        cA[wtrk]->cd(2);
        hLeadJTanyP->Draw();

        cB[wtrk]->cd(2);
        hLeadJTchgP->Draw();

        cC[wtrk]->cd(2);
        hLeadJTany2dP->Draw("COLZ");

        cD[wtrk]->cd(2);
        hLeadJTchg2dP->Draw("COLZ");



    }
}


