#define MyClass_cxx

#include "include/MyClass.h"
#include "include/coordinateTools.h"

#include <iostream>
#include <iomanip>

#include <vector>
#include "math.h"
#include <numeric>
#include <string>
#include <fstream>
#include <algorithm>
#include <iterator>

#include <TStyle.h>
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TProfile.h"
#include "TGraph.h"
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TCutG.h"
#include "TRandom3.h"

using TMath::ATan;
using TMath::Exp;

const int   trackbin                    = 3;
//const int   ptbin                     = 15;
const int   ptbin                       = 4;
const float ptmin                       = 0.0;
const float ptmax                       = 4.0;
const int   trackbinbounds[trackbin]    = {85,87,89};
//const float ptbinbounds_lo[ptbin]       = {0.0, 0.3, 0.4, 0.5, 0.6, 1.0, 0.0, 0.0, 0.0, 0.3, 0.3, 0.5, 0.5, 1.0, 1.0}; 
const float ptbinbounds_lo[ptbin]       = {0.4, 0.5, 0.6, 1.0};
//const float ptbinbounds_hi[ptbin]       = {3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 1.0, 2.0, 4.0, 2.0, 4.0, 2.0, 4.0, 2.0, 4.0}; 
const float ptbinbounds_hi[ptbin]       = {3.0, 3.0, 3.0, 3.0};
const float ptbinmin                    = 0.4;
const float ptbinmax                    = 3.0;
const double        PI = 3.14159265359;
const float         EtaBW = 0.3;
const float         PhiBW = TMath::Pi()/16;
const int bin_WRTJ_Eta        = 150;
const int low_WRTJ_Eta_Bin  =  0;
const int high_WRTJ_Eta_Bin = 10;
const int bin_Eta             = 60;
const int low_Eta_Bin       = -3;
const int high_Eta_Bin      = 3;

const int EPD_xb  = 150;
const int EPD_yb  = 120;
const int EPD_xhi = 10;
const int EPD_xlo = 0;
const int EPD_yhi = 4;
const int EPD_ylo = -4;

const int bin5 = 5;
const int bin1 = 30;
const int bin3 = 50;
const int bin60 = 60;
const int bin2 = 100;
const int bin120 = 120;
const int bin150 = 150;
const int bin100 = 100;
const int bin200 = 200;
const int bin2000 = 2000;
const int bin1600 = 1600;

const int HB1 = 1;
const int HB2 = 2;
const int HB4 = 4;
const int HB15 = 15;
const int HB20 = 20;
const int HB110 = 110;
const int HB120 = 120;
const int HB4000 = 4000;

const double HB2pi = 2*PI;
const double HB15pi = 1.5*PI;

const int LB0 = 0;
const int LB1 = 1;
const int LB50 = 50;
const int LB470 = 470;
const int LB500 = 500;
const int LB2 = -15;
const double LB2pi = -2*PI;

//---------------------------------------------------------------------CUTS
const float EtaCut      = 0.0;
const float jetEtaCut   = 1.6;
const float jetPtCut    = 500.0;


//Bool folds

bool F_eventpass(std::vector< float > *genJetPt, int jetN, float jetPtCut){
    if(jetN > 200)              return false;
    if((*genJetPt)[0] < jetPtCut) return false;

    return true;
}

bool F_jetpass(std::vector< float > * genJetEta, std::vector< float > * genJetPt, int     ijet, float   jetPtCut){
    if(fabs( (*genJetEta)[ijet])   >jetEtaCut)   return false;
    if((*genJetPt)[ijet]           <jetPtCut)    return false;
    return true;
}

std::vector<int> F_genNix(std::vector< float > * i1, std::vector< int > * i2, std::vector< float > * i3, float i4, float i5) {
    std::vector<int> result(i1->size());
    std::iota(result.begin(), result.end(), 0);
    std::sort(result.begin(), result.end(),
            [&](int A, int B) -> bool {
            return (*i2)[A] > (*i2)[B];
            });
    result.erase(std::remove_if(
                result.begin(), result.end(),
                [&](const int x) {
                return (*i3)[x] > i4;
                }), result.end());
    result.erase(std::remove_if(
                result.begin(), result.end(),
                [&](const int x) {
                return (*i1)[x] < i5;
                }), result.end());
    return result;
}

std::vector<int> F_jetNix(std::vector< float > * i1, std::vector< int > * i2, std::vector< float > * i3, float i4, float i5, int i6) {
    std::vector<int> result(i1->size());
    std::iota(result.begin(), result.end(), 0);
    std::sort(result.begin(), result.end(),
            [&](int A, int B) -> bool {
            return (*i2)[A] > (*i2)[B];
            });
    result.erase(std::remove_if(
                result.begin(), result.end(),
                [&](const int x) {
                return (*i3)[x] > i4;
                }), result.end());
    result.erase(std::remove_if(
                result.begin(), result.end(),
                [&](const int x) {
                return (*i1)[x] < i5;
                }), result.end());
    result.erase(std::remove_if(
                result.begin(), result.end(),
                [&](const int x) {
                return (*i2)[x] < i6;
                }), result.end());
    return result;
}



void MyClass::Loop(int job){

    TH1::SetDefaultSumw2(kTRUE);
    TH2::SetDefaultSumw2(kTRUE);

    //Initializing Histograms
    TH2D* hPairs   = new TH2D("hPairs","hPairs",  trackbin, LB0,trackbin, ptbin, LB0, ptbin);

    TH1D* hTotal_ijet   = new TH1D("hTotal_ijet","hTotal_ijet", trackbin,LB0,trackbin);
    TH1D* hExtra_ijet   = new TH1D("hExtra_ijet","hExtra_ijet", HB4,LB0,HB4);

    TH1D* hEvent_Pass   = new TH1D("hEvent_Pass","hEvent_Pass", trackbin,LB0,trackbin);
    TH1D* hEvent_Fail   = new TH1D("hEvent_Fail","hEvent_Fail", trackbin,LB0,trackbin);
    TH1D* hJet_Pass     = new TH1D("hJet_Pass"  ,"hJet_Pass"  , trackbin,LB0,trackbin);
    TH1D* hJet_Fail     = new TH1D("hJet_Fail"  ,"hJet_Fail"  , bin120,LB0,bin120);
    TH1D* hMult_Fail    = new TH1D("hMult_Fail" ,"hMult_Fail" , bin120,LB0,bin120);

    TH1D* hLeadJTany      = new TH1D("hLeadJTany" ,"hLeadJTany"     , bin2 , LB0   , HB20);
    TH2D* hLeadJTany2d    = new TH2D("hLeadJTany2d" ,"hLeadJTany2d" , bin2 , LB0   , HB20, bin150, LB0, bin150);
    TH1D* hLeadJTchg      = new TH1D("hLeadJTchg" ,"hLeadJTchg"     , bin2 , LB0   , HB20);
    TH2D* hLeadJTchg2d    = new TH2D("hLeadJTchg2d" ,"hLeadJTchg2d" , bin2 , LB0   , HB20, bin120, LB0, bin120);
    TH2D* h2GenN_JetN         = new TH2D("h2GenN_JetN"       ,"h2GenN_JetN"       , bin5,LB0,bin5, bin5,LB0,bin5);
    //TH2D* h2GenN_JetN_W       = new TH2D("h2GenN_JetN_W"       ,"h2GenN_JetN_W"       , bin5,LB0,bin5, bin5,LB0,bin5);
    TH1D* h1JetN_d_GenN       = new TH1D("h1JetN_d_GenN"    ,"h1JetN_d_GenN"    , bin5, -LB1,bin5);
    //TH1D* h1JetN_d_GenN_W     = new TH1D("h1JetN_d_GenN_W"    ,"h1JetN_d_GenN_W"    , bin5,-LB1,bin5);
    TH1D* h1GenN              = new TH1D("h1GenN"    ,"h1GenN"    , bin5,LB0,bin5);
    //TH1D* h1GenN_W            = new TH1D("h1GenN_W"    ,"h1GenN_W"    , bin5,LB0,bin5);
    TH1D* h1JetN              = new TH1D("h1JetN"    ,"h1JetN"    , bin5,LB0,bin5);
    //TH1D* h1JetN_W            = new TH1D("h1JetN_W"    ,"h1JetN_W"    , bin5,LB0,bin5);

    TH2D* hGenExtra_pt_ChgMult      = new TH2D("hGenExtra_pt_ChgMult"   ,"hGenExtra_pt_ChgMult"   , bin120,LB470,bin1600, bin60,LB50,HB110);
    TH2D* hGenExtra_pt_Eta          = new TH2D("hGenExtra_pt_Eta"       ,"hGenExtra_pt_Eta"       , bin120,LB470,bin1600, bin3,-HB4,HB4);
    TH2D* hGenExtra_Eta_ChgMult     = new TH2D("hGenExtra_Eta_ChgMult"  ,"hGenExtra_Eta_ChgMult"  , bin3,-HB4,HB4,        bin60,LB50,HB110);
    TH2D* hRecExtra_pt_ChgMult      = new TH2D("hRecExtra_pt_ChgMult"   ,"hRecExtra_pt_ChgMult"   , bin120,LB470,bin1600, bin60,LB50,HB110);
    TH2D* hRecExtra_pt_Eta          = new TH2D("hRecExtra_pt_Eta"       ,"hRecExtra_pt_Eta"       , bin120,LB470,bin1600, bin3,-HB4,HB4);
    TH2D* hRecExtra_Eta_ChgMult     = new TH2D("hRecExtra_Eta_ChgMult"  ,"hRecExtra_Eta_ChgMult"  , bin3,-HB4,HB4,        bin60,LB50,HB110);

    TH2D* h2ChgMult_v_branch = new TH2D("h2ChgMult_v_branch", "h2ChgMult_v_branch" , bin120, LB0, bin120, bin120, LB0, bin120);
    TH1D* h1ChgMult_d_branch = new TH1D("h1ChgMult_d_branch", "h1ChgMult_d_branch" , bin1, -HB15, HB15);


    TH1D* hJet_Kin_Eta[trackbin];
    TH1D* hJet_Kin_Phi[trackbin];

    TH2D* hdNdEvsPT[trackbin][ptbin];
    TH1D* hdNdE[trackbin][ptbin];

    TH2D* hEPDraw[trackbin][ptbin];

    TH2D* hSignalShifted[trackbin][ptbin];
    TH2D* hBckrndShifted[trackbin][ptbin];
    TH2D* hSignalCentral[trackbin][ptbin];
    TH2D* hBckrndCentral[trackbin][ptbin];


    TH2D* hEPDrawMinus[trackbin][ptbin];
    TH2D* hEPDrawPluss[trackbin][ptbin];

    TH1D* hDau_Kin_WRTJ_Phi[trackbin][ptbin];
    TH1D* hDau_Kin_WRTJ_Pt[trackbin][ptbin];
    TH1D* hDau_Kin_WRTJ_Theta[trackbin][ptbin];
    TH1D* hDau_Kin_WRTJ_Eta[trackbin][ptbin];

    TH1D* hDau_Kin_Phi[trackbin][ptbin];
    TH1D* hDau_Kin_Pt[trackbin][ptbin];
    TH1D* hDau_Kin_Theta[trackbin][ptbin];
    TH1D* hDau_Kin_Eta[trackbin][ptbin];

    TH1D* hRawDPhi[trackbin][ptbin];
    TH1D* hRawDEta[trackbin][ptbin];

    TH1D* hN_ChgDAU[trackbin];
    TH2D* hN_ChgDAUvsJPT[trackbin];

    for(int wtrk = 1; wtrk<trackbin+1; wtrk++){

        hJet_Kin_Eta[wtrk-1]    = new TH1D(Form("hJet_Kin_Eta_%d",wtrk),Form("hJet_Kin_Eta_%d",wtrk)              ,bin1  , LB2   , HB2);
        hJet_Kin_Phi[wtrk-1]    = new TH1D(Form("hJet_Kin_Phi_%d",wtrk),Form("hJet_Kin_Phi_%d",wtrk)              ,bin1  , LB2pi , HB2pi);

        hN_ChgDAU[wtrk-1]           = new TH1D( Form("hN_ChgDAU_%d",wtrk)     ,  Form("hN_ChgDAU_%d",wtrk)     , bin120 , LB0 , HB120);
        hN_ChgDAUvsJPT[wtrk-1]      = new TH2D( Form("hN_ChgDAUvsJPT_%d",wtrk),  Form("hN_ChgDAUvsJPT_%d",wtrk), bin120 , LB0 , HB120 , bin100 , LB500 , HB4000);

        for(int wppt = 1; wppt<ptbin+1; wppt++){

            hBckrndShifted[wtrk-1][wppt-1]             = new TH2D(Form("hBckrndS_trk_%d_ppt_%d",wtrk,wppt) ,Form("hBckrndS_trk_%d_ppt_%d",wtrk,wppt) ,41,-(20*EtaBW)-(0.5*EtaBW),(20*EtaBW)+(0.5*EtaBW),33,-(8*PhiBW)-0.5*PhiBW,(24*PhiBW)+0.5*PhiBW);
            hSignalShifted[wtrk-1][wppt-1]             = new TH2D(Form("hSignalS_trk_%d_ppt_%d",wtrk,wppt) ,Form("hSignalS_trk_%d_ppt_%d",wtrk,wppt) ,41,-(20*EtaBW)-(0.5*EtaBW),(20*EtaBW)+(0.5*EtaBW),33,-(8*PhiBW)-0.5*PhiBW,(24*PhiBW)+0.5*PhiBW);

            hBckrndCentral[wtrk-1][wppt-1]             = new TH2D(Form("hBckrndC_trk_%d_ppt_%d",wtrk,wppt) ,Form("hBckrndC_trk_%d_ppt_%d",wtrk,wppt) ,40,-(20*EtaBW),(20*EtaBW),32,-(8*PhiBW),(24*PhiBW));
            hSignalCentral[wtrk-1][wppt-1]             = new TH2D(Form("hSignalC_trk_%d_ppt_%d",wtrk,wppt) ,Form("hSignalC_trk_%d_ppt_%d",wtrk,wppt) ,40,-(20*EtaBW),(20*EtaBW),32,-(8*PhiBW),(24*PhiBW));

            hEPDraw[wtrk-1][wppt-1]             = new TH2D(Form("hEPDraw_trk_%d_ppt_%d",wtrk,wppt) ,Form("hEPDraw_trk_%d_ppt_%d",wtrk,wppt) , EPD_xb   , EPD_xlo, EPD_xhi , EPD_yb      , EPD_ylo    , EPD_yhi);
            hEPDrawMinus[wtrk-1][wppt-1]        = new TH2D(Form("hEPDrawMinus_trk_%d_ppt_%d",wtrk,wppt) ,Form("hEPDrawMinus_trk_%d_ppt_%d",wtrk,wppt) , EPD_xb   , EPD_xlo, EPD_xhi , EPD_yb      , EPD_ylo    , EPD_yhi);
            hEPDrawPluss[wtrk-1][wppt-1]        = new TH2D(Form("hEPDrawPlus_trk_%d_ppt_%d",wtrk,wppt)  ,Form("hEPDrawPlus_trk_%d_ppt_%d",wtrk,wppt)  , EPD_xb   , EPD_xlo, EPD_xhi , EPD_yb      , EPD_ylo    , EPD_yhi);

            hRawDPhi[wtrk-1][wppt-1]            = new TH1D(Form("hRawDPhi_trk_%d_ppt_%d",wtrk,wppt) ,Form("hRawDPhi_trk_%d_ppt_%d",wtrk,wppt) , 150,  -0.6*TMath::Pi(), 1.6*TMath::Pi());
            hRawDEta[wtrk-1][wppt-1]            = new TH1D(Form("hRawDEta_trk_%d_ppt_%d",wtrk,wppt) ,Form("hRawDEta_trk_%d_ppt_%d",wtrk,wppt) , 180,  -30*EtaBW, 30*EtaBW);

            hDau_Kin_WRTJ_Phi[wtrk-1][wppt-1]   = new TH1D( Form("hRotPhi_trk_%d_ppt_%d",wtrk,wppt)   , Form("hRotPhi_trk_%d_ppt_%d",wtrk,wppt)    , bin2 , LB2pi , HB2pi);
            hDau_Kin_WRTJ_Pt[wtrk-1][wppt-1]    = new TH1D( Form("hRotPt_trk_%d_ppt_%d",wtrk,wppt)    , Form("hRotPt_trk_%d_ppt_%d",wtrk,wppt)     , bin2 , LB0   , HB20);
            hDau_Kin_WRTJ_Theta[wtrk-1][wppt-1] = new TH1D( Form("hRotTheta_trk_%d_ppt_%d",wtrk,wppt) , Form("hRotTheta_trk_%d_ppt_%d",wtrk,wppt)  , bin2 , LB0   , HB1);
            hDau_Kin_WRTJ_Eta[wtrk-1][wppt-1]   = new TH1D( Form("hRotEta_trk_%d_ppt_%d",wtrk,wppt)   , Form("hRotEta_trk_%d_ppt_%d",wtrk,wppt)    , bin_WRTJ_Eta , low_WRTJ_Eta_Bin   , high_WRTJ_Eta_Bin);

            hDau_Kin_Phi[wtrk-1][wppt-1]        = new TH1D( Form("hPhi_trk_%d_ppt_%d",wtrk,wppt)      ,  Form("hPhi_trk_%d_ppt_%d",wtrk,wppt)      , bin2 , LB2pi , HB2pi);
            hDau_Kin_Pt[wtrk-1][wppt-1]         = new TH1D( Form("hPt_trk_%d_ppt_%d",wtrk,wppt)       ,  Form("hPt_trk_%d_ppt_%d",wtrk,wppt)       , bin2 , LB0   , HB20);
            hDau_Kin_Theta[wtrk-1][wppt-1]      = new TH1D( Form("hTheta_trk_%d_ppt_%d",wtrk,wppt)    ,  Form("hTheta_trk_%d_ppt_%d",wtrk,wppt)    , bin2 , LB0   , HB15pi);
            hDau_Kin_Eta[wtrk-1][wppt-1]        = new TH1D( Form("hEta_trk_%d_ppt_%d",wtrk,wppt)      ,  Form("hEta_trk_%d_ppt_%d",wtrk,wppt)      , bin_Eta , low_Eta_Bin   , high_Eta_Bin);

            hdNdEvsPT[wtrk-1][wppt-1]      = new TH2D(Form("hdNdEvsPT_trk_%d_ppt_%d",wtrk,wppt) ,Form("dNdEvsPT_trk_%d_ppt_%d",wtrk,wppt)              , bin_WRTJ_Eta , low_WRTJ_Eta_Bin   , high_WRTJ_Eta_Bin, bin100 , LB0 , EPD_yhi);
            hdNdE[wtrk-1][wppt-1]          = new TH1D(Form("hdNdE_trk_%d_ppt_%d",wtrk,wppt) ,Form("hdNdE_trk_%d_ppt_%d",wtrk,wppt)                         , bin_WRTJ_Eta , low_WRTJ_Eta_Bin   , high_WRTJ_Eta_Bin);

        }
    }

    std::cout << "Starting event loop" << std::endl;
    std::cout << "Total Number of Files in this Job: " << fileList.size() << std::endl;
    for(int f = 0; f<fileList.size(); f++){

        fFile = TFile::Open(fileList.at(f).c_str(),"read");
        TTree *tree = (TTree*)fFile->Get("analyzer/trackTree");
        Init(tree);

        std::cout << "File " << f+1 << " out of " << fileList.size() << std::endl;
        Long64_t nbytes = 0, nb = 0;
        Long64_t nentries = fChain->GetEntriesFast();
        cout<<"Total Entries is:"<<endl;
        cout<< nentries <<endl;

        // ENTERING EVENT LOOP
        for (Long64_t ievent=0; ievent <nentries; ievent ++){

            Long64_t jevent = LoadTree(ievent);
            nb = fChain->GetEntry(ievent);   nbytes += nb;
            /*
               hpthat->Fill(genQScale, xs[qHatHist->FindBin(genQScale) - 1 ] / qHatHist->GetBinContent(qHatHist->FindBin(genQScale)) );
               if(genJetPt->size()==0) continue;
               if(genJetChargedMultiplicity->size()==0) continue;
               hleadingJetPt->Fill(genJetPt->at(0), xs[qHatHist->FindBin(genQScale) - 1 ] / qHatHist->GetBinContent(qHatHist->FindBin(genQScale)) );
               float genWeightPy = xs[qHatHist->FindBin(genQScale) - 1 ] / qHatHist->GetBinContent(qHatHist->FindBin(genQScale)); 
               */

            float genWeightPy = 1;
            //if(jetN > 200) continue;
            //if(jetN < 1) continue;

            std::vector<int> indicesG = F_genNix(genJetPt, genJetChargedMultiplicity, genJetEta, jetEtaCut, jetPtCut);

            int genNix = indicesG.size();
            int genN = genNix;

            //COMPARING GEN AND RECO JET NUMBERS AND DEMOGRAPHICS

            hEvent_Pass->Fill(1);
            //h2GenN_JetN->Fill(genNix, jetNix);
            //h1JetN_d_GenN->Fill(jetNix - genNix);
            h1GenN->Fill(genNix);
            //h1JetN->Fill(jetNix);

            //ENTERING JET LOOP


            for(int kjet=0; kjet<genN; kjet++){

                int ijet = indicesG[kjet];
                long int NNtrk = (genDau_chg->at(ijet)).size();
                float maxanypt = 0;
                float maxchgpt = 0;
                long int ChgMult  =0;

                hTotal_ijet->Fill(1);
                //SKIPPING JETS WHERE THE JET IS COMPLETELY NOT ACCOUNTED FOR IN GEN... A JET THAT SHOULDNT PASS THE CUTS
                if(genN <= kjet){
                    hExtra_ijet->Fill(1);
                    continue;
                }

                //hjetPt->Fill(genJetPt->at(ijet), xs[qHatHist->FindBin(genQScale) - 1 ] / qHatHist->GetBinContent(qHatHist->FindBin(genQScale)) );

                //leading any pt

                for(long int  XXtrk=0; XXtrk < NNtrk; XXtrk++ ){

                    //nothing to do with charge, just checking to make sure the size is not zero, there ARE daughters.

                    if(genN <= ijet){
                        hExtra_ijet->Fill(2);
                        continue;
                        }
                    if((genDau_chg->at(ijet)).size() == 0 ){
                        hExtra_ijet->Fill(3);
                        continue;
                        }

                    double LeadingJT = ptWRTJet( (double) (*genJetPt)[ijet] , (double) (*genJetEta)[ijet] , (double) (*genJetPhi)[ijet] , (double)(*genDau_pt)[ijet][XXtrk] , (double)(*genDau_eta)[ijet][XXtrk], (double)(*genDau_phi)[ijet][XXtrk]);

                    if(maxanypt < LeadingJT ){maxanypt = LeadingJT;}

                    //in case we got to the end and this one happens to be non charged, fill

                    if(XXtrk == NNtrk - 1){

                        hLeadJTany  ->Fill(maxanypt                 , genWeightPy);

                        hLeadJTany2d->Fill(maxanypt,NNtrk           , genWeightPy);

                    }

                    if((*genDau_chg)[ijet][XXtrk] == 0) continue;

                    ChgMult=ChgMult+1;

                }

                //leading charged pt
                for(long int  XXtrk=0; XXtrk < NNtrk; XXtrk++ ){
                    double LeadingJT = ptWRTJet( (double) (*genJetPt)[ijet] , (double) (*genJetEta)[ijet] , (double) (*genJetPhi)[ijet] , (double)(*genDau_pt)[ijet][XXtrk] , (double)(*genDau_eta)[ijet][XXtrk], (double)(*genDau_phi)[ijet][XXtrk]);
                    //if the reco charge is zero but we got to the end we cant skip the fill step!
                    if((*genDau_chg)[ijet][XXtrk] == 0 && XXtrk == NNtrk - 1){
                        hLeadJTchg  ->Fill(maxchgpt                 , genWeightPy);
                        hLeadJTchg2d->Fill(maxchgpt,ChgMult         , genWeightPy);
                    }
                    //ok now that we fill an entry no matter what, we can skip without fear. 
                    if((*genDau_chg)[ijet][XXtrk] == 0) continue;
                    if(maxchgpt < LeadingJT ){
                        maxchgpt = LeadingJT;
                    }
                    //if the reco is charged, non zero, and we get to the end, fill it!
                    if(XXtrk == NNtrk - 1){
                        hLeadJTchg  ->Fill(maxchgpt                 , genWeightPy);
                        hLeadJTchg2d->Fill(maxchgpt, ChgMult        , genWeightPy);
                    }
                }

                //curious if defining ChgMult is redundant 
                h2ChgMult_v_branch->Fill(ChgMult, (*genJetChargedMultiplicity)[ijet]); 
                h1ChgMult_d_branch->Fill(ChgMult - (*genJetChargedMultiplicity)[ijet]); 

                //this is how we know if we should proceed for this track bin loop or if we should skip it.
                //this is not going to over fill the above every time i loop through track bins. the track bin loop starts lower.
                // this makes sense, obviously advantages to skip events lower than the lowest limit.
                if(ChgMult < trackbinbounds[0]){
                    hMult_Fail->Fill(ChgMult, genWeightPy);
                    continue;
                }

                // filling distrivutions within track bins
                // ALSO VERY IMPORTANLTY changing the tkBool to 1 for this particular jet. This will be usefull later wen I create conditons for filling other historgams.
                int tkBool[trackbin] = {0};
                for(int i = 0; i < trackbin; i++){
                    if(ChgMult >= trackbinbounds[i]){
                        tkBool[i] = 1;
                        hJet_Pass           ->Fill(i);
                        hJet_Kin_Eta[i]     ->Fill((  *genJetEta)[ijet], genWeightPy);
                        hJet_Kin_Phi[i]     ->Fill((  *genJetPhi)[ijet], genWeightPy);
                        hN_ChgDAUvsJPT[i]   ->Fill(ChgMult, (  *genJetPt)[ijet]    , genWeightPy);
                        hN_ChgDAU[i]        ->Fill(ChgMult                      , genWeightPy);
                    }
                }


                int Ntrig[trackbin][ptbin] = {0};
                int A_ptBool[NNtrk][ptbin] = {0};
                // VERY IMPORTANT calculating daughter pt wrt to jet axis. so this needs to be 2d vector, for pt bin and for daughter index. in each case I create a true falsse for that daughter falling in to the specific pt bin. not this is NOT jet x daughter. its pt bin x daughtr
                for(int  A_trk=0; A_trk < NNtrk; A_trk++ ){
                    if((*genDau_chg)[ijet][A_trk] == 0) continue;

                    double jet_dau_pt    =  ptWRTJet( (double) (*genJetPt)[ijet] , (double) (*genJetEta)[ijet] , (double) (*genJetPhi)[ijet] , (double)(*genDau_pt)[ijet][A_trk] , (double)(*genDau_eta)[ijet][A_trk], (double)(*genDau_phi)[ijet][A_trk]);
                    //excluding outside outermost limits.
                    if( jet_dau_pt < ptbinmin || jet_dau_pt > ptbinmax) continue;
                    //loop through pt bins and fill the boolean array
                    for(int i = 0; i < ptbin; i++){
                        if(jet_dau_pt >= ptbinbounds_lo[i] && jet_dau_pt < ptbinbounds_hi[i]){
                            A_ptBool[A_trk][i] = 1;
                        }
                    }
                    //in the case where the total jet mult and the individual daughter pt is acceptable for this track bin and pt bin, we increase the Ntrig count.
                    for(int i = 0; i < trackbin; i++){
                        for(int j = 0; j < ptbin; j++){
                            if(tkBool[i] + A_ptBool[A_trk][j] == 2){
                                Ntrig[i][j] += 1;
                            }
                        }
                    }
                }//here ends the boolean array creations.

                //continuation of main loops. here is where the 2D Corr plots are created using the above booleans and 
                for(int  A_trk=0; A_trk < (NNtrk-1); A_trk++ ){
                    //this should be redundant if it passes the bools above? i guess it helps skip daughters faster. maybe i can reindex and run through the daughters quickly by aranging all the charged dauhghter sat the front.
                    if((*genDau_chg)[ijet][A_trk] == 0) continue;

                    double jet_dau_pt    =  ptWRTJet( (double) (*genJetPt)[ijet] , (double) (*genJetEta)[ijet] , (double) (*genJetPhi)[ijet] , (double)(*genDau_pt)[ijet][A_trk] , (double)(*genDau_eta)[ijet][A_trk], (double)(*genDau_phi)[ijet][A_trk]);
                    if( jet_dau_pt < ptbinmin || jet_dau_pt > ptbinmax) continue;

                    int c1 = 0;
                    if((*genDau_chg)[ijet][A_trk] < 0) c1 = -1;
                    if((*genDau_chg)[ijet][A_trk] > 0) c1 =  1;

                    double jet_dau_theta = 2*ATan(Exp(-(etaWRTJet((double)(*genJetPt)[ijet],(double)(*genJetEta)[ijet],(double)(*genJetPhi)[ijet],(double)(*genDau_pt)[ijet][A_trk],(double)(*genDau_eta)[ijet][A_trk],(double)(*genDau_phi)[ijet][A_trk]))));
                    double jet_dau_eta   = etaWRTJet( (double) (*genJetPt)[ijet] , (double) (*genJetEta)[ijet] , (double) (*genJetPhi)[ijet] , (double)(*genDau_pt)[ijet][A_trk] , (double)(*genDau_eta)[ijet][A_trk], (double)(*genDau_phi)[ijet][A_trk]);
                    double jet_dau_phi   = phiWRTJet( (double) (*genJetPt)[ijet] , (double) (*genJetEta)[ijet] , (double) (*genJetPhi)[ijet] , (double)(*genDau_pt)[ijet][A_trk] , (double)(*genDau_eta)[ijet][A_trk], (double)(*genDau_phi)[ijet][A_trk]);

                    //                    if(jet_dau_eta <  0.86) continue;
                    //                    if(jet_dau_eta >  3.90) continue;

                    double beam_dau_eta   = (*genDau_eta)[ijet][A_trk];
                    double beam_dau_phi   = (*genDau_phi)[ijet][A_trk];
                    double beam_dau_pt    = (*genDau_pt)[ijet][A_trk];

                    //A_ptBool[ptbin]    = {0};
                    //for(int i = 0; i < ptbin; i++){
                    //if(jet_dau_pt >= ptbinbounds_lo[i] && jet_dau_pt < ptbinbounds_hi[i]){
                    //A_ptBool[i] = 1;
                    //}
                    //}       

                    for(int i = 0; i < trackbin; i++){
                        for(int j = 0; j < ptbin; j++){
                            if(tkBool[i] + A_ptBool[A_trk][j] == 2){

                                hEPDraw[i][j]->Fill(jet_dau_eta, jet_dau_phi, genWeightPy/Ntrig[i][j]);

                                if(c1 ==  1){ hEPDrawPluss[i][j]->Fill(jet_dau_eta, jet_dau_phi, genWeightPy/Ntrig[i][j]);}
                                if(c1 == -1){ hEPDrawMinus[i][j]->Fill(jet_dau_eta, jet_dau_phi, genWeightPy/Ntrig[i][j]);}

                                hdNdEvsPT               [i][j]->Fill(jet_dau_eta, jet_dau_pt, genWeightPy);
                                hdNdE                   [i][j]->Fill(jet_dau_eta, genWeightPy);

                                hDau_Kin_WRTJ_Phi       [i][j]->Fill(jet_dau_phi    , genWeightPy);
                                hDau_Kin_WRTJ_Pt        [i][j]->Fill(jet_dau_pt     , genWeightPy);
                                hDau_Kin_WRTJ_Theta     [i][j]->Fill(jet_dau_theta  , genWeightPy);
                                hDau_Kin_WRTJ_Eta       [i][j]->Fill(jet_dau_eta    , genWeightPy);
                                hDau_Kin_Phi            [i][j]->Fill(beam_dau_phi   , genWeightPy);
                                hDau_Kin_Pt             [i][j]->Fill(beam_dau_pt    , genWeightPy);
                                hDau_Kin_Eta            [i][j]->Fill(beam_dau_eta   , genWeightPy);
                            }
                        }
                    }

                    for(long int T_trk=A_trk+1; T_trk< NNtrk; T_trk++ ){

                        if((*genDau_chg)[ijet][T_trk] == 0) continue;

                        double T_jet_dau_pt = ptWRTJet(  (double) (*genJetPt)[ijet] , (double) (*genJetEta)[ijet] , (double) (*genJetPhi)[ijet] , (double)(*genDau_pt)[ijet][T_trk] , (double)(*genDau_eta)[ijet][T_trk], (double)(*genDau_phi)[ijet][T_trk] );
                        if( T_jet_dau_pt < ptbinmin || T_jet_dau_pt > ptbinmax) continue;

                        //int c2 = 0;
                        //if((*dau_chg)[ijet][T_trk] < 0) c2 = -1;
                        //if((*dau_chg)[ijet][T_trk] > 0) c2 =  1;              

                        double T_jet_dau_eta = etaWRTJet( (double) (*genJetPt)[ijet] , (double) (*genJetEta)[ijet] , (double) (*genJetPhi)[ijet] , (double)(*genDau_pt)[ijet][T_trk] , (double)(*genDau_eta)[ijet][T_trk], (double)(*genDau_phi)[ijet][T_trk]);
                        double T_jet_dau_phi = phiWRTJet( (double) (*genJetPt)[ijet] , (double) (*genJetEta)[ijet] , (double) (*genJetPhi)[ijet] , (double)(*genDau_pt)[ijet][T_trk] , (double)(*genDau_eta)[ijet][T_trk], (double)(*genDau_phi)[ijet][T_trk]);

                        double deltaEta = (jet_dau_eta - T_jet_dau_eta);
                        double deltaPhi = (TMath::ACos(TMath::Cos(jet_dau_phi - T_jet_dau_phi)));

                        for(        int i = 0; i < trackbin; i++){
                            for(    int j = 0; j < ptbin;    j++){
                                if(tkBool[i] + A_ptBool[A_trk][j] + A_ptBool[T_trk][j] == 3){

                                    hPairs->Fill(i,j);

                                    hRawDPhi[i][j]->Fill(deltaPhi, 1./Ntrig[i][j]);
                                    hRawDEta[i][j]->Fill(deltaEta, 1./Ntrig[i][j]);

                                    hSignalShifted[i][j]->Fill(deltaEta, deltaPhi, 1./Ntrig[i][j]);
                                    hSignalShifted[i][j]->Fill(-deltaEta, deltaPhi, 1./Ntrig[i][j]);
                                    hSignalShifted[i][j]->Fill(deltaEta, -deltaPhi, 1./Ntrig[i][j]);
                                    hSignalShifted[i][j]->Fill(-deltaEta, -deltaPhi, 1./Ntrig[i][j]); 

                                    hSignalShifted[i][j]->Fill( deltaEta,2*TMath::Pi() - deltaPhi, 1./Ntrig[i][j]);
                                    hSignalShifted[i][j]->Fill(-deltaEta,2*TMath::Pi() - deltaPhi, 1./Ntrig[i][j]);

                                    hSignalCentral[i][j]->Fill(deltaEta, deltaPhi, 1./Ntrig[i][j]);
                                    hSignalCentral[i][j]->Fill(-deltaEta, deltaPhi, 1./Ntrig[i][j]);
                                    hSignalCentral[i][j]->Fill(deltaEta, -deltaPhi, 1./Ntrig[i][j]);
                                    hSignalCentral[i][j]->Fill(-deltaEta, -deltaPhi, 1./Ntrig[i][j]); 

                                    hSignalCentral[i][j]->Fill( deltaEta,2*TMath::Pi() - deltaPhi, 1./Ntrig[i][j]);
                                    hSignalCentral[i][j]->Fill(-deltaEta,2*TMath::Pi() - deltaPhi, 1./Ntrig[i][j]);
                                }
                            }
                        }
                    }
                }
            }
        }
        fFile->Close();
    }
    int backMult =10;
    for(int wtrk = 1; wtrk < trackbin+1; wtrk++){
        for(int wppt = 1; wppt < ptbin+1; wppt++){

            long int NENT =  hPairs->GetBinContent(wtrk, wppt);
            long int XENT =  ((1+floor(sqrt(1+(4*2*backMult*NENT))))/2) ;
            //Nent is the number of pairs in the signal which we will try to 10x
            //Xent is the number of pseudoparticles requried such that when we build the pairs nCp = Xent CHOOSE 2 will give us 10 times as many pairs as we have in the signal histogrm.

            float A_ETA[XENT] = {0};
            float A_PHI[XENT] = {0};

            for(int x = 0; x<XENT; x++){
                double WEta1, WPhi1;//making the pseudoparticles
                hEPDraw[wtrk-1][wppt-1]->GetRandom2(WEta1, WPhi1);
                A_ETA[x] = WEta1;
                A_PHI[x] = WPhi1;
            }

            for(long int i = 0; i < (XENT-1); i++){
                for(long int j = (i+1); j < XENT; j++){

                    double WdeltaEta = (A_ETA[i]-A_ETA[j]);
                    double WdeltaPhi = (TMath::ACos(TMath::Cos(A_PHI[i]-A_PHI[j])));

                    hBckrndShifted[wtrk-1][wppt-1]->Fill(WdeltaEta, WdeltaPhi, 1./XENT);
                    hBckrndShifted[wtrk-1][wppt-1]->Fill(-WdeltaEta, WdeltaPhi, 1./XENT);
                    hBckrndShifted[wtrk-1][wppt-1]->Fill(WdeltaEta, -WdeltaPhi, 1./XENT);
                    hBckrndShifted[wtrk-1][wppt-1]->Fill(-WdeltaEta, -WdeltaPhi, 1./XENT);

                    hBckrndShifted[wtrk-1][wppt-1]->Fill(WdeltaEta, 2*TMath::Pi() - WdeltaPhi, 1./XENT);
                    hBckrndShifted[wtrk-1][wppt-1]->Fill(-WdeltaEta,2*TMath::Pi() - WdeltaPhi, 1./XENT);

                    hBckrndCentral[wtrk-1][wppt-1]->Fill(WdeltaEta, WdeltaPhi, 1./XENT);
                    hBckrndCentral[wtrk-1][wppt-1]->Fill(-WdeltaEta, WdeltaPhi, 1./XENT);
                    hBckrndCentral[wtrk-1][wppt-1]->Fill(WdeltaEta, -WdeltaPhi, 1./XENT);
                    hBckrndCentral[wtrk-1][wppt-1]->Fill(-WdeltaEta, -WdeltaPhi, 1./XENT);

                    hBckrndCentral[wtrk-1][wppt-1]->Fill(WdeltaEta, 2*TMath::Pi() - WdeltaPhi, 1./XENT);
                    hBckrndCentral[wtrk-1][wppt-1]->Fill(-WdeltaEta,2*TMath::Pi() - WdeltaPhi, 1./XENT);

                }
            }
        }
    }

    TFile* fS_tempA = new TFile(Form("src/Dists/MC_Full470_F/generator/Primary/Primary_job%d.root",job), "recreate");

    for(int wtrk =1; wtrk <trackbin+1; wtrk++){
        for(int wppt =1; wppt <ptbin+1; wppt++){

            hSignalShifted             [wtrk-1][wppt-1]->Write(Form("hSigS_%d_%d_to_%d",trackbinbounds[wtrk-1],(int)(10*ptbinbounds_lo[wppt-1]),(int)(10*ptbinbounds_hi[wppt-1])));
            hBckrndShifted             [wtrk-1][wppt-1]->Write(Form("hBckS_%d_%d_to_%d",trackbinbounds[wtrk-1],(int)(10*ptbinbounds_lo[wppt-1]),(int)(10*ptbinbounds_hi[wppt-1])));

            hSignalCentral             [wtrk-1][wppt-1]->Write(Form("hSigC_%d_%d_to_%d",trackbinbounds[wtrk-1],(int)(10*ptbinbounds_lo[wppt-1]),(int)(10*ptbinbounds_hi[wppt-1])));
            hBckrndCentral             [wtrk-1][wppt-1]->Write(Form("hBckC_%d_%d_to_%d",trackbinbounds[wtrk-1],(int)(10*ptbinbounds_lo[wppt-1]),(int)(10*ptbinbounds_hi[wppt-1])));

            hEPDraw                    [wtrk-1][wppt-1]->Write(Form("hEPD_%d_%d_to_%d",trackbinbounds[wtrk-1],(int)(10*ptbinbounds_lo[wppt-1]),(int)(10*ptbinbounds_hi[wppt-1])));
        }
    }
    hGenExtra_pt_ChgMult      ->Write();
    hGenExtra_pt_Eta          ->Write();
    hGenExtra_Eta_ChgMult     ->Write();
    hRecExtra_pt_ChgMult      ->Write();
    hRecExtra_pt_Eta          ->Write();
    hRecExtra_Eta_ChgMult     ->Write();

    hEvent_Pass   ->Write();
    hJet_Pass     ->Write();

    hEvent_Fail   ->Write();
    hJet_Fail     ->Write();
    hMult_Fail    ->Write();

    hLeadJTany2d  ->Write();
    hLeadJTany    ->Write();

    hLeadJTchg2d  ->Write();
    hLeadJTchg    ->Write();

    hTotal_ijet         ->Write();
    hExtra_ijet         ->Write();

    fS_tempA->Close();

    TFile* fS_tempB = new TFile(Form("src/Dists/MC_Full470_F/generator/Secondary/Secondary_job%d.root",job), "recreate");
    for(int wtrk =1; wtrk <trackbin+1; wtrk++){

        hJet_Kin_Eta            [wtrk-1]->Write(Form("hJet_Kin_Eta_%d",wtrk));
        hJet_Kin_Phi            [wtrk-1]->Write(Form("hJet_Kin_Phi_%d",wtrk));

        hN_ChgDAU               [wtrk-1]->Write(Form("hN_ChgDAU_%d",wtrk));
        hN_ChgDAUvsJPT          [wtrk-1]->Write(Form("hN_ChgDAUvsJPT_%d",wtrk));

        for(int wppt =1; wppt <ptbin+1; wppt++){
            hEPDrawMinus        [wtrk-1][wppt-1]->Write(Form("hEPDrawMinus_%d_%d",wtrk,wppt));
            hEPDrawPluss        [wtrk-1][wppt-1]->Write(Form("hEPDrawPluss_%d_%d",wtrk,wppt));

            hRawDEta            [wtrk-1][wppt-1]->Write(Form("hRawDEta_%d_%d",wtrk,wppt));
            hRawDPhi            [wtrk-1][wppt-1]->Write(Form("hRawDPhi_%d_%d",wtrk,wppt));

            hDau_Kin_WRTJ_Theta [wtrk-1][wppt-1]->Write(Form("hDau_Kin_WRTJ_Theta_%d_%d",wtrk,wppt));
            hDau_Kin_WRTJ_Eta   [wtrk-1][wppt-1]->Write(Form("hDau_Kin_WRTJ_Eta_%d_%d",wtrk,wppt));
            hDau_Kin_WRTJ_Phi   [wtrk-1][wppt-1]->Write(Form("hDau_Kin_WRTJ_Phi_%d_%d",wtrk,wppt));
            hDau_Kin_WRTJ_Pt    [wtrk-1][wppt-1]->Write(Form("hDau_Kin_WRTJ_Pt_%d_%d",wtrk,wppt));

            hDau_Kin_Theta      [wtrk-1][wppt-1]->Write(Form("hDau_Kin_Theta_%d_%d",wtrk,wppt));
            hDau_Kin_Eta        [wtrk-1][wppt-1]->Write(Form("hDau_Kin_Eta_%d_%d",wtrk,wppt));
            hDau_Kin_Phi        [wtrk-1][wppt-1]->Write(Form("hDau_Kin_Phi_%d_%d",wtrk,wppt));
            hDau_Kin_Pt         [wtrk-1][wppt-1]->Write(Form("hDau_Kin_Pt_%d_%d",wtrk,wppt));

            hdNdEvsPT           [wtrk-1][wppt-1]->Write(Form("hdNdEvsPT_%d_%d",wtrk, wppt));
            hdNdE               [wtrk-1][wppt-1]->Write(Form("hdNdE_%d_%d",wtrk, wppt));

        }
    }

    fS_tempB->Close();


}


//Code enters execution here
int main(int argc, const char* argv[])
{
    if(argc != 4)
    {
        std::cout << "Usage: Z_mumu_Channel <fileList> <jobNumber> <nJobs>" << std::endl;
        return 1;
    }


    //read input parameters
    std::string fList = argv[1];
    std::string buffer;
    std::vector<std::string> listOfFiles;
    std::ifstream inFile(fList.data());

    int job = (int)std::atoi(argv[2]);
    int nJobs = (int)std::atoi(argv[3]);


    //read the file list and spit it into a vector of strings based on how the parallelization is to be done
    //each vector is a separate subset of the fileList based on the job number
    if(!inFile.is_open())
    {
        std::cout << "Error opening jet file. Exiting." <<std::endl;
        return 1;
    }
    else
    {
        int line = 0;
        while(true)
        {
            inFile >> buffer;
            if(inFile.eof()) break;
            if( line%nJobs == job) listOfFiles.push_back(buffer);
            line++;
        }
    }

    //create the MyClass Object
    MyClass m = MyClass(listOfFiles);
    m.Loop(job);

    return 0;
}

