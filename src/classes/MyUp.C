
#define MyClass_cxx

#include "include/MyClass.h"
#include "include/Timer.h"
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

const int   trackbin                    = 7; 
//const int   ptbin		                = 15;
const int   ptbin		                = 5;
const float ptmin                       = 0.0;
const float ptmax                       = 4.0;
const int   trackbinbounds[trackbin]    = {0,85,86,87,88,89,90}; 
//const float ptbinbounds_lo[ptbin]       = {0.0, 0.3, 0.4, 0.5, 0.6, 1.0, 0.0, 0.0, 0.0, 0.3, 0.3, 0.5, 0.5, 1.0, 1.0}; 
const float ptbinbounds_lo[ptbin]       = {0.0,  0.4, 0.5, 0.6, 1.0}; 
//const float ptbinbounds_hi[ptbin]       = {3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 1.0, 2.0, 4.0, 2.0, 4.0, 2.0, 4.0, 2.0, 4.0}; 
const float ptbinbounds_hi[ptbin]       = {20.0, 3.0, 3.0, 3.0, 3.0}; 
const float ptbinmin                    = 0.39;
const float ptbinmax                    = 3.01;


const double        PI = 3.14159265359;
const float         EtaBW = 0.3;
const float         PhiBW = TMath::Pi()/16;
//const int           DeltaPBin =32;
//const float         etabound  =6;
//const int           DeltaEBin =40;//(2*etabound)/(0.3);
const int bin_WRTJ_Eta        = 150;
const int low_WRTJ_Eta_Bin  =  0;
const int high_WRTJ_Eta_Bin = 10;
const int bin_Eta             = 60; 
const int low_Eta_Bin       = -3; 
const int high_Eta_Bin      = 3;

const int bin1 = 30;
const int bin2 = 100;
const int bin3 = 50;
const int LB0 = 0;
const double LB2pi = -2*PI;
const int LB2 = -15;
const int LB1 = -1;
const double HB2pi = 2*PI;
const int HB2 = 2;
const int HB1 = 1;
const double HB15pi = 1.5*PI;
const int HB20 = 20;

const int EPD_xb  = 150;
const int EPD_yb  = 120;
const int EPD_xhi = 10;
const int EPD_xlo = 0;
const int EPD_yhi = 4;
const int EPD_ylo = -4;

//---------------------------------------------------------------------CUTS
const float EtaCut      = 0.0;
const float jetEtaCut   = 1.6;
const float jetPtCut    = 500.0;


//Bool folds

bool F_eventpass(
        std::vector< float > *jetPt,
        int jetN,
        float jetPtCut
        ){
    float max_element =0;
    if(jetN > 200) return false;
    if(jetN < 1) return false;
    for(int i=0; i<jetN; i++){
        float this_element = (*jetPt)[i];
        if(this_element > max_element){
            max_element=(*jetPt)[i];
        }
    }
    if(max_element < jetPtCut)        return false;
    return true;
}

bool F_jetpass(
        std::vector< float > * jetEta,
        std::vector< float > * jetPt,
        int     ijet,
        float   jetPtCut
        ){
    if(fabs( (*jetEta)[ijet]  ) >jetEtaCut)   return false;
    if((*jetPt)[ijet]          <jetPtCut)    return false;
    return true;
}

void MyClass::Loop(int job){
    //timing information
    Timer timer = Timer();
    timer.Start();
    timer.StartSplit("Start Up");

    std::cout << "Made it 1" << std::endl;


    TH1::SetDefaultSumw2(kTRUE);
    TH2::SetDefaultSumw2(kTRUE);

    //Output and bin variable init
    const int bin120 = 120;
    const int bin150 = 150;
    const int bin100 = 100;
    const int HB120 = 120;
    const int HB4000 = 4000;
    const int LB500 = 500;
    const int HB15 = 15;
    const int HB4 = 4;
    const int bin1600 = 1600;
    const int bin5 = 5;
    const int LB1 = 1;
    const int HB110 = 110;
    const int LB50 = 50;
    const int LB470 = 470;
    const int bin200 = 2000;
    const int bin2000 = 2000;
    const int bin60 = 60;

    //Initializing Histograms

    TH1D* hTotal_ijet   = new TH1D("hTotal_ijet","hTotal_ijet", 1,LB0,1);
    TH1D* hExtra_ijet   = new TH1D("hExtra_ijet","hExtra_ijet", 1,LB0,1);


    TH1D* hEvent_Pass   = new TH1D("hEvent_Pass","hEvent_Pass", 1,LB0,1);
    TH1D* hEvent_Fail   = new TH1D("hEvent_Fail","hEvent_Fail", 1,LB0,1);
    TH1D* hJet_Pass     = new TH1D("hJet_Pass"  ,"hJet_Pass"  , trackbin,LB0,trackbin);
    TH1D* hJet_Fail     = new TH1D("hJet_Fail"  ,"hJet_Fail"  , bin120,LB0,bin120);
    TH1D* hMult_Fail    = new TH1D("hMult_Fail" ,"hMult_Fail" , bin120,LB0,bin120);


    TH1D* hLeadJTany      = new TH1D("hLeadJTany" ,"hLeadJTany"     , bin2 , LB0   , HB20);
    TH2D* hLeadJTany2d    = new TH2D("hLeadJTany2d" ,"hLeadJTany2d" , bin2 , LB0   , HB20, bin150, LB0, bin150);

    TH1D* hLeadJTchg      = new TH1D("hLeadJTchg" ,"hLeadJTchg"     , bin2 , LB0   , HB20);
    TH2D* hLeadJTchg2d    = new TH2D("hLeadJTchg2d" ,"hLeadJTchg2d" , bin2 , LB0   , HB20, bin120, LB0, bin120);

    TH2D* h2GenN_JetN         = new TH2D("h2GenN_JetN"       ,"h2GenN_JetN"       , bin5,LB0,bin5, bin5,LB0,bin5);
    TH2D* h2GenN_JetN_W       = new TH2D("h2GenN_JetN_W"       ,"h2GenN_JetN_W"       , bin5,LB0,bin5, bin5,LB0,bin5);

    TH1D* h1JetN_d_GenN       = new TH1D("h1JetN_d_GenN"    ,"h1JetN_d_GenN"    , bin5, -LB1,bin5);
    TH1D* h1JetN_d_GenN_W     = new TH1D("h1JetN_d_GenN_W"    ,"h1JetN_d_GenN_W"    , bin5,-LB1,bin5);

    TH1D* h1GenN              = new TH1D("h1GenN"    ,"h1GenN"    , bin5,LB0,bin5);
    TH1D* h1GenN_W            = new TH1D("h1GenN_W"    ,"h1GenN_W"    , bin5,LB0,bin5);

    TH1D* h1JetN              = new TH1D("h1JetN"    ,"h1JetN"    , bin5,LB0,bin5);
    TH1D* h1JetN_W            = new TH1D("h1JetN_W"    ,"h1JetN_W"    , bin5,LB0,bin5);

    TH2D* h2_J1M_J2M_any      = new TH2D("h2_J1M_J2M_any"       ,"h2_J1M_J2M_any"       , bin120,LB0,bin120, bin120,LB0,bin120);
    TH2D* h2_J1M_J2M_chg      = new TH2D("h2_J1M_J2M_chg"       ,"h2_J1M_J2M_chg"       , bin120,LB0,bin120, bin120,LB0,bin120);
    TH2D* h2_gen_J1M_J2M_any  = new TH2D("h2_gen_J1M_J2M_any"   ,"h2_gen_J1M_J2M_any"   , bin120,LB0,bin120, bin120,LB0,bin120);
    TH2D* h2_gen_J1M_J2M_chg  = new TH2D("h2_gen_J1M_J2M_chg"   ,"h2_gen_J1M_J2M_chg"   , bin120,LB0,bin120, bin120,LB0,bin120);
    TH1D* h1_single_J1M_any   = new TH1D("h1_single_J1M_any"    ,"h1_single_J1M_any"    , bin120,LB0,bin120);
    TH1D* h1_single_J1M_chg   = new TH1D("h1_single_J1M_chg"    ,"h1_single_J1M_chg"    , bin120,LB0,bin120);
    TH1D* h1_single_J1pt      = new TH1D("h1_single_J1pt"       ,"h1_single_J1pt"       , bin200 , LB0, bin2000);
    TH1D* h1_double_J1M_any   = new TH1D("h1_double_J1M_any"    ,"h1_double_J1M_any"    , bin120,LB0,bin120);
    TH1D* h1_double_J1M_chg   = new TH1D("h1_double_J1M_chg"    ,"h1_double_J1M_chg"    , bin120,LB0,bin120);
    TH1D* h1_double_J1pt      = new TH1D("h1_double_J1pt"       ,"h1_double_J1pt"       , bin200 , LB0, bin2000);
    TH1D* h1_double_J2M_any   = new TH1D("h1_double_J2M_any"    ,"h1_double_J2M_any"    , bin120,LB0,bin120);
    TH1D* h1_double_J2M_chg   = new TH1D("h1_double_J2M_chg"    ,"h1_double_J2M_chg"    , bin120,LB0,bin120);
    TH1D* h1_double_J2pt      = new TH1D("h1_double_J2pt"       ,"h1_double_J2pt"       , bin200 , LB0, bin2000);

    TH1D* h1_J1E_d_J2E        = new TH1D("h1_J1E_d_J2E"         ,"h1_J1E_d_J2E"         , bin100 , -HB2, HB2);
    TH1D* h1_gJ1E_d_gJ2E      = new TH1D("h1_gJ1E_d_gJ2E"       ,"h1_gJ1E_d_gJ2E"       , bin100 , -HB2, HB2);

    TH1D* h1_J1P_d_J2P        = new TH1D("h1_J1P_d_J2P"         ,"h1_J1P_d_J2P"         , bin100 , -HB2pi, HB2pi);
    TH1D* h1_gJ1P_d_gJ2P      = new TH1D("h1_gJ1P_d_gJ2P"       ,"h1_gJ1P_d_gJ2P"       , bin100 , -HB2pi, HB2pi);

    TH2D* hGenExtra_pt_ChgMult      = new TH2D("hGenExtra_pt_ChgMult"   ,"hGenExtra_pt_ChgMult"   , bin120,LB470,bin1600, bin60,LB50,HB110);
    TH2D* hGenExtra_pt_Eta          = new TH2D("hGenExtra_pt_Eta"       ,"hGenExtra_pt_Eta"       , bin120,LB470,bin1600, bin3,-HB4,HB4);
    TH2D* hGenExtra_Eta_ChgMult     = new TH2D("hGenExtra_Eta_ChgMult"  ,"hGenExtra_Eta_ChgMult"  , bin3,-HB4,HB4,        bin60,LB50,HB110);
    TH2D* hRecExtra_pt_ChgMult      = new TH2D("hRecExtra_pt_ChgMult"   ,"hRecExtra_pt_ChgMult"   , bin120,LB470,bin1600, bin60,LB50,HB110);
    TH2D* hRecExtra_pt_Eta          = new TH2D("hRecExtra_pt_Eta"       ,"hRecExtra_pt_Eta"       , bin120,LB470,bin1600, bin3,-HB4,HB4);
    TH2D* hRecExtra_Eta_ChgMult     = new TH2D("hRecExtra_Eta_ChgMult"  ,"hRecExtra_Eta_ChgMult"  , bin3,-HB4,HB4,        bin60,LB50,HB110);

    TH1D* hJet_Kin_Eta[trackbin];  
    TH1D* hJet_Kin_Phi[trackbin];  

    TH2D* hdNdEvsPT[trackbin][ptbin];
    TH1D* hdNdE[trackbin][ptbin];

    TH2D* hEPDraw[trackbin][ptbin];

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

    //Defining Histograms

    for(int wtrk = 1; wtrk<trackbin+1; wtrk++){

        hJet_Kin_Eta[wtrk-1]    = new TH1D(Form("hJet_Kin_Eta_%d",wtrk),Form("hJet_Kin_Eta_%d",wtrk)              ,bin1  , LB2   , HB2);
        hJet_Kin_Phi[wtrk-1]    = new TH1D(Form("hJet_Kin_Phi_%d",wtrk),Form("hJet_Kin_Phi_%d",wtrk)              ,bin1  , LB2pi , HB2pi);

        hN_ChgDAU[wtrk-1]           = new TH1D( Form("hN_ChgDAU_%d",wtrk)     ,  Form("hN_ChgDAU_%d",wtrk)     , bin120 , LB0 , HB120);
        hN_ChgDAUvsJPT[wtrk-1]      = new TH2D( Form("hN_ChgDAUvsJPT_%d",wtrk),  Form("hN_ChgDAUvsJPT_%d",wtrk), bin120 , LB0 , HB120 , bin100 , LB500 , HB4000);

        for(int wppt = 1; wppt<ptbin+1; wppt++){

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



    std::cout << "Made it 2" << std::endl;
    /*
    // uncomment for pthat*****************
    std::cout << "Made it 2.1a" << std::endl;

    const int nqHats = 8;
    std::cout << "Made it 2.1b" << std::endl;

    const float qHatBoundaries[nqHats+1] = {470, 600, 800, 1000, 1400, 1800, 2400, 3200, 13000};
    std::cout << "Made it 2.1c" << std::endl;

    //units in pb

    float xs[nqHats] = {552.1, 156.5, 26.28, 7.47, 0.6484, 0.08743, 0.005236, 0.0001357};
    std::cout << "Made it 2.1d" << std::endl;

    TH1D * qHatHist = new TH1D("qhatHit",";;qHat",nqHats,qHatBoundaries);

    std::cout << "Made it 2.1" << std::endl;

    for(unsigned int ff = 0; ff < fileList.size(); ff++){  
    for(unsigned int ff = 0; ff < fileList.size(); ff++){  
std::cout << "Made it 2.2" << std::endl;
    fFile = TFile::Open(fileList.at(ff).c_str(),"read");

    std::cout << "Made it 2.3" << std::endl;

    TTree *tree = (TTree*)fFile->Get("analyzer/trackTree");
    Init(tree);

    for( int i = 0; i < tree->GetEntries(); i++){
    tree->GetEntry(i);
    qHatHist->Fill(genQScale);    
    }

    std::cout << "Made it 2.4" << std::endl;
    fFile->Close();
    } 

    const int i40 = 40;
    const int i470 = 470;
    const int i3500 = 3500;
    const int i50 = 50;
    const int i500 = 500;
    const int i1500 = 1500;
    const int i100 = 100;

    TH1D * hpthat = new TH1D("pthat",";#hat{q};#sigma (pb)",i40,i470,i3500);
    TH1D * hleadingJetPt = new TH1D("leadingJetPt",";Leading p_{T}^{gen};#sigma (pb)",i50,i500,i1500);
    TH1D * hjetPt = new TH1D("JetPt",";p_{T}^{gen};#sigma (pb)",i50,i100,i1500);



    // uncomment for pthat*****************
     */

    std::cout << "Made it 3" << std::endl;

    // Main track loops
    std::cout << "Starting event loop" << std::endl; 
    std::cout << "Total Number of Files in this Job: " << fileList.size() << std::endl;
    for(int f = 0; f<fileList.size(); f++){


        timer.StartSplit("Opening Files");
        fFile = TFile::Open(fileList.at(f).c_str(),"read");
        TTree *tree = (TTree*)fFile->Get("analyzer/trackTree"); 
        Init(tree);


        std::cout << "File " << f+1 << " out of " << fileList.size() << std::endl;
        Long64_t nbytes = 0, nb = 0;
        Long64_t nentries = fChain->GetEntriesFast();
        cout<<"Total Entries is:"<<endl;
        cout<< nentries <<endl;


        std::cout << "Made it 4" << std::endl;


        //Main loops
        for (Long64_t ievent=0; ievent <nentries; ievent ++){


            timer.StartSplit("Loading Events");

            Long64_t jevent = LoadTree(ievent);
            nb = fChain->GetEntry(ievent);   nbytes += nb;
            //if(ievent%1000==0) cout<< ievent << "/" << nentries  <<endl;

            timer.StartSplit("Event Selection");
            /*
               hpthat->Fill(genQScale, xs[qHatHist->FindBin(genQScale) - 1 ] / qHatHist->GetBinContent(qHatHist->FindBin(genQScale)) );
               if(genJetPt->size()==0) continue;
               if(genJetChargedMultiplicity->size()==0) continue;
               hleadingJetPt->Fill(genJetPt->at(0), xs[qHatHist->FindBin(genQScale) - 1 ] / qHatHist->GetBinContent(qHatHist->FindBin(genQScale)) );
               float genWeightPy = xs[qHatHist->FindBin(genQScale) - 1 ] / qHatHist->GetBinContent(qHatHist->FindBin(genQScale)); 
             */
            float genWeightPy = 1;
            if(jetN > 200) continue;
            if(jetN < 1) continue;

            // genNix
            std::vector<int> indicesG(genJetPt->size());
            std::iota(indicesG.begin(), indicesG.end(), 0);

            std::sort(indicesG.begin(), indicesG.end(),
                    [&](int A, int B) -> bool {
                    return (*genJetChargedMultiplicity)[A] > (*genJetChargedMultiplicity)[B];
                    });

            indicesG.erase(std::remove_if(
                        indicesG.begin(), indicesG.end(),
                        [&](const int x) { 
                        return (*genJetEta)[x] > 1.6; // put your condition here
                        }), indicesG.end());

            indicesG.erase(std::remove_if(
                        indicesG.begin(), indicesG.end(),
                        [&](const int x) { 
                        return (*genJetPt)[x] < 500; // put your condition here
                        }), indicesG.end());

            // jetNix
            std::vector<int> indicesR(jetPt->size());
            std::iota(indicesR.begin(), indicesR.end(), 0);

            std::sort(indicesR.begin(), indicesR.end(),
                    [&](int A, int B) -> bool {
                    return (*chargedMultiplicity)[A] > (*chargedMultiplicity)[B];
                    });

            indicesR.erase(std::remove_if(
                        indicesR.begin(), indicesR.end(),
                        [&](const int x) { 
                        return (*jetEta)[x] > 1.6; // put your condition here
                        }), indicesR.end());

            indicesR.erase(std::remove_if(
                        indicesR.begin(), indicesR.end(),
                        [&](const int x) { 
                        return (*jetPt)[x] < 500; // put your condition here
                        }), indicesR.end());

            indicesR.erase(std::remove_if(
                        indicesR.begin(), indicesR.end(),
                        [&](const int x) { 
                        return (*chargedMultiplicity)[x] < 60; // put your condition here
                        }), indicesR.end());

            int genNix = indicesG.size();
            int jetNix = indicesR.size();
            int genN = genNix;

            if(genNix != jetNix){
                if(genNix > jetNix){
                hGenExtra_pt_ChgMult ->Fill((*genJetPt)[indicesG[-1]]     , (*genJetChargedMultiplicity)[indicesG[-1]] );
                hGenExtra_pt_Eta     ->Fill((*genJetPt)[indicesG[-1]]     , (*genJetEta)[indicesG[-1]] );
                hGenExtra_Eta_ChgMult->Fill((*genJetEta)[indicesG[-1]]    , (*genJetChargedMultiplicity)[indicesG[-1]] );
                }
                if(jetNix > genNix){
                hRecExtra_pt_ChgMult ->Fill((*jetPt)[indicesR[-1]]     , (*chargedMultiplicity)[indicesR[-1]] );
                hRecExtra_pt_Eta     ->Fill((*jetPt)[indicesR[-1]]     , (*jetEta)[indicesR[-1]] );
                hRecExtra_Eta_ChgMult->Fill((*jetEta)[indicesR[-1]]    , (*chargedMultiplicity)[indicesR[-1]] ); 
                }
            }

            hEvent_Pass->Fill(1);

            h2GenN_JetN->Fill(genNix, jetNix);
            h2GenN_JetN_W->Fill(genNix, jetNix                                                  , genWeightPy);

            h1JetN_d_GenN->Fill(jetNix - genNix);
            h1JetN_d_GenN_W->Fill(jetNix - genNix                                               , genWeightPy);

            h1GenN->Fill(genNix);
            h1GenN_W->Fill(genNix                                                             , genWeightPy);

            h1JetN->Fill(jetNix);
            h1JetN_W->Fill(jetNix                                                             , genWeightPy);

            if(jetNix > 1){
                h2_J1M_J2M_any->Fill((*jetNumDaughters)[indicesR[0]],(*jetNumDaughters)[indicesR[1]]            , genWeightPy);
                h2_J1M_J2M_chg->Fill((*chargedMultiplicity)[indicesR[0]],(*chargedMultiplicity)[indicesR[1]]    , genWeightPy);

                h1_J1E_d_J2E->Fill((*jetEta)[indicesR[0]]-(*jetEta)[indicesR[1]]                                , genWeightPy);
                h1_J1P_d_J2P->Fill((*jetPhi)[indicesR[0]]-(*jetPhi)[indicesR[1]]                                , genWeightPy);
            }
            if(genNix > 1){ 
                h2_gen_J1M_J2M_any->Fill(genDau_chg->at(indicesG[0]).size(),genDau_chg->at(indicesG[1]).size()  , genWeightPy);
                h2_gen_J1M_J2M_chg->Fill((*genJetChargedMultiplicity)[indicesG[0]], (*genJetChargedMultiplicity)[indicesG[1]], genWeightPy);

                h1_gJ1E_d_gJ2E->Fill((*genJetEta)[indicesG[0]]-(*genJetEta)[indicesG[1]]                        , genWeightPy);
                h1_gJ1P_d_gJ2P->Fill((*genJetPhi)[indicesG[0]]-(*genJetPhi)[indicesG[1]]                        , genWeightPy);
            }
            if(jetNix == 1){
                h1_single_J1M_any->Fill((*jetNumDaughters)[indicesR[0]]                               , genWeightPy);
                h1_single_J1M_chg->Fill((*chargedMultiplicity)[indicesR[0]]                           , genWeightPy);
                h1_single_J1pt->Fill((*jetPt)[indicesR[0]]                                            , genWeightPy);
            }
            if(jetNix >= 2){
                h1_double_J1M_any->Fill((*jetNumDaughters)[indicesR[0]]                               , genWeightPy);
                h1_double_J1M_chg->Fill((*chargedMultiplicity)[indicesR[0]]                           , genWeightPy);
                h1_double_J1pt->Fill((*jetPt)[indicesR[0]]                                            , genWeightPy);

                h1_double_J2M_any->Fill((*jetNumDaughters)[indicesR[1]]                               , genWeightPy);
                h1_double_J2M_chg->Fill((*chargedMultiplicity)[indicesR[1]]                           , genWeightPy);
                h1_double_J2pt->Fill((*jetPt)[indicesR[1]]                                            , genWeightPy);
            }

            int ChgMult = 1;
            timer.StartSplit("Jet Selection");
            for(int kjet=0; kjet<jetNix; kjet++){

                int ijet = indicesR[kjet];

                long int NNtrk = (*jetNumDaughters)[ijet];
                float maxanypt = 0;
                float maxchgpt = 0;
                long int ChgMult  =0;


                hTotal_ijet->Fill(1);
                if(genN <= ijet){
                    hExtra_ijet->Fill(1);
                    continue;
                }

                //hjetPt->Fill(genJetPt->at(ijet), xs[qHatHist->FindBin(genQScale) - 1 ] / qHatHist->GetBinContent(qHatHist->FindBin(genQScale)) );

                for(long int  XXtrk=0; XXtrk < NNtrk; XXtrk++ ){
                    if( (genDau_chg->at(ijet)).size() == 0 ) continue;
                    double LeadingJT = ptWRTJet( (double) (*jetPt)[ijet] , (double) (*jetEta)[ijet] , (double) (*jetPhi)[ijet] , (double)(*dau_pt)[ijet][XXtrk] , (double)(*dau_eta)[ijet][XXtrk], (double)(*dau_phi)[ijet][XXtrk]);
                    if(maxanypt < LeadingJT ){maxanypt = LeadingJT;}
                    if(XXtrk == NNtrk - 1){ 
                        hLeadJTany  ->Fill(maxanypt                 , genWeightPy);
                        hLeadJTany2d->Fill(maxanypt,NNtrk           , genWeightPy);
                    }
                    if((*dau_chg)[ijet][XXtrk] == 0) continue;
                    ChgMult=ChgMult+1;
                }

                for(long int  XXtrk=0; XXtrk < NNtrk; XXtrk++ ){
                    double LeadingJT = ptWRTJet( (double) (*jetPt)[ijet] , (double) (*jetEta)[ijet] , (double) (*jetPhi)[ijet] , (double)(*dau_pt)[ijet][XXtrk] , (double)(*dau_eta)[ijet][XXtrk], (double)(*dau_phi)[ijet][XXtrk]);
                    if((*dau_chg)[ijet][XXtrk] == 0 && XXtrk == NNtrk - 1){
                        hLeadJTchg  ->Fill(maxchgpt                 , genWeightPy);
                        hLeadJTchg2d->Fill(maxchgpt,ChgMult         , genWeightPy);
                    }
                    if((*dau_chg)[ijet][XXtrk] == 0) continue;
                    if(maxchgpt < LeadingJT ){
                        maxchgpt = LeadingJT;
                    }
                    if(XXtrk == NNtrk - 1){ 
                        hLeadJTchg  ->Fill(maxchgpt                 , genWeightPy);
                        hLeadJTchg2d->Fill(maxchgpt, ChgMult        , genWeightPy);
                    }
                }

                if(ChgMult < trackbinbounds[0]){ 
                    hMult_Fail->Fill(ChgMult, genWeightPy);
                    continue;
                }

                int tkBool[trackbin] = {0};
                for(int i = 0; i < trackbin; i++){
                    if(ChgMult >= trackbinbounds[i]){
                        tkBool[i] = 1;
                        hJet_Pass           ->Fill(i);
                        hJet_Kin_Eta[i]     ->Fill((  *jetEta)[ijet], genWeightPy);
                        hJet_Kin_Phi[i]     ->Fill((  *jetPhi)[ijet], genWeightPy);
                        hN_ChgDAUvsJPT[i]   ->Fill(ChgMult, (  *jetPt)[ijet]    , genWeightPy);
                        hN_ChgDAU[i]        ->Fill(ChgMult                      , genWeightPy);
                    }
                    //cout << tkBool[i] << ", ";
                }
                //cout << " " << endl;

                timer.StartSplit("Track Loops");

                int Ntrig[trackbin][ptbin] = {0};
                int A_ptBool[NNtrk][ptbin] = {0};
                for(int  A_trk=0; A_trk < NNtrk; A_trk++ ){ 
                    if((*dau_chg)[ijet][A_trk] == 0) continue;

                    double jet_dau_pt    =  ptWRTJet( (double) (*jetPt)[ijet] , (double) (*jetEta)[ijet] , (double) (*jetPhi)[ijet] , (double)(*dau_pt)[ijet][A_trk] , (double)(*dau_eta)[ijet][A_trk], (double)(*dau_phi)[ijet][A_trk]);

                    if( jet_dau_pt < ptbinmin || jet_dau_pt > ptbinmax) continue; 

                    //cout << "New Daughter Bool Grid" << endl;
                    for(int i = 0; i < ptbin; i++){
                        if(jet_dau_pt >= ptbinbounds_lo[i] && jet_dau_pt < ptbinbounds_hi[i]){
                            A_ptBool[A_trk][i] = 1;
                        }    
                        //cout << A_ptBool[A_trk][i] << ", ";
                    }
                    //cout << "End" << endl;

                    for(int i = 0; i < trackbin; i++){
                        for(int j = 0; j < ptbin; j++){
                            if(tkBool[i] + A_ptBool[A_trk][j] == 2){
                                Ntrig[i][j] += 1;
                            }
                        }
                    }
                }

                for(int  A_trk=0; A_trk < (NNtrk-1); A_trk++ ){

                    if((*dau_chg)[ijet][A_trk] == 0) continue;

                    double jet_dau_pt    =  ptWRTJet( (double) (*jetPt)[ijet] , (double) (*jetEta)[ijet] , (double) (*jetPhi)[ijet] , (double)(*dau_pt)[ijet][A_trk] , (double)(*dau_eta)[ijet][A_trk], (double)(*dau_phi)[ijet][A_trk]);
                    if( jet_dau_pt < ptbinmin || jet_dau_pt > ptbinmax) continue; 

                    int c1 = 0;
                    if((*dau_chg)[ijet][A_trk] < 0) c1 = -1;
                    if((*dau_chg)[ijet][A_trk] > 0) c1 =  1;              

                    double jet_dau_theta = 2*ATan(Exp(-(etaWRTJet((double)(*jetPt)[ijet],(double)(*jetEta)[ijet],(double)(*jetPhi)[ijet],(double)(*dau_pt)[ijet][A_trk],(double)(*dau_eta)[ijet][A_trk],(double)(*dau_phi)[ijet][A_trk]))));
                    double jet_dau_eta   = etaWRTJet( (double) (*jetPt)[ijet] , (double) (*jetEta)[ijet] , (double) (*jetPhi)[ijet] , (double)(*dau_pt)[ijet][A_trk] , (double)(*dau_eta)[ijet][A_trk], (double)(*dau_phi)[ijet][A_trk]);
                    double jet_dau_phi   = phiWRTJet( (double) (*jetPt)[ijet] , (double) (*jetEta)[ijet] , (double) (*jetPhi)[ijet] , (double)(*dau_pt)[ijet][A_trk] , (double)(*dau_eta)[ijet][A_trk], (double)(*dau_phi)[ijet][A_trk]);

                    //                    if(jet_dau_eta <  0.86) continue;
                    //                    if(jet_dau_eta >  3.90) continue;

                    double beam_dau_theta = (*dau_theta)[ijet][A_trk];
                    double beam_dau_eta   = (*dau_eta)[ijet][A_trk];
                    double beam_dau_phi   = (*dau_phi)[ijet][A_trk];
                    double beam_dau_pt    = (*dau_pt)[ijet][A_trk];

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
                                hDau_Kin_Theta          [i][j]->Fill(beam_dau_theta , genWeightPy);
                                hDau_Kin_Eta            [i][j]->Fill(beam_dau_eta   , genWeightPy);
                            }
                        }
                    }
                }
            }
        }
        fFile->Close();
    }

    timer.StartSplit("Write Files");

    //WRITING FILES FOLD

    TFile* fS_tempA = new TFile(Form("src/Dists/MC_Full470_F/Primary/Primary_job%d.root",job), "recreate");


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

    h2GenN_JetN         ->Write();
    h2GenN_JetN_W       ->Write();

    h1JetN_d_GenN       ->Write();
    h1JetN_d_GenN_W     ->Write();

    h1GenN              ->Write();
    h1GenN_W            ->Write();

    h1JetN              ->Write();
    h1JetN_W            ->Write();

    h2_J1M_J2M_any      ->Write();
    h2_J1M_J2M_chg      ->Write();
    h2_gen_J1M_J2M_any  ->Write();
    h2_gen_J1M_J2M_chg  ->Write();
    h1_single_J1M_any   ->Write();
    h1_single_J1M_chg   ->Write();
    h1_single_J1pt      ->Write();
    h1_double_J1M_any   ->Write();
    h1_double_J1M_chg   ->Write();
    h1_double_J1pt      ->Write();
    h1_double_J2M_any   ->Write();
    h1_double_J2M_chg   ->Write();
    h1_double_J2pt      ->Write();    

    h1_J1E_d_J2E        ->Write();
    h1_J1P_d_J2P        ->Write();
    h1_gJ1E_d_gJ2E      ->Write();
    h1_gJ1P_d_gJ2P      ->Write();

    hTotal_ijet         ->Write();
    hExtra_ijet         ->Write();


    fS_tempA->Close();

    TFile* fS_tempB = new TFile(Form("src/Dists/MC_Full470_F/Secondary/Secondary_job%d.root",job), "recreate");
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

    timer.Stop();
    timer.Report();

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


