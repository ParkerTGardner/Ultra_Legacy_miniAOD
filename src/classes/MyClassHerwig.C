// This vim file has folds built in to it. 
// The .vimrc file for user Parker Gardner on Bonner cluster has been configured to save the fold configuration upon every :wq

// ---for action on all folds:
// to open ALL folds: hit esc and type zR
// to close ALL folds: hit esc and type zM

// ---for action on single fold: go to the specific fold.
// to open a fold: esc and zo 
// to close a fold: esc and zc

//this a change #2

#define MyClass_cxx
#include "include/MyClass.h"
#include "include/Timer.h"
#include <TStyle.h>
#include "TH1D.h"
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
#include <iomanip>
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include "include/coordinateTools.h"

using TMath::ATan;
using TMath::Exp;

//Quick change:

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


//Const variable folds


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

    //Initializing Histograms


    TH1D* hEvent_Pass   = new TH1D("hEvent_Pass","hEvent_Pass", 1,LB0,1);
    TH1D* hEvent_Fail   = new TH1D("hEvent_Fail","hEvent_Fail", 1,LB0,1);
    TH1D* hJet_Pass     = new TH1D("hJet_Pass"  ,"hJet_Pass"  , trackbin,LB0,trackbin);
    TH1D* hJet_Fail     = new TH1D("hJet_Fail"  ,"hJet_Fail"  , bin120,LB0,bin120);
    TH1D* hMult_Fail    = new TH1D("hMult_Fail" ,"hMult_Fail" , bin120,LB0,bin120);


    TH1D* hLeadJTany      = new TH1D("hLeadJTany" ,"hLeadJTany"     , bin2 , LB0   , HB20);
    TH2D* hLeadJTany2d    = new TH2D("hLeadJTany2d" ,"hLeadJTany2d" , bin2 , LB0   , HB20, bin150, LB0, bin150);

    TH1D* hLeadJTchg      = new TH1D("hLeadJTchg" ,"hLeadJTchg"     , bin2 , LB0   , HB20);
    TH2D* hLeadJTchg2d    = new TH2D("hLeadJTchg2d" ,"hLeadJTchg2d" , bin2 , LB0   , HB20, bin120, LB0, bin120);


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


        //Main loops
        for (Long64_t ievent=0; ievent <nentries; ievent ++){


            timer.StartSplit("Loading Events");

            Long64_t jevent = LoadTree(ievent);
            nb = fChain->GetEntry(ievent);   nbytes += nb;
            //if(ievent%1000==0) cout<< ievent << "/" << nentries  <<endl;

            timer.StartSplit("Event Selection");

            if(!F_eventpass(jetPt, jetN, jetPtCut)){
                hEvent_Fail->Fill(1);
                continue;
            }
            hEvent_Pass->Fill(1);



            timer.StartSplit("Jet Selection");
            for(int ijet=0; ijet<jetN; ijet++){

                long int NNtrk = (*jetNumDaughters)[ijet];
                float maxanypt = 0;
                float maxchgpt = 0;
                long int Nmult  =0;

                for(long int  XXtrk=0; XXtrk < NNtrk; XXtrk++ ){
                    double LeadingJT = ptWRTJet( (double) (*jetPt)[ijet] , (double) (*jetEta)[ijet] , (double) (*jetPhi)[ijet] , (double)(*dau_pt)[ijet][XXtrk] , (double)(*dau_eta)[ijet][XXtrk], (double)(*dau_phi)[ijet][XXtrk]);
                    if(maxanypt < LeadingJT ){maxanypt = LeadingJT;}
                    if(XXtrk == NNtrk - 1){ 
                        hLeadJTany  ->Fill(LeadingJT, genWeight);
                        hLeadJTany2d->Fill(LeadingJT,NNtrk, genWeight);
                    }
                    if((*dau_chg)[ijet][XXtrk] == 0) continue;
                    Nmult=Nmult+1;
                }

                for(long int  XXtrk=0; XXtrk < NNtrk; XXtrk++ ){
                    double LeadingJT = ptWRTJet( (double) (*jetPt)[ijet] , (double) (*jetEta)[ijet] , (double) (*jetPhi)[ijet] , (double)(*dau_pt)[ijet][XXtrk] , (double)(*dau_eta)[ijet][XXtrk], (double)(*dau_phi)[ijet][XXtrk]);
                    if((*dau_chg)[ijet][XXtrk] == 0 && XXtrk == NNtrk - 1){
                        hLeadJTchg  ->Fill(LeadingJT, genWeight);
                        hLeadJTchg2d->Fill(LeadingJT,Nmult, genWeight);
                    }
                    if((*dau_chg)[ijet][XXtrk] == 0) continue;
                    if(maxchgpt < LeadingJT ){
                        maxchgpt = LeadingJT;
                    }
                    if(XXtrk == NNtrk - 1){ 
                        hLeadJTchg  ->Fill(LeadingJT, genWeight);
                        hLeadJTchg2d->Fill(LeadingJT,Nmult, genWeight);
                    }
                }

                if(Nmult < trackbinbounds[0]){ 
                    hMult_Fail->Fill(Nmult, genWeight);
                    continue;
                }
                if( !F_jetpass(jetEta, jetPt, ijet, jetPtCut)){
                    hJet_Fail->Fill(Nmult, genWeight);
                    continue;
                }

                int tkBool[trackbin] = {0};
                for(int i = 0; i < trackbin; i++){
                    if(Nmult >= trackbinbounds[i]){
                        tkBool[i] = 1;
                        hJet_Pass           ->Fill(i);
                        hJet_Kin_Eta[i]     ->Fill((  *jetEta)[ijet], genWeight);
                        hJet_Kin_Phi[i]     ->Fill((  *jetPhi)[ijet], genWeight);
                        hN_ChgDAUvsJPT[i]   ->Fill(Nmult, (  *jetPt)[ijet], genWeight);
                        hN_ChgDAU[i]        ->Fill(Nmult, genWeight);
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

                                hEPDraw[i][j]->Fill(jet_dau_eta, jet_dau_phi, genWeight/Ntrig[i][j]);

                                if(c1 ==  1){ hEPDrawPluss[i][j]->Fill(jet_dau_eta, jet_dau_phi, genWeight/Ntrig[i][j]);}
                                if(c1 == -1){ hEPDrawMinus[i][j]->Fill(jet_dau_eta, jet_dau_phi, genWeight/Ntrig[i][j]);}

                                hdNdEvsPT               [i][j]->Fill(jet_dau_eta, jet_dau_pti, genWeight);
                                hdNdE                   [i][j]->Fill(jet_dau_eta, genWeight);

                                hDau_Kin_WRTJ_Phi       [i][j]->Fill(jet_dau_phi    , genWeight);
                                hDau_Kin_WRTJ_Pt        [i][j]->Fill(jet_dau_pt     , genWeight);
                                hDau_Kin_WRTJ_Theta     [i][j]->Fill(jet_dau_theta  , genWeight);
                                hDau_Kin_WRTJ_Eta       [i][j]->Fill(jet_dau_eta    , genWeight);
                                hDau_Kin_Phi            [i][j]->Fill(beam_dau_phi   , genWeight);
                                hDau_Kin_Pt             [i][j]->Fill(beam_dau_pt    , genWeight);
                                hDau_Kin_Theta          [i][j]->Fill(beam_dau_theta , genWeight);
                                hDau_Kin_Eta            [i][j]->Fill(beam_dau_eta   , genWeight);
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

    TFile* fS_tempA = new TFile(Form("src/Dists/Full2016/Primary/Primary_job%d.root",job), "recreate");

    hEvent_Pass   ->Write();
    hJet_Pass     ->Write();

    hEvent_Fail   ->Write(); 
    hJet_Fail     ->Write();
    hMult_Fail    ->Write();

    hLeadJTany2d  ->Write();
    hLeadJTany    ->Write();

    hLeadJTchg2d  ->Write();
    hLeadJTchg    ->Write();


    fS_tempA->Close();

    TFile* fS_tempB = new TFile(Form("src/Dists/Full2016/Secondary/Secondary_job%d.root",job), "recreate");
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


