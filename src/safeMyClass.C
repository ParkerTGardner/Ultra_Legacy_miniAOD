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

//---------------------------------------------------------------------Global Vars

const int           asize=70;
const long int      bsize=1000;
const double        PI = 3.14159265359;
const int           DeltaPBin =32;
const float         etabound  =6;
const int           DeltaEBin =40;//(2*etabound)/(0.3);

const int trackbin 		= 6; 
const int ptbin		    = 13;
const float ptmin       = 0.0;
const float ptmax       = 4.0;
const int trackbinbounds[trackbin] = {85,86,87,88,89,90}; 
const float ptbinbounds_lo[ptbin] = {0.0, 0.3, 0.5, 1.0, 0.0, 0.0, 0.0, 0.3, 0.3, 0.5, 0.5, 1.0, 1.0}; 
const float ptbinbounds_hi[ptbin] = {3.0, 3.0, 3.0, 3.0, 1.0, 2.0, 4.0, 2.0, 4.0, 2.0, 4.0, 2.0, 4.0}; 

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
const int HB2 = 15;
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

//---------------------------------------------------------------------BOOL FUNCTIONS

bool F_eventpass(
        std::vector< float > *jetPt,
        int jetN,
        float jetPtCut
        ){
    float max_element =0;
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

    TH1::SetDefaultSumw2(kFALSE);
    TH2::SetDefaultSumw2(kFALSE);

    long int skipA_lo = 0;
    long int skipA_hi = 0;
    long int skipB = 0;
    long int tripA = 0;
    long int tripB = 0;

    long int skipper =0;
    long int tripper =0;

    const int bin120 = 120;
    const int bin100 = 100;
    const int HB120 = 120;
    const int HB4000 = 4000;
    const int LB500 = 500;
    const int HB15 = 15;

    //-----------------------------------------------------------------Canvas, Print  and Hist
    TH1D* hError_catch   = new TH1D("hError_catch","hError_catch", 1,LB0,1);

    TH1D* hWhatA   = new TH1D("hWhatA","hWhatA", 10,-2,8);
    TH1D* hWhenA   = new TH1D("hWhenA","hWhenA", 10,-2,8);
    TH1D* hEvent_Pass   = new TH1D("hEvent_Pass","hEvent_Pass", 1,LB0,1);

    TH1D* hPass_Trig_Jet = new TH1D("hPass_Trig_Jet","hPass_Trig_Jet"   , trackbin, LB0, trackbin);
    TH1D* hFail_Trig_Jet = new TH1D("hFail_Trig_Jet","hFail_Trig_Jet"   , trackbin, LB0, trackbin);

    TH1D* hJet_Kin_Eta    = new TH1D("hJet_Kin_Eta","hJet_Kin_Eta"              ,bin1  , LB2   , HB2);
    TH1D* hJet_Kin_Phi    = new TH1D("hJet_Kin_Phi","hJet_Kin_Phi"              ,bin1  , LB2pi , HB2pi);
    TH1D* hJet_Kin_Theta  = new TH1D("hJet_Kin_Theta" ,"hJet_Kin_Theta"         ,bin1  , LB2pi , HB2pi);

    TH1D* hDau_Kin_WRTJ_Phi   = new TH1D( "hRotPhi"   ,"hRotPhi"                , bin2 , LB2pi , HB2pi);
    TH1D* hDau_Kin_WRTJ_Pt    = new TH1D( "hRotPt"    ,"hRotPt"                 , bin2 , LB0   , HB20);
    TH1D* hDau_Kin_WRTJ_Theta = new TH1D( "hRotTheta" ,"hRotTheta"              , bin2 , LB0   , HB1);
    TH1D* hDau_Kin_WRTJ_Eta   = new TH1D( "hRotEta"   ,"hRotEta"                , bin_WRTJ_Eta , low_WRTJ_Eta_Bin   , high_WRTJ_Eta_Bin);

    TH1D* hDau_Kin_Phi   = new TH1D( "hPhi"   ,"hPhi"                           , bin2 , LB2pi , HB2pi);
    TH1D* hDau_Kin_Pt    = new TH1D( "hPt"    ,"hPt"                            , bin2 , LB0   , HB20);
    TH1D* hDau_Kin_Theta = new TH1D( "hTheta" ,"hTheta"                         , bin2 , LB0   , HB15pi);
    TH1D* hDau_Kin_Eta   = new TH1D( "hEta"   ,"hEta"                           , bin_Eta , low_Eta_Bin   , high_Eta_Bin);

    TH2D* hEtaPhi_dd_Dau_Jet = new TH2D("hEtaPhi_dd_Dau_Jet","hEtaPhi_dd_Dau_Jet" , bin3, LB1 , HB1 , bin3 , LB1  , HB1);

    TH2D* hdNdEvsPT[trackbin];       //ptl level
    TH1D* hdNdE[trackbin];          //ptl level

    TH2D* hSignal[trackbin][ptbin];
    TH2D* hBckrnd[trackbin][ptbin];
    TH2D* hEPDraw[trackbin][ptbin];

    TH1D* hRawDPhi[trackbin][ptbin];

    TH1D* hN_ChgDAU      = new TH1D("hN_ChgDAU"     ,"hN_ChgDAU"       , bin120 , LB0 , HB120);
    TH2D* hN_ChgDAUvsJPT = new TH2D("hN_ChgDAUvsJPT","hN_ChgDAUvsJPT"  , bin120 , LB0 , HB120 , bin100 , LB500 , HB4000);

    for(int wtrk = 1; wtrk<trackbin+1; wtrk++){

        hdNdEvsPT[wtrk-1]      = new TH2D(Form("hdNdEvsPT_trk_%d",wtrk) ,Form("dNdEvsPT_trk_trk_%d",wtrk)              , bin_WRTJ_Eta , low_WRTJ_Eta_Bin   , high_WRTJ_Eta_Bin, bin100 , LB0 , HB15);
        hdNdE[wtrk-1]          = new TH1D(Form("hdNdE_trk_%d",wtrk) ,Form("hdNdE_trk_%d",wtrk)                         , bin_WRTJ_Eta , low_WRTJ_Eta_Bin   , high_WRTJ_Eta_Bin);

        for(int wppt = 1; wppt<ptbin+1; wppt++){
            hBckrnd[wtrk-1][wppt-1] = new TH2D(Form("hBckrnd_trk_%d_ppt_%d",wtrk,wppt) ,Form("hBckrnd_trk_%d_ppt_%d",wtrk,wppt) , DeltaEBin, -etabound   , etabound    , DeltaPBin, -0.5*PI  , HB15pi);
            hSignal[wtrk-1][wppt-1] = new TH2D(Form("hSignal_trk_%d_ppt_%d",wtrk,wppt) ,Form("hSignal_trk_%d_ppt_%d",wtrk,wppt) , DeltaEBin, -etabound   , etabound    , DeltaPBin, -0.5*PI  , HB15pi);
            hEPDraw[wtrk-1][wppt-1] = new TH2D(Form("hEPDraw_trk_%d_ppt_%d",wtrk,wppt) ,Form("hEPDraw_trk_%d_ppt_%d",wtrk,wppt) , EPD_xb   , EPD_xlo, EPD_xhi , EPD_yb      , EPD_ylo    , EPD_yhi);

            hRawDPhi[wtrk-1][wppt-1] = new TH1D(Form("hRawDPhi_trk_%d_ppt_%d",wtrk,wppt) ,Form("hRawDPhi_trk_%d_ppt_%d",wtrk,wppt) , 150,  -2.5*PI, 2.5*PI);
            hRawDEta[wtrk-1][wppt-1] = new TH1D(Form("hRawDEta_trk_%d_ppt_%d",wtrk,wppt) ,Form("hRawDEta_trk_%d_ppt_%d",wtrk,wppt) , 180,  -1.5*etabound, 1.5*etabound);
        }
    }

    //-------------------------------------------------------------------MAIN_LOOPS
    //-------------------------------------------------------------------MAIN_LOOPS

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

        for (Long64_t ievent=0; ievent <nentries; ievent ++){
            timer.StartSplit("Loading Events");
            Long64_t jevent = LoadTree(ievent);
            nb = fChain->GetEntry(ievent);   nbytes += nb;

            if(ievent%1000==0) cout<< ievent << "/" << nentries  <<endl;

            timer.StartSplit("Event Selection");
            if(!F_eventpass(jetPt, jetN, jetPtCut)) continue;

            hEvent_Pass->Fill(0);

            timer.StartSplit("Jet Selection");
            for(int ijet=0; ijet<jetN; ijet++){

                long int NNtrk = (*jetNumDaughters)[ijet];
                long int Ntrig  =0;

                for(long int  XXtrk=0; XXtrk < NNtrk; XXtrk++ ){
                    if((*dau_chg)[ijet][XXtrk] == 0) continue;
                    Ntrig=Ntrig+1;
                }

                if(Ntrig < trackbinbounds[0]) continue;

                int a = -1;

                hWhenA->Fill(1);

                for(int i = 0; i < trackbin; i++){
                    if (Ntrig >= trackbinbounds[i] && Ntrig < trackbinbounds[i+1] ){
                        a = i;
                        break;
                    }
                }

                hWhatA->Fill(a);

                if(a == -1){
                    hError_catch->Fill(1);
                    continue;
                }


                //if(a > (trackbin-1)){
                //    a = (trackbin-1);
                //}
                //if( a < trackbin-2) continue;

                if( !F_jetpass(jetEta, jetPt, ijet, jetPtCut)){
                    hFail_Trig_Jet->Fill(a);
                    continue;
                }
                hPass_Trig_Jet->Fill(a);

                hJet_Kin_Eta  ->Fill((  *jetEta)[ijet]);
                hJet_Kin_Theta->Fill((*jetTheta)[ijet]);
                hJet_Kin_Phi  ->Fill((  *jetPhi)[ijet]);

                hN_ChgDAU     ->Fill(Ntrig);
                hN_ChgDAUvsJPT->Fill(Ntrig,(*jetPt)[ijet]);

                timer.StartSplit("Track Loops");
                for(long int  XXtrk=0;XXtrk < NNtrk; XXtrk++ ){
                    if((*dau_chg)[ijet][XXtrk] == 0) continue;

                    float THISeventTheta = 2*ATan(Exp(-(etaWRTJet((double)(*jetPt)[ijet],(double)(*jetEta)[ijet],(double)(*jetPhi)[ijet],(double)(*dau_pt)[ijet][XXtrk],(double)(*dau_eta)[ijet][XXtrk],(double)(*dau_phi)[ijet][XXtrk]))));
                    float THISeventEta   = etaWRTJet( (double) (*jetPt)[ijet] , (double) (*jetEta)[ijet] , (double) (*jetPhi)[ijet] , (double)(*dau_pt)[ijet][XXtrk] , (double)(*dau_eta)[ijet][XXtrk], (double)(*dau_phi)[ijet][XXtrk]);
                    float THISeventPhi   = phiWRTJet( (double) (*jetPt)[ijet] , (double) (*jetEta)[ijet] , (double) (*jetPhi)[ijet] , (double)(*dau_pt)[ijet][XXtrk] , (double)(*dau_eta)[ijet][XXtrk], (double)(*dau_phi)[ijet][XXtrk]);
                    float THISeventPt    =  ptWRTJet( (double) (*jetPt)[ijet] , (double) (*jetEta)[ijet] , (double) (*jetPhi)[ijet] , (double)(*dau_pt)[ijet][XXtrk] , (double)(*dau_eta)[ijet][XXtrk], (double)(*dau_phi)[ijet][XXtrk]);

                    float Eta_d_Dau_Jet = (*dau_eta)[ijet][XXtrk] - (*jetEta)[ijet];
                    float Phi_d_Dau_Jet = (*dau_phi)[ijet][XXtrk] - (*jetPhi)[ijet];

                    hDau_Kin_WRTJ_Theta->Fill(THISeventTheta);
                    hDau_Kin_WRTJ_Eta->Fill(THISeventEta);
                    hDau_Kin_WRTJ_Phi->Fill(THISeventPhi);
                    hDau_Kin_WRTJ_Pt->Fill(THISeventPt);

                    hdNdEvsPT[a]->Fill(THISeventEta,THISeventPt);
                    hdNdE[a]    ->Fill(THISeventEta);

                    float beam_dau_theta = (*dau_theta)[ijet][XXtrk];
                    float beam_dau_eta   = (*dau_eta)[ijet][XXtrk];
                    float beam_dau_phi   = (*dau_phi)[ijet][XXtrk];
                    float beam_dau_pt    = (*dau_pt)[ijet][XXtrk];

                    hDau_Kin_Theta->Fill(beam_dau_theta);
                    hDau_Kin_Eta->Fill(beam_dau_eta);
                    hDau_Kin_Phi->Fill(beam_dau_phi);
                    hDau_Kin_Pt->Fill(beam_dau_pt);

                    hEtaPhi_dd_Dau_Jet->Fill(Eta_d_Dau_Jet,Phi_d_Dau_Jet);

                    double Db_1 = ptWRTJet(  (double) (*jetPt)[ijet] , (double) (*jetEta)[ijet] , (double) (*jetPhi)[ijet] , (double)(*dau_pt)[ijet][XXtrk] , (double)(*dau_eta)[ijet][XXtrk], (double)(*dau_phi)[ijet][XXtrk] );

                    int b_1 = -99;

                    //int ptbin = 4;
                    //float ptbin_bounds[ptbin] = {0, .3, .9, 1.7};

                    for(int i = 0; i < ptbin; i++){

                        if(Db_1 < ptbinbounds[i])   skipA_lo +=1;
                        if(Db_1 > ptbinbounds[i+1]) skipA_hi +=1;

                        if (Db_1 >= ptbinbounds[i] && Db_1 < ptbinbounds[i+1] ){
                            b_1 = i;
                            tripA +=1;
                            break;
                        }
                    }

                    if(b_1 == -99) continue;

                    hEPDraw[a][b_1]->Fill(THISeventEta, THISeventPhi, 1./Ntrig);

                    for(long int  YYtrk=XXtrk+1; YYtrk< NNtrk; YYtrk++ ){

                        if((*dau_chg)[ijet][YYtrk] == 0) continue;
                        if(YYtrk == XXtrk) continue;

                        double Db_3 = ptWRTJet(  (double) (*jetPt)[ijet] , (double) (*jetEta)[ijet] , (double) (*jetPhi)[ijet] , (double)(*dau_pt)[ijet][YYtrk] , (double)(*dau_eta)[ijet][YYtrk], (double)(*dau_phi)[ijet][YYtrk] );
                        int b_3 = 99;

                        for(int i = 0; i < ptbin; i++){
                            if (Db_3 >= ptbinbounds[i] && Db_3 < ptbinbounds[i+1] ){
                                b_3 = i;
                                tripB +=1;
                                break;
                            }
                            skipB +=1;
                        }

                        if(b_3 == 99) continue;

                        if(b_3 != b_1){ 
                            skipper +=1;
                            continue;
                        }
                        tripper +=1;

                        float YYTHISeventEta = etaWRTJet( (double) (*jetPt)[ijet] , (double) (*jetEta)[ijet] , (double) (*jetPhi)[ijet] , (double)(*dau_pt)[ijet][YYtrk] , (double)(*dau_eta)[ijet][YYtrk], (double)(*dau_phi)[ijet][YYtrk]);
                        float YYTHISeventPhi = phiWRTJet( (double) (*jetPt)[ijet] , (double) (*jetEta)[ijet] , (double) (*jetPhi)[ijet] , (double)(*dau_pt)[ijet][YYtrk] , (double)(*dau_eta)[ijet][YYtrk], (double)(*dau_phi)[ijet][YYtrk]);

                        float deltaEta = THISeventEta - YYTHISeventEta;
                        float deltaPhi = THISeventPhi - YYTHISeventPhi;

                        hRawDPhi[a][b_1]->Fill(deltaPhi);

                        if(deltaPhi<0) deltaPhi+=2*PI;
                        if(deltaPhi>1.5*PI) deltaPhi -= 2*PI;

                        hSignal[a][b_1]->Fill(deltaEta, deltaPhi, 1./Ntrig);
                        hSignal[a][b_1]->Fill(-deltaEta, deltaPhi, 1./Ntrig);
                        hSignal[a][b_1]->Fill(deltaEta, -deltaPhi, 1./Ntrig);
                        hSignal[a][b_1]->Fill(-deltaEta, -deltaPhi, 1./Ntrig);

                        //if(fabs(deltaEta) <etabound && fabs(deltaEta) >1.1) hSignal_sub[a][b_1]->Fill(fabs(deltaEta), deltaPhi, 1./Ntrig);

                    }
                }
            }
        }

        fFile->Close();
    }

    timer.StartSplit("EndJob");
    cout << "skipper is: " << skipper <<endl;
    cout << "tripper is: " << tripper <<endl;

    cout << "skipA lo and hi are: " << skipA_lo << " " << skipA_hi <<endl;
    cout << "skipB is: " << skipB <<endl;

    cout << "tripA is: " << tripA <<endl;
    cout << "tripB is: " << tripB <<endl;



    //auto rng = new TRandom3();
    int backMult =5;
    for(int wtrk = 1; wtrk < trackbin+1; wtrk++){
        for(int wppt = 1; wppt < ptbin+1; wppt++){
            cout << "ppt is " << wppt << " and trk is " << wtrk << endl;

            long int NENT =  hSignal[wtrk-1][wppt-1]->GetEntries();
            long int XENT =  1 + floor((1+sqrt(1+(backMult*8*NENT)))/2);

            float A_ETA[XENT] = {0};
            float A_PHI[XENT] = {0};

            for(int x = 0; x<XENT; x++){

                double WEta1, WPhi1;
                hEPDraw[wtrk-1][wppt-1]->GetRandom2(WEta1, WPhi1);

                A_ETA[x] = WEta1;
                A_PHI[x] = WPhi1;

            }

            for(long int i = 0; i < (XENT-1); i++){

                for(long int j = (i+1); j < XENT; j++){

                    double WdeltaEta = A_ETA[i]-A_ETA[j];
                    double WdeltaPhi = A_PHI[i]-A_PHI[j];

                    if(WdeltaPhi<0) WdeltaPhi+=2*PI;
                    if(WdeltaPhi>1.5*PI) WdeltaPhi -= 2*PI;

                    hBckrnd[wtrk-1][wppt-1]->Fill(WdeltaEta, WdeltaPhi, 1./XENT);
                    hBckrnd[wtrk-1][wppt-1]->Fill(-WdeltaEta, WdeltaPhi, 1./XENT);
                    hBckrnd[wtrk-1][wppt-1]->Fill(WdeltaEta, -WdeltaPhi, 1./XENT);
                    hBckrnd[wtrk-1][wppt-1]->Fill(-WdeltaEta, -WdeltaPhi, 1./XENT);

                }
            }
            cout << "***" << endl;
            cout << "hSig ent are: " << NENT << " ... and hBack ent are: " << hBckrnd[wtrk-1][wppt-1]->GetEntries() << endl;
            cout << "***" << endl;
        }
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    TFile* fS_temp = new TFile(Form("src/testNew/job%d.root",job), "recreate");

    //TFile* fS_temp = new TFile(Form("src/2018_out/job%d.root",job), "recreate");
    //TFile* fS_temp = new TFile(Form("src/2017_out/job%d.root",job), "recreate");
    //TFile* fS_temp = new TFile(Form("src/2016_out/job%d.root",job), "recreate");
    //TFile* fS_temp = new TFile(Form("src/pythia_600_800/job%d.root",job), "recreate");

    for(int wtrk =1; wtrk <trackbin+1; wtrk++){

        //hdNdEvsPT[wtrk-1]->Scale(1./(NtrkClass_Jtrig[wtrk-1]));
        //hdNdE[wtrk-1]    ->Scale(1./(NtrkClass_Jtrig[wtrk-1]));

        //hdNdEvsPT[wtrk-1]->Scale(1./( hdNdEvsPT[wtrk-1]->GetXaxis()->GetBinWidth(3)  ));
        //hdNdE[wtrk-1]    ->Scale(1./( hdNdE[wtrk-1]    ->GetXaxis()->GetBinWidth(3)  ));

        hdNdEvsPT[wtrk-1]       ->Write(Form("hdNdEvsPT_%d",wtrk));
        hdNdE[wtrk-1]           ->Write(Form("hdNdE_%d",wtrk));

        for(int wppt =1; wppt <ptbin+1; wppt++){
            hSignal[wtrk-1][wppt-1]->Write(Form("hSig_%d_%d",wtrk,wppt));
            hBckrnd[wtrk-1][wppt-1]->Write(Form("hBck_%d_%d",wtrk,wppt));
            hEPDraw[wtrk-1][wppt-1]->Write(Form("hEPD_%d_%d",wtrk,wppt));
            hRawDPhi[wtrk-1][wppt-1]->Write(Form("hRawDPhi_%d_%d",wtrk,wppt));
        }
    } 
    hN_ChgDAU       ->Write();
    hN_ChgDAUvsJPT  ->Write();

    hWhenA->Write();
    hWhatA->Write();

    hError_catch  ->Write();

    hEvent_Pass   ->Write();
    hPass_Trig_Jet->Write();
    hFail_Trig_Jet->Write();

    hDau_Kin_WRTJ_Theta->Write();
    hDau_Kin_WRTJ_Eta  ->Write();
    hDau_Kin_WRTJ_Phi  ->Write();
    hDau_Kin_WRTJ_Pt   ->Write();

    hDau_Kin_Theta->Write();
    hDau_Kin_Eta  ->Write();
    hDau_Kin_Phi  ->Write();
    hDau_Kin_Pt   ->Write();

    hEtaPhi_dd_Dau_Jet ->Write();

    fS_temp->Close();

    /*
       TCanvas* c1 = new TCanvas("c1","c1" , 800, 800);
       hDau_Kin_Eta->Draw();
       TCanvas* c2 = new TCanvas("c2","c2" , 800, 800);
       hDau_Kin_WRTJ_Eta->Draw();

       TCanvas* c3          = new TCanvas("c3","c3" , 800, 800);
       c3->Divide(ptbin,trackbin);
       TCanvas* c4          = new TCanvas("c4","c4", 800, 800);
       c4->Divide(ptbin,trackbin);

       TFile *f = new TFile("Narrow_trk_binning_full_save.root");

       for(int wtrk =1; wtrk <trackbin+1; wtrk++){
       for(int wppt =1; wppt <ptbin+1; wppt++){
       int where = wppt + (wtrk-1)*(ptbin);
       c3->cd(where);

       TH2D* h1 = (TH2D*)f->Get(Form("hSig_%d_%d",wtrk,wppt));
       TH2D* h2 = (TH2D*)f->Get(Form("hBck_%d_%d",wtrk,wppt));

       h1->Divide(h2);

       TH1D *histfit1 = (TH1D*) h1->ProjectionY("",30,80)->Clone();

       std::string function = "[0]/(TMath::Pi()*2)*(1+2*([1]*TMath::Cos(x)+[2]*TMath::Cos(2*x)+[3]*TMath::Cos(3*x)+[4]*TMath::Cos(4*x)+[5]*TMath::Cos(5*x)))";

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
       histfit1->Fit(&func1, "m q E");
       auto fitResult1 = histfit1->Fit(&func1, "m S E q");
    //fitResult->Print();            
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

    c4->cd(where);
    histfit1->DrawCopy("PE");

    delete fv1;
    delete fv2;
    delete fv3;
    delete histfit1;
    delete h2;
    delete h1;

    }
}*/

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
