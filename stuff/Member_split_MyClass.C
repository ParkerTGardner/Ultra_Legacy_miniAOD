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
const int           Resolution = 40;
const long int      RUNLENGTH = 4357342;
const int           p_inc = 2;
const int           gridinc = 2;
const int           DeltaPBin =60;
const int           DeltaEBin =80;
const int           DeltaPBin2 =20;
const int           DeltaEBin2 =60;
const float         etabound  =4;

//---------------------------------------------------------------------CUTS
const float PtTrihi     = 1000;
const float PtTrilow    = 0.0;
const float Thetahigh   = 0.55;
const float Thetalow    = 0.0;
const float EtaCut      = 0.0;
const float jetEtaCut   = 2.1;
const float jetPtCut    = 500.0;

//---------------------------------------------------------------------BOOL FUNCTIONS

/*
 void Modified_GetInt()
 {
    Int_t nbinsx = GetNbinsX();
    Int_t nbinsy = GetNbinsY();
    Int_t nbins  = nbinsx*nbinsy;
    //Double_t integral;
    // compute integral checking that all bins have positive content (see ROOT-5894)
    if (fIntegral) {
       if (fIntegral[nbins+1] != fEntries) integral = ComputeIntegral(true);
       else integral = fIntegral[nbins];
    } else {
       integral = ComputeIntegral(true);
    }
    //if (integral == 0 ) { x = 0; y = 0; return;}
    // case histogram has negative bins
    if (integral == TMath::QuietNaN() ) { x = TMath::QuietNaN(); y = TMath::QuietNaN(); return;}
  
    //if (!rng) rng = gRandom;
 }

 void Modified_GetRandom2(Double_t &x, Double_t &y, Double_t integral, TRandom * rng, nbins, nbinsx)
 {
    Double_t r1 = rng->Rndm();
    Int_t ibin = TMath::BinarySearch(nbins,integral,(Double_t) r1);
    Int_t biny = ibin/nbinsx;
    Int_t binx = ibin - nbinsx*biny;
    x = fXaxis.GetBinLowEdge(binx+1);
    if (r1 > integral[ibin]) x +=
       fXaxis.GetBinWidth(binx+1)*(r1-integral[ibin])/(integral[ibin+1] - integral[ibin]);
    y = fYaxis.GetBinLowEdge(biny+1) + fYaxis.GetBinWidth(biny+1)*rng->Rndm();
 }
*/

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

bool F_trigpassROTATED(
		std::vector< float > * jetPt,
		std::vector< float > * jetEta,
		std::vector< float > * jetPhi,
		vector<vector<int> > *dau_chg,
		vector<vector<float> > *dau_pt,
		vector<vector<float> > *dau_eta,
		vector<vector<float> > *dau_phi,
		long int            XXtrk,
		int                 ijet,
		float PtTrihi,
		float PtTrilow
		){
	if(                                                                                            (*dau_pt)[ijet][XXtrk] < 0.3)             return false;
	if(                                                                                     fabs((*dau_eta)[ijet][XXtrk] ) > 2.6)             return false;
	//double THETANOW1 = 2*ATan(Exp(-(etaWRTJet((*jetPt)[ijet],(*jetEta)[ijet],(*jetPhi)[ijet],(*dau_pt)[ijet][XXtrk],(*dau_Eta)[ijet][XXtrk],(*dau_Phi)[ijet][XXtrk]))));
	//if(                                                                THETANOW1 > Thetahigh && THETANOW1 < (PI-0*Thetahigh ))       return false;
	//if(                                                                THETANOW1 < Thetalow  && THETANOW1 < (PI-0*Thetalow  ))       return false;
	//if(           ptWRTJet((*jetPt)[ijet],(*jetEta)[ijet],(*jetPhi)[ijet],(*dau_pt)[ijet][XXtrk],(*dau_eta)[ijet][XXtrk],(*dau_phi)[ijet][XXtrk]) > PtTrihi)             return false;
	//if(           ptWRTJet((*jetPt)[ijet],(*jetEta)[ijet],(*jetPhi)[ijet],(*dau_pt)[ijet][XXtrk],(*dau_eta)[ijet][XXtrk],(*dau_phi)[ijet][XXtrk]) < PtTrilow)            return false;
	if(                                                                                             (*dau_chg)[ijet][XXtrk] == 0)             return false;
	return true;
}

void MyClass::Loop(int job){
    //timing information
    Timer timer = Timer();
    timer.Start();
    timer.StartSplit("Start Up");


    TH1::SetDefaultSumw2(kFALSE);
    TH2::SetDefaultSumw2(kFALSE);

	int Pass_Trig_Jet =0; long int skipper =0; long int tripper =0;
    int bin_WRTJ_Eta        = 150; float low_WRTJ_Eta_Bin  =  0; float high_WRTJ_Eta_Bin = 15;
    int bin_Eta             = 60; float low_Eta_Bin       = -3; float high_Eta_Bin      = 3;
	int trackbin 		= 7; int ptbin		    = 3; int tracksplitter 	= 15;
	//-----------------------------------------------------------------Canvas, Print  and Hist

	TH1D* hJet_Kin_Eta    = new TH1D("hJet_Kin_Eta","hJet_Kin_Eta" ,30  ,-15    ,15);
    TH1D* hJet_Kin_Phi    = new TH1D("hJet_Kin_Phi","hJet_Kin_Phi" ,30  ,-2*PI ,2*PI);
    TH1D* hJet_Kin_Theta  = new TH1D("hJet_Kin_Theta" ,"hJet_Kin_Theta"  ,30  ,-2*PI ,2*PI);

    TH1D* hDau_Kin_WRTJ_Phi   = new TH1D( "hRotPhi"   ,"hRotPhi"   , 100 , -2*PI, 2*PI);
    TH1D* hDau_Kin_WRTJ_Pt    = new TH1D( "hRotPt"    ,"hRotPt"    , 100 ,  0   , 20);
	TH1D* hDau_Kin_WRTJ_Theta = new TH1D( "hRotTheta" ,"hRotTheta" , 100 , 0, 1);
	TH1D* hDau_Kin_WRTJ_Eta   = new TH1D( "hRotEta"   ,"hRotEta"   , bin_WRTJ_Eta , low_WRTJ_Eta_Bin   , high_WRTJ_Eta_Bin);

	TH1D* hDau_Kin_Phi   = new TH1D( "hPhi"   ,"hPhi"   , 100 , -2*PI, 2*PI);
	TH1D* hDau_Kin_Pt    = new TH1D( "hPt"    ,"hPt"    , 100 ,  0   , 20);
    TH1D* hDau_Kin_Theta = new TH1D( "hTheta" ,"hTheta" , 100 , 0, 1.5*PI);
	TH1D* hDau_Kin_Eta   = new TH1D( "hEta"   ,"hEta"   , bin_Eta , low_Eta_Bin   , high_Eta_Bin);

	TH2D* hEtaPhi_dd_Dau_Jet = new TH2D("hEtaPhi_dd_Dau_Jet","hEtaPhi_dd_Dau_Jet", 50, -1 , 1 , 50 , -1  , 1);


/*

    TH1D* hDau_Kin_WRTJ_Phi   
    TH1D* hDau_Kin_WRTJ_Pt    
	TH1D* hDau_Kin_WRTJ_Theta 
	TH1D* hDau_Kin_WRTJ_Eta   

	TH1D* hDau_Kin_Phi   
	TH1D* hDau_Kin_Pt    
    TH1D* hDau_Kin_Theta 
	TH1D* hDau_Kin_Eta   



= new TH1D( "hRotPhi"   ,"hRotPhi"   , 100 , -2*PI, 2*PI);
= new TH1D( "hRotPt"    ,"hRotPt"    , 100 ,  0   , 20);
= new TH1D( "hRotTheta" ,"hRotTheta" , 100 , 0, 1);
= new TH1D( "hRotEta"   ,"hRotEta"   , bin_WRTJ_Eta , low_WRTJ_Eta_Bin   , high_WRTJ_Eta_Bin);

= new TH1D( "hPhi"   ,"hPhi"   , 100 , -2*PI, 2*PI);
= new TH1D( "hPt"    ,"hPt"    , 100 ,  0   , 20);
= new TH1D( "hTheta" ,"hTheta" , 100 , 0, 1.5*PI);
= new TH1D( "hEta"   ,"hEta"   , bin_Eta , low_Eta_Bin   , high_Eta_Bin);

*/
	TH2D* hSignal[trackbin][ptbin];
	TH2D* hBckrnd[trackbin][ptbin];
	TH2D* hEPDraw;

	hEPDraw                 = new TH2D(     "hEPDraw_trk_all_ppt_all"          ,"Eta Phi Drawing" , 50       ,-1.5*etabound, 1.5*etabound,                             100       , -2*PI  , 2*PI  );
	for(int wtrk = 1; wtrk<trackbin+1; wtrk++){
		for(int wppt = 1; wppt<ptbin+1; wppt++){
			hBckrnd[wtrk-1][wppt-1] = new TH2D(Form("hBckrnd_trk_%d_ppt_%d",wtrk,wppt) ,Form("hBckrnd_trk_%d_ppt_%d",wtrk,wppt) , DeltaEBin, 0    , etabound    , DeltaPBin, -0.5*PI  , 1.5*PI);
			hSignal[wtrk-1][wppt-1] = new TH2D(Form("hSignal_trk_%d_ppt_%d",wtrk,wppt) ,Form("hSignal_trk_%d_ppt_%d",wtrk,wppt) , DeltaEBin, 0    , etabound    , DeltaPBin, -0.5*PI  , 1.5*PI);

			//hBckrnd_sub[wtrk-1][wppt-1] = new TH2D(Form("hBckrnd_sub_trk_%d_ppt_%d",wtrk,wppt) ,Form("hBckrnd_sub_trk_%d_ppt_%d",wtrk,wppt) , DeltaEBin2, 1    , 4    , DeltaPBin, -0.5*PI  , 1.5*PI);
			//hSignal_sub[wtrk-1][wppt-1] = new TH2D(Form("hSignal_sub_trk_%d_ppt_%d",wtrk,wppt) ,Form("hSignal_sub_trk_%d_ppt_%d",wtrk,wppt) , DeltaEBin2, 1    , 4    , DeltaPBin, -0.5*PI  , 1.5*PI);
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
	
    //cout<<"*****************"<<endl;
    //cout<<"*****************"<<endl;
    cout<<"Total Entries is:"<<endl;
    cout<< nentries <<endl;
    //cout<<""<<endl;
    //cout<<"Run Length is:"<<endl;
    //cout<< RUNLENGTH << " ( " << floor(100*RUNLENGTH/nentries) << "% )" <<endl;
    //long int progressmeter1 = floor(RUNLENGTH/400);


	for (Long64_t ievent=0; ievent <nentries; ievent ++){
	//for (Long64_t ievent=0; ievent <RUNLENGTH; ievent ++){
	        timer.StartSplit("Loading Events");
		Long64_t jevent = LoadTree(ievent);
		//if (jevent < 0) break;
		nb = fChain->GetEntry(ievent);   nbytes += nb;

		//if (Cut(ievent) < 0) continue;

		//float percentdone =(float)100*(float)ievent/(float)RUNLENGTH;
		//if(ievent%progressmeter1==0) cout<< percentdone  << setprecision(2) << fixed  <<endl;
		if(ievent%1000==0) cout<< ievent << "/" << nentries  <<endl;

                timer.StartSplit("Event Selection");
		if(!F_eventpass(jetPt, jetN, jetPtCut)) continue;

                timer.StartSplit("Jet Selection");
		for(int ijet=0; ijet<jetN; ijet++){
			if( !F_jetpass(jetEta, jetPt, ijet, jetPtCut)) continue;
			Pass_Trig_Jet +=1;
			long int NNtrk = (*jetNumDaughters)[ijet];
			long int Ntrig  =0;

			hJet_Kin_Eta  ->Fill((  *jetEta)[ijet]);
            hJet_Kin_Theta->Fill((*jetTheta)[ijet]);
            hJet_Kin_Phi  ->Fill((  *jetPhi)[ijet]);

			for(long int  XXtrk=0; XXtrk < NNtrk; XXtrk++ ){
				if((*dau_chg)[ijet][XXtrk] == 0) continue;

				Ntrig=Ntrig+1;
			}

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

				float beam_dau_theta = (*dau_theta)[ijet][XXtrk];
				float beam_dau_eta   = (*dau_eta)[ijet][XXtrk];
				float beam_dau_phi   = (*dau_phi)[ijet][XXtrk];
				float beam_dau_pt    = (*dau_pt)[ijet][XXtrk];
                 
                hDau_Kin_Theta->Fill(beam_dau_theta);
		 		hDau_Kin_Eta->Fill(beam_dau_eta);
		 		hDau_Kin_Phi->Fill(beam_dau_phi);
		 		hDau_Kin_Pt->Fill(beam_dau_pt);
                 
		 		hEtaPhi_dd_Dau_Jet->Fill(Eta_d_Dau_Jet,Phi_d_Dau_Jet);
                 
		 		hEPDraw->Fill(THISeventEta, THISeventPhi, 1./Ntrig);

		 		int a = floor(Ntrig/tracksplitter);
		 		if(a > (trackbin-1)) a = (trackbin-1);

		 		int b_1 = floor(ptWRTJet(  (double) (*jetPt)[ijet] , (double) (*jetEta)[ijet] , (double) (*jetPhi)[ijet] , (double)(*dau_pt)[ijet][XXtrk] , (double)(*dau_eta)[ijet][XXtrk], (double)(*dau_phi)[ijet][XXtrk] ));
		 		int b_2 = b_1+1;
		 		if(b_1 > (ptbin-1)){
		 			b_1 = (ptbin-1);
		 			b_2 = 10000;	
		 		}
				
		 		for(long int  YYtrk=0; YYtrk< NNtrk; YYtrk++ ){
		 		if(YYtrk == XXtrk) continue;
                	if((*dau_chg)[ijet][YYtrk] == 0) continue;
		
                    int b_3 = floor(ptWRTJet(  (double) (*jetPt)[ijet] , (double) (*jetEta)[ijet] , (double) (*jetPhi)[ijet] , (double)(*dau_pt)[ijet][YYtrk] , (double)(*dau_eta)[ijet][YYtrk], (double)(*dau_phi)[ijet][YYtrk] ));

					if(b_3 > (ptbin-1)){
                    b_3 = (ptbin-1);
                    }
					if(b_3 != b_1){ 
					skipper +=1;
					continue;
					}
					tripper +=1;
					
		 			float YYTHISeventEta = etaWRTJet( (double) (*jetPt)[ijet] , (double) (*jetEta)[ijet] , (double) (*jetPhi)[ijet] , (double)(*dau_pt)[ijet][YYtrk] , (double)(*dau_eta)[ijet][YYtrk], (double)(*dau_phi)[ijet][YYtrk]);
		 			float YYTHISeventPhi = phiWRTJet( (double) (*jetPt)[ijet] , (double) (*jetEta)[ijet] , (double) (*jetPhi)[ijet] , (double)(*dau_pt)[ijet][YYtrk] , (double)(*dau_eta)[ijet][YYtrk], (double)(*dau_phi)[ijet][YYtrk]);
                 
		 			float deltaEta = THISeventEta - YYTHISeventEta;
		 			float deltaPhi = THISeventPhi - YYTHISeventPhi;

		 			if(deltaPhi<0) deltaPhi+=2*PI;
		 			if(deltaPhi>1.5*PI) deltaPhi -= 2*PI;
                 
					hSignal[a][b_1]->Fill(deltaEta, deltaPhi, 1./Ntrig);
					hSignal[a][b_1]->Fill(-deltaEta, deltaPhi, 1./Ntrig);

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

	//int backMult =1;
	int backMult =8;
	for(int wtrk = 1; wtrk < trackbin+1; wtrk++){
		for(int wppt = 1; wppt < ptbin+1; wppt++){

        cout << "ppt is " << wppt << " and trk is " << wtrk << endl;

    //TRandom * rng =3;
    double integral;
    Int_t nbinsx = hEPDraw->GetNbinsX();
    Int_t nbinsy = hEPDraw->GetNbinsY();
    Int_t nbins  = nbinsx*nbinsy;
    integral = hEPDraw->ComputeIntegral(true);
    

			for(long int ix=0; ix<backMult*(hSignal[wtrk-1][wppt-1]->GetEntries()); ix++){

				double WEta1, WPhi1;

    Int_t ibin1 = TMath::BinarySearch(nbins,integral,(Double_t) 3);
    Int_t biny1 = ibin1/nbinsx;
    Int_t binx1 = ibin1 - nbinsx*biny;
    x = fXaxis.GetBinLowEdge(binx1+1);
    if (r1 > integral[ibin1]) x +=
       fXaxis.GetBinWidth(binx1+1)*(r1-integral[ibin1])/(integral[ibin1+1] - integral[ibin1]);
    y = fYaxis.GetBinLowEdge(biny1+1) + fYaxis.GetBinWidth(biny1+1)*rng->Rndm();

				double WEta2, WPhi2; 

    Int_t ibin2 = TMath::BinarySearch(nbins,integral,(Double_t) 3);
    Int_t biny2 = ibin2/nbinsx;
    Int_t binx2 = ibin2 - nbinsx*biny;
    x = fXaxis.GetBinLowEdge(binx+1);
    if (r1 > integral[ibin]) x +=
       fXaxis.GetBinWidth(binx+1)*(r1-integral[ibin])/(integral[ibin+1] - integral[ibin]);
    y = fYaxis.GetBinLowEdge(biny+1) + fYaxis.GetBinWidth(biny+1)*rng->Rndm();





 }
				double WdeltaEta = WEta1 - WEta2;
				double WdeltaPhi = WPhi1 - WPhi2;

                //EXhBckrnd[wtrk-1][wppt-1]->Fill(WdeltaEta, WdeltaPhi, 1./backMult);
                //EXhBckrnd[wtrk-1][wppt-1]->Fill(-WdeltaEta, WdeltaPhi, 1./backMult);

				if(WdeltaPhi<0) WdeltaPhi+=2*PI;
				if(WdeltaPhi>1.5*PI) WdeltaPhi -= 2*PI;

				hBckrnd[wtrk-1][wppt-1]->Fill(WdeltaEta, WdeltaPhi, 1./backMult);
                hBckrnd[wtrk-1][wppt-1]->Fill(-WdeltaEta, WdeltaPhi, 1./backMult);

				//if(fabs(WdeltaEta) <etabound && fabs(WdeltaEta) >1.1) hBckrnd_sub[wtrk-1][wppt-1]->Fill(fabs(WdeltaEta), WdeltaPhi, 1./backMult);
				//hBckrnd[wtrk-1][wppt-1]->Scale(1./backMult);
    		}
		}
	}

    hDau_Kin_WRTJ_Eta->Scale(1./Pass_Trig_Jet);
    hDau_Kin_Eta->Scale(1./Pass_Trig_Jet);
    
    hDau_Kin_WRTJ_Eta->Scale(1./( hDau_Kin_WRTJ_Eta->GetXaxis()->GetBinWidth(3)  ));
    hDau_Kin_Eta->Scale(1./(           hDau_Kin_Eta->GetXaxis()->GetBinWidth(3)  ));

	TFile* fS_temp = new TFile(Form("unmergedOutputs/narrow_ak8_job%d.root",job), "recreate");
	for(int wtrk =1; wtrk <trackbin+1; wtrk++){
		for(int wppt =1; wppt <ptbin+1; wppt++){
			hSignal[wtrk-1][wppt-1]->Write(Form("hSig_%d_%d",wtrk,wppt));
			hBckrnd[wtrk-1][wppt-1]->Write(Form("hBck_%d_%d",wtrk,wppt));

			//hSignal_sub[wtrk-1][wppt-1]->Write(Form("hSig_sub_%d_%d",wtrk,wppt));
			//hBckrnd_sub[wtrk-1][wppt-1]->Write(Form("hBck_sub_%d_%d",wtrk,wppt));

            //hSignal[wtrk-1][wppt-1]->ProjectionY("",20,80)->Write(Form("hSig_profile_%d_%d",wtrk,wppt));

            //hBckrnd[wtrk-1][wppt-1]->ProjectionY("",20,80)->Write(Form("hBck_profile_%d_%d",wtrk,wppt));

		}
	}
  
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
