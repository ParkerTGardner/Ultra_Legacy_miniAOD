#define My2Class_cxx
#include "My2Class.h"

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
#include "coordinateTools.h"

using TMath::ATan;
using TMath::Exp;

//---------------------------------------------------------------------Global Vars

const int           asize=70;
const long int      bsize=1000;
const double        PI = 3.14159265359;
const int           Resolution = 40;
const long int      RUNLENGTH = 20000;
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
//const float PtAsshi     = 1000; 
//const float PtAsslow    = 0.0;
const float Thetahigh   = 0.55;
const float Thetalow    = 0.0;
const float EtaCut      = 0.0;
const float jetEtaCut   = 2.1;
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

void My2Class::Loop(){
	//void MyClass(){
	//std::cout << "made it 1" <<std::endl;
	//   In a ROOT session, you can do:
	//      root> .L MyClass.C
	//      root> MyClass t
	//      root> t.GetEntry(12); // Fill t data members with entry number 12
	//      root> t.Show();       // Show values of entry 12
	//      root> t.Show(16);     // Read and show values of entry 16
	//      root> t.Loop();       // Loop on all entries
	//

	//     This is the loop skeleton where:
	//    jentry is the global entry number in the chain
	//    ientry is the entry number in the current Tree
	//  Note that the argument to GetEntry must be:
	//    jentry for TChain::GetEntry
	//    ientry for TTree::GetEntry and TBranch::GetEntry
	//
	//       To read only selected branches, Insert statements like:
	// METHOD1:
	//    fChain->SetBranchStatus("*",0);  // disable all branches
	//    fChain->SetBranchStatus("branchname",1);  // activate branchname
	// METHOD2: replace line
	//    fChain->GetEntry(jentry);       //read all branches
	//by  b_branchname->GetEntry(ientry); //read only this branch
	//if (fChain == 0) return;
	//std::cout << "made it 1" <<std::endl;



	//-------------------------------------------------------------------MAIN_CODE-------------------------------------------------------------------------------
	//-------------------------------------------------------------------MAIN_CODE-------------------------------------------------------------------------------

	/*
	   TH2D* hRandomBack  = new TH2D("hRandomBack"  ,"Background (long vect)"              ,        DeltaEBin, -etabound, etabound, DeltaPBin, -0.5*PI, 1.5*PI);
	   TH2D* hWRandomBack = new TH2D("hWRandomBack" ,"Background (digitized hist draw)"    ,        DeltaEBin, -etabound, etabound, DeltaPBin, -0.5*PI, 1.5*PI);
	   TH2D* hEtaPhiDraw  = new TH2D("hEtaPhiDraw"  ,"2D EtaPhi drawing"                   ,        50       , -1.5*etabound, 1.5*etabound, 50,-2.0*PI, 2.0*PI);
	   TH2D* hSigBEAM     = new TH2D("hSigBEAM"     ,"2D Signal BEAM"                      ,        DeltaEBin, -etabound, etabound, DeltaPBin, -0.5*PI, 1.5*PI);
	   TH2D* hREVWRandomBack = new TH2D("hREVWRandomBack" ,"REV Background (digitized hist draw)"    ,        DeltaEBin, -etabound, etabound, DeltaPBin, -0.5*PI, 1.5*PI);
	   TH2D* hREVEtaPhiDraw  = new TH2D("hREVEtaPhiDraw"  ,"2D REV EtaPhi drawing"                   ,        50       , -1.5*etabound, 1.5*etabound, 50,-2.0*PI, 2.0*PI);
	   TH2D* hREVSigBEAM     = new TH2D("hREVSigBEAM"     ,"2D REV Signal BEAM"                      ,        DeltaEBin, -etabound, etabound, DeltaPBin, -0.5*PI, 1.5*PI);
	   */

	int NMinus = 0;
	int NPlus  = 0;
	int NEq    = 0;
	int NMinusB = 0;
	int NPlusB  = 0;
	int NEqB    = 0;

	int Pass_Trig_Event =0;
	int Pass_Trig_Jet   =0;
	int Pass_Back_Event =0;
	int Pass_Back_Jet   =0;

	int skipper =0;
	int tripper =0;


	//-----------------------------------------------------------------Canvas, Print  and Hist
	//const long int ALLEVENTS = trackTree->GetEntries();
	//cout << "All events  =  " << ALLEVENTS <<endl;
	//long int nevents =(ALLEVENTS);
	//cout << "Run events  =  " << ALLEVENTS << " ( " << floor(100*ALLEVENTS/ALLEVENTS) << " % )" <<endl;
	//long int progressmeter1 = floor(nevents/10);

	TH1D* hJet_Kin_Eta    = new TH1D("hJet_Kin_Eta","hJet_Kin_Eta" ,30  ,-15    ,15);
    TH1D* hJet_Kin_Phi    = new TH1D("hJet_Kin_Phi","hJet_Kin_Phi" ,30  ,-2*PI ,2*PI);
    TH1D* hJet_Kin_Theta  = new TH1D("hJet_Kin_Theta" ,"hJet_Kin_Theta"  ,30  ,-2*PI ,2*PI);

	TH1D* hDau_Kin_WRTJ_Theta = new TH1D( "hRotTheta" ,"hRotTheta" , 100 , 0, 1);
	TH1D* hDau_Kin_WRTJ_Eta   = new TH1D( "hRotEta"   ,"hRotEta"   , 20 , -5   , 15);
	TH1D* hDau_Kin_WRTJ_Phi   = new TH1D( "hRotPhi"   ,"hRotPhi"   , 100 , -2*PI, 2*PI);
    TH1D* hDau_Kin_WRTJ_Pt    = new TH1D( "hRotPt"    ,"hRotPt"    , 100 ,  0   , 20);

    TH1D* hDau_Kin_Theta = new TH1D( "hTheta" ,"hTheta" , 100 , 0, 1.5*PI);
	TH1D* hDau_Kin_Eta   = new TH1D( "hEta"   ,"hEta"   , 20 , -5   , 15);
	TH1D* hDau_Kin_Phi   = new TH1D( "hPhi"   ,"hPhi"   , 100 , -2*PI, 2*PI);
	TH1D* hDau_Kin_Pt    = new TH1D( "hPt"    ,"hPt"    , 100 ,  0   , 20);


    //h1= new TH2D(Form("hSignal_trk_%d_ppt_%d",wtrk,wppt) ,Form("hSignal_trk_%d_ppt_%d",wtrk,wppt) , DeltaEBin, -etabound    , etabound    , DeltaPBin, -0.5*PI  , 1.5*PI);



	TH2D* hEtaPhi_dd_Dau_Jet = new TH2D("hEtaPhi_dd_Dau_Jet","hEtaPhi_dd_Dau_Jet", 50, -1 , 1 , 50 , -1  , 1);

	int trackbin 		= 4;
	int ptbin		    = 3;
	int tracksplitter 	= 20;

	TH2D* hBckrnd[trackbin][ptbin];
	TH2D* hBckrnd_sub[trackbin][ptbin];
	TH2D* hEPDraw;
	TH2D* hSignal[trackbin][ptbin];
	TH2D* hSignal_sub[trackbin][ptbin];

	hEPDraw                 = new TH2D(     "hEPDraw_trk_all_ppt_all"          ,"Eta Phi Drawing" , 50       , -1.5*etabound, 1.5*etabound,                             100       , -2*PI  , 2*PI  );
	for(int wtrk = 1; wtrk<trackbin+1; wtrk++){
		for(int wppt = 1; wppt<ptbin+1; wppt++){
			hBckrnd[wtrk-1][wppt-1] = new TH2D(Form("hBckrnd_trk_%d_ppt_%d",wtrk,wppt) ,Form("hBckrnd_trk_%d_ppt_%d",wtrk,wppt) , DeltaEBin, 0    , etabound    , DeltaPBin, -0.5*PI  , 1.5*PI);
			hSignal[wtrk-1][wppt-1] = new TH2D(Form("hSignal_trk_%d_ppt_%d",wtrk,wppt) ,Form("hSignal_trk_%d_ppt_%d",wtrk,wppt) , DeltaEBin, 0    , etabound    , DeltaPBin, -0.5*PI  , 1.5*PI);

			hBckrnd_sub[wtrk-1][wppt-1] = new TH2D(Form("hBckrnd_sub_trk_%d_ppt_%d",wtrk,wppt) ,Form("hBckrnd_sub_trk_%d_ppt_%d",wtrk,wppt) , DeltaEBin2, 1    , 4    , DeltaPBin2, 0  , 2);
			hSignal_sub[wtrk-1][wppt-1] = new TH2D(Form("hSignal_sub_trk_%d_ppt_%d",wtrk,wppt) ,Form("hSignal_sub_trk_%d_ppt_%d",wtrk,wppt) , DeltaEBin2, 1    , 4    , DeltaPBin2, 0  , 2);
		}
	}

/*


    TH2D* EXhBckrnd[trackbin][ptbin];
    TH2D* EXhBckrnd_sub[trackbin][ptbin];
    TH2D* EXhEPDraw;
    TH2D* EXhSignal[trackbin][ptbin];
    TH2D* EXhSignal_sub[trackbin][ptbin];

    EXhEPDraw                 = new TH2D(     "EXhEPDraw_trk_all_ppt_all"          ,"EXEta EXPhi EXDrawing" ,                         50          , -1.5*etabound, 1.5*etabound, 100       , -3*PI  , 3*PI  );
    for(int wtrk = 1; wtrk<trackbin+1; wtrk++){
        for(int wppt = 1; wppt<ptbin+1; wppt++){
            EXhBckrnd[wtrk-1][wppt-1] = new TH2D(Form("EXhBckrnd_trk_%d_ppt_%d",wtrk,wppt) ,Form("EXhBckrnd_trk_%d_ppt_%d",wtrk,wppt) , DeltaEBin, -etabound    , etabound    , DeltaPBin, -3*PI  , 3*PI);
            EXhSignal[wtrk-1][wppt-1] = new TH2D(Form("EXhSignal_trk_%d_ppt_%d",wtrk,wppt) ,Form("EXhSignal_trk_%d_ppt_%d",wtrk,wppt) , DeltaEBin, -etabound    , etabound    , DeltaPBin, -3*PI  , 3*PI);
        }
    }

*/







	/*

	   TH1D* hRotEta                 = new TH1D(     "hRotEta"          ,"hRotEta" , 10       , -5, 5);
	   TH1D* hRotPhi                 = new TH1D(     "hRotPhi"          ,"hRotPhi" , 50       , -PI,PI);
	   TH1D* hRotPt                  = new TH1D(     "hRotPt"          ,"hRotPt"   , 100       , 0,20);

	   TH1D* hEta                 = new TH1D(     "hEta"          ,"hEta" , 10       , -5, 5);
	   TH1D* hPhi                 = new TH1D(     "hPhi"          ,"hPhi" , 50       , -PI,PI);
	   TH1D* hPt                  = new TH1D(     "hPt"          ,"hPt"   , 100       , 0,20);
	   */


	//-------------------------------------------------------------------MAIN_LOOPS
	//-------------------------------------------------------------------MAIN_LOOPS

	Long64_t nentries = fChain->GetEntriesFast();
    cout<<"*****************"<<endl;
    cout<<"*****************"<<endl;
    cout<<"Total Entries is:"<<endl;
    cout<< nentries <<endl;
    cout<<""<<endl;
    cout<<"Run Length is:"<<endl;
    cout<< RUNLENGTH << " ( " << floor(100*RUNLENGTH/nentries) << "% )" <<endl;
    long int progressmeter1 = floor(RUNLENGTH/400);

    
	Long64_t nbytes = 0, nb = 0;
	for (Long64_t ievent=0; ievent <RUNLENGTH; ievent ++){
		Long64_t jevent = LoadTree(ievent);
		//if (jevent < 0) break;
		nb = fChain->GetEntry(ievent);   nbytes += nb;

		//if (Cut(ievent) < 0) continue;

		float percentdone =(float)100*(float)ievent/(float)RUNLENGTH;
		if(ievent%progressmeter1==0) cout<< percentdone  << setprecision(2) << fixed  <<endl;

		if(!F_eventpass(jetPt, jetN, jetPtCut)) continue;

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

					if(fabs(deltaEta) <etabound && fabs(deltaEta) >1.1) hSignal_sub[a][b_1]->Fill(fabs(deltaEta), deltaPhi, 1./Ntrig);

				}
			}
		}
	}

	cout << "skipper is: " << skipper <<endl;
	cout << "tripper is: " << tripper <<endl;

	int backMult =10;
	for(int wtrk = 1; wtrk < trackbin+1; wtrk++){
		for(int wppt = 1; wppt < ptbin+1; wppt++){

            cout<< wtrk  << "and " << wppt  <<endl;

			for(long int ix=0; ix<backMult*(hSignal[wtrk-1][wppt-1]->GetEntries()); ix++){

				double WEta1, WPhi1;
				hEPDraw->GetRandom2(WEta1, WPhi1);
				double WEta2, WPhi2;
				hEPDraw->GetRandom2(WEta2, WPhi2);

				double WdeltaEta = WEta1 - WEta2;
				double WdeltaPhi = WPhi1 - WPhi2;

                //EXhBckrnd[wtrk-1][wppt-1]->Fill(WdeltaEta, WdeltaPhi, 1./backMult);
                //EXhBckrnd[wtrk-1][wppt-1]->Fill(-WdeltaEta, WdeltaPhi, 1./backMult);

				if(WdeltaPhi<0) WdeltaPhi+=2*PI;
				if(WdeltaPhi>1.5*PI) WdeltaPhi -= 2*PI;

				hBckrnd[wtrk-1][wppt-1]->Fill(WdeltaEta, WdeltaPhi, 1./backMult);
                hBckrnd[wtrk-1][wppt-1]->Fill(-WdeltaEta, WdeltaPhi, 1./backMult);
            
				if(fabs(WdeltaEta) <etabound && fabs(WdeltaEta) >1.1) hBckrnd_sub[wtrk-1][wppt-1]->Fill(fabs(WdeltaEta), WdeltaPhi, 1./backMult);

				//hBckrnd[wtrk-1][wppt-1]->Scale(1./backMult);
			}
		}
	}

    cout<< "NOW SAVING" <<endl;

	TFile* fS_temp = new TFile("test2k.root", "recreate");
	for(int wtrk =1; wtrk <trackbin+1; wtrk++){
		for(int wppt =1; wppt <ptbin+1; wppt++){
			hSignal[wtrk-1][wppt-1]->Write(Form("hSig_%d_%d",wtrk,wppt));
			hBckrnd[wtrk-1][wppt-1]->Write(Form("hBck_%d_%d",wtrk,wppt));

			hSignal_sub[wtrk-1][wppt-1]->Write(Form("hSig_sub_%d_%d",wtrk,wppt));
			hBckrnd_sub[wtrk-1][wppt-1]->Write(Form("hBck_sub_%d_%d",wtrk,wppt));

            hSignal[wtrk-1][wppt-1]->ProjectionY("",20,80)->Write(Form("hSig_profile_%d_%d",wtrk,wppt));

            hBckrnd[wtrk-1][wppt-1]->ProjectionY("",20,80)->Write(Form("hBck_profile_%d_%d",wtrk,wppt));
            cout<< wtrk  << "and " << wppt  <<endl;
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

    TCanvas* c1          = new TCanvas("c1","c1" , 800, 800);
    c1->Divide(ptbin,trackbin);
    TCanvas* c2          = new TCanvas("c2","c2" , 800, 800);
    c2->Divide(ptbin,trackbin);
    TCanvas* c3          = new TCanvas("c3","c3" , 800, 800);
    c3->Divide(ptbin,trackbin);
    TCanvas* c4          = new TCanvas("c4","c4" , 800, 800);
    c4->Divide(ptbin,trackbin);


    for(int wtrk =1; wtrk <trackbin+1; wtrk++){
        for(int wppt =1; wppt <ptbin+1; wppt++){
            int where = wppt + (wtrk-1)*(ptbin);
            c1->cd(where);
            hSignal[wtrk-1][wppt-1]->Divide(hBckrnd[wtrk-1][wppt-1]);
            hSignal[wtrk-1][wppt-1]->DrawCopy("COLZ");
            c2->cd(where);
            hSignal[wtrk-1][wppt-1]->ProjectionY("",20,80)->DrawCopy();


            c3->cd(where);
            hSignal_sub[wtrk-1][wppt-1]->Divide(hBckrnd_sub[wtrk-1][wppt-1]);
            hSignal_sub[wtrk-1][wppt-1]->DrawCopy("COLZ");
            c4->cd(where);
            hSignal_sub[wtrk-1][wppt-1]->ProjectionY()->DrawCopy();





        }
    }
*/
/*
    TCanvas* c2          = new TCanvas("c2","c2" , 800, 800);
    c2->Divide(ptbin,trackbin);
    //TCanvas* c2          = new TCanvas("c2","c2" , 800, 800);
    //c2->Divide(ptbin,trackbin);

    for(int wtrk =1; wtrk <trackbin+1; wtrk++){
        for(int wppt =1; wppt <ptbin+1; wppt++){
            int where = wppt + (wtrk-1)*(ptbin);
	        c2->cd(where);
            EXhSignal[wtrk-1][wppt-1]->Divide(EXhBckrnd[wtrk-1][wppt-1]);
            EXhSignal[wtrk-1][wppt-1]->DrawCopy("COLZ");
	    }
	}
*/
/*
	TCanvas* c1          = new TCanvas("c1","c1" , 800, 800);
	c1->Divide(2);
	TCanvas* cB          = new TCanvas("cB","cB" , 800, 800);
	cB->Divide(2);
	TCanvas* c3          = new TCanvas("c3","c3" , 800, 800);
	c3->Divide(2);
        TCanvas* c4          = new TCanvas("c4","c4" , 800, 800);
        c4->Divide(2);

        TCanvas* c5          = new TCanvas("c5","c5" , 800, 800);
	hEtaPhi_dd_Dau_Jet->Draw("COLZ");

	TCanvas* cA	     = new TCanvas("cA","cA" , 800, 800);
	cA->Divide(3);	
	cA->cd(1); hJet_Kin_Eta->Draw();
        cA->cd(2); hJet_Kin_Phi->Draw();
        cA->cd(3); hJet_Kin_Theta->Draw();

	c1->cd(1);
	hDau_Kin_WRTJ_Eta->Scale(1./Pass_Trig_Jet);
	hDau_Kin_WRTJ_Eta->Draw();

        TCanvas* cB          = new TCanvas("cB","cB" , 800, 800);
        cB->Divide(2);

	cB->cd(1);
	hDau_Kin_WRTJ_Theta->Draw();

	c3->cd(1);
	hDau_Kin_WRTJ_Phi->Draw();

	c4->cd(1);
	hDau_Kin_WRTJ_Pt->Draw();

	c1->cd(2);
	hDau_Kin_Eta->Scale(1./Pass_Trig_Jet);
	hDau_Kin_Eta->Draw();

	cB->cd(2);
	hDau_Kin_Theta->Draw();

	c3->cd(2);
	hDau_Kin_Phi->Draw();

	c4->cd(2);
	hDau_Kin_Pt->Draw();*/



}
