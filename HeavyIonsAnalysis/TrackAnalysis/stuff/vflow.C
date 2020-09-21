fndef __MY_INCLUDE__
#define __MY_INCLUDE__
R__ADD_INCLUDE_PATH(../include);
R__ADD_INCLUDE_PATH(../src);
#include "myAnaConsts.h"
#include "functions.cxx"
#endif

using namespace std;

typedef TH1* TH1ptr;

double lw = 0.;
double up = 0.;

void proj1D_Ref_longrange(TH2* h2DSignal, TH2D* h2DBackground, TH1ptr& h, const char* name)
{
   const double deltaEtaBound = 1.;
   int negBinMin = 0;
   //int negBinMax = h2DSignal->GetXaxis()->FindBin(-1.* deltaEtaBound);
   //int posBinMin = h2DSignal->GetXaxis()->FindBin(1.* deltaEtaBound);
   int negBinMax = h2DSignal->GetXaxis()->FindBin(-1.* deltaEtaBound)-1;
   int posBinMin = h2DSignal->GetXaxis()->FindBin(1.* deltaEtaBound)+1;
   int posBinMax = h2DSignal->GetXaxis()->GetNbins()+1;
   TH1D* hNeg = h2DSignal->ProjectionY("hneg", negBinMin, negBinMax);
   TH1D* hPos = h2DSignal->ProjectionY("hpos", posBinMin, posBinMax);
   hNeg->Add(hPos);

   TH1D* temp_neg = h2DBackground->ProjectionY("hneg_bkg", negBinMin, negBinMax);
   TH1D* temp_pos = h2DBackground->ProjectionY("hpos_bkg", posBinMin, posBinMax);
   temp_neg->Add(temp_pos);

   int center = h2DBackground->FindBin(0., 0.);
   hNeg->Divide(temp_neg);
   hNeg->Scale(h2DBackground->GetBinContent(center) / temp_neg->GetBinWidth(1)/h2DBackground->GetXaxis()->GetBinWidth(1));

   delete hPos;
   delete temp_neg;
   delete temp_pos;

   hNeg->SetName(name);
   h = hNeg;
}

void proj1D_ref_sr(TH2* h2DSignal, TH2D* h2DBackground, TH1ptr& h, const char* name)
{
   int negBin = h2DSignal->GetXaxis()->FindBin(-1.);
   int posBin = h2DSignal->GetXaxis()->FindBin(1.);
   h = h2DSignal->ProjectionY("hsig", negBin, posBin);

   TH1D* temp = h2DBackground->ProjectionY("hbkg", negBin, posBin);

   int center = h2DBackground->FindBin(0., 0.);
   h->Divide(temp);
   h->Scale(h2DBackground->GetBinContent(center) / temp->GetBinWidth(1)/h2DBackground->GetXaxis()->GetBinWidth(1));
   h->SetName(name);
}

void proj1D_ref_lr(TH2* h2DSignal, TH2D* h2DBackground, TH1ptr& h, const char* name)
{
   int negBinMin = 0;
   int negBinMax = h2DSignal->GetXaxis()->FindBin(-2.);
   int posBinMin = h2DSignal->GetXaxis()->FindBin(2.);
   int posBinMax = 33+1;

   TH1D* hNeg = h2DSignal->ProjectionY("hneg", negBinMin, negBinMax);
   TH1D* hPos = h2DSignal->ProjectionY("hpos", posBinMin, posBinMax);
   hNeg->Add(hPos);
   TH1D* hsig = new TH1D(*hNeg);
   hsig->SetName("hsig");
   delete hNeg;
   delete hPos;

   TH1D* temp_neg = h2DBackground->ProjectionY("hneg_bkg", negBinMin, negBinMax);
   TH1D* temp_pos = h2DBackground->ProjectionY("hpos_bkg", posBinMin, posBinMax);
   temp_neg->Add(temp_pos);

   TH1D* temp = new TH1D(*temp_neg);
   temp->SetName("temp");
   delete temp_neg;
   delete temp_pos;

   int center = h2DBackground->FindBin(0., 0.);
   hsig->Divide(temp);
   hsig->Scale(h2DBackground->GetBinContent(center) / temp->GetBinWidth(1)/h2DBackground->GetXaxis()->GetBinWidth(1));

   h = (TH1D*) hsig->Clone();
   h->SetName(name);

   delete hsig;
}

void drawRef_proj1D(string dataset = "", string input = "", string output = "",
      int itrk=-1, 
      const vector<unsigned int>* trkbin=nullptr, 
      unordered_map<string, TGraphErrors*>* graph_map = nullptr)
{
   //gErrorIgnoreLevel = kWarning;
   gStyle->SetOptStat(0);
   //gStyle->SetOptFit(1111);

   TFile* fin = TFile::Open(input.c_str());
   if(isFailed(fin, "input files")) return;

   TH2D* h2DSignal_Ref;
   TH2D* h2DBackground_Ref;
   TH1D* hMultRef;

   fin->GetObject(Form("hSignal_Ref_trk%d", itrk), h2DSignal_Ref);
   if(isFailed(h2DSignal_Ref, "hsignal")) return;
   fin->GetObject(Form("hBackground_Ref_trk%d", itrk), h2DBackground_Ref);
   if(isFailed(h2DBackground_Ref, "hbkg")) return;
   fin->GetObject(Form("hMult_trk%d",  itrk), hMultRef);
   if(isFailed(hMultRef, "hmult")) return;

   // scaled by event number
   auto nevents = hMultRef->Integral(3, 100000);
   h2DSignal_Ref->Scale(1./nevents);

   // projection
   TH1ptr hDeltaPhi_Ref;
   proj1D_Ref_longrange(h2DSignal_Ref, h2DBackground_Ref, hDeltaPhi_Ref, 
            Form("h1D_trk%d_%s", itrk, dataset.c_str()));

   double Vn[1];
   double Vn_err[1];

   double nass[1];
   double nass_err[1];

   TCanvas c("c", "", 550, 450);
   TLatex ltx;
   if(true){
      std::string function = "[0]/(TMath::Pi()*2)*(1+2*([1]*TMath::Cos(x)+[2]*TMath::Cos(2*x)+[3]*TMath::Cos(3*x)))";

      gPad->SetBottomMargin(0.14);
      TF1 func("deltaPhi", function.c_str(), -3.14159*0.5, 3.14159*1.5);
      func.SetParameter(0, hDeltaPhi_Ref->GetMaximum());
      func.SetParameter(1, 0.1);
      func.SetParameter(2, 0.1);
      func.SetParameter(3, 0.1);

      hDeltaPhi_Ref->SetMarkerStyle(20);

      hDeltaPhi_Ref->Fit(&func, "q");
      hDeltaPhi_Ref->Fit(&func, "q");
      hDeltaPhi_Ref->Fit(&func, "m q");
      hDeltaPhi_Ref->Fit(&func, "m q E");
      auto fitResult = hDeltaPhi_Ref->Fit(&func, "m S E q");

      hDeltaPhi_Ref->SetTitle(";#Delta#phi;");
      hDeltaPhi_Ref->GetXaxis()->CenterTitle();
      hDeltaPhi_Ref->GetXaxis()->SetTitleSize(0.05);

      auto h = hDeltaPhi_Ref->DrawCopy();
      /*
      gPad->Update();
      TPaveStats* pave = (TPaveStats*) h->FindObject("stats");
      pave->SetX1NDC(0.16);
      pave->SetX2NDC(0.56);
      gPad->Modified();
      gPad->Update();
      */

      TString ntrk_str;
      if(trkbin->at(itrk+1) == numeric_limits<unsigned int>::max()) ntrk_str = Form("N_{trk}^{offline} #geq %u", trkbin->at(itrk));
      else ntrk_str = Form("%u #leq N_{trk}^{offline} < %u", trkbin->at(itrk), trkbin->at(itrk+1));
      ltx.DrawLatexNDC(0.65, 0.9, ntrk_str.Data());
      ltx.DrawLatexNDC(0.65, 0.38, "|#Delta#eta|>1");
      //ltx.DrawLatexNDC(0.65, 0.24, "|y|<1");

      Vn[0] = func.GetParameter(2);
      Vn_err[0] = func.GetParError(2);

      nass[0] = func.GetParameter(0);
      nass_err[0] = func.GetParError(0);
   }
   c.Print(Form("%s_ref.pdf", output.c_str()));
   delete hDeltaPhi_Ref;

   // jet calculation
   // short range
   TH1ptr hsr;
   double min_sr;
   double minX_sr;
   double yields_sr;
   double yields_sr_err;
   proj1D_ref_sr(h2DSignal_Ref, h2DBackground_Ref, hsr,
         Form("h1Dsr_ref_trk%d_%s", itrk, dataset.c_str()));
   TCanvas csr("csr", "", 550, 450);
   hsr->SetMarkerStyle(20);
   if(true){
      TF1 func("func", "[0]+[1]*x+[2]*x*x", -2.0, 2.0);
      hsr->Fit(&func, "Q R", "", 0.6, 2.0);
      hsr->Fit(&func, "Q R", "", 0.6, 2.0);
      hsr->Fit(&func, "E Q R", "", 0.6, 2.0);
      auto fitResult = hsr->Fit(&func, "ESQR", "", 0.6, 2.0);
      cout << "sr fit result: " << fitResult->Status() << endl;

      hsr->DrawCopy();
      
      TString ntrk_str;
      if(trkbin->at(itrk+1) == numeric_limits<unsigned int>::max()) ntrk_str = Form("N_{trk}^{offline} #geq %u", trkbin->at(itrk));
      else ntrk_str = Form("%u #leq N_{trk}^{offline} < %u", trkbin->at(itrk), trkbin->at(itrk+1));
      ltx.DrawLatexNDC(0.65, 0.9, ntrk_str.Data());
      ltx.DrawLatexNDC(0.65, 0.38, "|#Delta#eta|<1");

      min_sr = func.GetMinimum(0.6, 2.0);
      minX_sr = func.GetMinimumX(0.6, 2.0);
   }
   csr.Print(Form("%s_ref_sr.pdf", output.c_str()));
   if(true){
      TF1 fmin("fmin", "[0]", -0.5*ana::PI, 1.5*ana::PI);
      fmin.SetParameter(0, min_sr);
      hsr->Add(&fmin, -1);
      int binlw = hsr->GetXaxis()->FindBin(0.);
      int binup = hsr->GetXaxis()->FindBin(minX_sr);
      yields_sr = hsr->IntegralAndError(binlw, binup, yields_sr_err, "width");
   }
   delete hsr;

   TH1ptr hlr;
   double min_lr;
   double minX_lr;
   double yields_lr;
   double yields_lr_err;
   proj1D_ref_lr(h2DSignal_Ref, h2DBackground_Ref, hlr,
         Form("h1Dlr_ref_trk%d_%s", itrk, dataset.c_str()));
   TCanvas clr("clr", "", 550, 450);
   hlr->SetMarkerStyle(20);
   if(true){
      TF1 func("func", "[0]+[1]*x+[2]*x*x", lw, up);
      hlr->Fit(&func, "Q R", "", lw, up);
      hlr->Fit(&func, "Q R", "", lw, up);
      hlr->Fit(&func, "E Q R", "", lw, up);
      auto fitResult = hlr->Fit(&func, "ESQR", "", lw, up);
      cout << "lr fit result: " << fitResult->Status() << endl;

      hlr->DrawCopy();

      TString ntrk_str;
      if(trkbin->at(itrk+1) == numeric_limits<unsigned int>::max()) ntrk_str = Form("N_{trk}^{offline} #geq %u", trkbin->at(itrk));
      else ntrk_str = Form("%u #leq N_{trk}^{offline} < %u", trkbin->at(itrk), trkbin->at(itrk+1));
      ltx.DrawLatexNDC(0.25, 0.8, ntrk_str.Data());
      ltx.DrawLatexNDC(0.65, 0.6, "|#Delta#eta|>2");

      min_lr = func.GetMinimum(lw, up);
      minX_lr = func.GetMinimumX(lw, up);
   }
   clr.Print(Form("%s_ref_lr.pdf", output.c_str()));
   if(true){
      TF1 fmin("fmin", "[0]", -0.5*ana::PI, 1.5*ana::PI);
      fmin.SetParameter(0, min_lr);
      hlr->Add(&fmin, -1);
      int binlw = hlr->GetXaxis()->FindBin(0.);
      int binup = hlr->GetXaxis()->FindBin(minX_lr);
      yields_lr = hlr->IntegralAndError(binlw, binup, yields_lr_err, "width");
   }
   delete hlr;

   double yields_jet[1];
   yields_jet[0] = yields_sr - yields_lr;
   double yields_jet_err[1];
   yields_jet_err[0] = sqrt(
            pow(yields_sr, 2) - pow(yields_lr, 2)
         );

   double e[1] = {1.};

   (*graph_map)["Vn"] = new TGraphErrors(1, e, Vn, e, Vn_err);
   (*graph_map)["jets"] = new TGraphErrors(1, e, yields_jet, e, yields_jet_err);
   (*graph_map)["nass"] = new TGraphErrors(1, e, nass, e, nass_err);
}

void drawRef()
{
   string dataset[] = {"PAMB", "PAHM1-6", "PAHM7"};
   string dataMult[] = {"PAMB0-185", "PAHM185-250", "PAHM250-inf"};
   string dataTrigger[] = {"PAMB", "PAHM", "PAHM"};

   string cmdDirPlots = 
   Form("if [ ! -d \"../plots\" ]; then\n"
         "    mkdir ../plots\n"
         "fi");
   string cmdDirPlotsCorr2D = 
   Form("if [ ! -d \"../plots/proj1D\" ]; then\n"
         "    mkdir ../plots/proj1D\n"
         "fi");
   gSystem->Exec(cmdDirPlots.c_str());
   gSystem->Exec(cmdDirPlotsCorr2D.c_str());

   double fit_lw[] = {0., 0.5, 0.6, 0.6, 0.7, 0.7};
   double fit_up[] = {2., 2., 2., 2., 2.1, 2.1};

   for(int iset=0; iset<3; iset++){

      cout << dataset[iset] << endl;

      auto trkbin = ana::get_Mult_Edges(dataset[iset]);
      int ntrk = trkbin.size() - 1;

      unordered_map<string, TGraphErrors*> graph_map[ntrk];

      for(int itrk=0; itrk<ntrk; itrk++){
         lw = fit_lw[itrk+ntrk*iset];
         up = fit_up[itrk+ntrk*iset];

         TString output;
         if(trkbin.at(itrk+1) == numeric_limits<unsigned int>::max()) 
            output= Form("../plots/proj1D/%s%u-inf", dataTrigger[iset].c_str(), 
                  trkbin[itrk]);
         else 
            output = Form("../plots/proj1D/%s%u-%u", dataTrigger[iset].c_str(), 
                  trkbin[itrk], trkbin[itrk+1]);

         drawRef_proj1D(dataset[iset], Form("../data/corr2D_trg_pd0_%s_ref.root", dataMult[iset].c_str()), 
                  output.Data(), itrk, &trkbin, &graph_map[itrk]);
      }
      TFile fout(Form("../data/%s_ref.root", dataMult[iset].c_str()), "recreate");
      for(int itrk=0; itrk<ntrk; itrk++){
         for(auto& graph : graph_map[itrk]) graph.second->Write(Form("%s_ref_trk%d", graph.first.c_str(), itrk));
      }
   }
}

