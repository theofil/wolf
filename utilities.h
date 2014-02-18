#include <iostream>
#include <string>
#include <vector>
#include <sstream>
#include <iomanip>
#include <fstream>
#include <algorithm>
#include <map>
#include <set>


#include "TFile.h"
#include "TLine.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TCut.h"
#include "TTree.h"
#include "TH1I.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TF1.h"
//#include "TMath.h"
#include "TGaxis.h"
#include "TLegend.h"
#include "THStack.h"
#include "TPaveText.h"
#include "TColor.h"
#include "TGraphAsymmErrors.h"
#include "TLorentzVector.h"
#include "TProfile.h"
#include "TBox.h"
#include "TGraphErrors.h"
#include "TPaveStats.h"
#include "TPolyLine.h"
#include "setTDRStyle.C" // README we need this in the same dir  "./" with the root macro or make c++ to know where to find it
#include "pnumber.h"
#include "TRandom3.h"
#include "TChain.h"

using namespace std;

template<class T> std::string any2string(T i){std::ostringstream buffer;    buffer << i;    return buffer.str();}
template<class T> std::string any2string2(T i){std::ostringstream buffer;    buffer << setprecision(3) << i;    return buffer.str();}
float ratioError(float a,float b);
float ratioError(float a,float b, float da , float db);
float productError(float a,float b);
float productError(float a,float b, float da , float db);
float histIntegral(TH1F *hist, float minMet);
float histIntegral(TH1F *hist, float minMet, float maxMet);
float histIntegralAndError(TH1F *hist, float minMet, float & error);
float histIntegralAndError(TH1F *hist, float minMet, float maxMet, float & error);
void errorBand(TH1F*);
void statErrorBand(TH1F*);
TH1F *squareRootHist(TH1F*);
float significance(float,float);
//TH1F  *effRatio(TH1F*h1,TH1F*h2,TH1F*h3);
float effRatioError(float A, float B);
float statErrorN(float x){return x - 0.5*TMath::ChisquareQuantile(0.3173/2,2*x);}
float statErrorP(float x){return 0.5*TMath::ChisquareQuantile(1-0.3173/2,2*(x+1))-x;}

float wjetsXS = 37509.0;
float ttbarXS = 227;
float ttbarFullLepXS = 22.4;
float ttbarSemiLepXS = 88.6;
float ttbarHadronicXS = 88.8;
float dym50XS = 3503.71;
float ZJetsCrossSection = 3503.71;
float LowZJetsCrossSection = 11050*0.069*1.15; // LO xs * filter eff * k-factor (same as ZJets) should be updated to 11050*0.069*1.18

// new since Dec 2013, 3*1177.3 = 3531.9

/*

PREP cross section (LO)
TTJets_FullLeptMGDecays_8TeV-madgraph-tauola 13.43
TTJets_SemiLeptMGDecays_8TeV-madgraph-tauola 53.2
TTJets_HadronicMGDecays_8TeV-madgraph        53.3
TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola 136.3

MCFM 225.197

CMS(2L) = 227 

DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball LO 2950.0 pb
DYJetsToLL_M-10To50filter_8TeV-madgraph LO 11050.0pb  (filter 0.069)
*/


float significance(float n, float b)
{
   return (n-b)/sqrt(n + b);
}


void errorBand(TH1F*hist)
{
 for(int i =0 ; i <= hist->GetNbinsX()+1;i++)
 {
   //hist->SetBinError(i,hist->GetBinContent(i)*0.1);
   float sysError  =  hist->GetBinContent(i)*0.1; // 10% sys uncertainty
   float statError =  sqrt(hist->GetBinContent(i));
   float totError  =  sqrt(sysError*sysError + statError*statError); 
   hist->SetBinError(i,totError);
 } 
}

void statErrorBand(TH1F*hist)
{
 for(int i =0 ; i <= hist->GetNbinsX()+1;i++)
 {
   //hist->SetBinError(i,hist->GetBinContent(i)*0.1);
   float sysError  =  0; // 0.0% sys uncertainty
   float statError =  sqrt(hist->GetBinContent(i));
   float totError  =  sqrt(sysError*sysError + statError*statError); 
   hist->SetBinError(i,totError);
 } 
}

TH1F  *squareRootHist(TH1F *hist)
{
 TH1F *newHist = (TH1F*)hist->Clone();
 for(int i =0 ; i <= newHist->GetNbinsX()+1;i++)
 {
   float binContent = newHist->GetBinContent(i);
   float binError  = newHist->GetBinError(i);
   
   float myBinContent=0; 
   float myBinError=0;
   if(binContent>0) 
   {
     myBinContent = sqrt(binContent);
     myBinError = 0.5*binError*(1/myBinContent);
   }
   newHist->SetBinError(i,myBinError);
   newHist->SetBinContent(i,myBinContent);
 } 
 return newHist;
}

void normHist(TH1F *hist)
{
 hist->Scale(1.0/hist->Integral());
 hist->GetYaxis()->SetTitle("PDF");
 hist->SetLineWidth(2);
}

TH1F  *normHistPointer(TH1F *hist)
{
 TH1F *newHist = (TH1F*)hist->Clone();
 newHist->Scale(1.0/hist->Integral());
 newHist->GetYaxis()->SetTitle("PDF");
 return newHist;
}

TH1F  *h1fclone(TH1F *hist)
{
 TH1F *newHist = (TH1F*)hist->Clone();
 return newHist;
}

void histStyle(TH1F*hist)
{
 hist->GetYaxis()->SetTitleSize(0.06);
 hist->GetYaxis()->SetTitleFont(42);
 hist->GetYaxis()->SetTitleOffset(1.05);
 hist->GetYaxis()->SetLabelSize(0.05);
 hist->GetYaxis()->SetLabelFont(42);
 hist->GetYaxis()->SetLabelOffset(0.007);
 hist->GetXaxis()->SetTitleSize(0.06);
 hist->GetXaxis()->SetTitleFont(42);
 hist->GetXaxis()->SetTitleOffset(1.05);
 hist->GetXaxis()->SetLabelSize(0.05);
 hist->GetXaxis()->SetLabelFont(42);
 hist->GetXaxis()->SetLabelOffset(0.007);
}
void histStyle(TH2F*hist)
{
 hist->GetYaxis()->SetTitleSize(0.06);
 hist->GetYaxis()->SetTitleFont(42);
 hist->GetYaxis()->SetTitleOffset(1.05);
 hist->GetYaxis()->SetLabelSize(0.05);
 hist->GetYaxis()->SetLabelFont(42);
 hist->GetYaxis()->SetLabelOffset(0.007);
 hist->GetXaxis()->SetTitleSize(0.06);
 hist->GetXaxis()->SetTitleFont(42);
 hist->GetXaxis()->SetTitleOffset(1.05);
 hist->GetXaxis()->SetLabelSize(0.05);
 hist->GetXaxis()->SetLabelFont(42);
 hist->GetXaxis()->SetLabelOffset(0.007);
}

float ratioError(float a,float b)
{
//  float res=0;
  return sqrt(a/(a*a) + b/(b*b))*(a/b);
}

float ratioError(float a,float b, float da , float db)
{
  return sqrt(pow(da/a,2) + pow(db/b,2))*(a/b);
}


float productError(float a,float b)
{
//  float res=0;
  return sqrt(a/(a*a) + b/(b*b))*(a*b);
}

float productError(float a,float b, float da , float db)
{
  return sqrt(pow(da/a,2) + pow(db/b,2))*(a*b);
}

float histIntegral(TH1F *hist, float minMet)
{
  //      bin = 0;       underflow bin
  //      bin = 1;       first bin with low-edge xlow INCLUDED
  //      bin = nbins;   last bin with upper-edge xup EXCLUDED
  //      bin = nbins+1; overflow bin

      float val=0.;
      val = hist->Integral(hist->FindBin(double(minMet)),hist->GetNbinsX()+1);
 //             if(hist->GetBinContent(hist->GetNbinsX()+1)!=0)cout<< " hist " << hist->GetTitle() << " overflowed " << endl;

                  return val;
}

float histIntegral(TH1F *hist, float minMet,float maxMet)
{
  //      bin = 0;       underflow bin
  //      bin = 1;       first bin with low-edge xlow INCLUDED
  //      bin = nbins;   last bin with upper-edge xup EXCLUDED
  //      bin = nbins+1; overflow bin

      float val=0.;
          val = hist->Integral(hist->FindBin(double(minMet)),hist->FindBin(double(maxMet)));
   //           if(hist->GetBinContent(hist->GetNbinsX()+1)!=0)cout<< " hist " << hist->GetTitle() << " overflowed " << endl;
                  return val;
}



float histIntegralAndError(TH1F *hist, float minMet, float &error)
{
  //      bin = 0;       underflow bin
  //      bin = 1;       first bin with low-edge xlow INCLUDED
  //      bin = nbins;   last bin with upper-edge xup EXCLUDED
  //      bin = nbins+1; overflow bin

      float val=0.;
          val = hist->Integral(hist->FindBin(double(minMet)),hist->GetNbinsX()+1);
              double valError=0.;
                  hist->IntegralAndError(hist->FindBin(double(minMet)),hist->GetNbinsX()+1,valError);
                      error = float(valError);
                          return val;
}

float effRatioError(float A, float B)
{
   float eff = A/B;
   float tmpError2 = eff*(1-eff)/B;
   return sqrt(tmpError2);
} 

/*
TH1F *effRatio(TH1F*h1,TH1F*h2)
{
 TH1F * h3 = (TH1F*) h1->Clone();
 for(int i =0 ; i <= h1->GetNbinsX()+1;i++)
 {
   float ratio=0;
   float ratioError=0;
   float A= h1->GetBinContent(i);
//   float AError = h1->GetBinError(i);
   float B= h2->GetBinContent(i);
//   float BError = h2->GetBinError(i);

   //float tmpError2 = pow(  pow((1-B/A),-2)*(B/(A*A))*AError, 2) + pow(  pow((1-B/A),-2)*(1/A)*BError,  2);
   float eff = A/B;
   float tmpError2 = eff*(1-eff)/B;

   if(B!=0)
   {
     ratio = A/B;
     ratioError = sqrt(tmpError2);
   }else{continue;}


   h3->SetBinContent(i,ratio);
   h3->SetBinError(i,ratioError);
 }
 return h3;
}
*/


float histIntegralAndError(TH1F *hist, float minMet, float maxMet, float &error)
{
  //      bin = 0;       underflow bin
  //      bin = 1;       first bin with low-edge xlow INCLUDED
  //      bin = nbins;   last bin with upper-edge xup EXCLUDED
  //      bin = nbins+1; overflow bin

      float val=0.;
          //val = hist->Integral(hist->FindBin(double(minMet)),hist->GetNbinsX()+1);
          val = hist->Integral  (hist->FindBin(double(minMet)),hist->FindBin(double(maxMet)));
          double valError=0.;
          hist->IntegralAndError(hist->FindBin(double(minMet)),hist->FindBin(double(maxMet)),valError);
          error = float(valError);
          return val;
}
