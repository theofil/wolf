// --- Class SuperHist
class SuperHist
{
  public:
  SuperHist(){};
  SuperHist(TH1F *hist)
  {
      hist_ = hist;
      X_Ndivisions = 509;
      Y_Ndivisions = 509;
      hist_->GetXaxis()->SetNdivisions(X_Ndivisions);
      hist_->GetYaxis()->SetNdivisions(Y_Ndivisions);
  }

  pnumber IntUp(float minMet);
  pnumber Int(float minMet,float maxMet);
  TH1F *Rebin(int nGroup);
  TH1F *RebinSimple(int nGroup);  // don't divide by bin width
  TH1F *Rebin(int nBins,float xMin, float xMax,float power=1);
  TH1F *Rebin(int nBins,float *binArray,float power=1);
  TH1F *CloneSimple();
  TH1F *Auto(int n=0);
  void Multiply(pnumber);
  void Multiply(TH1F *scaleHist);

  TH1F *hist_;
  int X_Ndivisions;
  int Y_Ndivisions;
};

void SuperHist::Multiply(pnumber scale)
{
  //      bin = 0;       underflow bin
  //      bin = 1;       first bin with low-edge xlow INCLUDED
  //      bin = nbins;   last bin with upper-edge xup EXCLUDED
  //      bin = nbins+1; overflow bin
  int nbins = hist_->GetNbinsX();
  for(int myBin=0 ; myBin<=nbins+1 ; ++myBin) 
  {
    pnumber binContent(hist_->GetBinContent(myBin), hist_->GetBinError(myBin));
   
    binContent = binContent * scale;
    hist_->SetBinContent(myBin, binContent.x );
    hist_->SetBinError(myBin, binContent.xE);
  }
}

void SuperHist::Multiply(TH1F *scaleHist)
{
  //      bin = 0;       underflow bin
  //      bin = 1;       first bin with low-edge xlow INCLUDED
  //      bin = nbins;   last bin with upper-edge xup EXCLUDED
  //      bin = nbins+1; overflow bin
  int nbins = hist_->GetNbinsX();
  for(int myBin=0 ; myBin<=nbins+1 ; ++myBin) 
  {
    float binCenter = hist_->GetBinCenter(myBin);
    pnumber scale(0,0);

    scale.x  = scaleHist->GetBinContent( scaleHist->FindBin(binCenter) );  
    scale.xE = scaleHist->GetBinError  ( scaleHist->FindBin(binCenter) );  
                                                     
    pnumber binContent(hist_->GetBinContent(myBin), hist_->GetBinError(myBin));

    binContent = binContent * scale;
    
    hist_->SetBinContent(myBin, binContent.x );
    hist_->SetBinError(myBin, binContent.xE);
  }
}

TH1F *SuperHist::Auto(int n)
{
  TH1F *tmp;
  if(n==0)
  {
    float myBins[] = {0,20,70,120,170,500};
    int nBins = sizeof(myBins)/4 -1;
    tmp = this->Rebin(nBins,myBins);
    tmp->GetXaxis()->SetRangeUser(20,500);
  }
  if(n==5)
  {
    tmp = this->RebinSimple(5);
    tmp->GetXaxis()->SetRangeUser(20,320);
    tmp->GetYaxis()->SetTitle("events / 5 GeV");
  }
  if(n==10)
  {
    tmp = this->RebinSimple(10);
    tmp->GetXaxis()->SetRangeUser(20,320);
    tmp->GetYaxis()->SetTitle("events / 10 GeV");
  }
  if(n==13)
  {
    float myBins[] = {0,25,50,75,100,125,150,175,200,225,500};
    int nBins = sizeof(myBins)/4 -1;
    tmp = this->Rebin(nBins,myBins);
//    tmp->GetXaxis()->SetRangeUser(20,500);
  }
  if(n==20)
  {
    tmp = this->RebinSimple(20);
    tmp->GetYaxis()->SetTitle("events / 20 GeV");
  }

  return tmp;
}

TH1F *SuperHist::CloneSimple()
{
   return (TH1F*)hist_->Clone();
}

TH1F *SuperHist::Rebin(int nGroup)
{
  TH1F *tmp_hist = (TH1F*)hist_->Clone();
  tmp_hist->Rebin(nGroup);
  tmp_hist->Scale(1/float(nGroup));
  return tmp_hist;
}

TH1F *SuperHist::RebinSimple(int nGroup)
{
  TH1F *tmp_hist = (TH1F*)hist_->Clone();
  tmp_hist->Rebin(nGroup);
  return tmp_hist;
}

TH1F *SuperHist::Rebin(int nBins, float xMin, float xMax, float power)
{
  float binArray[nBins+1];
  float binStep = (xMax - xMin)/float(nBins);

  for(int myBin=0; myBin<nBins+1; ++myBin)
  {
    binArray[myBin] = xMin + binStep*myBin;
  }

  return Rebin(nBins, binArray, power);
}

TH1F *SuperHist::Rebin(int nBins,float *binArray, float power)
{
  TH1F *new_hist = new TH1F( "" , hist_->GetTitle(), nBins, binArray);

  new_hist->GetXaxis()->SetTitle(hist_->GetXaxis()->GetTitle());
  new_hist->GetYaxis()->SetTitle(hist_->GetYaxis()->GetTitle());
  new_hist->GetXaxis()->SetNdivisions(this->X_Ndivisions);  
  new_hist->GetYaxis()->SetNdivisions(this->Y_Ndivisions);  

  TAxis *newAxis  = new_hist->GetXaxis();

  bool isGood = true;
  if( hist_->GetBinLowEdge(hist_->GetXaxis()->GetFirst()) != new_hist->GetBinLowEdge(new_hist->GetXaxis()->GetFirst())) isGood = false;
  
  int new_histLastBinHighEdge =  new_hist->GetBinLowEdge(new_hist->GetXaxis()->GetLast()) + new_hist->GetBinWidth(new_hist->GetXaxis()->GetLast());
  int hist_LastBinHighEdge =  hist_->GetBinLowEdge(hist_->GetXaxis()->GetLast()) + hist_->GetBinWidth(hist_->GetXaxis()->GetLast());
  if(new_histLastBinHighEdge!=hist_LastBinHighEdge) isGood = false;

  
  float binContent;
  float binError;

  // loop over the bins of the new histogram
  for(int binIt =  newAxis->GetFirst() ; binIt<= newAxis->GetLast() ; ++binIt) 
  {
    float lowEdge = new_hist->GetBinLowEdge(binIt);
    float highEdge = new_hist->GetBinLowEdge(binIt) + new_hist->GetBinWidth(binIt);
    //cout << "new: " << lowEdge << " " << highEdge << endl;
   
    //float origLowEdge  = hist_->GetBinLowEdge(hist_->FindBin(lowEdge));
    //float origHighEdge = hist_->GetBinLowEdge(hist_->FindBin(highEdge));
    //cout << "org: " << origLowEdge << " " << origHighEdge << endl;

    if(isGood)
    {
      binContent      = this->Int(lowEdge,highEdge).x;
      binError = this->Int(lowEdge,highEdge).xE;
   
      pnumber binValue(binContent, binError);

      int numbersOfEffectiveBins = hist_->FindBin(highEdge) - hist_->FindBin(lowEdge);
      pnumber nOEB(numbersOfEffectiveBins,0);     
      binValue = binValue/nOEB;

    //  cout << "before power : " << binValue << endl;      
      binValue = binValue.pow(power); // take power of this
    //  cout << "after power : " << binValue << endl;      

      new_hist->SetBinContent(binIt,binValue.x);
      new_hist->SetBinError(binIt,binValue.xE);
    }
  }

  // set underflow/overflow
  int myBinIt = 0;
  if(hist_->IsBinUnderflow(myBinIt) && new_hist->IsBinUnderflow(myBinIt))
  {
    binContent = hist_->GetBinContent(myBinIt);
    binError = hist_->GetBinError(myBinIt);
    pnumber binValue(binContent, binError);
 
    binValue = binValue.pow(power);
    new_hist->SetBinContent(myBinIt, binValue.x);
    new_hist->SetBinError(myBinIt, binValue.xE);
  } 

  int overflow_bin_hist_ = hist_->GetXaxis()->GetLast() + 1;
  int overflow_bin_new_hist = new_hist->GetXaxis()->GetLast() + 1;
  if(hist_->IsBinUnderflow(myBinIt) && new_hist->IsBinUnderflow(myBinIt))
  {
    binContent = hist_->GetBinContent(overflow_bin_hist_);
    binError = hist_->GetBinError(overflow_bin_hist_);
    pnumber binValue(binContent, binError);

    binValue = binValue.pow(power);
    new_hist->SetBinContent(overflow_bin_new_hist,binValue.x);
    new_hist->SetBinError(overflow_bin_new_hist,binValue.xE);
  } 
  
  if(!isGood) new_hist->Reset();
  return new_hist;
}


pnumber SuperHist::IntUp(float minMet)
{
  //      bin = 0;       underflow bin
  //      bin = 1;       first bin with low-edge xlow INCLUDED
  //      bin = nbins;   last bin with upper-edge xup EXCLUDED
  //      bin = nbins+1; overflow bin
  bool exactBorder = true;

  if(  hist_->GetBinLowEdge( hist_->FindBin(double(minMet)) ) != minMet )exactBorder=false; 

  pnumber val(0,0);
  if(exactBorder)
  {
    val.x = hist_->Integral(hist_->FindBin(double(minMet)),hist_->GetNbinsX()+1);
    Double_t error =0 ; hist_->IntegralAndError(hist_->FindBin(double(minMet)),hist_->GetNbinsX()+1,error);
    val.xE = float(error);
  } else
  {
    cout << setprecision(6) <<  "### Warning :: requested integral [" << minMet<< ", inf) of" << hist_->GetName();
    cout << " can't be performed due to segmentation incombatibility" << endl;
  } 
  return val;
}

pnumber SuperHist::Int(float minMet, float maxMet)
{
  // calulates the integral in [minMet,maxMet)
  // calculate if minMet and maxMet correspond to bin edges
  bool exactBorder = true;

  if(  hist_->GetBinLowEdge( hist_->FindBin(double(minMet)) ) != minMet )exactBorder=false; 
  if(  hist_->GetBinLowEdge( hist_->FindBin(double(maxMet)) ) != maxMet )exactBorder=false; 

  //cout << "calculating integral for [" << minMet << "," << maxMet << "]" << endl;

  float epsilon=0.0001;
  pnumber val(0,0);
  if(exactBorder) 
  {
    Double_t error =0 ; 
    val.x = hist_->Integral(hist_->FindBin(double(minMet)), hist_->FindBin(double(maxMet-epsilon)));
    hist_->IntegralAndError(hist_->FindBin(double(minMet)), hist_->FindBin(double(maxMet-epsilon)),error);
    val.xE = float(error);
       
//    cout << "hist_->FindBin(double(minMet) = " << hist_->FindBin(double(minMet)) << endl;
//    cout << "hist_->FindBin(double(maxMet) = " << hist_->FindBin(double(maxMet-epsilon)) << endl;
//    cout << "integral = " << val << endl;
  }
  else
  {
    cout << setprecision(6) <<  "### Warning :: requested integral [" << minMet<< "," << maxMet << ") of "<< hist_->GetName();
    cout << " can't be performed due to segmentation incombatibility" << endl;
  }

  cout << setprecision(3);
  return val;
}

// --- end: SuperHist

// --- Class SuperCanvas
class SuperCanvas 
{
  public:
  SuperCanvas(string,TH1F*,TH1F*, int mode=0, float x1ndc = 0.55, float x2ndc = 0.70,  float y1ndc = 0.91, float y2ndc =0.93);
  SuperCanvas(string,SuperHist*,SuperHist*, int mode=0, float x1ndc = 0.55, float x2ndc = 0.70,  float y1ndc = 0.91, float y2ndc =0.93);
  TH1F *MakeRatioError(TH1F*);
  void Update();
  TLegend *leg_;
  TCanvas *can_;
  TPad* can_up_;
  TPad* can_dw_;
  TPaveText *title_;
  TString legHeader_;
  TH1F *h1_;
  TH1F *h2_;
  TH1F *h3_;
  TH1F *h4_;
  TH1F *ratio_;
  TH1F *ratioError_;
  TLine *line_;
  string canName_;
  string suffix_;
  int mode_;
};

SuperCanvas::SuperCanvas(string name,SuperHist* aSH, SuperHist*bSH, int mode, float x1ndc, float x2ndc,  float y1ndc, float y2ndc)
{
  SuperCanvas(name, aSH->hist_, bSH->hist_, mode, x1ndc, x2ndc, y1ndc, y2ndc);
}

void SuperCanvas::Update()
{
  if(mode_ == 1) // Update with ratio
  {
    can_up_->cd();
    h1_->Draw("e1 same");
    h2_->Draw("hist same");
    h1_->Draw("same");
    h1_->Draw("axis same");
    leg_->Clear();
    TString leg_entry1 = "Data";
    TString leg_entry2 = "MC";
    if(h1_->GetTitle()!=TString("")) leg_entry1 = h1_->GetTitle();
    if(h2_->GetTitle()!=TString("")) leg_entry2 = h2_->GetTitle();
    leg_->SetHeader(legHeader_);
    leg_->AddEntry(h1_,leg_entry1,"p");
    leg_->AddEntry(h2_,leg_entry2,"f");
    leg_->Draw("same");
    can_dw_->cd();
    ratio_->Draw("e1 same");
    can_->SaveAs(("~/fig/"+canName_+suffix_+".pdf").c_str());
    can_->SaveAs(("~/fig/"+canName_+suffix_+".png").c_str());
    can_->SaveAs(("~/fig/"+canName_+suffix_+".C").c_str());
  }

  if(mode_ ==0 ) // Update without ratio
  {
    can_->cd();
    h1_->Draw("e1");
    h2_->Draw("e2 same");
    h1_->Draw("same"); 
    h1_->Draw("axis same");
    leg_->Clear();
    TString leg_entry1 = "Data";
    TString leg_entry2 = "MC";
    if(h1_->GetTitle()!=TString("")) leg_entry1 = h1_->GetTitle();
    if(h2_->GetTitle()!=TString("")) leg_entry2 = h2_->GetTitle();
    leg_->AddEntry(h1_,leg_entry1,"p");
    leg_->AddEntry(h2_,leg_entry2,"f");
    leg_->SetHeader(legHeader_);
    leg_->Draw("same");
    can_->SaveAs(("~/fig/"+canName_+suffix_+".pdf").c_str());
    can_->SaveAs(("~/fig/"+canName_+suffix_+".png").c_str());
    can_->SaveAs(("~/fig/"+canName_+suffix_+".C").c_str());
  }
}


SuperCanvas::SuperCanvas(string canName, TH1F *h1, TH1F *h2, int mode, float x1ndc, float x2ndc,  float y1ndc, float y2ndc)
{
  h1_ = (TH1F*) h1->Clone();
  h2_ = (TH1F*) h2->Clone();
  canName_ = canName;
  mode_ = mode;

  title_ = new TPaveText(0.20, 0.95, 0.80, 1.0,"blNDC");
  title_->SetFillStyle(0);
  title_->SetFillColor(kWhite);
  title_->SetBorderSize(0);
  title_->AddText("CMS Preliminary 2013, #sqrt{s} = 8 TeV, L=9.2 fb-1"); // README put correct header

  legHeader_ = TString("central rapidity");
  //legHeader_ = TString("forward rapidity");
  suffix_ = "";


  if(mode_ == 1)
  {
    can_ = new TCanvas(canName_.c_str(), canName_.c_str(), 500,650);
    can_->Divide(1, 2);
    can_up_ = (TPad*) can_->GetListOfPrimitives()->FindObject((canName_+"_1").c_str());
    can_dw_ = (TPad*) can_->GetListOfPrimitives()->FindObject((canName_+"_2").c_str());
    can_up_->SetPad(0., 0.35, 1., 1.);
    can_dw_->SetPad(0., 0., 1., 0.35);
    can_up_->SetFrameFillColor(0);
    can_up_->SetFillColor(0);
    can_dw_->SetFillColor(0);
    can_dw_->SetFrameFillColor(0);
    can_dw_->SetTopMargin(0);
    can_up_->SetBottomMargin(0.03);
    can_dw_->SetBottomMargin(0.3);
    can_up_->cd();
    
    h1_->GetXaxis()->SetTitleFont(43);
    h1_->GetXaxis()->SetTitleSize(28);
    h1_->GetXaxis()->SetLabelFont(43);
    h1_->GetXaxis()->SetLabelSize(26);
    h1_->GetXaxis()->SetTitleOffset(3);
    h1_->GetYaxis()->SetTitleFont(43);
    h1_->GetYaxis()->SetTitleSize(28);
    h1_->GetYaxis()->SetLabelFont(43);
    h1_->GetYaxis()->SetLabelSize(26);
    h1_->GetYaxis()->SetTitleOffset(1.4);
    leg_ = new TLegend(x1ndc,x2ndc,y1ndc,y2ndc);
    leg_->SetTextSize(28);
    leg_->SetTextFont(43);
    leg_->SetFillColor(kWhite);
    leg_->SetFillStyle(0);
    leg_->SetBorderSize(0);
    TString leg_entry1 = "Data";
    TString leg_entry2 = "MC";
    if(h1->GetTitle()!=TString("")) leg_entry1 = h1->GetTitle();
    if(h2->GetTitle()!=TString("")) leg_entry1 = h2->GetTitle();
    leg_->AddEntry(h1_,"Data","p");
    leg_->AddEntry(h2_,"MC","f");
    leg_->SetHeader(legHeader_);
    can_dw_->SetGridy();
    can_dw_->SetGridx();
    
  
    can_up_->cd();
    h4_ = (TH1F*)h1_->Clone();
    h4_->GetXaxis()->SetLabelSize(0);
    h4_->Draw("e1");
    h2_->SetFillColor(kWhite);
    h2_->SetFillStyle(0);
    h2_->SetLineWidth(2);
    h2_->Draw("hist same");
    h4_->Draw("same");
    h4_->Draw("axis same");
    leg_->Draw("same");
    title_->Draw("same");
    can_dw_->cd();
    ratio_ = (TH1F*)h1_->Clone();
    ratio_->Divide(h2_);
    ratio_->GetYaxis()->SetNdivisions(504);
    ratio_->GetYaxis()->SetTitle("ratio");
    if(h1->GetTitle()!=TString(""))
    if(h2->GetTitle()!=TString(""))
    ratio_->GetYaxis()->SetTitle(h1->GetTitle() + TString(" / ") + h2->GetTitle());
    ratio_->GetYaxis()->SetTitle("ratio");
    ratio_->GetYaxis()->SetRangeUser(0.3,1.7);
    ratio_->GetYaxis()->CenterTitle();
    ratioError_ = MakeRatioError(h2);
    ratio_->Draw("e1");
    ratioError_->Draw("e2 same");
    float minX = ratio_->GetBinLowEdge(ratio_->GetXaxis()->GetFirst());
    float maxX = ratio_->GetBinLowEdge(ratio_->GetXaxis()->GetLast()) + ratio_->GetBinWidth(ratio_->GetXaxis()->GetLast());
    line_ = new TLine(minX , 1 , maxX , 1 ); 
    line_->SetLineColor(kGray+2);
    line_->SetLineStyle(3);
    line_->Draw("same");
    ratio_->Draw("axis same");
    can_->SaveAs(("~/fig/"+canName_+suffix_+".pdf").c_str());
    can_->SaveAs(("~/fig/"+canName_+suffix_+".png").c_str());
    can_->SaveAs(("~/fig/"+canName_+suffix_+".C").c_str());
  }

  if(mode_ == 0)
  {
    can_ = new TCanvas(canName_.c_str(), canName_.c_str());

    h1_->Draw("e1");
    h2_->Draw("e2 same");
    h1_->Draw("same");
    h1_->Draw("axis same");

    float compresionFactorX = 0.07 ;
    float compresionFactorY = 0.08 ;
    leg_ = new TLegend(x1ndc + compresionFactorX,x2ndc + compresionFactorY,y1ndc,y2ndc);
    leg_->SetTextSize(28);
    leg_->SetTextFont(43);
    leg_->SetFillColor(kWhite);
    leg_->SetFillStyle(0);
    leg_->SetBorderSize(0);
    TString leg_entry1 = "Data";
    TString leg_entry2 = "MC";
    if(h1->GetTitle()!=TString("")) leg_entry1 = h1->GetTitle();
    if(h2->GetTitle()!=TString("")) leg_entry1 = h2->GetTitle();
    leg_->AddEntry(h1_,"Data","p");
    leg_->AddEntry(h2_,"MC","f");
    leg_->SetHeader(legHeader_);
    leg_->Draw("same");
    title_->Draw("same");
    can_->SaveAs(("~/fig/"+canName_+suffix_+".pdf").c_str());
    can_->SaveAs(("~/fig/"+canName_+suffix_+".png").c_str());
    can_->SaveAs(("~/fig/"+canName_+suffix_+".C").c_str());
  }

}

TH1F * SuperCanvas::MakeRatioError(TH1F* ratioError)
{

  //      bin = 0;       underflow bin
  //      bin = 1;       first bin with low-edge xlow INCLUDED
  //      bin = nbins;   last bin with upper-edge xup EXCLUDED
  //      bin = nbins+1; overflow bin
  TH1F *ratioErrorClone = (TH1F*)ratioError->Clone();

  int nbins = ratioErrorClone->GetNbinsX();
  for(int myBin=0 ; myBin<=nbins+1 ; ++myBin) 
  {
    pnumber binContent(ratioErrorClone->GetBinContent(myBin), ratioErrorClone->GetBinError(myBin));
    binContent.x=1;
    //float frac_error = fabs(binContent.xE/binContent.x);
   // binContent.xE= frac_error;
    binContent.xE= 0.05;
    
    ratioErrorClone->SetBinContent(myBin, binContent.x);
    ratioErrorClone->SetBinError(myBin, binContent.xE);
  }

  ratioErrorClone->SetMarkerColor(kGray);
  ratioErrorClone->SetLineColor(kGray);
  ratioErrorClone->SetFillColor(kGray);
  ratioErrorClone->SetMarkerStyle(1);
   
  return ratioErrorClone;
}
// --- end: Class SuperCanvas

// --- Class DataSamples
class DataSample
{
  public:
  DataSample(string fileName_,float CrossSection_, int UserProcessID_, bool goFast_)
  {
    fp = new TFile(fileName_.c_str(),"OPEN");
    events = (TTree*)fp->Get("events");
    maxN = events->GetEntries();
    if(goFast_)maxN = ULong64_t(0.01*maxN);

    xsection = CrossSection_;
    Lumi=9200;
    LumiNorm=Lumi*1.0;
    eventWeight = (xsection*(LumiNorm))/(float(maxN)); // I know we don't have such precission
    UID = UserProcessID_;
    if(CrossSection_==-1) // -1 for DATA
    {
       xsection = -1;
       eventWeight = 1.0;
    }
    events->SetAlias("DR","sqrt(abs(eta1-eta2)**2 + dphi**2)");
    //string myEveW = "PUweight*"+any2string(eventWeight);
    string myEveW = "(PUweight<100 ? PUweight:0)*"+any2string(eventWeight);
    if(CrossSection_==-1) {myEveW = "1";}
    eveW = TCut(myEveW.c_str());
    cout << "opening file: " <<fileName_.c_str() << "with eveW = "<< myEveW << endl;
  }

  TFile *fp;
  TTree *events;
  TH1I *hist;
  float xsection;
  float eventWeight;
  float Lumi;
  float LumiNorm;
  float percentage;
  ULong64_t maxN;
  int  FillColor;
  int  LineStyle;
  int  LineColor;
  int  LineWidth;
  int UID;
  string Title;
  TCut eveW;
};
