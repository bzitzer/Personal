
#include "TMinuit.h"

void fcn(int &npar, double *gin, double &f, double *par, int iflag);

double* NumTot;
double* NumBkg;
double* NumSig;

const int nBins = 100;

void calcFractionCrabTest(string inResultsFile,string inSimFile,const double Elow = 0, const double Ehigh = 3e4)
{
  
  TH1D* hGammaEbin1 = new TH1D("hGammaEbin1","hGammaEbin1",nBins,-1.0,1);
  TH1D* hBackgroundEbin1 = new TH1D("hBackgroundEbin1","hBackgroundEbin1",nBins,-1.0,1);
  TH1D* hDataEbin1 = new TH1D("hDataEbin1","hDataEbin1",nBins,-1.0,1);
  TH1D* hTotal = new TH1D("hTotal","hTotal",nBins,-1.0,1);
  TH1D* hDiff  = new TH1D("hDiff","hDiff",nBins,-1.0,1);
  
  hGammaEbin1->SetStats(0);
  hBackgroundEbin1->SetStats(0);
  hDataEbin1->SetStats(0);

  hDataEbin1->Sumw2(1);
  hGammaEbin1->Sumw2(1);
  hBackgroundEbin1->Sumw2(1);

  fillBDTHistStg6(inResultsFile,hDataEbin1,true,Elow,Ehigh);
  fillBDTHistStg6(inResultsFile,hBackgroundEbin1,false,Elow,Ehigh);
  fillBDTHistStg5(inSimFile,hGammaEbin1,Elow,Ehigh);

  TCanvas* c1 = new TCanvas("c1","c1",60,60,600,500);	   
  hBackgroundEbin1->Draw();

  TCanvas* c2 = new TCanvas("c2","c2",70,70,600,500);		   
  hGammaEbin1->Draw();

  TCanvas* c3 = new TCanvas("c3","c3",50,50,600,500);			   
  hDataEbin1->Draw();

  double nData = hDataEbin1->GetSumOfWeights();
  double nBackground = hBackgroundEbin1->GetSumOfWeights();
  double nGamma = hGammaEbin1->GetSumOfWeights();

  cout << "Total data events: " << nData << " Total Bg events: " << nBackground << " Total sim events: " << nGamma << endl;
  
  hGammaEbin1->SetLineColor(kRed);
  hBackgroundEbin1->SetLineColor(kBlue);
  hDataEbin1->SetLineColor(kBlack);
  
  hBackgroundEbin1->Scale(nData/(nBackground+DBL_EPSILON));
  hGammaEbin1->Scale(nData/(nGamma+DBL_EPSILON));

  NumSig = hGammaEbin1->GetArray();
  NumBkg = hBackgroundEbin1->GetArray();
  NumTot = hDataEbin1->GetArray();

  TMinuit* t = new TMinuit(1);
  t->SetFCN(fcn);
  int ierflg;

  t->mnparm(0,"Fraction ",0.19,1e-3,0.0,1.0,ierflg);

  double arglist[10];
  arglist[0] = 500;
  //arglist[1] = 1.;

  t->mnexcm("MIGRAD",arglist,2,ierflg);
  
  double f = calcMinElecFraction(hDataEbin1,hBackgroundEbin1,hGammaEbin1);

  TGraph* g = new TGraph();
  double chi2;
  double* gin;
  double par[1];
  for(int i=0; i<1000; i++)
    {
      par[0] = (double)i/1000;
      fcn(1,gin,chi2,par,ierflg);
      g->SetPoint(i,par[0],chi2/(nBins-1));
      cout << " fraction: " << par[0] << " chi^2: " << chi2 << endl;
    }

  double f_tmin,f_tmin_err;
  t->GetParameter(0,f_tmin,f_tmin_err);
  cout << "Minimized fraction: " << f << endl;
  cout << "Min. fraction using TMinuit: " << f_tmin << " +/- " << f_tmin_err << endl;
  
  hBackgroundEbin1->Scale(1.0 - f);
  hGammaEbin1->Scale(f);
  hBackgroundEbin1->Draw("same");
  hGammaEbin1->Draw("same");

  double Ng,Nb,Nt,Nd;
  
  for(int i=1; i<=hTotal->GetNbinsX(); i++)
    {
      Ng = hGammaEbin1->GetBinContent(i);
      Nb = hBackgroundEbin1->GetBinContent(i);
      Nt = Nb + Ng;
      hTotal->SetBinContent(i,Nt);

      Nd = hDataEbin1->GetBinContent(i);
      hDiff->SetBinContent(i,(Nd - Nt)/(Nt + DBL_EPSILON));
      hDiff->SetBinError(i,sqrt(Nd + Nd*Nd/(Nt+DBL_EPSILON))/(Nt+DBL_EPSILON));
      //cout << hTotal->GetBinCenter(i) << " " << Nt << endl;
      
      
    }

  hTotal->SetLineColor(kGreen+2);
  hTotal->Draw("same");

  cout << "Min fraction: " << f << endl;
  cout << "Number of Electron-like events: " << f*hDataEbin1->GetSumOfWeights() << endl;
  cout << "Number of Proton-like events: " << (1 - f)*hDataEbin1->GetSumOfWeights() << endl;
  cout << "Total Number of events: " << hDataEbin1->GetSumOfWeights() << endl;
  cout << "Total Simulated events (scaled): " << hTotal->GetSumOfWeights() << endl;

  cout << "KS test results: " << hTotal->KolmogorovTest(hDataEbin1) << endl;
  TCanvas* c4 = new TCanvas("c4","c4",80,80,600,500);	
  //hDiff->Draw("E");
  g->Draw("A*");


}

void fcn(int &npar, double *gin, double &f, double *par, int iflag)
{

  double ratio = par[0];

  f = 0;
  for(int i=0; i<nBins; i++)
    {
      if(NumTot[i] > 0)
	f += pow(NumTot[i] - (1 - ratio)*NumBkg[i] - ratio*NumSig[i],2)/(NumTot[i] + DBL_EPSILON);
    }

  
}

double calcMinElecFraction(TH1D* hData,TH1D* hPr,TH1D* hEl)
{

  double Nd,Np,Ne;
  double a = 0;
  double b = 0;

  for(int i=1; i<=hData->GetNbinsX(); i++)
    {
      Nd = hData->GetBinContent(i);
      Np = hPr->GetBinContent(i);
      Ne = hEl->GetBinContent(i);

      b += pow(Np - Ne,2);
      a += (Np - Ne)*(Nd - Np);
    }

  double f = TMath::Abs(a/(b+DBL_EPSILON));
  return(f);
}


void fillBDTHistStg6(string inFile,TH1D* &h,bool isOnHist = true,const double Elow = 0, const double Ehigh = 3e4)
{

 if(h == NULL)
    {
      cerr << "no hist!" << endl;
      return;
    }
  h->SetDirectory(0);
  
  TFile* f = new TFile(inFile.c_str(),"READ");

  if(!f->IsOpen())
    {
      cerr << "File not open!" << endl;
      return;
    }


  TTree* t = (TTree*)f->Get("EventStatsTree");
  if(t == NULL)
    {
      cerr << "Events Tree missing for: " << inFile << endl;
      return;
    }

  double bdt;
  float energy;
  bool isOn,isOff;
  
  t->SetBranchAddress("BDTScore",&bdt);
  t->SetBranchAddress("EnergyGeV",&energy);
  t->SetBranchAddress("OnEvent",&isOn);
  t->SetBranchAddress("OffEvent",&isOff);

  for(int i=0; i<t->GetEntries(); i++)
    {
      t->GetEntry(i);
      if(i%100000 == 0)
	cout << "On File: " << inFile << " On Event: " << i << " of: " << t->GetEntries() << endl;

      if( energy < Elow || Ehigh < energy ){ continue; } 

      if(isOnHist && isOn)
	h->Fill(bdt);

      if(!isOnHist && isOff)
	h->Fill(bdt);
    }
}

void fillBDTHistStg5(string inFile,TH1D* &h,const double Elow = 0,const double Ehigh = 3e4)
{
  if(h == NULL)
    {
      cerr << "no hist!" << endl;
      return;
    }
  h->SetDirectory(0);
  
  TFile* f = new TFile(inFile.c_str(),"READ");

  if(!f->IsOpen())
    {
      cerr << "File not open!" << endl;
      return;
    }


  TTree* t = (TTree*)f->Get("SelectedEvents/CombinedEventsTree");
  if(t == NULL)
    {
      cerr << "Combined Tree missing for: " << inFile << endl;
      return;
    }
  VAShowerData* S = new VAShowerData();
  double bdt,xdeg,ydeg;

  t->SetBranchAddress("BDTScore",&bdt);
  t->SetBranchAddress("S",&S);
  t->SetBranchAddress("XDeg",&xdeg);
  t->SetBranchAddress("YDeg",&ydeg);

  
  for(int i=0; i<t->GetEntries(); i++)
    {
      t->GetEntry(i);
      if(i%100000 == 0)
	cout << "On File: " << inFile << " On Event: " << i << " of: " << t->GetEntries() << endl;

      //if( pow(xdeg,2) + pow(ydeg,2) > 1.0 ){ continue; }
      if( S->fTheta2_Deg2 > 0.03 ){ continue; }
      
      if( S->fEnergy_GeV < Elow || Ehigh < S->fEnergy_GeV ){ continue; } 

      h->Fill(bdt);
    }

  f->Close();
}
