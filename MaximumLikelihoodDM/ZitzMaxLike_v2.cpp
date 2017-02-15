
#include <vector>
#include <sstream>
#include <iostream>
#include <fstream>

#include <TH1D.h>
#include <TH2D.h>
#include <TGraph.h>
#include <TGraph2D.h>
#include <TGraphAsymmErrors.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TAxis.h>
#include <TF1.h>
#include <TF2.h>
#include <TLegend.h>
#include <TLine.h>
#include <TFile.h>
#include <TTree.h>
#include <TROOT.h>
#include "VASlalib.h"
#include "VASlalib_abr.cpp"
//#include "ZitzEventWeighting.cpp"
#include "ZitzMaxLike_v2.h"
using namespace std;


void ZitzMaxLike_v2(string inEventFile,string inEAFile,string dispPath)
{
  //  gROOT->ProcessLine(".L ZitzEventWeighting.cpp+");
  double expTot;
  double E;
  bool isOn;
  int runNum;
  int runNum_old = -1;
  double NonTot,NoffTot;
  double NbgTot;
  double w;
  double w_avg;
  double dE;
  TH1D* hBg = new TH1D("hBg","hBg",1000,1.5,5); // units of log10 GeV
  TH1D* hNon = new TH1D("hNon","hNon",1000,1.5,5); // units of log10 GeV

  TH2D* hDisp[100];
  int j = 0;
  double val;
  cout << "Reading event file.. ";
  TTree* tEv = (TTree*)readEventFile(inEventFile.c_str(),crab);
  if( tEv == NULL )
    {
      cout << "problem with events tree!" << endl;
      return;
    }
  cout << "Done!" << endl;
  cout << "Getting effective area... ";
  TGraphAsymmErrors* gEA = getEffectiveArea(inEAFile,0,expTot);
  if( gEA == NULL )
    {
      cout << "problem with effective areas!" << endl;
      return;
    }
  cout << "Done!" << endl;
  cout << "Total Exposure: " << expTot/3600 << endl;
  TCanvas* c1 = new TCanvas();
  //c1->SetLogx(1);
  gEA->Draw("A*");
  tEv->SetBranchAddress("isOn",&isOn);
  tEv->SetBranchAddress("energy",&E);
  tEv->SetBranchAddress("runNum",&runNum);
  tEv->SetBranchAddress("w",&w);
  for(int i=0; i<tEv->GetEntries(); i++)
    {
      tEv->GetEntry(i);
      if( runNum_old != runNum )
	{
	  hDisp[j] = getEnergyDispersion(dispPath,runNum);
	  j++;
	}
      if(!isOn)
	{
	  hBg->Fill(TMath::Log10(E));
	  w_avg += w;
	  NoffTot++;
	}
      else
	{
	  hNon->Fill(TMath::Log10(E));
	  NonTot++;
	}
      runNum_old = runNum;
    }
  w_avg = w_avg/NoffTot;
  TH2D* hDispAvg = renormalizeDisp(hDisp,j);
  double norm = 0;
  hBg->Smooth(4);
  hNon->Smooth(4);
  for(int i=0; i<=hBg->GetNbinsX(); i++)
    {
      dE = TMath::Power(10,hBg->GetBinWidth(i) + hBg->GetBinLowEdge(i)) - TMath::Power(10,hBg->GetBinLowEdge(i));
      norm += hBg->GetBinContent(i)*dE;
    }
  for(int i=0; i<=hBg->GetNbinsX(); i++)
    {
      val = hBg->GetBinContent(i);
      hBg->SetBinContent(i,val/norm);
    }
   
  vector<double> p_off(tEv->GetEntries(),0);
  vector<double> p_on(tEv->GetEntries(),0);
  //  TCanvas* c1 = new TCanvas();
  //hDispAvg->Draw("colz");
  
  TH1D* hBgConv = convolveEnergyDisp(hBg,hDispAvg);
  /*
  norm = hBg->GetSumOfWeights();
  hBg->Scale(1.0/norm);
  norm = hBgConv->GetSumOfWeights();
  hBgConv->Scale(1.0/norm);
  norm = hNon->GetSumOfWeights();
  hNon->Scale(1.0/norm);
  */
  TCanvas* c2 = new TCanvas("c2","c2",50,50,600,500);
  hBg->Draw();
  hNon->SetLineColor(kRed);
  hNon->Draw("same");

  double norm_check = 0;
  for(int i=0; i<=hBg->GetNbinsX(); i++)
    norm_check += hBg->GetBinContent(i)*hBg->GetBinWidth(i);

  cout << "normalization check: " << norm_check << endl;
  norm_check = 0;
  for(int i=0; i<=hBgConv->GetNbinsX(); i++)
    norm_check += hBgConv->GetBinContent(i)*hBgConv->GetBinWidth(i);

  cout << "second norm. check: " << norm_check << endl;

  //return;
  vector<double> F0;
  double F0_min = -9;
  double F0_max = -5;
  const int nStepsF0 = 20;
  //double dF0 = (TMath::Log10(F0_max) - TMath::Log10(F0_min))/nStepsF0;
  double dF0 = (F0_max - F0_min)/nStepsF0;
  for(int i=0; i<=nStepsF0; i++)
    F0.push_back(F0_min + i*dF0);
  //F0.push_back(F0_min*TMath::Power(10,i*dF0));

  vector<double> I;
  double I_min = 1.5;
  double I_max = 3.5;
  const int nStepsI = 20;
  double dI = (I_max - I_min)/nStepsI;
  for(int i=0; i<=nStepsI; i++)
    I.push_back(I_min + i*dI);

  //vector<double> p_off,p_on;
  double LL_null = 0;
  double LL_tmp;
  double f_i;
  double g_norm;
  vector<double> g_i(tEv->GetEntries(),0);
  // ----------------------
  //  NULL hypothesis loop:
  // ----------------------
  double b_null = (NonTot + NoffTot)/(1+w_avg);
  LL_null += NonTot*TMath::Log(w_avg*b_null) - w_avg*b_null;
  for(int i=0; i<tEv->GetEntries(); i++)
    {
      tEv->GetEntry(i);
      //p_off.at(i) = hBgConv->Interpolate(TMath::Log10(E));
      p_off.at(i) = hBg->Interpolate(TMath::Log10(E));
      if(p_off.at(i) == 0 )
	{
	  cout << "Warning p = 0 ! E = " << E << endl;
	  continue;
	}
      if( E < 150.0 ){ continue; }
      LL_null += TMath::Log(p_off.at(i));
    }
  cout << "LL null: " << LL_null << endl;
  //return;
  // -----------
  //  Main loop:
  // -----------
  TH1D* hPL  = new TH1D();
  TH1D* hPL0 = new TH1D();
  TH2D* hLL  = new TH2D("hLL","Fit Results",nStepsI,I_min,I_max,nStepsF0,F0_min,F0_max);
  TH2D* hg   = new TH2D("hg","number expected signal events",nStepsI,I_min,I_max,nStepsF0,F0_min,F0_max);
  TH2D* hb_min = new TH2D("hb_min","minimized b",nStepsI,I_min,I_max,nStepsF0,F0_min,F0_max);
  TGraph* gb_min = new TGraph();
  double b;
  for(int n =0; n<=nStepsI; n++) // loop over index
    {
      hPL0 = calcPLOnDist(gEA,expTot,F0_min,I.at(n));
      LL_tmp = 0;
      
      for(int m = 0; m<=nStepsF0; m++) // loop over flux amplitude
	{
	  gb_min->Delete();
	  TGraph* gb_min = new TGraph();
	  //hPL = hPL0;
	  //hPL->Scale(F0.at(m)/F0_min);
	  hPL->Delete();
	  TH1D* hPL = calcPLOnDist(gEA,expTot,TMath::Power(10,F0.at(m)),I.at(n));
	  g_norm = 0;
	  for(int i=0; i<hPL->GetNbinsX(); i++)
	    {
	      dE = TMath::Power(10,hPL->GetBinWidth(i) + hPL->GetBinLowEdge(i)) - TMath::Power(10,hPL->GetBinLowEdge(i));
	      g_norm += hPL->GetBinContent(i)*dE;
	    }
	  // Poisson part:
	  
	  for(int i=0; i<tEv->GetEntries(); i++)
	    {
	      tEv->GetEntry(i);
	      if(isOn)
		{
		  g_i.at(i) = hPL->Interpolate(TMath::Log10(E));
		  /*
		  if(g_i.at(i) == 0.0)
		    cout << "g = 0 E = " << E << endl;
		    */
		}
	      else
		g_i.at(i) = 0.0;
	    }
	  
	  b = min_b(NoffTot,NonTot,g_norm,w_avg,p_off,g_i,gb_min);
		  //b = NoffTot*w_avg;
	  LL_tmp += NonTot*TMath::Log(g_norm + b*w_avg) - g_norm - b*(1.0 + w_avg);
	  cout << "g_norm: " << g_norm << " b_norm: " << b << " " << w_avg << endl;
	  for(int i=0; i<tEv->GetEntries(); i++)
	    {
	      tEv->GetEntry(i);
	      if(p_off.at(i) == 0){ continue; }
	      if( E < 150.0 ){ continue; }
	      if(isOn)
		{
		  p_on.at(i) = b*w_avg*p_off.at(i) + g_i.at(i);
		  LL_tmp += TMath::Log(p_on.at(i));
		}
	      else
		LL_tmp += TMath::Log(p_off.at(i));
	    }
	  cout << "Index: " << I.at(n) << " F0: " << F0.at(m) << " LL: " << LL_tmp << endl;
	  hLL->SetBinContent(n,m,-2.0*(LL_tmp - LL_null));
	  hb_min->SetBinContent(n,m,TMath::Log10(b));
	  hg->SetBinContent(n,m,TMath::Log10(g_norm));
	}
 
    }
  cout << "Total ON counts: " << NonTot << " Total OFF counts: " << NoffTot << " alpha: " << w_avg << endl;
  hLL->SetStats(0);
  TCanvas* c3 = new TCanvas();
  //c3->SetLogy(1);
  hLL->Draw("colz");
  hg->SetStats(0);
  TCanvas* c4 = new TCanvas();
  //c4->SetLogy(1);
  hg->Draw("colz");
  hb_min->SetStats(0);
  TCanvas* c5 = new TCanvas();
  //c5->SetLogy(1);
  hb_min->Draw("colz");
  //gb_min->Draw("A*");

}

TH1D* calcPLOnDist(TGraphAsymmErrors* g,double exp,double k,double index)
{
  TH1D* h = new TH1D("hOn","hOn",1000,1.5,5); // units of log 10 GeV
  double E,Elog10TeV;
  double F;
  for(int i=0; i<h->GetNbinsX(); i++)
    {
      E = TMath::Power(10.0,h->GetBinCenter(i));
      Elog10TeV = TMath::Log10(E/1000.0);
      F = exp*k*g->Eval(Elog10TeV)*TMath::Power(E/1000.0,-1.0*index);
      h->SetBinContent(i,F);
    }
  return(h);

} 

TH1D* convolveEnergyDisp(TH1D* h,TH2D* hDisp)
{
  //hDisp->Smooth(1);

  double norm = 0;
  double D;
  double val;
  TH1D* hOut =(TH1D*)h->Clone("hBgConv");
  const int nBins = hOut->GetNbinsX();
  const int nDispBinsX = hDisp->GetNbinsX();
  const int nDispBinsY = hDisp->GetNbinsY();
  
  double Erec,E,dE;
  for(int i=1; i<=h->GetNbinsX(); i++)
    {
      hOut->SetBinContent(i,0.0);
      Erec = h->GetBinCenter(i); // GeV
      val = 0;
      D = 0;
      norm = 0;
      cout << "Rec. Energy: " << Erec << endl;
      for(int j=0; j<h->GetNbinsX(); j++)
	{
	  dE = hDisp->GetXaxis()->GetBinWidth(j);
	  E = h->GetBinCenter(j); // GeV
	  //cout << "Rec. Energy: " << Erec << " MC energy: " << E << endl;
	  if( E < 0.0 ){ continue; }
	  D = h->GetBinContent(j)*hDisp->Interpolate(E,Erec)*dE;
	  norm += hDisp->Interpolate(E,Erec)*dE;
	  if( E    < hDisp->GetXaxis()->GetBinCenter(0) || hDisp->GetXaxis()->GetBinCenter(nDispBinsX) < E)
	    cout << "Rec. Energy: " << Erec << " MC energy: " << E << endl;
	  if( Erec < hDisp->GetYaxis()->GetBinCenter(0) || hDisp->GetYaxis()->GetBinCenter(nDispBinsY) < Erec)
	    cout << "Rec. Energy: " << Erec << " MC energy: " << E << endl;
	  
	  val += D;
	  
	}
      
      //cout << "normalization for energy disp: " << norm << endl;
      if(norm != 0)
	hOut->SetBinContent(i,val/norm);
      
    }
  norm = 0;
  for(int i=0; i<=hOut->GetNbinsX(); i++)
    {
      val = hOut->GetBinContent(i)*hOut->GetBinWidth(i);
      norm += val; 
    }
  //cout << "Energy Conv. normalization: " << norm << endl;
  hOut->Scale(1.0/norm);
  cout << "Convolved normalization check: " << hOut->Integral(0,nBins,"width") << endl;
  return(hOut);
}

TH1D* convolveEDisp_varBinWidth(TH1D* h,TH2D* hDisp)
{
  //hDisp->Smooth(1);

  double norm = 0;
  double D;
  double val;
  TH1D* hOut =(TH1D*)h->Clone("hBgConv");
  const int nBins = hOut->GetNbinsX();
  const int nDispBinsX = hDisp->GetNbinsX();
  const int nDispBinsY = hDisp->GetNbinsY();
  
  double Erec,E,dE;
  for(int i=1; i<=h->GetNbinsX(); i++)
    {
      hOut->SetBinContent(i,0.0);
      Erec = h->GetBinCenter(i); // GeV
      val = 0;
      D = 0;
      norm = 0;
      //cout << "Rec. Energy: " << Erec << endl;
      for(int j=1; j<=h->GetNbinsX(); j++)
	{
	  dE = hDisp->GetXaxis()->GetBinWidth(j);
	  E = h->GetBinCenter(j); // GeV
	  //cout << "Rec. Energy: " << Erec << " MC energy: " << E << endl;
	  if( E < 0.0 ){ continue; }
	  D = h->GetBinContent(j)*h->GetBinWidth(j)*hDisp->Interpolate(E,Erec)*dE;
	  norm += hDisp->Interpolate(E,Erec)*dE;
	  if( E    < hDisp->GetXaxis()->GetBinCenter(0) || hDisp->GetXaxis()->GetBinCenter(nDispBinsX) < E)
	    cout << "Rec. Energy: " << Erec << " MC energy: " << E << endl;
	  if( Erec < hDisp->GetYaxis()->GetBinCenter(0) || hDisp->GetYaxis()->GetBinCenter(nDispBinsY) < Erec)
	    cout << "Rec. Energy: " << Erec << " MC energy: " << E << endl;
	  
	  val += D;
	  
	}
      
      //cout << "normalization for energy disp: " << norm << endl;
      if(norm != 0)
	hOut->SetBinContent(i,val/norm);
      
    }
  norm = 0;
  for(int i=0; i<=hOut->GetNbinsX(); i++)
    {
      dE  = hOut->GetBinWidth(i);
      val = hOut->GetBinContent(i);
      hOut->SetBinContent(i,val/dE); 
    }
  //cout << "Energy Conv. normalization: " << norm << endl;
  //hOut->Scale(1.0/norm);
  //cout << "Convolved normalization check: " << hOut->Integral(0,nBins,"width") << endl;
  return(hOut);
}

TH2D* renormalizeDisp(TH2D* h[],int numFilled)
{
  TH2D* hNorm = (TH2D*)h[0]->Clone("hEDispAvg");
  for(int i=0; i<=hNorm->GetNbinsX(); i++)
    for(int j=0; j<=hNorm->GetNbinsY(); j++)
      hNorm->SetBinContent(i,j,0.0);

  double Erec,Etr,val,norm,dE;
  // adding all histograms:
  for(int k=0; k<numFilled; k++)
    {
      if(h[k] == NULL){ continue; }
      cout << "filling dispersion histogram # " << k << endl;
      for(int i=0; i<=hNorm->GetNbinsX(); i++)
	for(int j=0; j<=hNorm->GetNbinsY(); j++)
	  {
	    val = h[k]->GetBinContent(i,j) + hNorm->GetBinContent(i,j);
	    hNorm->SetBinContent(i,j,val);
	    //	    h[k]->Delete();
	  }
    }
  // renomralize:
  for(int j=0; j<=hNorm->GetNbinsY(); j++)
    {
      // cout << " Normalizing row: " << j << endl;
      norm = 0;
      for(int i=0; i<=hNorm->GetNbinsX(); i++)
	norm += hNorm->GetBinContent(i,j);
      for(int i=0; i<=hNorm->GetNbinsX(); i++)
	{
	  dE = hNorm->GetXaxis()->GetBinWidth(i);
	  val = hNorm->GetBinContent(i,j);
	  if(norm != 0)
	    hNorm->SetBinContent(i,j,val/(dE*norm));
	}
    }
  for(int j=0; j<=hNorm->GetNbinsY(); j++)
    {
      norm = 0;
      for(int i=0; i<=hNorm->GetNbinsX(); i++)
	norm += hNorm->GetBinContent(i,j)*hNorm->GetXaxis()->GetBinWidth(i); 

      cout << "Dispersion norm check for Erec = " << hNorm->GetYaxis()->GetBinCenter(j) << " norm = " << norm << endl;
    }

  return(hNorm);
}


TH2D* renormalizeDispWeighted(TH2D* h[],vector<double> w)
{
  int numFilled = w.size();
  TH2D* hNorm = (TH2D*)h[0]->Clone("hEDispAvg");
  for(int i=0; i<=hNorm->GetNbinsX(); i++)
    for(int j=0; j<=hNorm->GetNbinsY(); j++)
      hNorm->SetBinContent(i,j,0.0);

  double Erec,Etr,val,norm,dE;
  // adding all histograms:
  for(int k=0; k<numFilled; k++)
    {
      if(h[k] == NULL){ continue; }
      cout << "filling dispersion histogram # " << k << endl;
      for(int i=0; i<=hNorm->GetNbinsX(); i++)
	for(int j=0; j<=hNorm->GetNbinsY(); j++)
	  {
	    val = h[k]->GetBinContent(i,j)*w.at(k) + hNorm->GetBinContent(i,j);
	    hNorm->SetBinContent(i,j,val);
	    //	    h[k]->Delete();
	  }
    }

  // smooth summed histogram
  //hNorm->Smooth(1);
  // renomralize:
  
  for(int j=0; j<=hNorm->GetNbinsY(); j++)
    {
      // cout << " Normalizing row: " << j << endl;
      norm = 0;
      for(int i=0; i<=hNorm->GetNbinsX(); i++)
	norm += hNorm->GetBinContent(i,j);
      for(int i=0; i<=hNorm->GetNbinsX(); i++)
	{
	  dE = hNorm->GetXaxis()->GetBinWidth(i);
	  val = hNorm->GetBinContent(i,j);
	  if(norm != 0)
	    hNorm->SetBinContent(i,j,val/(dE*norm));
	}
    }

  // check:

  for(int j=0; j<=hNorm->GetNbinsY(); j++)
    {
      norm = 0;
      for(int i=0; i<=hNorm->GetNbinsX(); i++)
	norm += hNorm->GetBinContent(i,j)*hNorm->GetXaxis()->GetBinWidth(i); 

      cout << "Dispersion norm check for Erec = " << hNorm->GetYaxis()->GetBinCenter(j) << " norm = " << norm << endl;
    }

  return(hNorm);
}


TH2D* getEnergyDispersion(string path,int runNum)
{
  ostringstream os;
  ostringstream os1;
  ostringstream os2;
  ostringstream os3;
  
  os  << path.c_str() << runNum << ".mat.txt";
  os1 << path.c_str() << runNum << ".ErecBin.txt";
  os2 << path.c_str() << runNum << ".EtrBin.txt";
  os3 << "eDisp" << runNum;
  //cout << "looking in path: " << path << endl;
  int nBinsMC  = 101;
  int nBinsRec = 101;
  double* eRecBin = readBinFile2(os1.str(),nBinsRec);
  double* eMCBin  = readBinFile2(os2.str(),nBinsMC);
  TH2D* hDisp = new TH2D(os3.str().c_str(),os3.str().c_str(),nBinsRec-1,eRecBin,nBinsMC-1,eMCBin);
  ifstream in(os.str().c_str());
  if(!in.is_open())
    {
      cout << "" << endl;
      cout << os.str().c_str() << " not found!" << endl;
      return(NULL);
    }
  char str[300];
  in.getline(str,300);
  in.getline(str,300);
  int i = 1; // Rec. Energies
  int j = 1; // MC Energies
  double val;
  while( in >> val )
    {
      //if( val != 0.0 ){ cout << val << endl; }
      if(i >= nBinsRec)
	{
	  i = 1;
	  j++;
	}
      i++;
      hDisp->SetBinContent(j,i,val);
    }
  hDisp->SetStats(0);
  hDisp->GetXaxis()->SetTitle("Reconstructed Energy [GeV]");
  hDisp->GetYaxis()->SetTitle("True Energy [GeV]");
  
  return(hDisp);
}

double* readBinFile2(string inFile,int &N,int nHeaderLines)
{
  ifstream in(inFile.c_str());
  vector<double> x2;
  double tmp;
  if(!in.is_open())
    {
      cout << "Missing file! All is woe!" << endl;
      return(0);
    }
  
  char str[300];
  for(int i=0; i<nHeaderLines; i++)
    in.getline(str,300);
  //in.getline(str,300);
  while( in >> tmp )
    x2.push_back(tmp);
    
  N = x2.size();
  //cout << "Number of bins: " << N << endl;
  double* x = new double[N];
  for(int i=0; i<x2.size(); i++)
    {
      x[i] = x2.at(i);
      //cout << x[i] << endl;
    }
  return(x);
}
TGraphAsymmErrors* getEffectiveArea(string inFile,int runNum,double &exp)
{
  //  cout << "Got here" << endl;
  ifstream ifs(inFile.c_str());
  if(!ifs.is_open())
    {
      cout << "Warning! No effective area file!" << endl;
      return(NULL);
    }

  TGraphAsymmErrors* g = new TGraphAsymmErrors();
  double LT,E,EA,dEAplus,dEAminus;
  int i = 0;
  int runNum_tmp;
  char tmpString[250];
  ifs.getline(tmpString,250); 
  //cout << tmpString << endl;
  ifs.getline(tmpString,250);
  //cout << tmpString << endl;
  while( ifs >> runNum_tmp >> LT >> E >> EA >> dEAminus >> dEAplus )
    {
      //      cout << LT << endl;
      
      if(runNum_tmp == runNum)
	{
	  //	  if( i == 0)
	  //  cout << "Found EA file for run " << runNum_tmp << endl;
	  if( EA < 0.0 )
	    EA = 0.0;
	  g->SetPoint(i,TMath::Log10(E/1000.0),EA);
	  g->SetPointEYlow(i,-1.0*dEAminus);
	  g->SetPointEYhigh(i,dEAplus);
	  i++;
	  exp = LT;
	}
    }
  //cout << "number of points: " << i << endl;
  return(g);

}

TTree* readEventFile(string inFile,dwarf_t dwarf)
{
  ifstream in(inFile.c_str());
  if(!in.is_open())
    {
      cout << "No event file!" << endl;
      return(NULL);

    }
  TTree* t = new TTree("eventTree","eventTree");
  bool isOn;
  int runNum;
  double runLT,time,ra,dec,ra_tr,dec_tr,E,w,el,az,no,offset;
  double theta;
  double ra_dw,dec_dw;
  if( dwarf == segue_1 )
    {
      ra_dw = 15.0*(10.0 + 7.0/60 + 4.0/3600)*TMath::DegToRad();
      dec_dw = (16.0 + 4.0/60 + 55.0/3600)*TMath::DegToRad();
    }
  char str[300];
  in.getline(str,300);
  in.getline(str,300);

  t->Branch("runNum",&runNum,"runNum/I");
  t->Branch("runLT",&runLT,"runLT/D");
  t->Branch("time",&time,"time/D");
  t->Branch("ra",&ra,"ra/D");
  t->Branch("dec",&dec,"dec/D");
  t->Branch("isOn",&isOn,"isOn/O");
  t->Branch("w",&w,"w/D");
  t->Branch("energy",&E,"energy/D");
  t->Branch("theta",&theta,"theta/D");
  int i = 0;
  while(in >> runNum >> runLT >> time >> 
	ra >> dec >> ra_tr >> dec_tr >> 
	E >> isOn >> w >> 
	el >> az >> no >> offset)
    {
      if(i%10000 == 0){ cout << "Filling from event list - On event: " << i << endl; }
      theta = slaDsep(ra*TMath::DegToRad(),dec*TMath::DegToRad(),ra_dw,dec_dw);
      theta *= TMath::RadToDeg();
      t->Fill();
      i++;
    }
  return(t);

}


double LLBg(double b,double g,double Non,double Noff,double alpha,vector<double> h,vector<double> p)
{ 
  //double LL = 0;
  double LL = -b*(1.0 + alpha) + Noff;
  
  for(int i=0; i<h.size(); i++)
    if(p.at(i) != 0 || h.at(i) != 0 )
      LL += b*alpha*h.at(i)/(alpha*b*h.at(i) + p.at(i)); 
  
  return(LL);

}

double min_b(double Noff,double Non,double g,double alpha,vector<double> h_E,vector<double> g_E,TGraph* &g1)
{
  double x1 = 0.01;
  double x2 = Noff*10; 
  //TGraph* g1 = new TGraph();
  double dx = x2;
  int s1,s2,s3;
  double y,y1,y2;
  int i=0;
  double x3,y3;
  double tmp;
  double dy = 10000;
  //cout << "Size of h: " << h_E.size() << " since of g: " << g_E.size() << endl;
  while( dy > 0.05 )
    {
      
      y1 = LLBg(x1,g,Non,Noff,alpha,h_E,g_E);
      y2 = LLBg(x2,g,Non,Noff,alpha,h_E,g_E);
      if(i == 0)
	{
	  g1->SetPoint(i,x1,y1);
	  i++;
	  g1->SetPoint(i,x2,y2);
	  i++;	
	}
      x3 = (x1 + x2)/2.0;
      y3 = LLBg(x3,g,Non,Noff,alpha,h_E,g_E);
      s1 = y1/TMath::Abs(y1);
      s2 = y2/TMath::Abs(y2);
      s3 = y3/TMath::Abs(y3);
      g1->SetPoint(i,x3,y3);
      i++;
      if(s1 == s3)
	{
	  x1 = x3;
	}
      if(s2 == s3)
	{
	  x2 = x3;
	}
      
      dx = x2 - x1;
      if(i > 2000 )
	{
	  cout << "warning! No minimium b found... using b = " << x3 << endl; 
	  cout << "low endpoint: " << x1 << " high endpoint: " << x2 << endl;
	  cout << "f(x1) = " << y1 << " f(x2) = " << y2 << endl;
	  cout << "Non: " << Non << " Noff: " << Noff << " alpha: " << alpha << endl;
	  break;
	}
      
      if(TMath::IsNaN(y1) || TMath::IsNaN(y2))
	{
	  cout << "Warning dlogL/db = NaN! Non = " << Non << " Noff: " << Noff << " alpha: " << alpha <<endl;
	  cout << "x1 = " << x1 << " x2 = " << x2 << " g = " << g << endl; 
	  /*
	  for(int i=0; i<h_E.size(); i++)
	    if(h_E.at(i) == 0 || g_E.at(i) == 0)
	      cout << "Warning! h(E) = " << h_E.at(i) << " g(E) = " << g_E.at(i) << endl;
	  */
	}
      
      dy = TMath::Abs(y3);
       
    }
  //cout << "Minimized b: " << x3 << " at y = " << y3 << endl;
  if( x3 <= 0.0 )
    cout << "Warning! b_min <= 0" << endl;
  return(x3);

}
