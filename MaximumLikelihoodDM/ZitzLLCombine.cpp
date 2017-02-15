
#include <vector>
#include <sstream>
#include <iostream>
#include <fstream>

#include <TH1D.h>
#include <TH2D.h>
#include <TGraph.h>
#include <TGraph2D.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TMultiGraph.h>
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
#include <TRandom3.h>

TGraph* findMinAndSubtract(TGraph* g0,double y_null);

double findLimit(TGraph* g,const double cl = 2.71);

void plotAllLimits();

TGraphAsymmErrors* plotLimitsAlexPass5f(string inFile);

TGraph* ZitzLLCombine(string inFileList,const double cl = 2.71)
{
  string inFile;
  ifstream in(inFileList.c_str());
  vector<string> fileVector;
  string str;
  while(in >> str)
    fileVector.push_back(str);

  const int nDSphs = fileVector.size();
  if(nDSphs == 0)
    {
      cout << "Warning! no files in list!" << endl;
      return(NULL);
    }
  
  TGraph* gLL[nDSphs];
  
  TGraph* gTmp = new TGraph();
  TFile* f = new TFile(fileVector.at(0).c_str(),"READ");
  TFile* fIn[nDSphs];
  gTmp = (TGraph*)f->Get("gLLnull");
  const int nMass = gTmp->GetN();
  const double *m = gTmp->GetX();
  TGraph* gLL_combined[nMass];
  TGraph* gLL_2deltaLL[nMass];
  TGraph* gLL_null[nDSphs];
  TGraph* gLL_null_combined = new TGraph();
  TGraph* gLL_min[nDSphs];
  TGraph* gLL_limit = new TGraph();

  ostringstream os[nMass];
  ostringstream os1[nMass];
  double tmpLL,signu,limit_val,signu_limit,tmpLLnull;

  for(int i=0; i<nMass; i++) // loop over every mass value
    {
      cout << "Combining LL results for M = " << m[i] << endl;
      gLL_combined[i] = new TGraph();
      //gLL_min[i] = new TGraph();
      os[i]  << "gLL_signu_m" << (int)m[i];
      os1[i] << "Mass = " << (int)m[i];

      limit_val = 100;
      tmpLLnull = 0;
      for(int j=0; j<nDSphs; j++) // loop over files
	{
	  fIn[j] = new TFile(fileVector.at(j).c_str(),"READ");
	  cout << " opening: " << fileVector.at(j) << endl;
	 
	  if(!fIn[j]->IsOpen())
	    {
	      cout << "Warning! file " << fileVector.at(j) << " not found!" << endl;
	      return(NULL);
	    }

	  gLL[j]      = (TGraph*)fIn[j]->Get(os[i].str().c_str());
	  gLL_min[j]  = (TGraph*)fIn[j]->Get("gLLmin");
	  gLL_null[j] = (TGraph*)fIn[j]->Get("gLLnull");

	  if(gLL[j] == NULL)
	    {
	      cout << "Warning! can't find : " << os[i].str() << " in file: " << fileVector.at(j) << endl;
	      return(NULL);
	    }
	  if( gLL_min[j] == NULL )
	    {
	      cout << "Warning! Can't find minimum TGraph in file: " << fileVector.at(j) << endl;
	      return(NULL);
	    }
	  fIn[j]->Close();
	  tmpLLnull += gLL_null[j]->GetY()[i];

	}

      gLL_null_combined->SetPoint(i,m[i],tmpLLnull);

      cout << " NULL value for M = " << m[i] << " " << gLL_null_combined->GetY()[i] << endl;
      // adding all the likelihoods:
      for(int k=0; k<gLL[0]->GetN(); k++)
	{

	  tmpLL = 0;
	  signu = gLL[0]->GetX()[k];

	  for(int j=0; j<nDSphs; j++)
	    {
	      tmpLL += 0.5*( gLL[j]->GetY()[k] ) + gLL_min[j]->GetY()[i] ;
	      //tmpLL += 0.5*gLL[j]->GetY()[k];
	      //cout << fileVector.at(j) << " Mass: " << m[i] << " sig_nu " << signu << " LL: " << gLL[j]->GetY()[k] << endl;
	      //cout << fileVector.at(j) << " " << gLL_min[j]->GetY()[i] << endl;
	    }
	
	  gLL_combined[i]->SetPoint(k,signu,tmpLL);
	
	}
      
      gLL_2deltaLL[i] = findMinAndSubtract(gLL_combined[i],tmpLLnull);
      signu_limit = findLimit(gLL_2deltaLL[i],cl);
      gLL_limit->SetPoint(i,m[i],signu_limit);
    }



  // plotting results:
  
  TLegend* l = new TLegend(0.2,0.4,0.5,0.9);
  TCanvas* c1 = new TCanvas("c1","c1",50,50,700,500);
  c1->SetLogx(1);
  gLL_2deltaLL[0]->Draw("AL");
  for(int i=0; i<nMass; i++)
    {
      cout << gLL_2deltaLL[i]->GetN() << endl;
      l->AddEntry(gLL_2deltaLL[i],os1[i].str().c_str(),"L");
      gLL_2deltaLL[i]->SetLineColor(i+1);
      //      gLL_2deltaLL[i]->SetLineStyle(i);
      if(i==0)
	{
	  gLL_2deltaLL[i]->Draw("AL");
	  gLL_2deltaLL[i]->GetYaxis()->SetRangeUser(-5,25);
	}
      else
	gLL_2deltaLL[i]->Draw("L");
    }
  l->Draw("same");
  
  TCanvas* c2 = new TCanvas("c2","c2",60,60,700,500);
  c2->SetLogy(1);
  c2->SetLogx(1);
  gLL_limit->SetLineWidth(2);
  gLL_limit->SetLineColor(2);
  gLL_limit->GetXaxis()->SetTitle("Mass [GeV]");
  gLL_limit->GetYaxis()->SetTitle("#LT#sigma#nu#GT [cm^{3}s^{-1}]");
  gLL_limit->GetYaxis()->SetRangeUser(1e-25,1e-22);
  gLL_limit->Draw("AL");
  
  return(gLL_limit);
}

TGraphAsymmErrors* plotLimitsAlexPass5f(string inFile)
{
  ifstream in(inFile.c_str());
  if(!in.is_open())
    {
      cout << "Warning! File not found!" << endl;
      return(NULL);
    }
  char str[300];
  in.getline(str,300);
  in.getline(str,300);
  double m,signu_med,signu_p1,signu_m1;
  int i = 0;
  TGraphAsymmErrors* g = new TGraphAsymmErrors();
  while ( in >> m >> signu_p1 >> signu_med >> signu_m1 )
    {
      g->SetPoint(i,m,signu_med);
      g->SetPointEYhigh(i,signu_p1 - signu_med);
      g->SetPointEYlow(i,signu_med - signu_m1);
      i++;
    }
  g->SetFillStyle(3001);
  g->SetFillColor(kGray);
  return(g);
}

void plotAllLimits()
{

  TGraph* gLimitAll =  ZitzLLCombine("combineList_all.txt");
  TGraph* gLimitSeg =  ZitzLLCombine("combineList_seg.txt");
  TGraph* gLimitUMi =  ZitzLLCombine("combineList_umi.txt");
  TGraph* gLimitDra =  ZitzLLCombine("combineList_dra.txt");
  TGraph* gLimitBoo =  ZitzLLCombine("combineList_boo.txt");

  TGraphAsymmErrors* gLimitAlex = plotLimitsAlexPass5f("../Pass5fSyst/limits_combined_tautau.txt");

  gLimitAll->SetLineColor(kBlue);
  gLimitSeg->SetLineColor(kRed);
  //  gLimitSeg->SetLineStyle(2);
  gLimitUMi->SetLineColor(kGreen+1);
  //gLimitUMi->SetLineStyle(3);
  gLimitDra->SetLineColor(kMagenta);
  //gLimitDra->SetLineStyle(4);
  gLimitBoo->SetLineColor(kCyan);
  //gLimitBoo->SetLineStyle(5);
  gLimitAll->SetLineWidth(1);
  gLimitSeg->SetLineWidth(1);
  gLimitUMi->SetLineWidth(1);
  gLimitDra->SetLineWidth(1);
  gLimitBoo->SetLineWidth(1);



  TLegend* l = new TLegend(0.6,0.2,0.9,0.5);
  TCanvas* c1 = new TCanvas("c1","c1",50,50,700,500);
  c1->SetLogy(1);
  c1->SetLogx(1);
  gLimitAll->Draw("AL");
  gLimitSeg->Draw("L");
  gLimitUMi->Draw("L");
  gLimitDra->Draw("L");
  gLimitBoo->Draw("L");
  gLimitAlex->Draw("L3");

  l->AddEntry(gLimitAlex,"Event Weight Combined","L");
  l->AddEntry(gLimitAll,"Max. Like. Combined","L");
  l->AddEntry(gLimitSeg,"Max. Like. Segue","L");
  l->AddEntry(gLimitUMi,"Max. Like. Ursa Minor","L");
  l->AddEntry(gLimitDra,"Max. Like. Draco","L");
  l->AddEntry(gLimitBoo,"Max. Like. Bootes","L");
  l->Draw("same");
}

void plotSegLimits()
{

  TGraph* gLimitAll  = ZitzLLCombine("combineList_all.txt");
  TGraph* gLimitSeg  = ZitzLLCombine("combineList_seg.txt");
  TGraph* gLimitNSeg = ZitzLLCombine("combineList_nseg.txt");

  TGraphAsymmErrors* gLimitAlex = plotLimitsAlexPass5f("../Pass5fSyst/limits_combined_tautau.txt");

  gLimitAll->SetLineColor(kBlue);
  gLimitSeg->SetLineColor(kRed);
  //  gLimitSeg->SetLineStyle(2);
  gLimitNSeg->SetLineColor(kGreen+1);
  //gLimitUMi->SetLineStyle(3);
  gLimitAll->SetLineWidth(1);
  gLimitSeg->SetLineWidth(1);
  gLimitNSeg->SetLineWidth(1);
  
  TLegend* l = new TLegend(0.6,0.2,0.9,0.5);
  TCanvas* c1 = new TCanvas("c1","c1",50,50,700,500);
  c1->SetLogy(1);
  c1->SetLogx(1);
  gLimitAll->Draw("AL");
  gLimitSeg->Draw("L");
  gLimitNSeg->Draw("L");
  gLimitAlex->Draw("L3");

  l->AddEntry(gLimitAlex,"Event Weight Combined","L");
  l->AddEntry(gLimitAll,"Max. Like. Combined","L");
  l->AddEntry(gLimitSeg,"Max. Like. Segue","L");
  l->AddEntry(gLimitNSeg,"Max. Like. no Segue","L");
  l->Draw("same");
}

TGraph* findMinAndSubtract(TGraph* g0,double y_null)
{
  TGraph* g = new TGraph();
  const int nSteps = 1e4;
  double xlow  = TMath::Log10(g0->GetX()[0]);
  const int n  = g0->GetN();
  double xhigh = TMath::Log10(g0->GetX()[n-1]);
  double dx = (xhigh-xlow)/nSteps;
  double y,y_min = 0;
  double x,x_min;
  for(int i=0; i<nSteps; i++)
    {
      if( i%1000 == 0 )
	cout << " x = " << x << " y = " << y << endl;
      x = TMath::Power(10.0,xlow + i*dx);
      y = g0->Eval(x);
      
      if( y < y_min )
	{
	  
	  y_min = y;
	  x_min = x;
	}
    }
  
  if( y_null < y_min )
    {
      y_min = y_null;
      x_min = 0;
    }
  
  cout << " minimum at: " << x_min << " with LL: " << y_min << endl; 
  for(int i=0; i<n; i++)
    {
      x = g0->GetX()[i];
      y = g0->GetY()[i];
      g->SetPoint(i,x,2.0*(y - y_min));
    }
  return(g);
}

double findLimit(TGraph* g,const double cl)
{
  const int nSteps = 1e4;
  double xlow  = TMath::Log10(g->GetX()[0]);
  const int n  = g->GetN();
  double xhigh = TMath::Log10(g->GetX()[n-1]);
  double dx = (xhigh-xlow)/nSteps;
  double y,y_min = 100;
  double x,x_min;
  double y2,x2;
  
  for(int i=0; i<nSteps; i++)
    {
      x  = TMath::Power(10.0,xlow + i*dx);
      y  = TMath::Abs(g->Eval(x) - cl);
      x2 = TMath::Power(10.0,xlow + (i+1)*dx);
      y2 = TMath::Abs(g->Eval(x2) - cl);
      
      if( y < y_min && (g->Eval(x) < g->Eval(x2)) )
	{
	  y_min = y;
	  x_min = x;
	}
    }
  if( y_min == 100.0)
    x_min = 1e-20;

  return(x_min);

}
