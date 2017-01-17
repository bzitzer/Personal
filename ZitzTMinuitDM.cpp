
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

//#include "VASlalib.h"
//#include "VASlalib_abr.cpp"

//#include "ZitzEventWeighting.cpp"
#include "ZitzMaxLike_v2.cpp"
//#include "plotLimitsAlex.cpp"
//#include "readAlexJFactor.cpp"
#include "TMinuit.h"

bool unfoldEnergy = 1;

bool testMode = 0;

bool expectedLimit = 0;

const int np = 30000;

const int nBins = 100;

double binEdges[nBins+1];

const double Elow  = 10.0;

const double Ehigh = 5e5;

vector<double> p_off,g_i,E_i;

string eaFile;

void fcn(int &npar, double *gin, double &f, double *par, int iflag);

TH1D* calc_signal(double *par);

TH2D* readPPP4DMFile(vector<double> &m_vals,decay_t decay = tt,string inSpecFile = "../AtProduction_gammas.dat/AtProduction_gammas.dat");

TMultiGraph* plotLimitsAlex(string inFile,decay_t chan = tt,string title = "");

TGraph* readIntJFactorFile(string inFile,dwarf_t dwarf = segue_1);

TH2D* hPPP4_spec = new TH2D();

TGraph* jProfile = new TGraph();

TH2D* hDispAvg = new TH2D();

TF1* fEnergyDisp = new TF1("fEnergyDisp","gaus",TMath::Log10(Elow),TMath::Log10(Ehigh));

TGraphAsymmErrors* plotLimitsAlexJSyst(string inFile);

TGraph* getAlexJFactor1D(string inFile, string enBinFile = "JPSFtables_alldwarfs/Eknots.txt",
			 string thBinFile = "JPSFtables_alldwarfs/thetaknots.txt", const double thCut = 0.17);

TH2D* readAlexJFactor(string inFile, string enBinFile, string thBinFile);

void runMultiTMinuitDM(string inFileList,string inEAFile,string inJFactorFile,string dispPath);

void ZitzTMinuitDM(string inEventFile = "/raid/reedbuck/bzitzer/Pass5f/segue1_eventList_pass5f_wZnCorr.txt",
		   string inEAFile = "/raid/reedbuck/bzitzer/Pass5f/segue1_eventList_pass5f_stg5noise_nZnCorr.txt.mce",
		   string inJFactorFile = "JPSFtables_alldwarfs/integratedJconvPSF_Seg1.txt",
		   string dispPath = "/raid/reedbuck/bzitzer/Pass5f/segue1_bias/",
		   string outRootFile = "segue_LLDM_outFile.root",dwarf_t dwarf = segue_1)
{
  // intial steps. Getting J-factor and DM spectra.
  eaFile = inEAFile;
  vector<double> m_vals;
  hPPP4_spec = readPPP4DMFile(m_vals);
  if(hPPP4_spec == NULL)
    {
      cout << "Warning! No DM Spectra!" << endl;
      return;
    }

  //string inJFactorFile;
  /*
  if(dwarf == segue_1)
    inJFactorFile = "JPSFtables_alldwarfs/integratedJconfPSF_Seg1.txt";
  if(dwarf == ursa_minor)
    inJFactorFile = "JPSFtables_alldwarfs/integratedJconfPSF_UMi.txt";
  if(dwarf == draco)
    inJFactorFile = "JPSFtables_alldwarfs/integratedJconfPSF_Draco.txt";
  if(dwarf == bootes)
    inJFactorFile = "JPSFtables_alldwarfs/integratedJconfPSF_BooI.txt";
  */
  
  TFile* f = new TFile(inJFactorFile.c_str(),"READ");
  jProfile = (TGraph*)f->Get("gConvJvE");
  //jProfile = readIntJFactorFile(inJFactorFile,dwarf);

  // J profile for GC projections:
  //jProfile = new TGraph(inJFactorFile.c_str());

  //this is the one to use for the most updated J factors for dSphs (2016/08/22)
  //jProfile = getAlexJFactor1D(inJFactorFile);

  if( jProfile == NULL )
    {
      cout << "warning! No J factor! " << endl;
      return;
    }
  f->Close();
  //jProfile->Draw("A*");
  //return;

  double expTot;
  double E;
  bool isOn;
  int runNum;
  int runNum_old = -1;
  double NonTot  = 0;
  double NoffTot = 0;
  double NbgTot;
  double w;
  double w_avg;
  double dE;
  // Range of mass values:
  // dwarf galaxies:
   
  double m_lower = 100; // GeV
  double m_upper = 1.0e5; // GeV
  const double m_num = 15;
  
  // Galactic center:
  /*  
  double m_lower = 1000; // GeV
  double m_upper = 1.0e5; // GeV
  const double m_num = 10;
  */

  m_vals.clear();
  for(int i=0; i<=m_num; i++)
    m_vals.push_back(m_lower*TMath::Power(10.0,TMath::Log10(m_upper/m_lower)*i/m_num));
  // bin edges for all histograms.
  double deltaE = (TMath::Log10(Ehigh) - TMath::Log10(Elow))/nBins;
  for(int i=0; i<=nBins; i++)
    binEdges[i] = Elow*TMath::Power(10.0,deltaE*i);

  TH1D* hBg0  = new TH1D("hBg0","hBg0",nBins,binEdges); // units of GeV
  TH1D* hBg   = new TH1D("hBg","hBg",nBins,binEdges); // units of GeV
  TH1D* hNon0 = new TH1D("hNon0","hNon0",nBins,binEdges); // units of GeV
  TH1D* hNon  = new TH1D("hNon","hNon",nBins,binEdges); // units of GeV
  
  TH2D* hDisp[400];
  int j = 0;
  double k = 0;
  TRandom3* r = new TRandom3(0);
  double val;
  double P,N,dE_p,runLT;
  int n,m;
  double x0,sigE,A;
  double Etmp;
  vector<double> LT;

  cout << "Reading event file.. ";
  TTree* tEv = (TTree*)readEventFile(inEventFile.c_str(),segue_1);
  TGraphAsymmErrors* gEA = getEffectiveArea(inEAFile.c_str(),
					    0,expTot);
  if( tEv == NULL )
    {
      cout << "problem with events tree!" << endl;
      return;
    }
  cout << "Done!" << endl;
  cout << "Getting effective area... ";
  
  if( gEA == NULL )
    {
      cout << "problem with effective areas!" << endl;
      return;
    }
  cout << "Done!" << endl;
  cout << "Total Exposure: " << expTot/3600 << endl;
  //TCanvas* c1 = new TCanvas();
  //c1->SetLogx(1);
  // gEA->Draw("A*");
  tEv->SetBranchAddress("isOn",&isOn);
  tEv->SetBranchAddress("energy",&E);
  tEv->SetBranchAddress("runNum",&runNum);
  tEv->SetBranchAddress("w",&w);
  tEv->SetBranchAddress("runLT",&runLT);
  //  const double eUpperCut = 5e3; // 5 TeV;
  const double eUpperCut = 1e5; // 100 TeV
  //const double eUpperCut = 3e4; // 30 TeV
  for(int i=0; i<tEv->GetEntries(); i++)
    {
      tEv->GetEntry(i);
      if( runNum_old != runNum )
	{
	  hDisp[j] = getEnergyDispersion(dispPath,runNum);
	  LT.push_back(runLT);
	  j++;
	}
      if(!expectedLimit)
	{
	  if(!isOn)
	    {	
	      hBg0->Fill(E,w);
	      if(unfoldEnergy)
		{
		  x0 = E;
		  sigE = TMath::Power(10.0,0.8*TMath::Log10(E));
		  A = 1.0/(sqrt(TMath::TwoPi())*sigE);

		  fEnergyDisp->SetParameter(0,A);
		  fEnergyDisp->SetParameter(1,x0);
		  fEnergyDisp->SetParameter(2,sigE);
		  for(int k=1; k<hBg->GetNbinsX(); k++)
		    {
		      Etmp = hBg->GetBinCenter(k);
		      dE = hBg->GetBinWidth(k);
		      hBg->Fill(Etmp,w*fEnergyDisp->Eval(Etmp)*dE);
		  /*
		  if(E > 1000)
		    cout << "Rec. Energy: " << E << " True E: " << Etmp << " weight: " << w << " P: " << fEnergyDisp->Eval(Etmp) << " dE: " << dE << endl;
		  */
		}
		}
	      if( E < eUpperCut)
		{
		  w_avg += w;
		  NoffTot++;
		}
	    }
	}
      if(expectedLimit)
	{
	  if(isOn)
	    {
	      w_avg += 1.0;
	      NoffTot++;
	      hBg0->Fill(E,1.0);
	    }

	}
       
      if(testMode) // test mode: using the energy of a fraction of the OFF events instead of the ON data. 
	{
	  if(!isOn)
	    {
	      k = r->Uniform(0,1);
	      if( k < w )
		{
		  hNon0->Fill(E,1.0);
		  NonTot++;
		  E_i.push_back(E);
		}
	    }
	}
      else
	{
	  if(isOn)
	    {
	      
	      if( E < eUpperCut )
	        {	      
		  hNon0->Fill(E,1.0);
		  NonTot++;
		  E_i.push_back(E);
		}
	    }
	}
    
      runNum_old = runNum;
    }
  w_avg = w_avg/NoffTot;
  
  //hDispAvg = renormalizeDisp(hDisp,j); // adds all energy dispersion matrixes and averages them.
  hDispAvg = renormalizeDispWeighted(hDisp,LT); // adds all energy dispersion matrixes and averages them.
  //hDispAvg->Smooth(1);
  //return;
  
  double norm = 0;
  //if(!unfoldEnergy)
  //hBg0->Smooth(4);
  
  for(int i=0; i<hBg0->GetNbinsX(); i++)
    {
      dE = hBg0->GetBinWidth(i);
      val = hBg0->GetBinContent(i);
      hBg0->SetBinContent(i,val/dE);
    }
  
  for(int i=0; i<hBg->GetNbinsX(); i++)
    {
      dE = hBg->GetBinWidth(i);
      val = hBg->GetBinContent(i);
      hBg->SetBinContent(i,val/dE);
    }
  
  /*
  if(unfoldEnergy)
    hBg = convolveEnergyDisp((TH1D*)hBg0->Clone(),hDispAvg);
  else
    hBg = hBg0;
  */
  if(!unfoldEnergy){ hBg = hBg0; }

  norm = hBg->Integral(0,nBins,"width");
  cout << "bg normalization: " << norm << endl;
  if(norm != 0)
    hBg->Scale(1.0/norm);
  
  norm = hBg0->Integral(0,nBins,"width");
  cout << "bg normalization: " << norm << endl;
  if(norm != 0)
    hBg0->Scale(1.0/norm);
  

  for(int i=0; i<hNon->GetNbinsX(); i++)
    {
      dE  = hNon0->GetBinWidth(i);
      val = hNon0->GetBinContent(i);
      hNon0->SetBinContent(i,val/dE);
    }
  
  norm = hNon0->Integral(0,nBins,"width");
  if(norm != 0)
    hNon0->Scale(NonTot/norm);
  
  for(unsigned int i=0; i<E_i.size(); i++)
    {
      
      if(unfoldEnergy)
	val = hBg0->Interpolate(E_i.at(i));
      else
        val = hBg0->Interpolate(E_i.at(i)); // dN/dE
      
      if( val == 0 )
	cout << "p_off = 0, E = " << E_i.at(i) << endl;
      p_off.push_back( val );

      //cout << p_off.at(i) << endl;
    }

  norm = hNon->Integral(0,nBins,"width");
  if(norm != 0)
    hNon->Scale(NonTot/norm);

  //return;

  TMinuit* t = new TMinuit(5);
  t->SetFCN(fcn);
  int ierflg;
  double b_null = (NonTot + NoffTot)/(1.0 + w_avg);
  //t->mnparm(0,"b    ",NoffTot,0.1 ,NoffTot-1,NoffTot+1,ierflg);
  t->mnparm(0,"b    ",b_null,0.1 ,0.0,1e7,ierflg);
  //t->mnparm(0,"b    ",NonTot/w_avg,0.1,NonTot/w_avg-1 ,NonTot/w_avg+1,ierflg);
  t->mnparm(1,"Mass ",m_vals.at(0),0.1,m_vals.at(0) , m_vals.at(0)+0.01,ierflg);
  //t->mnparm(2,"signu",-25.0  ,1e-4 ,-30.0  ,-20.0       ,ierflg);
  t->mnparm(2,"signu",1e-23  ,1e-30 ,1e-30  ,1e-20       ,ierflg);

  t->mnparm(3,"alpha",w_avg  ,0.1  ,w_avg  ,w_avg+0.01  ,ierflg);
  t->mnparm(4,"Noff ",NoffTot,0.1  ,NoffTot,NoffTot+0.01,ierflg);
 
  t->FixParameter(1);
  //t->FixParameter(2);
  t->FixParameter(3);
  t->FixParameter(4);
  t->SetPrintLevel(1);

  double arglist[10];
  arglist[0] = 500;
  arglist[1] = 1.;
  
  double par_null[5];
  par_null[0] = b_null;
  par_null[1] = 100.0;
  par_null[2] = 0.0;
  par_null[3] = w_avg;
  par_null[4] = NoffTot;

  double logl0;
  double *gin;
  int npar;
  fcn(npar,gin,logl0,par_null,ierflg);
  cout << "NULL logL: " << logl0 << endl;
  double loglmin;
  double b_min;
  double b_err;
  
  double signu_min = 0;
  double signu_err = 0;
 
  double par_min[5];
  double TS,dLLcheck;

  TGraph* gTS = new TGraph();
  TGraphErrors* gSigNu_min = new TGraphErrors();
  TGraph* gContour[m_vals.size()];
  TGraph* gLimit = new TGraph();
  TGraph* gLLcheck = new TGraph();
  TGraph* gLLmin = new TGraph();
  TGraph* gLLnull = new TGraph();
  TGraphErrors* gbmin = new TGraphErrors();
  // BJZ plotting shape of signal counts:
  
  par_min[0] = b_min;
  par_min[1] = 1000;
  par_min[2] = 1e-23;
  par_min[3] = w_avg;
  par_min[4] = NoffTot;
  /*
  TH1D* hDMSig1e3 = calc_signal(par_min);
  par_min[1] = 1e4; // 10 TeV
  TH1D* hDMSig1e4 = calc_signal(par_min);
  par_min[1] = 1e5; // 100 TeV
  TH1D* hDMSig1e5 = calc_signal(par_min);
 
  par_min[1] = 5e2; // 500 GeV
  TH1D* hDMSig5e2 = calc_signal(par_min);
  par_min[1] = 5e3; // 5 TeV
  TH1D* hDMSig5e3 = calc_signal(par_min);
  par_min[1] = 5e4; // 50 TeV
  TH1D* hDMSig5e4 = calc_signal(par_min);
 
  hDMSig1e3->SetLineColor(kBlue);
  hDMSig1e4->SetLineColor(kGreen+1);
  hDMSig1e5->SetLineColor(kRed);
  
  hDMSig5e2->SetLineColor(kBlack);
  hDMSig5e3->SetLineColor(kCyan);
  hDMSig5e4->SetLineColor(kGray+2);
  
  TCanvas* c1 = new TCanvas();
  c1->SetLogx(1);
  hDMSig1e3->Draw();
  hDMSig1e3->GetYaxis()->SetTitle("dN/dE [counts/GeV]");
  hDMSig1e3->GetXaxis()->SetTitle("Reconstructed Energy [GeV]");
  hDMSig1e4->Draw("same");
  hDMSig1e5->Draw("same");

  hDMSig5e2->Draw("same");
  hDMSig5e3->Draw("same");
  hDMSig5e4->Draw("same");

  
  TLegend* l0 = new TLegend(0.7,0.7,0.95,0.95);
  l0->AddEntry(hDMSig1e3,"M = 1 TeV","L");
  l0->AddEntry(hDMSig1e4,"M = 10 TeV","L");
  l0->AddEntry(hDMSig1e5,"M = 100 TeV","L");
  l0->AddEntry(hDMSig5e2,"M = 0.5 TeV","L");
  l0->AddEntry(hDMSig5e3,"M = 5 TeV","L");
  l0->AddEntry(hDMSig5e4,"M = 50 TeV","L");

  l0->Draw("same");
  return;
  */

  for(int i=0; i<m_vals.size(); i++)
    {
      t->mnparm(0,"b    ",b_null,0.1 ,0.0,1e7,ierflg);
      t->mnparm(1,"Mass ",m_vals.at(i),1e-3 ,
		m_vals.at(i)*0.9,m_vals.at(i)*1.1,ierflg);
      t->mnparm(2,"signu",1e-23  ,1e-30 ,1e-30  ,1e-20       ,ierflg);
      //    t->mnparm(2,"signu",-23.0  ,1e-4 ,-30.0  ,-20.0       ,ierflg);
      //t->FixParameter(0);
      t->FixParameter(1);
      //t->SetErrorDef(3.84/2);
      t->SetErrorDef(2.71/2);
      t->mnexcm("MIGRAD", arglist ,2,ierflg);
      cout << "Error flag: " << ierflg << endl;

      t->GetParameter(0,b_min,b_err);
      //t->GetParameter(1,index_min,index_err);
      t->GetParameter(2,signu_min,signu_err);
      gContour[i] = (TGraph*)t->Contour(40,0,2);
      //  if(gContour[i] == NULL){ return; }
      par_min[0] = b_min;
      par_min[1] = m_vals.at(i);
      par_min[2] = signu_min;
      par_min[3] = w_avg;
      par_min[4] = NoffTot;
      fcn(npar,gin,loglmin,par_min,ierflg);
      cout << "NULL logL: " << logl0 << " Minimized logL: " << loglmin << endl;
  
      TS = sqrt(-2*(loglmin - logl0));
      if(TMath::IsNaN(TS)){ TS = 0.0; }
      cout << "----------------" << endl;
      cout << "TS value: " << TS << endl;
      cout << "----------------" << endl;

      gTS->SetPoint(i,m_vals.at(i),TS);
      gLLnull->SetPoint(i,m_vals.at(i),logl0); // doesn't change with mass... but whatever.
      gSigNu_min->SetPoint(i,m_vals.at(i),signu_min);
      gSigNu_min->SetPointError(i,0,signu_err); 
      gLLmin->SetPoint(i,m_vals.at(i),loglmin);
      gbmin->SetPoint(i,m_vals.at(i),b_min);
      gbmin->SetPointError(i,0,b_err);
    }

  
  double par_check[5];
  double loglcheck;
  TGraph* gLL[m_vals.size()];
  m = 0;
  cout << "---------------" << endl;
  cout << "Final results: " << endl;
  cout << "---------------" << endl;
  for(int i=0; i<m_vals.size(); i++)
    {
      gLL[i] = new TGraph();
      if(gContour[i] != NULL)
	{
	  cout << TMath::MaxElement(gContour[i]->GetN(),gContour[i]->GetY()) << endl;
	  gLimit->SetPoint(m,m_vals.at(i),TMath::MaxElement(gContour[i]->GetN(),gContour[i]->GetY()));
	  m++;
	}
      else
	{
	  cout << "Warning! Contour for M = " << m_vals.at(i) << " is missing! " <<  endl;
	}
      //par_check[0] = gbmin->GetY()[i];
      par_check[0] = NoffTot;
      par_check[1] = m_vals.at(i);
      //par_check[2] = gLimit->GetY()[i];
      par_check[3] = w_avg;
      par_check[4] = NoffTot;
      //fcn(npar,gin,loglcheck,par_check,ierflg);
      //gLLcheck->SetPoint(i,m_vals.at(i),2.0*(loglcheck-logl0));
      
      for(int j=0; j<=70; j++)
	{
	  par_check[2] = TMath::Power(10.0,-27.0 + 10.0*j/100);
	  //par_check[2] = -27.0 + 10.0*j/100.0;
	  fcn(npar,gin,loglcheck,par_check,ierflg);
	  
	  cout << "Mass: " << m_vals.at(i) << " sig_nu:" << par_check[2] << " delta -2*logL :" << 2*(loglcheck - logl0) << endl;
	  gLL[i]->SetPoint(j,par_check[2],2*(loglcheck - gLLmin->GetY()[i]));
	}
     
      //      cout << "Mass: " << m_vals.at(i) << " minimized sig_nu: " << gSigNu_min->GetY()[i] << " +/- " << gSigNu_min->GetErrorY(i) << endl;
      if(gContour[i] != NULL)
	cout << " TS: " << gTS->GetY()[i] << " 95% CL limit: " << TMath::MaxElement(gContour[i]->GetN(),gContour[i]->GetY()) <<endl; 
      //cout << " LL check: " << gLLcheck->GetY()[i] << endl;
    }

  
  TCanvas* c0 = new TCanvas("c0","c0",40,40,700,500);
  c0->SetLogx(1);
  gbmin->GetXaxis()->SetTitle("Mass [GeV]");
  gbmin->GetYaxis()->SetTitle("Bg Normalization");
  gbmin->Draw("AE*");
  TLine* lNoff = new TLine(m_lower,b_null,m_upper,b_null);
  lNoff->SetLineColor(kRed);
  lNoff->SetLineStyle(2);
  lNoff->Draw("same");
  
  TCanvas* c1 = new TCanvas("c1","c1",50,50,700,500);
  c1->SetLogx(1);
  gTS->GetXaxis()->SetTitle("Mass [GeV]");
  gTS->GetYaxis()->SetTitle("TS^{1/2}");
  gTS->Draw("AL");
  

  TCanvas* c2 = new TCanvas("c2","c2",60,60,700,500);
  c2->SetLogx(1);
  gLL[0]->Draw("AL");
  gLL[0]->GetYaxis()->SetRangeUser(-10,20);
  gLL[0]->GetXaxis()->SetRangeUser(1e-25,1e-21);
  gLL[0]->GetXaxis()->SetTitle("#LT#sigma#nu#GT [cm^{3}s^{-1}]");
  gLL[0]->GetYaxis()->SetTitle("-2*#Delta log L");

  TLegend* l = new TLegend(0.2,0.7,0.5,0.9);
  ostringstream os[m_vals.size()];
  os[0] << "Mass = " << m_vals.at(0);
  l->AddEntry(gLL[0],os[0].str().c_str(),"L");
  for(int i=1; i<m_vals.size(); i++)
    {
      os[i] << "Mass = " << m_vals.at(i);
      l->AddEntry(gLL[i],os[i].str().c_str(),"L");
      gLL[i]->SetLineColor(i+1);
      gLL[i]->SetLineStyle(i);
      gLL[i]->Draw("L");
    }
  l->Draw("same");
  
  TCanvas* c3 = new TCanvas("c3","c3",70,70,800,600);
  c3->SetLogx(1);
  c3->SetLogy(1);
  gLimit->GetXaxis()->SetTitle("Mass [GeV]");
  gLimit->GetYaxis()->SetTitle("#LT#sigma#nu#GT [cm^{3}s^{-1}]");
  gLimit->SetLineColor(kBlue);
  gLimit->SetLineWidth(2);
  gLimit->Draw("AL");
  cout << "Writing results to file: " << outRootFile.c_str() << endl;
  ostringstream os1[m_vals.size()];
  ostringstream os2[m_vals.size()];

  TFile* fOut = new TFile(outRootFile.c_str(),"RECREATE");
  for( int i=0; i<m_vals.size(); i++ )
    {
      os1[i] << "gLL_signu_m" << (int)m_vals.at(i);
      gLL[i]->Write(os1[i].str().c_str());

      os2[i] << "gContour_m" << (int)m_vals.at(i);
      if(gContour[i] != NULL)
	gContour[i]->Write(os2[i].str().c_str());
    }
  gLLmin->Write("gLLmin");
  gLLnull->Write("gLLnull");
  gLimit->Write("gLimit");
  hBg->Write("hBg");
  hBg0->Write("hBg0");
  hNon0->Write("hNon");
  fOut->Close();
  /*
  TMultiGraph* gAlex = plotLimitsAlex("../Pass5dExpected/limits_expected_tautau.txt");
  
  gAlex->GetXaxis()->SetTitle("Mass [GeV]");
  gAlex->GetYaxis()->SetTitle("#LT#sigma#nu#GT [cm^{3}s^{-1}]");
  gAlex->Draw("AL3");
  gLimit->Draw("L");
  */
  cout << "NonTot: " << NonTot << " NoffTot: " << NoffTot << " average alpha: " << w_avg << endl;
  
}

void fcn(int &npar, double *gin, double &f, double *par, int iflag)
{

  double logl  = 0;
  double b     = par[0];
  double M     = par[1];
  //double signu = par[2];
  
  double signu = TMath::Power(10.0,par[2]);
  //cout << "signu: " << signu << endl;
  
  double alpha = par[3];
  double Noff  = par[4];
 
  double g = 0;
  double p_on;
  double dE,E,dE_p;
  double val;
  double norm;
  double P;
  double N;
  double sigE,x0,A;
  int nBin_g;
  int n,m;
  TH1D*  hg0 = calc_signal(par);
  TH1D*  hg = new TH1D();
  //if(unfoldEnergy)
  //  hg = convolveEDisp_varBinWidth(hg0,hDispAvg);
  //else
  //  hg = hg0;

  if(hg0 == NULL)
    {
      f = logl;
      return;
    }

  g = hg0->Integral(1,nBins,"width");
  if(TMath::IsNaN(g))
    {
      cout << "warning! expected signal events is NaN! " << endl;
      for(int i=0; i<hg0->GetNbinsX(); i++)
	cout << hg0->GetBinCenter(i) << " " << hg0->GetBinContent(i) << endl;
      
      f = logl;
      return;
    }
  
  logl = Noff*TMath::Log(b) - g - (alpha + 1)*b;
  for(int i=0; i<p_off.size(); i++)
    {
      if(unfoldEnergy)
	p_on = hg0->Interpolate( E_i.at(i) );
      else
	p_on = hg0->Interpolate( E_i.at(i) );

      val = alpha*b*p_off.at(i) + p_on;
      
      if(val > 0)
	logl += TMath::Log(val);

    }

  hg->Delete();
  hg0->Delete();
  //  cout << "logl: " << -logl << " Number of expected events: " << g << endl;
  f = -logl;
}

TH1D* calc_signal(double *par)
{
  double b = par[0];
  double M = par[1];
  double sig_nu = par[2];
  /*
  if( par[2] == 0.0 )
    sig_nu = 0.0;
  else
    sig_nu = TMath::Power(10.0,par[2]);
  */
  double expTot;
  double E,Elog10TeV,A,dE;
  double N;
  double dNdx;
  double dNdlog10x;
  double x;
  double j_E;
  double Etr,sum,dEtr,D,norm;
  int n,m;
  int k = 0;
  //cout << "Getting effective area from: " << eaFile << endl;
  TGraphAsymmErrors* gEA = getEffectiveArea(eaFile.c_str(),0,expTot);
  
  if(gEA == NULL)
    {
      cout << "Warning! No Effective area!" << endl;
      return(NULL);
    }
  //  TH1D* hg = new TH1D("hg","hg",nBins,Elow,Ehigh);
  double phi0 = sig_nu*expTot/(8*TMath::Pi()*pow(M,2));
  double s;

  TH1D* hg = new TH1D("hg","hg",nBins,binEdges);

  for(int i=1; i<=hg->GetNbinsX(); i++)
    {
      E = hg->GetBinCenter(i);
      dE = hg->GetBinWidth(i);
     
      if(!unfoldEnergy)
	{
	  x = E/M;
	  if( x >= 1.0 )
	    {
	      hg->SetBinContent(i,0);
	      continue;
	    }
	  Elog10TeV = TMath::Log10(E/1000.0);
	  A = gEA->Eval(Elog10TeV);
	  A *= 1e4; // m^2 -> cm^2

	  dNdlog10x = hPPP4_spec->Interpolate(TMath::Log10(x),E);
	  dNdx = dNdlog10x*TMath::Log10(TMath::Exp(1.0))/x; // dN/dlog_{10}x -> dN/dx
	  j_E = jProfile->Eval(E);
	  //E = TMath::Power(10.0,hg->GetBinCenter(i));
      
	  s = phi0*j_E*dNdx*A/M;
	}
      else
	{
	  sum = 0;
	  norm = 0;
	  for(int j=1; j<=hDispAvg->GetXaxis()->GetNbins(); j++)
	    {
	      Etr  = hDispAvg->GetXaxis()->GetBinCenter(j);
	      dEtr = hDispAvg->GetXaxis()->GetBinWidth(j);
	      //k = hDispAvg->GetYaxis()->FindBin(E);
	      /*
	      Etr = hg->GetBinCenter(j);
	      k = hDispAvg->GetXaxis()->FindBin(Etr);
	      dEtr = hDispAvg->GetXaxis()->GetBinWidth(k);
	      */
	      x = Etr/M;	    
	      if(x >= 1.0)
		dNdlog10x = 0.0;
	      else
		dNdlog10x = hPPP4_spec->Interpolate(TMath::Log10(x),Etr);	      
	      dNdx = dNdlog10x*TMath::Log10(TMath::Exp(1.0))/x; // dN/dlog_{10}x -> dN/dx
	       
	      Elog10TeV = TMath::Log10(Etr/1000.0);
	      A = gEA->Eval(Elog10TeV);
	      A *= 1e4; // m^2 -> cm^2
	      if(Etr > 7e4)
		j_E = jProfile->Eval(7e4);
	      else
		j_E = jProfile->Eval(Etr);
	      
	      //D = hDispAvg->GetBinContent(j,k);
	      D = hDispAvg->Interpolate(Etr,E);
	      
	      		

	      norm += D*dEtr; 
	      sum  += A*j_E*dNdx*D*dEtr/M;
	    }
	  //cout << E << " " << norm << endl;
	  if(norm != 0)
	    s = phi0*sum/norm;
	  else
	    s = 0;
	}
      hg->SetBinContent(i,s);
    }
  /*
  if(unfoldEnergy)
    {
      for(int i=1; i<=hg->GetNbinsX(); i++)
	{
	  dE  = hg->GetBinWidth(i);
	  s = hg->GetBinContent(i);
	  hg->SetBinContent(i,s/dE); 
	}
    }
  */
  //cout << "Total exposure: " << expTot/3600 << endl;
  gEA->Delete();
  
  return(hg);

}


TH2D* readPPP4DMFile(vector<double> &m_vals,decay_t decay,string inSpecFile)
{
  double m,x;
  double eL,eR,e,muL,muR,mu,tauL,tauR,tau,q,c,b,t,WL,WT,W,ZL,ZT,Z,g,gamma,h,NuE,NuMu,NuTau,Ve,VMu,VTau;
  ifstream in(inSpecFile.c_str());
  if(!in.is_open())
    {
      cout << "File missing!" << endl;
      return(NULL);
    }
  char str[1000];
  in.getline(str,1000);
  //  in.getline(str,1000);

  int i = 0;
  bool foundLowerMass = 0;
  bool foundUpperMass = 0;

  vector<double> x_vals;
  //vector<double> m_vals;
  TGraph2D* gSpec = new TGraph2D();
  x = -8.9;
  while( x <= 0 )
    {
      x_vals.push_back(x);
      x += 0.05;
    }

  vector< vector<double> > spec(x_vals.size());
  
  double m_old = -1;
  int j=0;

  while( in >> m >> x >>
	 eL >> eR >> e >> muL >> muR >> mu >> tauL >> tauR >> tau >>
	 q >> c >> b >> t >> WL >> WT >> W >> ZL >> ZT >> Z >> g >> gamma >>
	 h >> NuE >> NuMu >> NuTau >> Ve >> VMu >> VTau
	 )
    {
      
      //cout << x << endl;
      if( m != m_old )
	{
	  m_vals.push_back(m);
	  i = 0;
	}

      if( TMath::Abs(x - x_vals.at(i)) > 0.01 )
	{
	  cout << "Warning! x values don't match! " << x << " " << x_vals.at(i)
	       << " Mass: " << m << endl;
	  cout << x - x_vals.at(i) << endl;
	}
      
      if( decay == tt )
	spec.at(i).push_back(tau);
      if( decay == bbar )
	spec.at(i).push_back(b);
      if( decay == ee )
	spec.at(i).push_back(e);
      if( decay == uu )
	spec.at(i).push_back(mu);
      if( decay == WW )
	spec.at(i).push_back(W);
      if( decay == ZZ )
	spec.at(i).push_back(Z);

      if( decay == tt )
	gSpec->SetPoint(j,x,m,tau);
      if( decay == bbar )
	gSpec->SetPoint(j,x,m,b);
      if( decay == ee )
	gSpec->SetPoint(j,x,m,e);
      if( decay == uu )
	gSpec->SetPoint(j,x,m,mu);
      if( decay == WW )
	gSpec->SetPoint(j,x,m,W);
      if( decay == ZZ )
	gSpec->SetPoint(j,x,m,Z);
	
      i++;
      j++;
      m_old = m;

      if( decay == WW && m < 90.0 && W > 0.0 )
	cout << m << " " << W << endl;
      if( decay == ZZ && m < 90.0 && Z > 0.0 )
	cout << m << " " << Z << endl;
    }
  
 
  double m_vals_array[m_vals.size()+1];
  double x_vals_array[x_vals.size()+1];
  for(int i=1; i<=m_vals.size(); i++)
    m_vals_array[i] = m_vals.at(i-1);
  for(int i=1; i<=x_vals.size(); i++)
    x_vals_array[i] = x_vals.at(i-1);

  m_vals_array[0] = 0;
  x_vals_array[0] = -9.0;
  /*
  TH1D* hSpec_m[m_vals.size()];
  ostringstream os[m_vals.size()];
  for(int i=0; i<m_vals.size(); i++)
    {
      os[i] << "Spectrum for M=" << m_vals.at(i);
      hSpec_m[i] = new TH1D(os[i].str().c_str(),os[i].str().c_str(),
			    x_vals.size(),-8.925,0.025);
      
      for(int j=1; j<=x_vals.size(); j++)
	{
	  //	  cout << "m = " << m_vals.at(i) << " x = "<< hSpec_m[i]->GetBinCenter(j) << endl;
	  hSpec_m[i]->SetBinContent(j,spec.at(j-1).at(i));
	}
      hSpec_m[i]->GetXaxis()->SetTitle("log_{10} E/M_{#chi}");
      hSpec_m[i]->GetXaxis()->SetTitle("dN/dx");
      
    }
  */
  TH2D* hSpec = new TH2D("hSpec","hSpec",
			 x_vals.size(),x_vals_array,m_vals.size(),m_vals_array);
  hSpec->GetXaxis()->SetTitle("log_{10} x");
  hSpec->GetYaxis()->SetTitle("M_{#chi} [GeV]");
  hSpec->SetStats(0);
  /*
  cout << spec.size() << " " << spec.at(0).size() << endl;
  cout << hSpec->GetNbinsX() << " " << hSpec->GetNbinsY() << endl;
  cout << x_vals.size() << " " << m_vals.size() << endl;
  */
  for(int i=0; i<spec.size(); i++) // for each x value
    for(int j=0; j<spec.at(0).size(); j++) // for each mass value
	{
	  //	  cout << spec.at(i).size() << endl;
	  if( TMath::Abs(hSpec->GetYaxis()->GetBinUpEdge(j) - m_vals_array[j]) > 0.01 )
	    cout << "Y values not equal! " << hSpec->GetYaxis()->GetBinUpEdge(j) << " " << m_vals_array[j] << endl;
	  if( TMath::Abs(hSpec->GetXaxis()->GetBinUpEdge(i) - x_vals_array[i]) > 0.01 )
	    cout << "X values not equal! " << hSpec->GetXaxis()->GetBinUpEdge(i) << " " << x_vals_array[i] << endl;
	 
	  hSpec->SetBinContent(i+1,j+1,spec.at(i).at(j));
	 	
	}

  return(hSpec);
}

TGraphAsymmErrors* plotLimitsAlexJSyst(string inFile)
{

  ifstream in(inFile.c_str());
  char str[300];
  in.getline(str,300);
  cout << str << endl;
  //string title = str;
  in.getline(str,300);
  in.getline(str,300);
  TGraphAsymmErrors* g = new TGraphAsymmErrors();
  int i = 0;
  double m,med,plus1sig,minus1sig;
  while( in >> m >> med >> minus1sig >> plus1sig )
    {
      g->SetPoint(i,m,med);
      g->SetPointEYlow(i,med - minus1sig);
      g->SetPointEYhigh(i,plus1sig - med);
      i++;
    } 

  g->SetFillStyle(3001);
  g->SetFillColor(kGray);
  g->SetLineColor(kBlack);

  g->GetXaxis()->SetTitle("M_{#chi} [GeV]");
  g->GetYaxis()->SetTitle("#LT#sigma#nu#RT [cm^{3}s^{-1}]");

  return(g);
}

TMultiGraph* plotLimitsAlex(string inFile, decay_t chan, string title)
{
  ifstream in(inFile.c_str());
  if(!in.is_open())
    {
      cout << "No file, Moron!" << endl;
      return(NULL);
    }
  char str[300];
  in.getline(str,300);
  cout << str << endl;
  in.getline(str,300);
  cout << str << endl;
  in.getline(str,300);
  cout << str << endl;
  //TGraph* gBayes = new TGraph("Bayesian/pass5f_trimmed_bayes.txt");
  //gBayes->SetLineColor(kBlue);
  TGraphAsymmErrors* g[3];
  for(int i=0; i<3; i++)
    g[i] = new TGraphAsymmErrors();

  TMultiGraph* mg = new TMultiGraph();

  double m,obsLimit,expLimit,minus1sigLimit,minus2sigLimit,plus1sigLimit,plus2sigLimit;

  int i=0;
  while( in >> m >> obsLimit >> minus2sigLimit >> minus1sigLimit >>
	 expLimit >> plus1sigLimit >> plus2sigLimit )
    {
      //cout << obsLimit << " << 
      g[0]->SetPoint(i,m,obsLimit);
      g[1]->SetPoint(i,m,expLimit);
      g[2]->SetPoint(i,m,expLimit);
      g[1]->SetPointEYlow(i,expLimit - minus1sigLimit);
      g[2]->SetPointEYlow(i,expLimit - minus2sigLimit);
      g[1]->SetPointEYhigh(i,plus1sigLimit - expLimit);
      g[2]->SetPointEYhigh(i,plus2sigLimit - expLimit);
      i++;
    }
  //mg->GetXaxis()->SetTitle("M_{#chi} [GeV]");
  //mg->GetYaxis()->SetTitle("#langle#sigma#nu#rangle [cm^{3}s^{-1}]");
  g[0]->SetLineColor(kBlack);
  
  for(int i=0; i<3; i++)
    {
      g[i]->SetLineWidth(2);
      g[i]->SetFillStyle(3001);
      g[i]->SetLineColor(i+1);
      g[i]->SetFillColor(i+1);
      mg->Add(g[2-i]);
    }

  return(mg);
}

TGraph* readIntJFactorFile(string inFile,dwarf_t dwarf)
{
  ifstream in(inFile.c_str());
  char str[300];
  in.getline(str,300);
  in.getline(str,300);
  double E,seg,umi,dra,boo1,J;
  TGraph* g = new TGraph();
  int i = 0;
  while( in >> E >> boo1 >> seg >> dra >> umi )
    {
      if(dwarf == segue_1){ J = seg; }

      if(dwarf == draco){ J = dra; }

      if(dwarf == bootes){ J = boo1; }

      if(dwarf == ursa_minor){ J = umi; }

      g->SetPoint(i,E,J);

      i++;
    }
  return(g);

}

TH2D* readAlexJFactor(string inFile, string enBinFile, string thBinFile)
{
  ifstream in(inFile.c_str());
  if(!in.is_open())
    {
      cout << "File not here... blah blah... " << endl;
      return(NULL);
    }
  char str[300];
  in.getline(str,300);
  in.getline(str,300);
  //vector<double> en,th;
  int enN,thN;
  double* th = readBinFile2(thBinFile,thN,1);
  double* en = readBinFile2(enBinFile,enN,1);
  double* enCen = new double[enN];
  for(int i=0; i<enN; i++)
    enCen[i] = pow(10.0,TMath::Log10(en[i]) - 0.05);
  TH2D* h = new TH2D("hJProfAlex","Integrate J vs. E",enN-1,enCen,thN-1,th);
  double tmp;
  int i = 1;
  int j = 1;
  //double c[] = {16.0, 16.25, 16.5, 16.75, 17.0, 17.25, 17.5, 17.75, 18.0, 18.25, 18.5, 18.75, 19.0, 19.25, 19.5};
  //double c[] = {17.0, 17.1, 17.2, 17.3, 17.4, 17.5, 17.6, 17.7, 17.8, 17.9, 18.0,
  //		18.0, 18.1, 18.2, 18.3, 18.4, 18.5, 18.6, 18.7, 18.8, 18.9, 19.0,
  //		19.0, 19.1, 19.2, 19.3, 19.4, 19.5};
  double c[] = { 18.0, 18.1, 18.2, 18.3, 18.4, 18.5, 18.6, 18.7, 18.8, 18.9, 19.0,
		 19.0, 19.1, 19.2, 19.3, 19.4, 19.5};
  cout << "Got here" << endl;
  while( in >> tmp )
    {
      if(tmp > 1e17)
	h->SetBinContent(j,i,TMath::Log10(tmp));
      //cout << i << " " << j << endl;
      //cout << h->GetXaxis()->GetBinCenter(i) << " " << h->GetYaxis()->GetBinCenter(j) << " " << tmp << endl;
      i++;
      if(i>thN)
	{
	  i = 1;
	  j++;
	}
    }
  h->SetMinimum(18.0);
  h->SetMaximum(19.5);
  h->SetContour(15,c);
  /*
  TLine* l = new TLine(60.0,0.17,2e5,0.17);
  l->SetLineStyle(2);
  l->SetLineWidth(2);
  TCanvas* c1 = new TCanvas();
  c1->SetLogy(1);
  c1->SetLogx(1);
  h->SetStats(0);
  h->Draw("colz");
  l->Draw("same");
  */
  h->GetXaxis()->SetTitle("Energy [GeV]");
  h->GetYaxis()->SetTitle("#theta [deg]");
  h->GetYaxis()->SetRangeUser(0.01,1.0);
  return(h);

}

TGraph* getAlexJFactor1D(string inFile, string enBinFile, 
			  string thBinFile,const double thCut)
{
  TH2D* h = readAlexJFactor(inFile,enBinFile,thBinFile);
  TGraph* g = new TGraph();
  
  int nBinY = h->GetYaxis()->FindBin(thCut);
  double E,J;
  int j = 0;
  for(int i=1; i<=h->GetNbinsX(); i++)
    {
      E = h->GetXaxis()->GetBinCenter(i);
      J = pow(10.0,h->Interpolate(E,thCut));
      if( E > 80.0 )
	//if(E > 0.0)
	{
	  g->SetPoint(j,E,J);
	  j++;
	}
            
      else
	{
	  g->SetPoint(j,E,0.0);
	  j++;
	}
      
      cout << "Energy: " << E << " Int. J Factor: " << J << endl;
    }
  g->GetXaxis()->SetTitle("Energy [GeV]");
  g->GetYaxis()->SetTitle("J Factor [GeV^{2}cm^{-5}]");
  return(g);
}

void runMultiTMinuitDM(string inFileList,string inEAFile,string inJFactorFile,string dispPath)
{
  gROOT->SetBatch(1);
  const int nFiles = 100;
  ostringstream os[nFiles];
  ifstream in(inFileList.c_str());
  string str;
  int i=0;
  testMode = 0;
  expectedLimit = 0;
  //while(i < nFiles)
  while( in >> str )
    {
      os[i] << "umaII/test_outDM" << i << "_tt_unfE_p50hrs_bdt_ext.root"; 
      ZitzTMinuitDM(str,inEAFile,inJFactorFile,dispPath,os[i].str());
      p_off.clear();
      g_i.clear();
      E_i.clear();

      i++;
    }



}
