
/*
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <functional>
#include <sstream>

#include "TCanvas.h"
#include "TGraph.h"
#include "TFile.h"
#include "TTree.h"
#include "TDirectory.h"
#include "TSystemDirectory.h"
#include "TList.h"
#include "TF1.h"
#include "TIterator.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TF2.h"
#include "TStyle.h"
#include "TProfile.h"
#include "TChain.h"
#include "TLine.h"
#include "TBox.h"
#include "TEllipse.h"
#include "TPaletteAxis.h"
#include "TAxis.h"
#include "TMarker.h"
#include "TMath.h"
#include "TPaveStats.h"
#include "TPaveText.h"
#include "TColor.h"
#include "TKey.h"
*/

//#include "aclicPreProcCommands.h"
//#include <Riostream.h>

//#include "VAUpperLimit.h"

enum decay_t {WW, ZZ, bbar, tt, ee, uu};
enum dwarf_t {segue_1, draco, ursa_minor, wilman_1, bootes, smith};
enum prof_t {einasto, NFW};

TGraph* plotDMCurve(string inFile = "upperlimit.combined.root",decay_t decay = WW,dwarf_t dwarf = segue_1, double theta2Cut = 0.03, double Emin_cut = 0, bool writeFile = true)
{
  string effAreaFile = "Pass3a/Segue_SoftCut_pass3a_mcE.txt.ea";
  double J = calcJFactor(dwarf,einasto,theta2Cut);
  //double J_old = 7.7e18; // units of GeV^2 cm^-5 sr
  double J_old = 9.0e18; // From Matthieu, no convolution
  //  double J = 7.7e18;

  double WWParam[] = {-1.36,0,0,0.5,-12.6,12.1,-9.86};
  double ZZParam[] = {-1.40,0,0,0.47,-14.4,16.3,-13.6};
  double bbarParam[] = {0,0,0,1.05,-17.8,12.3,-1.86};
  TF1* fDMSpectra = new TF1("fDMSpectra","pow(x,[0] + [1]*TMath::Log(x) + [2]*pow(TMath::Log(x),2))*exp([3] + [4]*x + [5]*(x**2.0) + [6]*(x**3.0))",0.01,1);
  
  if(decay == WW){ fDMSpectra->SetParameters(WWParam); }
  if(decay == ZZ){ fDMSpectra->SetParameters(ZZParam); }
  if(decay == bbar){ fDMSpectra->SetParameters(bbarParam); }
  if(decay == tt)
    {
      double ttParam[] = {-1.29,7.83,-2.70,-9e-3,-5.06};
      TF1* fDMSpectra_tau = new TF1("fDMSpectra_tau","pow(x,[0])*([1]*x + [2]*(x**2.0) + [3]*(x**3.0))*exp([4]*x)",1e-4,1);
      fDMSpectra_tau->SetParameters(ttParam);
      fDMSpectra = fDMSpectra_tau;
    }
  /*
  if( decay == ee || decay == uu )
    {
      TF1* fDMSpectra_lep = new TF1("fDMSpectra_lep","[0]/(x*TMath::Pi())*((1+(1+x)**2.0 - ([1]/[2])**2.0)*TMath::Log((1+[3])/(1-[3])) - 2*(1-x)*[3])",0.01,1.0);
      double m_l;
      if( decay == ee )

      fDMSpectra = fDMSpectra_lep; 
    }
  */
  //fDMSpectra->Draw();
  //c1->SetLogy(1);

  TFile* f = new TFile(inFile.c_str(),"READ");
  if( !f->IsOpen() )
    {
      cout << "Problem with UL! " << endl;
      return;
    }
 
  gDirectory->cd("UpperLimit");
  VAUpperLimit* vaUL = (VAUpperLimit*)gDirectory->Get("VAUpperLimit");
  if( vaUL == NULL )
    {
      cout << "No UL Object present!" << endl;
      return;
    }
  double T_obs = vaUL->GetLiveTime();
  double N_ul = vaUL->GetUL();
  // adding code to use the Helene method instead:
  /*if( N_ul < 10.0) // if UL counts is below 10, uee Helene instead of Rolke. This is an arbitary cut. 
    {
      vaUL->SetMethod("Helene");
      vaUL->CalculateUpperLimit();
      vaUL->PrintUpperLimit();
      
    }
  */
  //double N_ul = vaUL->GetUL();
  //if(vaUL->GetUL() < 10.0)
  //  {
  cout << "Using Bound Rolke Method..." << endl;
  TRolke* t = new TRolke(0.95);
  t->SetBounding(false);
  t->SetPoissonBkgKnownEff(vaUL->GetTotON(),vaUL->GetTotOFF(),1.0/vaUL->GetAlpha(),1.0);
  N_ul = t->GetUpperLimit();
      //   } 

  TGraphAsymmErrors* gEffArea = vaUL->GetEffectiveArea();
  if( gEffArea == NULL )
    {
      cout << "Problem with effective areas!" << endl;
      return;
    }
  // TGraphAsymmErrors* gEffArea = ReadEATextFile();
  //double EThresh = vaUL->GetEThreshold();

  //TF1* fAnnilSpectrum = new TF1("A_eff(E)dN/dE",EThresh,10);
  TGraphAsymmErrors* gAnnilSpectrum = new TGraphAsymmErrors();
  double E;
  for(int i=0; i<gEffArea->GetN(); i++)
    {
      E = pow(10.0,gEffArea->GetX()[i])*1000; // GeV
      gAnnilSpectrum->SetPoint(i,E,gEffArea->GetY()[i]*fDMSpectra->Eval(E/1e3));
    }
  
  // <sigma v> = 8pi/J * N_95*m_chi^2/(T_obs*int(A_eff*dN/dx dx)
  
  TGraph* gMchiSigmaV = new TGraph();
  TGraph* gMchiSigmaVOld = new TGraph();
  gMchiSigmaV->SetLineColor(kRed);
  gMchiSigmaVOld->SetLineColor(kBlue);
  
  double m_chi_min = 100;
  //double m_chi_max = 2.5e4;
  double m_chi_max = 1e5
  double m_chi = m_chi_min;
  int nSteps = 500;
  double sigma_v;
  double pi = TMath::Pi();
  double intDMFlux;
  int j = 0;
  double min_m_chi = 0;
  double sigma_v_min = 1;
  // Mathiew's result:
  double sigma_v_old;
  double T_obs_old = 47.524833*3600;
  double N_ul_old = 135.9;
  //  double s = 0.6;
  double s = 1.0;

  for(int i=0; i<nSteps; i++)
    {
    
      intDMFlux = integrateAnnilSpec2(gEffArea,fDMSpectra,m_chi,decay);
      //intDMFlux = integrateAnnilSpec(gEffArea,fDMSpectra,m_chi,decay,1e-4,Emin_cut);
      if(intDMFlux > 0.0 )
	{
	  if(i%10 == 0 )
	    cout << "m_chi=" << m_chi << " sigma_v="<<sigma_v << " sigma_v_old=" << sigma_v_old << endl;

	  sigma_v = s*(8*pi/J)*(N_ul*pow(m_chi,2)/(T_obs*intDMFlux));
	  sigma_v_old = s*(8*pi/J_old)*(N_ul_old*pow(m_chi,2)/(T_obs_old*intDMFlux));

	 
	  if(!TMath::IsNaN(sigma_v))
	    {
	      gMchiSigmaV->SetPoint(j,m_chi,sigma_v);
	      gMchiSigmaVOld->SetPoint(j,m_chi,sigma_v_old);
	  
	      j++;
	      if(sigma_v < sigma_v_min)
		{
		  sigma_v_min = sigma_v;
		  min_m_chi = m_chi;
		} 
	    }
	}	    
      m_chi = m_chi*pow(10.0,TMath::Log10(m_chi_max/m_chi_min)/nSteps);
    }
  TCanvas* c1 = new TCanvas();
  gMchiSigmaV->Draw("AL");
  //gMchiSigmaVOld->Draw("L");
  gMchiSigmaV->GetYaxis()->SetRangeUser(1e-27,1e-22);
  // c1->SetRangeUser(10e-26,10e-21);
  TMarker* pt = new TMarker();
  //pt->SetMarkerColor(kGreen);
  pt->SetMarkerStyle(kMultiply);
  pt->SetMarkerSize(1.5);

  TMarker* pt_mag = new TMarker();
  pt_mag->SetMarkerColor(kGreen+2);
  pt_mag->SetMarkerStyle(3);
  pt_mag->SetMarkerSize(1.5);

  TLine* th_limit = new TLine(m_chi_min,3e-26,m_chi_max,3e-26);
  th_limit->SetLineColor(kBlue);
  th_limit->SetLineStyle(2);
  th_limit->SetLineWidth(3);
  th_limit->Draw("same");
  if(decay == bbar)
    {
      pt->DrawMarker(1625.4,8.90e-24);
      gMchiSigmaV->SetTitle("95% Confidence DM Exclusion curves for b#bar{b}");
    }
  if(decay == tt)
    {
      pt->DrawMarker(382.56,1.92e-24);
      pt_mag->DrawMarker(500.0,1.24e-24);
      gMchiSigmaV->SetTitle("95% Confidence DM Exclusion curves for #tau^{+}#tau^{-}");
    }
  if(decay == WW)
    {
      pt->DrawMarker(1033.54,7.95e-24);
      gMchiSigmaV->SetTitle("95% Confidence DM Exclusion curves for W^{+}W^{-}");
    }
  if(decay == ZZ)
    gMchiSigmaV->SetTitle("95% Confidence DM Exclusion curves for Z^{0}Z^{0}");

  c1->SetLogy(1);
  c1->SetLogx(1);
  gMchiSigmaV->GetXaxis()->SetTitle("m_{#chi} (GeV)");
  gMchiSigmaV->GetYaxis()->SetTitle("<#sigma v>^{95% CL} (cm^{3} s^{-1})");

  TLegend* leg = new TLegend(0.6,0.6,0.9,0.9);
  leg->AddEntry(gMchiSigmaV,"94 hr limit","L");
  if(decay == tt)
    leg->AddEntry(pt_mag,"160 hr limit from MAGIC","P");
 
  //leg->AddEntry(gMchiSigmaVOld,"48 hr limit","L");
  leg->AddEntry(pt,"Minimum from Aliu et. al. (2012), 48hr limit","P");
  leg->AddEntry(th_limit,"Natural Thermal cross-section limit","L");
  leg->Draw("same");
  cout << " Number of UL counts: " << N_ul << endl;
  cout << " Observation Time: " << T_obs/3600 << " hrs" << endl;
  cout << " Number of UL counts (previously): " << N_ul_old << endl;
  cout << " Observation Time (previously): " << T_obs_old/3600 << " hrs" << endl;
  
  cout << " Minimal <sigma_v>^95% :" << sigma_v_min << " at m_chi: " << min_m_chi << " GeV" << endl;
  //TCanvas* c1 = new TCanvas();
  //vaUL->Draw();
  /*  TCanvas* c1 = new TCanvas();
  c1->SetLogy(1);
  gEffArea->Draw("AE*");
  */
  cout << "<ov>^{95%} at 170 GeV: " << gMchiSigmaV->Eval(170.0) << endl;
  cout << "<ov>^{95%} at 300 GeV: " << gMchiSigmaV->Eval(300.0) << endl;
  cout << "<ov>^{95%} at 500 GeV: " << gMchiSigmaV->Eval(500.0) << endl;
  cout << "<ov>^{95%} at 1 TeV: " << gMchiSigmaV->Eval(1000.0) << endl;
  cout << "<ov>^{95%} at 10 TeV: " << gMchiSigmaV->Eval(1e4) << endl;
  if(writeFile)
    {
      cout << "Writting to file...";
      ofstream os("DMLimits.txt");
      os << "E       Eff_Area (cm^2)" << endl;
      os << "----------------------" << endl;
      for(int i=0; i<gEffArea->GetN(); i++)
	os << pow(10.0, gEffArea->GetX()[i] + 3.0) << "   " << gEffArea->GetY()[i]*1e4 << endl;

      os << "m_chi    <ov>^{95%}" << endl;
      os << "-------------------" << endl;
      for(int i=0; i<gMchiSigmaV->GetN(); i++)
	os << gMchiSigmaV->GetX()[i] << "   " << gMchiSigmaV->GetY()[i] << endl;
    
      cout << "Done!" << endl;
    }
  return(gMchiSigmaV);
}

double integrateAnnilSpec2(TGraphAsymmErrors* gEffArea,TF1* fDMSpectra,double m_chi,decay_t decay_type)
{
  if( decay_type == bbar ) // for bbar
    {
      double a = -1.46 + 4.26e-2*TMath::Log(m_chi) - 4.37e-3*(pow(TMath::Log(m_chi),2.0));
      //double a = -1.46 + 4.26e-2*TMath::Log10(m_chi) - 4.37e-3*(pow(TMath::Log10(m_chi),2.0));

      fDMSpectra->SetParameter(0,a);
    }
  if( decay_type == tt )
    {
      double b = 6.08 + 0.27*TMath::Log(m_chi);
      double c = -6.9 + 5.2e-4*m_chi - 4.2e-8*(m_chi**2.0) + 1.6e-12*(m_chi**3.0);
      double e = -5.25 + 0.36*TMath::Log(m_chi) - 0.04*pow(TMath::Log(m_chi),2.0);

      fDMSpectra->SetParameter(1,b);
      fDMSpectra->SetParameter(2,c);
      fDMSpectra->SetParameter(4,e);
    }
  if( decay_type == ee || decay_type == uu )
    {
      TF1* fDMSpectra_lep = new TF1("fDMSpectra_lep","([0]/(x*TMath::Pi()))*((1.0 + ((1.0-x)**2.0) + (([1]/[2])**2.0))*TMath::Log((1+sqrt(1-(([1]/[2])**2.0)/(1-x)))/(1-sqrt(1-(([1]/[2])**2.0)/(1-x)))) - 2.0*(1.0 - x)*sqrt(1-(([1]/[2])**2.0)/(1-x)))",0.01,1.0);
      double m_l;
      if(decay_type == ee){ m_l = 5.11e-4; } // mass of electron in GeV
      if(decay_type == uu){ m_l = 0.1057; } // mass of muon in GeV
      fDMSpectra_lep->SetParameter(0,1.0/137.0);
      fDMSpectra_lep->SetParameter(1,m_l);
      fDMSpectra_lep->SetParameter(2,m_chi);
      //TF1* fDMSpectra_lep = new TF1("fDMSpectra_lep",leptonSpectra,0.01,1.0,2);
      //fDMSpectra_lep->SetParameter(0,m_l);
      //fDMSpectra_lep->SetParameter(1,m_chi);
      fDMSpectra = (TF1*)fDMSpectra_lep->Clone();
    }
  //cout << "from int. spec funtion: " << fDMSpectra->Eval(0.5) << endl;
  if(m_chi > 200.0)
    cout << "Integral check for m=" << m_chi << " N(E>200GeV): " << fDMSpectra->Integral(200.0/m_chi,1) << endl;
  double sum = 0;
  double E1_log10_TeV,E1_GeV;
  double E2_log10_TeV,E2_GeV;
  
  double x1,x2,dx;
 
  double EffArea1,EffArea2,DMSpec;
  double SpecPoint;
  int N = gEffArea->GetN();
  for(int i=0; i<N; i++)
    {
      
      gEffArea->GetPoint(i,E1_log10_TeV,EffArea1);
      if(i != N-1)
	gEffArea->GetPoint(i+1,E2_log10_TeV,EffArea2);
      else
	gEffArea->GetPoint(i,E2_log10_TeV,EffArea2);
      EffArea1 *= 1e4; // convert to cm
      EffArea2 *= 1e4; // convert to cm
      if(EffArea1 < 0.0){ EffArea1 = 0.0; }
      if(EffArea2 < 0.0){ EffArea2 = 0.0; }
      E1_GeV = 10**(E1_log10_TeV + 3.0);
      E2_GeV = 10**(E2_log10_TeV + 3.0);
    
      x1= E1_GeV/m_chi;
      x2= E2_GeV/m_chi;
      dx = x2 - x1;
      cout << "x1: " << x1 << " x2: " << x2 << endl;
      cout << "Spec1: " << fDMSpectra->Eval(x1) << " Spec2: " << fDMSpectra->Eval(x2) << endl;
      if(x2 < 1.0)
	sum += 0.5*(EffArea2*fDMSpectra->Eval(x2) + EffArea1*fDMSpectra->Eval(x1))*dx;
      else
	{
	  sum += 0.5*(EffArea2*fDMSpectra->Eval(x2) + EffArea1*fDMSpectra->Eval(x1))*(1 - x1);
	  break;
	}
    }
  return(sum);
}

Double_t leptonSpectra(double* x,double* par)
{
  double xx= x[0];
  //cout << "xx=" << xx << endl;
  double fs = 1.0/137.0;
  double m_l = par[0];
  double m_c = par[1];
  //cout << "m_c=" << m_c << endl;
  double b_l = sqrt(1.0 - ((m_l/m_c)**2.0)*(1 - xx));
  //cout << "b_l=" << b_l << endl;
  Double_t spec = fs/(xx*TMath::Pi())*((1.0 + ((1.0-xx)**2.0) - ((m_l/m_c)**2.0))*TMath::Log((1.0+b_l)/(1.0-b_l)) - 2*(1.0 - xx)*b_l);
  //cout << "spec=" << spec << endl;
  //  if(TMath::IsNan(spec))
  //cout << "m_c=" << m_c << " m_l=" << m_l << " b_l=" << b_l << endl;
  //cout << 1.0 - b_l << endl;
  //cout << "From lepton spec function: " << spec << endl;
  return(spec);

}
double integrateAnnilSpec(TGraphAsymmErrors* gEffArea, TF1* fDMSpectra,double m_chi, decay_t decay_type, double dx = 1e-4, double Emin = 0)
{
  // gEffArea is the averaged effective area curve
  // gEffArea has X-axis units of m^2, Y-axis units of log10(E/1TeV)
  // fDMSpectra is dN/dx , where x = E/m_chi (GeV/GeV).
  // we will say that the mass of the m_chi is in units of GeV.
  // x = E/m_chi (GeV/GeV)
  if( decay_type == bbar ) // for bbar
    {
      double a = -1.46 + 4.26e-2*TMath::Log(m_chi) - 4.37e-3*(pow(TMath::Log(m_chi),2.0));
      fDMSpectra->SetParameter(0,a);
    }

  if( decay_type == tt )
    {
      double b = 6.08 + 0.27*TMath::Log(m_chi);
      double c = -6.9 + 5.2e-4*m_chi - 4.2e-8*(m_chi**2.0) + 1.6e-12*(m_chi**3.0);
      double e = -5.25 + 0.36*TMath::Log(m_chi) - 0.04*pow(TMath::Log(m_chi),2.0);
      
      fDMSpectra->SetParameter(1,b);
      fDMSpectra->SetParameter(2,c);
      fDMSpectra->SetParameter(4,e);
    }
  if( decay_type == ee || decay_type == uu )
    {
      TF1* fDMSpectra_lep = new TF1("fDMSpectra_lep","([0]/(x*TMath::Pi()))*((1.0 + ((1.0-x)**2.0) - (([1]/[2])**2.0))*TMath::Log((1+sqrt(1-(([1]/[2])**2.0)/(1-x)))/(1-sqrt(1-(([1]/[2])**2.0)/(1-x)))) - 2.0*(1.0 - x)*sqrt(1-(([1]/[2])**2.0)/(1-x)))",0.01,1.0);
      double m_l;
      if(decay_type == ee){ m_l = 5.11e-4; } // mass of electron in GeV
      if(decay_type == uu){ m_l = 0.1057; } // mass of muon in GeV
      fDMSpectra_lep->SetParameter(0,1.0/137.0);
      fDMSpectra_lep->SetParameter(1,m_l);
      fDMSpectra_lep->SetParameter(2,m_chi);
      //TF1* fDMSpectra_lep = new TF1("fDMSpectra_lep",leptonSpectra,0.01,1.0,2);
      //fDMSpectra_lep->SetParameter(0,m_l);
      //fDMSpectra_lep->SetParameter(1,m_chi);
      if(decay_type == ee)
	fDMSpectra = (TF1*)fDMSpectra_lep->Clone();
      if(decay_type == uu)
	{
	  double* muonPar = new double[3];
	  muonPar[0] = m_l;
	  muonPar[1] = m_chi;
	  muonPar[2] = 0;
	  fDMSpectra->Delete();
	  //fDMSpectra_lep->Delete();
	  TF1* fDMSpectra = new TF1("fDMSpectra_muon",getMuonSpectra,0.01,1.0,3);
	  fDMSpectra->SetParameters(muonPar);
	}
    }
  if(m_chi > 200.0)
    cout << "Integral check for m=" << m_chi << " N(E>200GeV): " << fDMSpectra->Integral(200.0/m_chi,1) << endl;

  //const int nSteps = 1e8;
  //const int nSteps = 1e6;
  double sum = 0;
  double EffArea;
  double E_GeV,E_log10_TeV;
  double x = 50/m_chi;
  //  double x = 0.01;
  //double x = Emin/m_chi;
  //double dx = (1 - x)/nSteps;
  while(x<1)
    {
      E_GeV = x*m_chi;
      E_log10_TeV = TMath::Log10(E_GeV/1e3);
      EffArea = gEffArea->Eval(E_log10_TeV);
      if(EffArea < 0.0){ EffArea = 0.0; }
      //if((int)(x*100)%10==0 && m_chi < 1e3 )
      //	cout << m_chi << " " << x << " " << fDMSpectra->Eval(x) << " " << EffArea << endl;
      sum += EffArea*fDMSpectra->Eval(x)*dx;
      x += dx;
    }
  if( decay_type != ee && decay_type != uu )
    sum += gEffArea->Eval(TMath::Log10(m_chi/1e3))*fDMSpectra->Eval(1.0)*dx/2;
  sum *= 1e4; // final output in cm^2
  return(sum);


}

double getMuonSpectra(double* x,double* par)
{

  double X = x[0];
  const double m_u = par[0];
  const double m_chi = par[1];
  int n = par[2];
  //cout << n << endl;
  const double r = pow((5.11e-4/m_u),2.0);
  if( n > 2 ){ return(0); }
  
  //cout << "r = " << r << endl;
  // Final state radiation (FSR) part:
  //TF1* fDMSpectra_lep = new TF1("fDMSpectra_lep","([0]/(x*TMath::Pi()))*((1.0 + ((1.0-x)**2.0) + (([1]/[2])**2.0))*TMath::Log((1+sqrt(1-(([1]/[2])**2.0)/(1-x)))/(1-sqrt(1-(([1]/[2])**2.0)/(1-x)))) - 2.0*(1.0 - x)*sqrt(1-(([1]/[2])**2.0)/(1-x)))",0.01,1.0);
  //fDMSpectra_lep->SetParameter(0,1.0/137.0);
  //fDMSpectra_lep->SetParameter(1,m_u); // lepton mass
  //fDMSpectra_lep->SetParameter(2,m_chi); // DM particle mass
  
  double spec_FSR = 0;
  double spec_decay = 0;
  if( n == 1 || n == 0)
    {
      if(X == 1)
	spec_FSR = 0.0;
      else{
	double B_l = sqrt(1.0 - (m_u/m_chi)**2.0/(1-X));
	if(TMath::IsNaN(B_l))
	  spec_FSR = 0.0;
	else
	  spec_FSR = ( (1.0/137.0)/(X*TMath::Pi()) )*(  ( 1.0 + ((1.0-X)**2.0) - ((m_u/m_chi)**2.0) )*TMath::Log( (1+B_l)/(1-B_l) ) - 2.0*(1.0 - X)*B_l  );
      }
      if(TMath::IsNaN(spec_FSR))
	cout << "spec_FSR is NaN at: " << X << endl;
    }
  //cout << X << " " << spec_FSR << endl;
  
  // Decay part:
  if( n == 2 || n == 0)
    {
  // Attempting to do this numerically, because life is pain:
      const double c0 = 2.0*(1.0/137)*( 1.0/(3.0*TMath::Pi()) );
      const double c1 = 3.0;
      const double c2 = -2.0;
      const double c3 = 4.0;
      const double c4 = -2.0;
      const double c5 = -17.0/2;
      const double c6 = 23.0/6;
      const double c7 = -101.0/12;
      const double c8 = 55.0/12;
      const int nSteps = 1e3;
      const double dy = (1.0-X/2)/nSteps;
      //      double y = X/2;  
      double y = X;

      while(y < 1.0)
	{
	  spec_decay += c0*( (1-y)/pow(y,2) )*( (c1 + c2*y + c3*pow(y,2) + c4*pow(y,3))*TMath::Log((1-y)/r) + c5 + c6*y + c7*pow(y,2) + c8*pow(y,3) )*dy;
	  y += dy;
	  if(TMath::IsNaN(spec_decay))
	    cout << X << " " << y << endl;
	}
      if(spec_decay < 0.0)
	spec_decay = 0.0;
    }
  //cout << spec_decay << endl;
  /*
  TF1* fDMSpec_lep_decay = new TF1("fDMSpec_lep_decay","(1.0/x)*( [0]/(3*TMath::Pi()) )*( (1-x)/x )*( (3 + 2*x + 4*(x**2.0) - 2*(x**3.0))*TMath::Log((1-x)/[1]) - 17.0/2 + (23.0/6)*x - (101.0/12)*(x**2.0) + (55.0/12)*(x**3.0) )",0.01,1.0); // dN/dy*1/y
  double* par_decay = new double[2];
  par_decay[0] = 1.0/137.0;
  par_decay[1] = r;
  fDMSpec_lep_decay->SetParameters(par_decay);
  //fDMSpec_lep_decay->SetParameter(1,r);
  double* z1 = new double[10];
  double* z2 = new double[10];
  double spec_decay = 2.0*(fDMSpec_lep_decay->Integral(X,1.0,par_decay,1e-6)); // 2*int(1/y*dN/dy) from x to 1
  //double spec_decay = 2.0*(fDMSpec_lep_decay->IntegralFast(10,z1,z2,X,1.0));
  //  double spec_FSR = fDMSpectra_lep->Eval(X);
  //fDMSpec_lep_decay->Delete();
  //fDMSpectra_lep->Delete();
  */
  //spec_decay = 0;
  // return sum of two:
  if(spec_FSR < 0.0 )
    cout << "Spec FSR is negative at x=" << X << endl;
  if(spec_decay < 0.0 )
    cout << "Spec Decay is negative at x=" << X << endl;
  return(spec_FSR + spec_decay);

}

void plot2DMCuves(string inFile1, string inFile2,double theta2cut1 = 0.03,double theta2cut2 = 0.03, decay_t decay = WW, dwarf_t dwarf = segue_1)
{
  TGraph* g1 = plotDMCurve(inFile1,decay,dwarf,theta2cut1);  
  TGraph* g2 = plotDMCurve(inFile2,decay,dwarf,theta2cut2);

  cout << "Writting output to file..." << endl;
  ofstream ifs("DMLimits.txt");
  ifs << "m_chi   <sigma_v>_1     <sigma_v>_2" << endl;
  ifs << "-----------------------------------" << endl;
  ifs.precision(5);
  for(int i=0; i<g1->GetN(); i++)
    ifs << g1->GetX()[i] << "     " << g1->GetY()[i] << "   " << g2->GetY()[i] << endl;
  ifs.close();


  TLegend* leg = new TLegend(0.7,0.7,0.9,0.9);

  g1->SetLineColor(kBlue);
  g2->SetLineColor(kRed);
  g1->Draw("AL");
  g2->Draw("L");

  leg->AddEntry(g1,"93hrs, Soft Cuts","L");
  leg->AddEntry(g2,"93hrs, Medium Cuts","L");
  

  TMarker* pt = new TMarker();
  pt->SetMarkerStyle(kPlus);
  if(decay == bbar)
    {
      pt->DrawMarker(1625.4,8.90e-24);
    }
  if(decay == tt)
    {
      pt->DrawMarker(382.56,1.92e-24);
    }
  if(decay == WW)
    {
      pt->DrawMarker(1033.54,7.95e-24);
    }
  leg->AddEntry(pt,"Event Display limit (Matthieu)","p");
  leg->Draw("same");

}

void plotMultiDMCurve(string inFile = "~/veritas/segue_1/config/results_ulcombined_gamma3.0_soft2.root",dwarf_t dwarf = segue_1,double theta2Cut = 0.03)
{
 
  TGraph* g[4];
  g[WW] = plotDMCurve(inFile,WW,segue_1,theta2Cut);
  g[ZZ] = plotDMCurve(inFile,ZZ,segue_1,theta2Cut);
  g[bbar] = plotDMCurve(inFile,bbar,segue_1,theta2Cut);
  g[tt] = plotDMCurve(inFile,tt,segue_1,theta2Cut);
 
  g[WW]->SetLineColor(kBlue);
  g[ZZ]->SetLineColor(kRed);
  g[bbar]->SetLineColor(kBlack);
  g[tt]->SetLineColor(kGreen);
    
  g[WW]->SetLineWidth(2);
  g[ZZ]->SetLineWidth(2);
  g[bbar]->SetLineWidth(2);
  g[tt]->SetLineWidth(2);

  g[WW]->SetLineStyle(2);
  g[bbar]->SetLineStyle(6);
  g[tt]->SetLineStyle(9);

  TLegend* leg = new TLegend(0.7,0.7,0.9,0.9);
  leg->AddEntry(g[WW],"#chi#chi -> W^{+}W^{-}","L");
  leg->AddEntry(g[ZZ],"#chi#chi ->Z^{0}Z^{0}","L");
  leg->AddEntry(g[bbar],"#chi#chi ->b #bar{b}","L");
  leg->AddEntry(g[tt],"#chi#chi ->#tau^{+} #tau^{-}","L");

  g[WW]->GetXaxis()->SetTitle("m_{#chi} (GeV)");
  g[WW]->GetYaxis()->SetTitle("<#sigma v>^{95% CL} (cm^{3} s^{-1})");
  		
  TCanvas* c1 = new TCanvas();
  c1->SetLogy(1);
  c1->SetLogx(1);

  g[WW]->SetTitle("95% Confidence DM Exclusion Curves");

  g[WW]->Draw("AL");
  g[ZZ]->Draw("L");
  g[bbar]->Draw("L");
  g[tt]->Draw("L");
  leg->Draw("same");
  g[WW]->GetYaxis()->SetRangeUser(2e-27,2e-22);

  TLine* l = new TLine(150,2.5e-26,2.5e4,2.5e-26);
  l->Draw("same");
  l->SetLineWidth(3);
  l->SetLineColor(kGray);
}

double calcJFactor(dwarf_t dwarf = segue_1, prof_t prof = einasto, double theta_sqr = 0.015,double psf_68 = 0.1)
{
  double CmToPc = 3.24e-19; // 1cm in pc
  double MsolToKg = 1.989e30; // 1 solar mass to kg
  double GevToKg = 1.0/(1.78e-27);
  double KgToGeV = 1.0/GevToKg;
  double d; // distance to dwarf in pc
  double theta_max = sqrt(theta_sqr)*TMath::DegToRad(); // theta cut in rad
  //  double p_s = 3.8; // GeV cm^-3
  double p_s = 4.2;
  double r_s = 150; // pc
  double r_t = 500; // pc
  //double n = 1.0/0.3;
  double n = 3.3333;
  // numbers from Magic collab, 2013
  //n = 0.3;
  //  p_s = 1.1e8;
  //p_s = 4.168;
  // numbers for smith cloud (arXiV:1405.1030v2)
  if(dwarf == smith)
    {
      r_s = 1040.0;
      //p_s = 9.2e6*MsolToKg*KgToGeV*pow(CmToPc*1e3,3.0);
      //p_s = 0.350; // GeV/cm^3
      //n = 1.0/0.17;
      p_s = 0.571; //Nicols et al. 2014
      n = 1.0/0.18;
      d = 1.24e4;
      //if(prof == NFW){ p_s = 1.41; }
      if(prof == einasto)
	{
	  r_s = 850.0;
	  p_s = p_s/4.0;
	}
    }
  if(dwarf == segue_1){ d = 2.3e4; }

  // All distances in cm
  r_s = r_s/CmToPc; 
  r_t = r_t/CmToPc;  
  d = d/CmToPc;
  psf_68 *= TMath::DegToRad(); // PSF in radians
    

 // Einasto profile:
  if(prof == einasto)
    { 
      TF1* fp_r = new TF1("Einastro Profile","[0]*exp(-2*[1]*(pow(x/[2],1/[3]) - 1))",0,d/1000);
      fp_r->SetParameter(0,p_s); 
      fp_r->SetParameter(1,n);
      fp_r->SetParameter(2,r_s);
      fp_r->SetParameter(3,n);
    }
  else if(prof == NFW)
    {
      TF1* fp_r = new TF1("NFW Profile","[0]/((x/[1])*(1.0 + x/[1])**2.0)",0,d/1000);
      fp_r->SetParameter(0,p_s);
      fp_r->SetParameter(1,r_s);

    }
  /*
  TCanvas* c1 = new TCanvas();
  fp_r->Draw();
  c1->SetLogy(1);
  c1->SetLogx(1);
  return(0);
  */
  static const double theta_int_max = 1.2*TMath::DegToRad();
  static const int NumInt = 1e2;

  TF1* fpsf = new TF1("Gaussian PSF","[0]*exp(-0.5*((x-[1])/[2])**2.0)",-theta_int_max,theta_int_max);
  fpsf->SetParameter(2,psf_68);
  fpsf->SetParameter(1,0.0); // just to intialize it
  double norm = 1.0/(psf_68*2*TMath::Pi());
  cout << norm << endl;
  fpsf->SetParameter(0,norm);
  //  fpsf->Draw();
  // return(0);
  //static const double theta_int_max = 1.2*TMath::DegToRad();
  //static const int NumInt = 1e3;
  double theta = 0;
  double dtheta = theta_int_max/NumInt;
  double s,s_min,s_max;
  double ds = dtheta;
  double L[NumInt + 1];
  vector<double> L_p(NumInt+1,0);
  double r;
  int i = 0;
  double J = 0;

  s_min = d*cos(theta_max) - sqrt(pow(r_t,2.0) - pow(d*sin(theta_max),2.0));
  s_max = d*cos(theta_max) + sqrt(pow(r_t,2.0) - pow(d*sin(theta_max),2.0));

  while(theta <= theta_int_max)
    {
      if(i%100 == 0)
	cout << "theta = " << theta*TMath::RadToDeg() << endl;
      
      if( r_t > d*sin(theta) )
	{
	  s_min = d*cos(theta) - sqrt(pow(r_t,2.0) - pow(d*sin(theta),2.0));
	  s_max = d*cos(theta) + sqrt(pow(r_t,2.0) - pow(d*sin(theta),2.0));	 
	}
      else
	{
	  // I'm not 100% sure on this part, but should only be an issue with extended sources
	  s_min = d*cos(theta) - sqrt( -1.0*pow(r_t,2.0) + pow(d*sin(theta),2.0) );
	  s_max = d*cos(theta) + sqrt( -1.0*pow(r_t,2.0) + pow(d*sin(theta),2.0) );
	  cout << "Warning! Tidal radius smaller than d*sin(theta)!" << endl;
	}

      s = s_min;
      L[i] = 0;     
      ds = (s_max - s_min)/NumInt;

      while(s <= s_max)
	{
	  r = sqrt(s**2.0 + d**2.0 - 2*s*d*cos(theta));
	  //J += pow(fp_r->Eval(r),2.0)*ds;	  
	  L[i]   += pow(fp_r->Eval(r),2.0)*ds;
	  L_p[i] += calcJprofConv(fp_r,fpsf,r,theta,s,d)*ds;
	  s += ds;
	}

      // L[i] *= 2*TMath::Pi()*sin(theta); // Luminosity profile at a particular value of theta
      theta += dtheta;
      i++;

    }
  //cout << "starting convolution... ";
  
  //convolveWithEXP(L,0,theta_int_max,NumInt,L_p);
  //convolveWithPSF(L,fpsf,0,1.2,NumInt,L_p);

  //cout << " done!" << endl;
  // Integrate Luminosity profile:
  plotLuminProfile(L,L_p,0,theta_int_max,NumInt);
  double J_nConv = 0;
  theta = 0;
  int j = 0;
  while(theta <= theta_max )
    { 
      
      if(TMath::IsNaN(L_p[j]))
	 cout << "Warning! Nan number in integral: " << j << endl; 
      else
	J += 2*TMath::Pi()*L_p[j]*sin(theta)*dtheta;
      
      J_nConv += 2*TMath::Pi()*L[j]*sin(theta)*dtheta;
      j++;
      theta += dtheta;
    }
  cout << "J factor is: " << J << endl;
  cout << "J without convolution is: " << J_nConv << endl;
  return J_nConv;
}

double calcJprofConv(TF1* fJprof,TF1* fPSF,double r0,double th0,double s,double d,const int nSteps = 1000,const double th_max = 1.2*TMath::DegToRad())
{
  //cout << "Convolving radius = " << r0 << endl;
  double dth = th_max/nSteps;
  double th = 0;
  double sum = 0;
  double r;
  fPSF->SetParameter(1,th0);
  
  while(th <= th_max)
    {
      r = sqrt(pow(s,2.0) + pow(d,2.0) - 2*d*s*cos(th)); 
      sum += pow(fJprof->Eval(r),2)*fPSF->Eval(th)*dth;
      th += dth;		 
    }
  return(sum);


}

void plotLuminProfile(double* L,vector<double> L_p,double x_min,double x_max,const int N)
{
  TGraph* gL = new TGraph(N);
  TGraph* gL_p = new TGraph(N);
  x_min *= TMath::RadToDeg();
  x_max *= TMath::RadToDeg();
  double dx = (x_max - x_min)/N;
  double x = x_min;
  for(int n=0; n<N; n++)
    {
      gL->SetPoint(n,x,2*TMath::Pi()*sin(x*TMath::DegToRad())*L[n]);
      gL_p->SetPoint(n,x,2*TMath::Pi()*sin(x*TMath::DegToRad())*L_p[n]);
      x += dx;
    }
  TCanvas* c1 = new TCanvas();
  TLegend* leg = new TLegend(0.7,0.7,0.9,0.9);
  gL->GetYaxis()->SetRangeUser(1e18,2e22);
  gL->SetLineColor(kRed);  
  gL->Draw("AL");
  gL_p->Draw("L");
  leg->AddEntry(gL,"L(#theta)","L");
  leg->AddEntry(gL_p,"L(#theta)*PSF(#theta)","L");
  leg->Draw("same");


  c1->SetLogy(1);

}

void convolveWithPSF(double* y, TF1* fPSF, double x_min, double x_max, const int N, vector<double> &y_p)
{
  // fPSF is the PSF function, x is the array to be convolved
  double dx = (x_max - x_min)/N;
  double x = x_min;
  double x_0 = x_min;
  //fPSF->SetParameter(0,1.0/0.1);
  //double y_p[N];
  for(int n=0; n<=N; n++)
    {
      y_p[n] = 0.0;
      x = x_min;
      //x *= TMath::DegToRad();
      fPSF->SetParameter(1,x_0);
      fPSF->SetParameter(1,sin(x_0));
      // cout << x_0 << endl;
      for(int m=0; m<=N; m++)
	{
	  //y_p[n] += TMath::TwoPi()*sin(x)*y[m] * fPSF->Eval(sin(x)) * dx;
	  y_p[n] += TMath::TwoPi()*x_0*y[m] * fPSF->Eval(x) * dx;
	  //y_p[n] += y[m] * fPSF->Eval(x) * dx;
	  //y_p[n] += y[m]*exp(-0.5*((x-x_0)/0.1)**2.0)*(x_0 - x)*dx;
	  x += dx;	  
	}
      
      x_0 += dx;
    }
  for(int i=0; i<=N; i++)
    { 
      if(TMath::IsNaN(y_p[i]))
	 cout << "Warning! Nan number: " << i << endl; 
    }

  // return(y_p);
}

void convolveWithEXP(double* y, double x_min, double x_max, const int N, vector<double> &y_p)
{

  double dx = (x_max - x_min)/N;
  double x = x_min;
  double x_0 = x_min;
  double psf = 0.1*TMath::DegToRad();

  for(int n=0; n<=N; n++)
    {
      y_p[n] = 0.0;
      x = x_0;
      for(int m=n; m<=N; m++)
	{
	  y_p[n] += y[m]*exp(-(x_0 - x)/psf)*dx;
	  x += dx;
			     
	}
      x_0 += dx;
    }
  for(int i=0; i<=N; i++)
    { 
      if(TMath::IsNaN(y_p[i]))
	 cout << "Warning! Nan number: " << i << endl; 
    }

}
 
TGraphAsymmErrors* ReadEATextFile(string inFile = "Pass3a/Segue1_SoftCut_pass3a_avgEA_mcE.txt")
{
  cout << "Reading Effective Areas from " << inFile << endl;
 ifstream in(inFile.c_str());
  if(!in.is_open())
    {
      cerr << inFile << "Not Found! " << endl;
      return;
    }
  double runID,LT,E,EA,EA_l,EA_h;
  TGraphAsymmErrors* g = new TGraphAsymmErrors();
  int i=0;
  while(in >> runID >> LT >> E >> EA >> EA_l >> EA_h)
    {
      cout << E << " " << EA << " + " << EA_l << " - " << EA_h << endl;
      if(runID == 0 )
	{
	  g->SetPoint(i,TMath::Log10(E/1000),EA);
	  g->SetPointEYlow(i,-1.0*EA_l);
	  g->SetPointEYhigh(i,EA_h);
	  i++;
	}
    }
  g->GetXaxis()->SetTitle("Energy (log_{10} TeV)");
  g->GetYaxis()->SetTitle("Effective Area (m^{2})");
  return(g);

  
}

TGraph* plotAlexCheckDMCurve(string inFile = "lims.txt")
{
  ifstream in(inFile.c_str());
  if(!in.is_open())
    {
      cerr << "No file!" <<endl;
      return;
    }
  char str[250];
  in.getline(str,250);
  cout << str << endl;
  in.getline(str,250);
  cout << str << endl;
  TGraph* g_tt_mat  = new TGraph();
  TGraph* g_tt_pppc = new TGraph();
  TGraph* g_bb_mat  = new TGraph();
  TGraph* g_bb_pppc = new TGraph();

  double m_chi,sigv_tt_mat,sigv_tt_pppc,sigv_bb_mat,sigv_bb_pppc;
  int i=0;
  
  while( in >> m_chi >> sigv_bb_mat >> sigv_bb_pppc >> sigv_tt_mat >> sigv_tt_pppc )
    {
      g_bb_mat->SetPoint(i,m_chi,sigv_bb_mat);
      g_bb_pppc->SetPoint(i,m_chi,sigv_bb_pppc);
      g_tt_mat->SetPoint(i,m_chi,sigv_tt_mat);
      g_tt_pppc->SetPoint(i,m_chi,sigv_tt_pppc);
      i++;
    }

  TGraph* g_tt_zit = plotVivDMCurve("ed_ea/z20_n200.txt",tt);
  TGraph* g_bb_zit = plotVivDMCurve("ed_ea/z20_n200.txt",bbar);

  TMarker* tMarkMin_tt = new TMarker(383.0,1.9e-24,22);
  TMarker* tMarkMin_bb = new TMarker(1625.0,8.9e-24,22);
  tMarkMin_tt->SetMarkerColor(kGreen+2);
  tMarkMin_bb->SetMarkerColor(kGreen+2);
 
  g_bb_mat->SetLineColor(kBlue);
  g_bb_pppc->SetLineColor(kRed);
  g_tt_mat->SetLineColor(kBlue);
  g_tt_pppc->SetLineColor(kRed);

  g_tt_zit->GetXaxis()->SetTitle("M_{#chi} [GeV]");
  g_tt_zit->GetYaxis()->SetTitle("<#sigma#nu>^{95%CL} [cm^{3}s^{-1}]");
  g_tt_zit->SetTitle("48 hrs. Limit, #tau#tau");

  TCanvas* c1 = new TCanvas("c1","c1",50,50,700,500);
  TLegend* leg = new TLegend(0.6,0.6,0.95,0.95);
  leg->AddEntry(g_tt_zit,"Ben, Matthieu Spec.","L");
  leg->AddEntry(g_tt_mat,"Alex, Matthieu Spec.","L");
  leg->AddEntry(g_tt_pppc,"Alex, PPPC4DM Spec.","L");
  leg->AddEntry(tMarkMin_tt,"Minimum Published value","p");
  c1->SetLogy(1);
  c1->SetLogx(1);
  g_tt_zit->GetYaxis()->SetRangeUser(1e-24,1e-21);
  g_tt_zit->Draw("AL");
  g_tt_mat->Draw("L");
  g_tt_pppc->Draw("L");
  leg->Draw("same");
  tMarkMin_tt->Draw("same");

  g_bb_zit->GetXaxis()->SetTitle("M_{#chi} [GeV]");
  g_bb_zit->GetYaxis()->SetTitle("<#sigma#nu>^{95%CL} [cm^{3}s^{-1}]");
  g_bb_zit->SetTitle("48 hrs. Limit, b#bar{b}");

  TCanvas* c2 = new TCanvas("c2","c2",60,60,700,500);
  c2->SetLogy(1);
  c2->SetLogx(1);
  g_bb_zit->GetYaxis()->SetRangeUser(1e-24,1e-21);
  g_bb_zit->Draw("AL");
  g_bb_mat->Draw("L");
  g_bb_pppc->Draw("L");
  tMarkMin_bb->Draw("same");
  leg->Draw("same");
}

TGraph* plotVivDMCurve(string eaFile,decay_t decay = tt,bool doShift = 1,double index = -2.2,double E_min = 0.3)
{
  if( eaFile.find(".root") != string::npos )
    {
      cout << "EA file is a ROOT file" << endl;
      TFile* f = new TFile(eaFile.c_str(),"READ");
      if(!f->IsOpen())
	{
	  cerr << "File not open!" << endl;
	  return;
	}
      TGraphAsymmErrors* gEffArea = new TGraphAsymmErrors();
      TTree* EffFit = (TTree*)gDirectory->Get("EffFit");
      if(EffFit==NULL)
	{
	  cerr << "No TTree object!" << endl;
	  return;
	}
      EffFit->SetBranchAddress("gEffArea",&gEffArea);
      EffFit->GetEntry(0);
      TF1* fEffFit = f->Get("fEff");
      f->Close();
    }
  elseif( eaFile.find(".txt") != string::npos )
    {
      cout << "EA File is an ascii file" << endl;
      TGraphAsymmErrors* gEffArea = readEDEffFile(eaFile); // m^2
      //TGraph* gEffArea0 = new TGraph("ed_ea/data.csv","%lg, %lg");
      //cout << gEffArea0->GetN();
      //TGraphAsymmErrors* gEffArea = new TGraphAsymmErrors(gEffArea0->GetN(),gEffArea0->GetX(),gEffArea0->GetY());
    }
  //gEffArea->Draw("A*");
  //fEff->Draw("same");
    
  // Do the shift in the EA energy scale (or not). Function also sets any EA < 0 = 0
  shiftEA(gEffArea,doShift);

  double N_ul = 135.9; // counts
  double J = 7.7e18; // GeV^2 cm^-5
  //double T_obs = 47.5*3600.0; // sec
  double T_obs = 171089.4; // sec
  double WWParam[] = {-1.36,0,0,0.5,-12.6,12.1,-9.86};
  double ZZParam[] = {-1.40,0,0,0.47,-14.4,16.3,-13.6};
  double bbarParam[] = {0,0,0,1.05,-17.8,12.3,-1.86};
  TF1* fDMSpectra = new TF1("fDMSpectra","pow(x,[0] + [1]*TMath::Log(x) + [2]*pow(TMath::Log(x),2))*exp([3] + [4]*x + [5]*(x**2.0) + [6]*(x**3.0))",0.01,1);
  
  if(decay == WW){ fDMSpectra->SetParameters(WWParam); }
  if(decay == ZZ){ fDMSpectra->SetParameters(ZZParam); }
  if(decay == bbar){ fDMSpectra->SetParameters(bbarParam); }
  if(decay == tt)
    {
      double ttParam[] = {-1.29,7.83,-2.70,-9e-3,-5.06};
      TF1* fDMSpectra_tau = new TF1("fDMSpectra_tau","pow(x,[0])*([1]*x + [2]*(x**2.0) + [3]*(x**3.0))*exp([4]*x)",1e-4,1);
      fDMSpectra_tau->SetParameters(ttParam);
      fDMSpectra = fDMSpectra_tau;
    }
  if(decay == uu || decay == ee )
    fDMSpectra = new TF1();
  //  checkDMSpecInt(fDMSpectra,decay);

  // Flux ul:
  double E_max = 100.0;
  TF1* fPLSpectra = new TF1("fPLSpectra","[0]*pow(x,[1])",E_min,E_max);
  fPLSpectra->SetParameter(0,1.0);
  fPLSpectra->SetParameter(1,index);
  double A = 0;
  double E_step = E_min;
  double dE = 1e-5;
  //double nSteps0 = (TMath::Log10(E_max) - TMath::Log10(E_min))/dE;
  double nSteps0 = (E_max - E_min)/dE;
  //cout << nSteps0 << endl;
  double intPl = 0;
  /*
  if(doShift)
    A = fPLSpectra->Eval(E_step)*fEffFit->Eval(TMath::Log10(E_step*2.3))*dE/2.0;
  else
    A = fPLSpectra->Eval(E_step)*fEffFit->Eval(TMath::Log10(E_Step))*dE/2.0;
  */
  A = fPLSpectra->Eval(E_step)*gEffArea->Eval(TMath::Log10(E_step))*dE/2.0;
  while(E_step <= E_max)
    {
      E_step += dE;
      //A += fPLSpectra->Eval(10.0**E_step)*gEffArea->Eval(E_step)*dE;
      int i = E_step/dE;
      if( i%10000 == 0 )
	cout << i << " " <<E_step << " " << gEffArea->Eval(TMath::Log10(E_step)) << " " << endl;// << fEffFit->Eval(E_step) << endl;
      A += fPLSpectra->Eval(E_step)*gEffArea->Eval(TMath::Log10(E_step))*dE;
      /*
      if(doShift)
	A += fPLSpectra->Eval(E_step)*fEffFit->Eval(TMath::Log10(E_step*2.3))*dE;
      else
	A += fPLSpectra->Eval(E_step)*fEffFit->Eval(TMath::Log10(E_step))*dE;
      */
      //intPl += fPLSpectra->Eval(E_step)*dE; 
    }
  //A = A/(fPLSpectra->Integral(E_min,E_max));
  //A = A/intPl;
  A = A/(pow(E_max,index+1)/(index+1) - pow(E_min,index+1)/(index+1));
  /*
  TCanvas* c1 = new TCanvas();
  c1->SetLogy(1);
  c1->SetLogx(1);
  fDMSpectra->Draw();
  TCanvas* c2 = new TCanvas();
  c2->SetLogy(1);
  c2->SetLogx(1);
  gEffArea->Draw("A*");
  */
  //A *= 1.0e4.0 // cm^2
  double F_ul = 102.5/(T_obs*A);
  
  cout << "Flux upper limit: " << F_ul << endl;
  //return;
  double m_chi_min = 100;
  double m_chi_max = 1e5;
  //double m_chi_min = 100;
  //double m_chi_max = 2.0e4;
  double m_chi = m_chi_min;
  if( decay == uu )
    int nSteps = 100;
  else
    int nSteps = 1000;

  TGraph* gMchiSigmaV = new TGraph();
  double sigma_v = 0;
  double sigma_v_min = 100;
  double min_m_chi = 0;
  double pi = TMath::Pi();
  double intDMFlux = 0;
  int j=0;
  double Emin_cut = 0;
  for(int i=0; i<nSteps; i++)
    {
      if(i%10 == 0 )
	cout << "m_chi=" << m_chi << " sigma_v="<<sigma_v << endl;
      //intDMFlux = integrateAnnilSpec2(gEffArea,fDMSpectra,m_chi,decay);
      
      //if( m_chi > 1e3 )
      intDMFlux = integrateAnnilSpec(gEffArea,fDMSpectra,m_chi,decay,1e-4,Emin_cut);
	//else
	//intDMFlux = integrateAnnilSpec(gEffArea,fDMSpectra,m_chi,decay,1e-5,Emin_cut);
      
      //if(intDMFlux <= 0.0 ){continue;}
      
      sigma_v = (8*pi/J)*(N_ul*pow(m_chi,2)/(T_obs*intDMFlux));
      //cout << J << " " << N_ul << " " << m_chi << " " intDMFlux << endl;
      //sigma_v_old = s*(8*pi/J_old)*(N_ul_old*pow(m_chi,2)/(T_obs_old*intDMFlux));
      if(!TMath::IsNaN(sigma_v) && sigma_v > 0.0 && intDMFlux > 0.0)
	{
	  //cout << sigma_v << endl;
	  gMchiSigmaV->SetPoint(j,m_chi,sigma_v);
	  // gMchiSigmaVOld->SetPoint(j,m_chi,sigma_v_old);
	  
	  j++;
	  if(sigma_v < sigma_v_min)
	    {
	      sigma_v_min = sigma_v;
	      min_m_chi = m_chi;
	    } 
	}
      m_chi = m_chi*pow(10.0,TMath::Log10(m_chi_max/m_chi_min)/nSteps);
      //m_chi += 100.0;
    }
  cout << "Minimum of curve: " << sigma_v_min << " at mass: " << min_m_chi << endl; 
    
  TCanvas* c1 = new TCanvas("c1","c1",50,50,700,500);
  c1->SetLogy(1);
  c1->SetLogx(1);
  gMchiSigmaV->Draw("AL");	    
  ofstream fout("DMLimsErratum.txt");
  fout << "Mass (GeV)    SigmaV 95% (cm^3s^-1)" << endl;
  fout << "-----------------------------------" << endl;

  for(int i=0; i<gMchiSigmaV->GetN(); i++)
    fout << gMchiSigmaV->GetX()[i] << " " << gMchiSigmaV->GetY()[i] << endl;
  return(gMchiSigmaV);
  
}

void shiftEA(TGraphAsymmErrors* &g,bool doShift = true)
{
  double dx = TMath::Log10(TMath::Log(10.0));
  //double dx = TMath::Log10(exp(1));
  double* x = g->GetX();
  double* y = g->GetY();
  for(int i=0; i<g->GetN(); i++)
    {
      //if(y[i] < 0.0)
	//	y[i] = 0.0;     
      if(doShift)
	g->SetPoint(i,x[i]+dx,y[i]);
      
    }
}
TGraphAsymmErrors* readEDEffFile(string inFile,double s = 1.0)
{

  ifstream in(inFile.c_str());
  if(!in.is_open())
    {
      cerr << "No ED EA file!" << endl;
      return;
    }
  double E,EffA,dEffA;
  int i=0;
  TGraphAsymmErrors* g = new TGraphAsymmErrors();
  while( in >> E >> EffA >> dEffA )
    {
      g->SetPoint(i,E,s*EffA); // m^2
      g->SetPointError(i,0,0,s*dEffA,s*dEffA); // m^2
      i++;
    }
  return(g);
}

