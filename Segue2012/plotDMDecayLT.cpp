
#include <TGraph.h>
#include <TGraphAsymmErrors.h>
#include <TFile.h>
#include <TTree.h>
#include <TF1.h>
#include <TMath.h>
#include <TCanvas.h>
#include <math.h>

enum decay_t {WW, ZZ, bbar, tt, ee, uu, We, Wu};
enum dwarf_t {segue_1, draco, ursa_minor, wilman_1, bootes};
enum prof_t {einasto, NFW};

double integrateDecaySpec(TGraphAsymmErrors* gEffArea,double m_chi, decay_t decay, double dx = 1e-4, double Emin = 0);

TF1* getDMSpectra(double m_chi,decay_t decay);

double getMuonSpectra(double* x,double* par);
double getWlSpectra(double* x,double* par);

void shiftEA(TGraphAsymmErrors* &g,bool doShift = true);

using namespace std;

TGraph* plotErratumDMDecayLT(string eaFile,decay_t decay,bool doShift = true)
{
  TGraphAsymmErrors* gEffArea = new TGraphAsymmErrors();
  TGraph* gDecayLT = new TGraph();
  if( eaFile.find(".root") != string::npos )
    {
      cout << "EA file is a ROOT file" << endl;
      TFile* f = new TFile(eaFile.c_str(),"READ");
      if(!f->IsOpen())
        {
          cerr << "File not open!" << endl;
          return(NULL);
        }
      
      TTree* EffFit = (TTree*)gDirectory->Get("EffFit");
      if(EffFit==NULL)
        {
          cerr << "No TTree object!" << endl;
          return(NULL);
        }
      EffFit->SetBranchAddress("gEffArea",&gEffArea);
      EffFit->GetEntry(0);
      TF1* fEffFit = (TF1*)f->Get("fEff");
      f->Close();
    }
  else if( eaFile.find(".txt") != string::npos )
    {
      cout << "EA File is an ascii file" << endl;
      //TGraphAsymmErrors* gEffArea = readEDEffFile(eaFile); // m^2
      //TGraph* gEffArea0 = new TGraph("ed_ea/data.csv","%lg, %lg");
      //cout << gEffArea0->GetN();
      //TGraphAsymmErrors* gEffArea = new TGraphAsymmErrors(gEffArea0->GetN(),gEffArea0->GetX(),gEffArea0->GetY());
    }


  shiftEA(gEffArea,doShift);

  double N_ul = 135.9; // counts
  double J = 3e17; // GeV^2 cm^-5
  //double T_obs = 47.5*3600.0; // sec
  double T_obs = 171089.4; // sec
  double m_chi_min = 150;
  double m_chi_max = 2e4;
  double m_chi = m_chi_min;
  const int nSteps = 100;
  double intDMFlux,decayLT;
  //int j=0;
  
  for(int i=0; i<nSteps; i++)
    {
      intDMFlux = integrateDecaySpec(gEffArea,m_chi,decay);
      decayLT = ( J/(4*TMath::Pi()) )*(T_obs*intDMFlux/(N_ul*m_chi));
      gDecayLT->SetPoint(i,m_chi,decayLT);
      if(i%10 == 0)
	cout << m_chi << " " << decayLT << endl;
      m_chi = m_chi*pow(10.0,TMath::Log10(m_chi_max/m_chi_min)/nSteps);
    }
  TCanvas* c1 = new TCanvas();
  c1->SetLogy(1);
  c1->SetLogx(1);
  gDecayLT->Draw("AL");

  return(gDecayLT);
}

double integrateDecaySpec(TGraphAsymmErrors* gEffArea,double m_chi, decay_t decay, double dx, double Emin)
{
  // gEffArea is the averaged effective area curve
  // gEffArea has X-axis units of m^2, Y-axis units of log10(E/1TeV)
  // fDMSpectra is dN/dx , where x = E/m_chi (GeV/GeV).
  // we will say that the mass of the m_chi is in units of GeV.
  // x = 2*E/m_chi (GeV/GeV)
  //m_chi = m_chi/2;

 TF1* fDMSpectra = getDMSpectra(m_chi,decay);
 double E_GeV,E_log10_TeV,EffArea;
 double sum = 0;
 double x = 2*50.0/m_chi;
 if(m_chi > 200.0)
    cout << "Integral check for m=" << m_chi << " N(E>200GeV): " << fDMSpectra->Integral(200.0/m_chi,1) << endl;

  while(x<1)
    {
      E_GeV = 0.5*x*m_chi;
      E_log10_TeV = TMath::Log10(E_GeV/1e3);
      EffArea = gEffArea->Eval(E_log10_TeV);
      if(EffArea < 0.0){ EffArea = 0.0; }
      //if((int)(x*100)%10==0 && m_chi < 1e3 )
      //        cout << m_chi << " " << x << " " << fDMSpectra->Eval(x) << " " << EffArea << endl;
      sum += EffArea*fDMSpectra->Eval(x)*dx;
      x += dx;
    }
  //if( decay != ee && decay != uu )
  //  sum += gEffArea->Eval(TMath::Log10(m_chi/1e3))*fDMSpectra->Eval(2.0)*dx/2;
  sum *= 1e4; // final output in cm^2
  return(sum);
}

TF1* getDMSpectra(double m_chi,decay_t decay)
{
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
 if( decay == bbar ) // for bbar
    {
      double a = -1.46 + 4.26e-2*TMath::Log(m_chi) - 4.37e-3*(pow(TMath::Log(m_chi),2.0));
      fDMSpectra->SetParameter(0,a);
    }

  if( decay == tt )
    {
      double b = 6.08 + 0.27*TMath::Log(m_chi);
      double c = -6.9 + 5.2e-4*m_chi - 4.2e-8*pow(m_chi,2.0) + 1.6e-12*pow(m_chi,3.0);
      double e = -5.25 + 0.36*TMath::Log(m_chi) - 0.04*pow(TMath::Log(m_chi),2.0);

      fDMSpectra->SetParameter(1,b);
      fDMSpectra->SetParameter(2,c);
      fDMSpectra->SetParameter(4,e);
    }
  if( decay == ee )
    {
      TF1* fDMSpectra_lep = new TF1("fDMSpectra_lep",
	 "( [0]/(x*TMath::Pi()) )*( (1.0 + pow(1.0-x,2.0) - pow([1]/[2],2.0))*TMath::Log( (1+sqrt(1-pow([1]/[2],2.0)/(1-x)))/(1-sqrt(1-(pow([1]/[2],2.0)/(1-x))) )) - 2.0*(1.0 - x)*sqrt(1-pow([1]/[2],2.0)/(1-x)))"
	 ,0.01,1.0);
      double m_l = 5.11e-4; // mass of electron in GeV
      //      if(decay == uu){ m_l = 0.1057; } // mass of muon in GeV
      fDMSpectra_lep->SetParameter(0,1.0/137.0);
      fDMSpectra_lep->SetParameter(1,m_l);
      fDMSpectra_lep->SetParameter(2,m_chi);
      //TF1* fDMSpectra_lep = new TF1("fDMSpectra_lep",leptonSpectra,0.01,1.0,2);
      //fDMSpectra_lep->SetParameter(0,m_l);
      //fDMSpectra_lep->SetParameter(1,m_chi);
      fDMSpectra = fDMSpectra_lep;
    }      
  if(decay == uu)
    {
      double* muonPar = new double[3];
      muonPar[0] = 0.1057;
      muonPar[1] = m_chi;
      muonPar[2] = 0;
      
      TF1* fDMSpectra_muon = new TF1("fDMSpectra_muon",getMuonSpectra,0.01,1.0,3);
      fDMSpectra_muon->SetParameters(muonPar);
      fDMSpectra = fDMSpectra_muon;
    }
  if(decay == We || decay == Wu)
    {
      TF1* fDMSpectra_Wl = new TF1("fDMSpectra_Wl",getWlSpectra,0.01,1.0,2);
      double* wlPar = new double[2];
      wlPar[0] = m_chi;
      if( decay == We )
	wlPar[1] = 5.11e-4;
      fDMSpectra_Wl->SetParameters(wlPar);
      fDMSpectra = fDMSpectra_Wl;
    }

  return(fDMSpectra);
}   

double getWlSpectra(double* x,double* par)
{
  double X = x[0];
  double m_chi = par[0];
  double m_l = par[1];
  const double M_W = 80.0;
  double E_W_eff = (pow(m_chi,2.0) + pow(M_W,2))/(2*m_chi);
  double E_e_eff = (pow(m_chi,2.0) - pow(M_W,2))/(2*m_chi);
  double E = 0.5*m_chi*X;
  double x_W = E/E_W_eff;
  double x_e = E/E_e_eff;
  //TF1* fW;
  double spec_tot,spec_W,spec_e = 0;
  
  spec_W = getDMSpectra(m_chi,WW)->Eval(x_W);
  if( m_l == 5.11e-4 )
    spec_e = getDMSpectra(m_chi,ee)->Eval(x_e);

  spec_tot = 0.5*spec_W*m_chi/E_W_eff + 0.5*spec_e*m_chi/E_e_eff;
  return(spec_tot);
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
        double B_l = sqrt(1.0 - pow(m_u/m_chi,2.0)/(1-X) );
        if(TMath::IsNaN(B_l))
          spec_FSR = 0.0;
        else
          spec_FSR = ( (1.0/137.0)/(X*TMath::Pi()) )*(  ( 1.0 + pow(1.0-X,2.0) - pow(m_u/m_chi,2.0) )*TMath::Log( (1+B_l)/(1-B_l) ) - 2.0*(1.0 - X)*B_l  );
      }
      if(TMath::IsNaN(spec_FSR))
        cout << "spec_FSR is NaN at: " << X << endl;
    }
  //cout << X << " " << spec_FSR << endl;

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

  //double spec_decay = 2.0*(fDMSpec_lep_decay->Integral(X,1.0,par_decay,1e-6)); // 2*int(1/y*dN/dy) from x to 1
  //double spec_decay = 2.0*(fDMSpec_lep_decay->IntegralFast(10,z1,z2,X,1.0));
  //  double spec_FSR = fDMSpectra_lep->Eval(X);
  //fDMSpec_lep_decay->Delete();
  //fDMSpectra_lep->Delete();
  
  //spec_decay = 0;
  // return sum of two:
  if(spec_FSR < 0.0 )
    cout << "Spec FSR is negative at x=" << X << endl;
  if(spec_decay < 0.0 )
    cout << "Spec Decay is negative at x=" << X << endl;
  return(spec_FSR + spec_decay);

}

void shiftEA(TGraphAsymmErrors* &g,bool doShift)
{
  double dx = TMath::Log10(TMath::Log(10.0));
  //double dx = TMath::Log10(exp(1));
  double* x = g->GetX();
  double* y = g->GetY();
  for(int i=0; i<g->GetN(); i++)
    {
      //if(y[i] < 0.0)
        //      y[i] = 0.0;     
      if(doShift)
        g->SetPoint(i,x[i]+dx,y[i]);
      
    }
}
