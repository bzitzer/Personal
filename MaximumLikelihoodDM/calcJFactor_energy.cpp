
#include <vector>
#include <sstream>
#include <iostream>
#include <fstream>

#include <TH1D.h>
#include <TH2D.h>
#include <TGraph.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TAxis.h>
#include <TF1.h>
#include <TF2.h>
#include <TLegend.h>
#include <TLine.h>
#include <TFile.h>

enum dwarf_t {segue_1,draco,ursa_minor,bootes,umaII};

double CmToPc   = 3.24e-19; // 1cm in pc
double MsolToKg = 1.989e30; // 1 solar mass to kg
double GeVToKg  = 1.78e-27;
double KgToGeV  = 1.0/GeVToKg;

using namespace std;

TGraph* getAlexJFactor1D(const double thCut = 0.17);

TF1* calcAlexDMProfile();

void readBinFile(vector<double> &x,string inFile,const int nBins = 100);

double* readBinFile2(string inFile,int &N);

TGraph* calcJProfile(dwarf_t dwarf = umaII, const double thCut = 0.17, const int nBinsTh = 81);

TGraph* getPSF1D(TH2D* h,double en = 500);

TH2D* calcConvJProf2D(TGraph* gJProf,TGraph* gPSF,double &J,const double thCut = 0.17);

TH2D* readPSFFile(string inFile = "/raid/reedbuck/bzitzer/Pass5f/PSF/No200_Zn30_ua_thHist_psf_overflow.txt",
		  bool doDraw = 1,const int nBinsX = 100,const int nBinsY = 500);

TH2D* readAlexJFactor(string inFile = "JPSFtables/integratedJconvPSF.txt",
		      string enBinFile = "JPSFtables/Eknots.txt",string thBinFile = "JPSFtables/thetaknots.txt");

TGraph* getRadialConvJProf(TH2D* h);

double calcJFactor_energy(string outFile, const double thCut = 0.17)
{

  TGraph* gJProf = calcJProfile(segue_1,thCut);
  TH2D* hTh_energy = readPSFFile("/raid/reedbuck/bzitzer/Pass5f/PSF/No200_Zn30_ua_thHist_psf_overflow.txt",0);

  vector<double> en;
  vector<double> J;
  readBinFile(en,"JPSFtables/Eknots.txt");

  TGraph* gPSF1DProf = getPSF1D(hTh_energy,1e4);
  double tmp;
  //calcConvJProf2D(gJProf,gPSF1DProf,tmp,thCut);
  //return(0);
  TGraph* convJvE = new TGraph();
  TH2D* hJProf2D[en.size()];
  TGraph* gRadJProf[en.size()];
  TFile* f = new TFile(outFile.c_str(),"RECREATE");
  ostringstream os[en.size()];
  ostringstream os1[en.size()];
  ostringstream os2[en.size()];

  TGraph* gAlex1DProf = getAlexJFactor1D();

  for(int i=0; i<en.size(); i++)
    {
      J.push_back(0.0);
      TGraph* gPSF1DProf = getPSF1D(hTh_energy,en.at(i));
      os[i] << "hE_" << (int)en.at(i);
      os1[i] << "gPSF_" << (int)en.at(i);
      os2[i] << "gRadJProf_" << (int)en.at(i);
      cout << "Working on energy: " << en.at(i) << endl;

      hJProf2D[i] = (TH2D*)calcConvJProf2D(gJProf,gPSF1DProf,J.at(i),thCut);
      
      convJvE->SetPoint(i,en.at(i),J.at(i));

      gRadJProf[i] = getRadialConvJProf(hJProf2D[i]);

      hJProf2D[i]->SetTitle(os[i].str().c_str());
      hJProf2D[i]->Write(os[i].str().c_str());
      gPSF1DProf->Write(os1[i].str().c_str());
      gRadJProf[i]->Write(os2[i].str().c_str());
    }
  gAlex1DProf->SetMarkerColor(kBlue);
  gAlex1DProf->SetMarkerStyle(22);
  gAlex1DProf->Write("gAlex1DProf");
  convJvE->Write("gConvJvE");
  f->Close();
}

TGraph* getRadialConvJProf(TH2D* h)
{
  const int N = 100;
  TH1D* hRadJProf = new TH1D("hRadJProf","hRadJProf",N,0,2);
  TH1D* hRadN     = new TH1D("hRadJProf","hRadJProf",N,0,2);
  TGraph* gRadJProf = new TGraph();
  double x,y,r;
  double J,n;
  int k = 0;
  for(int i=0; i<h->GetNbinsX(); i++)
    for(int j=0; j<h->GetNbinsY(); j++)
      {
	x = h->GetXaxis()->GetBinCenter(i);
	y = h->GetYaxis()->GetBinCenter(j);
	r = sqrt(x*x + y*y);
	hRadJProf->Fill(r,h->GetBinContent(i,j));
	hRadN->Fill(r,1.0);
      }
  for(int i=0; i<N; i++)
    {
      n = hRadN->GetBinContent(i);
      r = hRadJProf->GetBinCenter(i);
      J = hRadJProf->GetBinContent(i);
      if(n == 0) { continue;}
      gRadJProf->SetPoint(k,r,J/n);
      k++;
    }
  return(gRadJProf);
}

TGraph* getAlexJFactor1D(const double thCut)
{
  TH2D* h = readAlexJFactor();
  TGraph* g = new TGraph();
  
  int nBinY = h->GetYaxis()->FindBin(thCut);
  double E,J;
  for(int i=1; i<=h->GetNbinsX(); i++)
    {
      E = h->GetXaxis()->GetBinCenter(i);
      J = pow(10.0,h->Interpolate(E,thCut));
      g->SetPoint(i-1,E,J);
      cout << "Energy: " << E << " Int. J Factor: " << J << endl;
    }
  g->GetXaxis()->SetTitle("Energy [GeV]");
  g->GetYaxis()->SetTitle("J [GeV^{2}cm^{-5}]");
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
  double* th = readBinFile2(thBinFile,thN);
  double* en = readBinFile2(enBinFile,enN);
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

TH2D* calcConvJProf2D(TGraph* gJProf,TGraph* gPSF,double &J,const double thCut)
{

  bool foundThetaMax = false;
  double theta_max;
  const int NPSF   = gPSF->GetN();
  const int NJProf = gJProf->GetN();
  const int NBins = 100;

  double x_min = -2.0;
  double x_max = -x_min;
  double y_min = x_min;
  double y_max = -x_min;
  
  double* th = gJProf->GetX();
  for(int i=0; i<NPSF; i++)
    {
      if( gJProf->Eval(th[i]) <= 0.0 && !foundThetaMax )
	{
	  foundThetaMax = true;
	  theta_max = th[i];
	}
   }
  theta_max = 1.0;
  cout << "Maximum theta: " << theta_max << endl;

  TH2D* hProfConv  = new TH2D("hProfConv","J*PSF Profile",NBins,-theta_max,theta_max,NBins,-theta_max,theta_max);
  TH2D* hProf2D    = new TH2D("hProf2D","J Profile",NBins,-theta_max,theta_max,NBins,-theta_max,theta_max);
  TH2D* hProfPSF2D = new TH2D("hProfPSF","PSF Profile",NBins,-theta_max,theta_max,NBins,-theta_max,theta_max);
  TH2D* hProfConv_l  = new TH2D("hProfConv","J*PSF Profile",NBins,-theta_max,theta_max,NBins,-theta_max,theta_max);
  hProfConv->SetStats(0);
  hProf2D->SetStats(0);
  hProfConv_l->SetStats(0);
  hProfPSF2D->SetStats(0);
  const int NBinsX = hProf2D->GetXaxis()->GetNbins();
  const int NBinsY = hProf2D->GetYaxis()->GetNbins();
  double x,y,x0,y0,dx,dy,L_r,F_r,r,r0;
  double rJ,rP;
  int k = 0;
  double r_smooth = 0.5;
  int n_min,n_max,m_min,m_max;
  for(int i=0; i<=NBinsX; i++)
    for(int j=0; j<=NBinsY; j++)
      {
	k++;
	if( k%1000 == 0 )
	  cout << "On bin: " << k << " of: " << NBinsX*NBinsY << endl;

	if(k==1)
	  {
	    for(int i=0; i<=NBinsX; i++)
	      for(int j=0; j<=NBinsY; j++)
		{
		  x0 = hProf2D->GetXaxis()->GetBinCenter(i);
		  y0 = hProf2D->GetYaxis()->GetBinCenter(j);
		  
		  dx = hProf2D->GetXaxis()->GetBinWidth(i);
		  dy = hProf2D->GetYaxis()->GetBinWidth(j);

		  r = sqrt(x0*x0 + y0*y0);
		  L_r = gJProf->Eval(r)/(TMath::TwoPi()*sin(r*TMath::DegToRad()));
		  if( L_r > 0.0 )
		    hProf2D->SetBinContent(i,j,TMath::Log10(L_r));
		  hProfPSF2D->SetBinContent(i,j,gPSF->Eval(r));
		}
	  }
	
	x0 = hProfConv->GetXaxis()->GetBinCenter(i);
	y0 = hProfConv->GetYaxis()->GetBinCenter(j);
	r0 = sqrt(x0*x0 + y0*y0);
	if(r0 > theta_max){ continue; }
	x_min = x0 - r_smooth;
	y_min = y0 - r_smooth;
	x_max = x0 + r_smooth;
	y_max = y0 + r_smooth;
	
	n_min = hProfConv->GetXaxis()->FindBin(x_min);
	m_min = hProfConv->GetYaxis()->FindBin(y_min);
	n_max = hProfConv->GetXaxis()->FindBin(x_max);
	m_max = hProfConv->GetYaxis()->FindBin(y_max);

	for(int n=n_min; n<=n_max; n++)
	  {
	    
	    dx = hProfConv->GetXaxis()->GetBinWidth(n);
	    x = hProfConv->GetXaxis()->GetBinCenter(n);

	    for(int m=m_min; m<=m_max; m++)
	      {	

		dy = hProfConv->GetYaxis()->GetBinWidth(m);
		y = hProfConv->GetYaxis()->GetBinCenter(m);

		rJ = sqrt(x*x + y*y);
		rP = sqrt(pow(x-x0,2.0) + pow(y-y0,2.0));
		L_r = gJProf->Eval(rJ)/(TMath::TwoPi()*TMath::Sin(rJ*TMath::DegToRad()));

		F_r = gPSF->Eval(rP)*L_r*dx*dy;
		
		hProfConv->Fill(x0,y0,F_r);	
	      }
	  }
      }

  J = 0.0;
  double J_nConv = 0.0;
  double PSF_norm = 0.0;
  for(int i=0; i<NBinsX; i++)
    for(int j=0; j<NBinsY; j++)
      {
	x0 = hProfConv->GetXaxis()->GetBinCenter(i);
	y0 = hProfConv->GetYaxis()->GetBinCenter(j);
	
	dx = hProfConv->GetXaxis()->GetBinWidth(i)*TMath::DegToRad(); // working in degrees, need final J factors in terms of sr
	dy = hProfConv->GetYaxis()->GetBinWidth(j)*TMath::DegToRad();
	
	r = sqrt(x0*x0 + y0*y0);
	if(r<thCut)
	  {
	    J       += hProfConv->GetBinContent(i,j)*dx*dy;
	    J_nConv += pow(10.0,hProf2D->GetBinContent(i,j))*dx*dy;
	  }
	PSF_norm += gPSF->Eval(r)*dx*dy*pow(TMath::RadToDeg(),2.0);
	if(hProfConv->GetBinContent(i,j) >= 0.0 )
	  hProfConv_l->SetBinContent(i,j,TMath::Log10(hProfConv->GetBinContent(i,j)));
      }
  cout << "J Factor (w/ convolution): " << J << endl;
  cout << "J Factor (no convolution): " << J_nConv << endl;
  cout << "PSF integration: " << PSF_norm << endl;
  /*
  double c[] = {19.0,20.0,21,22,23,24,25};
  hProf2D->SetMinimum(19.0);
  hProf2D->SetMaximum(25.0);
  hProfConv_l->SetMinimum(19.0);
  hProfConv_l->SetMaximum(25.0);
  hProf2D->SetContour(7,c);
  hProfConv_l->SetContour(7,c);
  TCanvas* c1 = new TCanvas("c1","c1",50,50,1200,400);
  c1->Divide(3,1);
  c1->cd(1);
  hProfPSF2D->Draw("colz");
  //gPSF->Draw("A*");
  c1->cd(2);
  hProf2D->Draw("colz");
  c1->cd(3);
  hProfConv_l->Draw("colz");
  */
  return(hProfConv);
  
}

TGraph* getPSF1D(TH2D* h,double en)
{
  // Adding this because of low statistics at the low/high ends.
  //h->RebinX(4);

  int nBinEn = h->GetXaxis()->FindBin(TMath::Log10(en/1000.0));
  TGraph* gTh1D    = new TGraph();
  TGraph* gTh1D_nm = new TGraph();
  
  double th,N,dth;
  double norm = 0;
  double Ntot = 0;
  double P = 0;
  for(int i=1; i<=h->GetNbinsY(); i++)
    {
      th  = h->GetYaxis()->GetBinCenter(i);
      dth = h->GetYaxis()->GetBinWidth(i);
      //if( th > 0.01 )
      //	N = h->Interpolate(TMath::Log10(en/1000.0),th);
      //      else
      N = h->Interpolate(TMath::Log10(en/1000.0),th); 
      //      cout << "theta: " << th << " N: " << N << endl;
      //    norm += N*th;
      norm += N*dth;
      Ntot += N; 
      gTh1D->SetPoint(i,th,N/th);
    }
  norm *= TMath::TwoPi();
  
  for(int i=0; i<gTh1D->GetN(); i++)
    {
      gTh1D->GetPoint(i,th,P);
      gTh1D_nm->SetPoint(i,th,P/norm);
    }
  TH2D* hCheckNorm = new TH2D("hCheckNorm","hCheckNorm",1000,-2,2,1000,-2,2);
  hCheckNorm->SetStats(0);
  double x,y,r,dx,dy,psf,sum;
  sum = 0;
  for(int i=0; i<hCheckNorm->GetNbinsX(); i++)
    for(int j=0; j<hCheckNorm->GetNbinsX(); j++)
      {
	x = hCheckNorm->GetXaxis()->GetBinCenter(i);
	y = hCheckNorm->GetYaxis()->GetBinCenter(j);
	dx = hCheckNorm->GetXaxis()->GetBinWidth(i);
	dy = hCheckNorm->GetYaxis()->GetBinWidth(j);
	r = sqrt(x*x + y*y);
	psf = gTh1D_nm->Eval(r);
	hCheckNorm->SetBinContent(i,j,psf);
	sum += psf*dx*dy;
      }
  //hCheckNorm->Draw("colz");
  cout << " PSF Integral from +/- 2 deg: " << sum << endl;
  //hCheckNorm->Delete();
  if( Ntot < 1000.0 )
    cout << "Warning! Low number of events! " << endl;
  
  cout << "Getting PSF with energy of: " << en << endl;
  cout << "PSF normalization of: " << hCheckNorm->Interpolate(0,0) << endl;
  return(gTh1D_nm);
}

TGraph* calcJProfile(dwarf_t dwarf, const double thCut, const int nBinsTh)
{
  // Get the DM profile from 2014 G.S,Koupishappias, and Walker paper
  TF1* fDM_prof = (TF1*)calcAlexDMProfile()->Clone();
  //TF1* fDM_prof = (TF1*)getDMProfile(dwarf)->Clone();
  double p_s = fDM_prof->GetParameter(0);
  p_s *= MsolToKg*KgToGeV*pow(CmToPc,3.0);
  cout << "KgToGeV " << KgToGeV << endl; 
  cout << "p_s " << p_s << endl;
  fDM_prof->SetParameter(0,p_s);

  // Get the VERITAS PSF bin edges:
  vector<double> th; 

  readBinFile(th,"JPSFtables/thetaknots.txt",nBinsTh);

  for(int i=0; i<nBinsTh; i++)
    th[i] *= TMath::DegToRad();
  double th_max = th[nBinsTh-1]; //radians
  double th_min = th[0]; // radians
  cout << "Theta^2 range between " << th_min*TMath::RadToDeg() << " and " << th_max*TMath::RadToDeg() << endl;


  // Need distances: Taken from McConnachie 2012.
  double d;
  if(dwarf == segue_1){ d = 2.3e4; } // pc }
  if(dwarf == draco){ d = 7.6e4; }
  if(dwarf == ursa_minor){ d = 7.6e4; }
  if(dwarf == bootes){ d = 6.6e4; }
  if(dwarf == umaII){ d = 3.2e4; }
  
  double r_t = fDM_prof->GetXmax();
  cout << "Truncation radius: " << r_t << endl;
  double s_min = d*cos(th_max) - sqrt(pow(r_t,2.0) - pow(d*sin(th_max),2.0));
  double s_max = d*cos(th_max) + sqrt(pow(r_t,2.0) - pow(d*sin(th_max),2.0));
  double L[nBinsTh];
  double L_p[nBinsTh];

  const int nInt = 1000;

  TGraph* gL = new TGraph();
  
  double r,s,ds;
  int j = 0;

  // Line of Sight integral:  
  for(int i=0; i<nBinsTh; i++)
    {
      if(i%50 == 0)
	cout << "th[i] = " << th[i] << endl;
      
      if( r_t > d*sin(th[i]) )
	{
	  s_min = d*cos(th[i]) - sqrt(pow(r_t,2.0) - pow(d*sin(th[i]),2.0)); // Distance in pc
	  s_max = d*cos(th[i]) + sqrt(pow(r_t,2.0) - pow(d*sin(th[i]),2.0));	 
	}
      else
	{
	  // I'm not 100% sure on this part, but should only be an issue with extended sources
	  s_min = d*cos(th[i]) - sqrt( -1.0*pow(r_t,2.0) + pow(d*sin(th[i]),2.0) );
	  s_max = d*cos(th[i]) + sqrt( -1.0*pow(r_t,2.0) + pow(d*sin(th[i]),2.0) );
	  cout << "Warning! Tidal radius smaller than d*sin(th[i])!" << endl;
	  //s_min = 0;
	  //s_max = 0;
	}

      s = s_min;
      L[i] = 0;     
      ds = (s_max - s_min)/nInt;

      while(s < s_max)
	{
	  r = sqrt(s*s + d*d - 2*s*d*cos(th[i])); // pc
	  //J += pow(fp_r->Eval(r),2.0)*ds;	  
	  L[i]   += pow(fDM_prof->Eval(r),2.0)*(ds/CmToPc);
	  //L_p[i] += calcJprofConv(fp_r,fpsf,r,th[i],s,d)*ds;
	  s += ds;
	}
      
      L[i] *= 2*TMath::Pi()*sin(th[i]); // Luminosity profile at a particular value of th[i]
       if(!TMath::IsNaN(L[i]))
	 {
	   gL->SetPoint(j,th[i]*TMath::RadToDeg(),L[i]);
	   j++;  
	 }
      //L[i] = TMath::TwoPi()*L[i]
       
    }
 
  return(gL);
}


TF1* calcAlexDMProfile()
{
  /*
  double rhos  = 1.77958937145;
  double rs    = 310.75663656;
  double alpha = 0.544110936127;
  double beta  = 4.35536608401;
  double gamma = 0.642320059857;
  double rt    = 138.58;
  */
  /*
  double rhos  = pow(10.0,-1.0561);
  double rs    = pow(10.0,3.2576);
  double alpha = 1.9126;
  double beta  = 6.3932;
  double gamma = 0.76629;
  double rt    = 138.58;
  */
  /*  
  double rhos  = pow(10.0,-1.1331);
  double rs    = pow(10.0,3.6317);
  double alpha = 1.8606;
  double beta  = 6.3682;
  double gamma = 0.5828;
  double rt    = 294.0;
  */
  double rhos  = 0.187; 
  double rs    = 280.0;
  double alpha = 1.0;
  double beta  = 3.0;
  double gamma = 1.0;
  double rt    = 280.0;

  TF1* fNFW_prof = new TF1("Generalized NFW Segue",
			"[0]/(pow(x/[1],[4])*pow(1 + pow(x/[1],[2]),([3] - [4])/[2]))",
			0,rt);
  fNFW_prof->SetParameter(0,rhos);
  fNFW_prof->SetParameter(1,rs);
  fNFW_prof->SetParameter(2,alpha);
  fNFW_prof->SetParameter(3,beta);
  fNFW_prof->SetParameter(4,gamma);

  return(fNFW_prof);


}

double* readBinFile2(string inFile,int &N)
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
  in.getline(str,300);

  while( in >> tmp )
    x2.push_back(tmp);
    
  N = x2.size();
  double* x = new double[N];
  for(int i=0; i<x2.size(); i++)
    {
      x[i] = x2.at(i);
      //cout << x[i] << endl;
    }
  return(x);
}

void readBinFile(vector<double> &x,string inFile,const int nBins)
{
  //double x[nBins] = {0};
  ifstream in(inFile.c_str());
  vector<double> x2;
  if(!in.is_open())
    {
      cout << "Missing file! All is woe!" << endl;
    }
  
  char str[300];
  in.getline(str,300);

  double tmp;
  int i=0;
  while( in >> tmp )
    {
      x.push_back(tmp);
      i++;
    }

}

TH2D* readPSFFile(string inFile, bool doDraw,const int nBinsX,const int nBinsY)
{
  ifstream in(inFile.c_str());
  char str[300];

  in.getline(str,300);
  cout << str << endl;
  in.getline(str,300);
  cout << str << endl;

  TH2D* h = new TH2D("h","EnergyThHist",nBinsX,-2,2,nBinsY,0,2);
  h->GetXaxis()->SetTitle("Energy [log_{10} TeV]");
  h->GetYaxis()->SetTitle("Theta [deg]");
  int i = 1;
  int j = 1;
  double val;
  while( in >> val )
    {

      if( i > nBinsY )
	{
	  i = 0;
	  j++;
	}
      h->SetBinContent(j,i,val);
      i++;
    }

  //h->RebinX(4);
  //h->RebinY(5);
  if(doDraw)
    h->Draw("colz");

  return(h);
}

