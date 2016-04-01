
enum decay_t {WW, ZZ, bbar, tt, ee, uu};
enum somm_t {somm0, somm1, somm2};

TGraph* plotDMCurve_wSomm1(string inFile)
{

  TFile* f = new TFile(inFile.c_str(),"READ");
  if(!f->IsOpen())
    {
      cout << "Missing file!" << endl;
      return;
    }
  //TGraph* gSigV0 = (TGraph*)gDirectory->Get("g_WW");
  TGraph* gSigV0 = (TGraph*)gDirectory->Get("g_ee");
  if(gSigV0 == NULL)
    {
      cout << "Missing TGraph Object!" << endl;
      return;
    }
  TGraph* gSigVsomm1 = new TGraph();
  double sigv0,m_chi,S;
  for(int i=0; i<gSigV0->GetN(); i++)
    {
      gSigV0->GetPoint(i,m_chi,sigv0);
      S = calcSommBoost1MaxwellInt(m_chi);
      gSigVsomm1->SetPoint(i,m_chi,sigv0/S);
    }
  gSigVsomm1->SetLineStyle(2);
  gSigVsomm1->Draw("AL");
  gSigV0->Draw("L");
}

void plotSommBoost1()
{
  double B[] = {0.1,0.01,0.001,1e-4,1e-5};
  double m_chi_min = 1000;
  double m_chi_max = 1e5;
  TGraph* g[5];
  const int nSteps = 1000;
  double m_chi = m_chi_min;
  double S;
  ostringstream os[5];
  for(int i=0; i<nSteps; i++)
    {
      
      for(int j=0; j<5; j++)
	{
	  if(i==0)
	    g[j] = new TGraph();

	  S = calcSommBoost1(m_chi,B[j]*TMath::C());
	  if(S > 10)
	    cout << "Large Sommerfield at Beta " << B[j] << " and mass: " << m_chi << endl;
	  g[j]->SetPoint(i,m_chi,S);
	}
      m_chi = m_chi*pow(10.0,TMath::Log10(m_chi_max/m_chi_min)/nSteps);
    }
  TCanvas* c1 = new TCanvas("c1","c1",50,50,700,600);
  g[0]->GetXaxis()->SetLimits(m_chi_min,m_chi_max);
  g[0]->GetYaxis()->SetRangeUser(0.9,1.1e6);
  c1->SetLogy(1);
  c1->SetLogx(1);
  //g[0]->Draw("AL");
  TLegend* l = new TLegend(0.7,0.7,0.9,0.9);
  for(int i=0; i<5; i++)
    {
      os[i] << "B = " << B[i] << endl;
      g[i]->SetLineColor(i+2);
      l->AddEntry(g[i],os[i].str().c_str(),"l");
      if(i==0)
	g[0]->Draw("AL");
      else
	g[i]->Draw("L");
    }
  l->Draw("same");
}

double calcSommBoost1(double m_chi,double v)
{
  double S = 1;
  //double m_v = 90.0; // W mass in GeV;
  double m_v = 0.250;
  double R = 1.0/m_v; // Yukawa radius;
  double B = v/TMath::C(); // v/c for mediator particle
  double a = 1.0/30.0; // alpha
  //B *= 2/sqrt(TMath::Pi());
  double K_out = m_chi*B;
  double K_in = sqrt(K_out**2.0 + a*m_chi*m_v);
  //double K_in = sqrt(a*m_chi*m_v);
  double S_v = 1;
  // integrating times maxwell dist: 
  //S = ((TMath::Cos(K_in*R))**2.0 + ((K_out/K_in)*TMath::Sin(K_in*R))**2.0)**-1.0;
  double e_v = B/a;
  double e_phi = m_v/(a*m_chi);
  double pi = TMath::Pi();
  S = (pi/e_v)*TMath::SinH(12.0*e_v/(pi*e_phi))/(TMath::CosH(12.0*e_v/(pi*e_phi)) - TMath::Cos(2*pi*sqrt(6.0/(pow(pi,2)*e_phi) - e_v**2.0/((pi**2.0)*e_phi/6.0)**2.0)) );
  //  cout << m_chi << " " << S << endl;
  return(S);
}

double calcSommBoost1MaxwellInt(double m_chi, double m_v = 90.0, double a = 1.0/30.0)
{
  //double S = 1;
  double R = 1.0/m_v;
  double v0 = 6.4*1e3;
  double b0 = v0/TMath::C();
  double kT0 = (1.0/3.0)*m_v*(b0**2.0);
  double vMax = 150.0*1e3;
  int numSteps = 1e5;
  double dv = vMax/numSteps;
  double v = 0;
  double f_v;
  double k = 4*TMath::Pi()*pow(m_v/(TMath::TwoPi()*kT0),3.0/2.0);
  double S = 0;
  while(v<=vMax)
    {     
      f_v = k*(v/TMath::C())**2.0*TMath::Exp(-1.0*m_v*(v/TMath::C())**2.0/(2.0*kT0));
      v += dv;
      S += f_v*calcSommBoost1(m_chi,v)*dv/TMath::C();
    }
  cout << m_chi << " " << S << endl;
  return(S);
}
