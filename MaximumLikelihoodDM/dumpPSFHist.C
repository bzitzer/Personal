TGraphAsymmErrors* dumpPSFHist(string inFile,string outFile,bool doNorm = false)
{

  VARootIO io(inFile,true);
  io.loadTheRootFile();
  if(!io.IsOpen())
    {
      cout << "No Root file here, bub!" << endl;
      return(NULL);
    }
  vector<double> Etr,Th2;
  TH2D* h = new TH2D();
  hTh2E = getSimDataCombTree(inFile,Etr,Th2);
  //getRatioONRegion(Etr,Th2);
  //return;
  double Elog10first = -2; // units log10(E_TeV)
  double Elog10last = 2;
  const int nBins = 100;
  double dElog10 = (Elog10last - Elog10first)/nBins;
  double E;
  vector<double> Elow,Ehigh,Emid;
  vector<int> NEventsEBin(nBins,0);
  TH1F* hTh2[nBins];
  
  string hId;
  string title;
  ofstream fout(outFile.c_str());
  ofstream fout1("GraphPSFvE.txt");
  ofstream fout2("ExampleTh2Hist.txt");
  ofstream fout3("ThBinEdgeSqrtSpacing.txt");
  ofstream fout4("EBinEdge.txt");
  gStyle->SetOptStat(0);

  // loop over Energy bins and make a TH1F object for each E bin:
  cout << "Seting up histograms...";
  const int nThBins = 500;
  const double thMax = 2.0;
  double thBin[nThBins];
  thBin[0] = 0.0;
  for(int i=1; i<=nThBins; i++)
    thBin[i] = (sqrt(thMax)*(i)/nThBins)**2.0;

  //for(int i=1; i<=nThBins; i++)
  //  cout << thBin[i] - thBin[i-1]<< endl;


  for(int i=0; i<nBins; i++)
    {
      sprintf(hId.c_str(),"h%d",i);
      hTh2[i] = new TH1F(hId.c_str(),"",nThBins,thBin);
      hTh2[i]->GetXaxis()->SetTitle("#theta [deg]");
      //Elow.push_back( pow(10.0,Elog10first + i*dElog10 + 3));
      //Ehigh.push_back( pow(10.0,Elog10first + (i+1)*dElog10 + 3));
      Elow.push_back(Elog10first + i*dElog10);
      Ehigh.push_back(Elog10first + (i+1)*dElog10);
      Emid.push_back(Elog10first + (i+0.5)*dElog10);
      sprintf(title.c_str(),"#theta Distribution between %d GeV and %d GeV",pow(10.0,Elow.at(i)+3),pow(10.0,Ehigh.at(i)+3));
      hTh2[i]->SetTitle(title.c_str());
      cout << pow(10.0,Elow.at(i)+3) << " " << pow(10.0,Ehigh.at(i)+3) << endl;
      fout4 << Elow.at(i) << endl;	    
    }
  fout4 << Elog10last << endl;
  cout << " done!" << endl;
  for(int i=1; i<=hTh2[0]->GetNbinsX(); i++)
    fout3 << hTh2[0]->GetBinLowEdge(i) << endl;

  fout3 << hTh2[0]->GetBinLowEdge(i) << endl;

  // loop over each event in the Stg5 combined tree, bin it in th^2/E space:
  cout << "Filling Histograms...";
  for(int j=0; j<Etr.size(); j++)
    {
      for(int i=0; i<nBins; i++)
	{
	  if( Elow.at(i) < Etr.at(j) && Etr.at(j) <= Ehigh.at(i) )
	    {
	      NEventsEBin.at(i)++;
	      hTh2[i]->Fill(sqrt(Th2.at(j)));
	      continue;
	      //cout << "Got Here!" << endl; 
	    }
	}
    }
  cout << endl;
  double PSF_68,PSF_95;
  double norm,intVal;
  bool foundPSF = false;
  TGraphAsymmErrors* gPSFvE = new TGraphAsymmErrors();
  TF1* fFit = new TF1("fFit","[0]*ROOT::Math::cauchy_pdf(sqrt(x),[1]) + [2]/TMath::CosH(-sqrt(x)/[3])",0,0.1);
  fFit->SetLineColor(kBlue);

  fout << "N(th | E_tr) For MC file: " << inFile << endl;
  fout << "--------------------------------------------------------------------" << endl;
  fout.width(5);
  fout.precision(5);
  fout << scientific;
  double sum,dsum,P,dP = 0;
  bool foundPSF_low,foundPSF_high = false;
  double PSF_68_low,PSF_68_high = 0;
  double N,dOmega;
  for(int i=0; i<nBins; i++)
    {
      if(NEventsEBin.at(i) != 0 && doNorm)
       	hTh2[i]->Scale(1.0/NEventsEBin.at(i));
      
      writeHistogram2(hTh2[i],fout);
      //norm = hTh2[i]->GetSumOfWeights();
      norm = NEventsEBin.at(i);
      sum = 0;
      dsum = 0;
      foundPSF = false;
      foundPSF_high = false;
      foundPSF_low = false;
      for(int j=0; j<hTh2[i]->GetNbinsX(); j++)
	{
	  //	  dOmega = TMath::TwoPi()*(TMath::Cos(TMath::DegToRad()*(hTh2[i]->GetBinLowEdge(j)+hTh2[i]->GetBinWidth(j))) - TMath::Cos(TMath::DegToRad()*hTh2[i]->GetBinWidth(j)));
	  //N = hTh2[i]->GetBinContent(j);
	  //if(dOmega != 0)
	  //  hTh2[i]->SetBinContent(j,N/dOmega);
	  sum += hTh2[i]->GetBinContent(j);
	  dsum += sqrt(hTh2[i]->GetBinContent(j));
	  if(doNorm)
	    {
	      P = sum/norm;
	      dP = dsum/norm;
	    }
	  else
	    {
	      P = sum;
	      dP = dsum;
	    }
	  if(P >= 0.68 & !foundPSF)
	    {
	      PSF_68 = hTh2[i]->GetBinCenter(j); 
	      foundPSF = true;
	    }
	  if(P-dP >= 0.68 & !foundPSF_high)
	    {
	      PSF_68_high = hTh2[i]->GetBinCenter(j); 
	      foundPSF_high = true;
	    }
	  if(P+dP >= 0.68 & !foundPSF_low)
	    {
	      PSF_68_low = hTh2[i]->GetBinCenter(j); 
	      foundPSF_low = true;	      
	    }
	}
      gPSFvE->SetPoint(i,Emid.at(i),PSF_68);
      gPSFvE->SetPointEYlow(i,TMath::Abs(PSF_68-PSF_68_low));
      gPSFvE->SetPointEYhigh(i,TMath::Abs(PSF_68_high-PSF_68));
      
      cout << "68% PSF at: " << PSF_68 << " - " << PSF_68_low << " + " << PSF_68_high << endl;
    }
  gPSFvE->GetXaxis()->SetTitle("Energy [log_{10} TeV]");
  gPSFvE->GetYaxis()->SetTitle("PSF (68% cont.) [deg]");
  gPSFvE->SetTitle("PSF vs. Energy");
  
  gPSFvE->Draw("A*");
  //ifstream out1("output1.txt");
  //  writePSFGraph(gPSFvE,fout1);
  
  return(gPSFvE);
}

void writeHistogram2(TH1F* f,const ostream &out = cout)
{
  out.precision(7);
  for(int i=1; i<=f->GetNbinsX()+1; i++)
    out << f->GetBinContent(i) << " ";

  out << endl;
}

void writePSFGraph(TGraphAsymmErrors* g,const ostream &out = cout)
{
  out << "Energy (GeV)    PSF (68%, deg)   +PSF Err (deg)   -PSF Err (deg)";
  out << "----------------------------------------------------------------";

  out.width(5);
  out.precision(5);
  out << scientific;

  double E,PSF,dPSFHigh,dPSFLow;
  for(int i=0; i<g->GetN(); i++)
    {
      g->GetPoint(i,E,PSF);
      E = pow(10.0 , E + 3.0);
      dPSFHigh = g->GetErrorYhigh(i);
      dPSFLow = g->GetErrorYlow(i);
      out << E << " " << PSF << " " << dPSFHigh << " " << -1.0*dPSFLow << endl;
    }
}

void writeHistogram(TH1F* h,TF1* fFit,const ostream &out = cout)
{
  out << "P(th^2 | E_tr) = " << fFit->GetExpFormula() << endl;
  out << "E_low (GeV) E_high (GeV)   [0]    [1]    [2]    [3]    chi^2     NDF" << endl; 
  out << "--------------------------------------------------------------------" << endl;
  out.width(5);
  out.precision(5);
  out << scientific;
  double* p = fFit->GetParameters();
  double* dp = fFit->GetParErrors();
  
  for(int i=0; i<fFit->GetNpar(); i++)
    {
      out << p[i] << " ";
      out << " +/- " << dp[i] << " ";
      out << endl;
    }
  out << "---------------------------------------------------------------------" << endl;
  out << "Bin Center (deg^2)  Number Counts" << endl;
  out << "---------------------------------------------------------------------" << endl;
  for(int i=0; i<h->GetNbinsX(); i++)
    {
      out << h->GetBinCenter(i) << " ";
      out << h->GetBinContent(i) << " ";
      out << endl;
    }
}

void writeFitFun(TF1* f1,double Elow,double Ehigh,double N,const ostream &out = cout)
{
  bool doNorm = true;
  //normalizePol2Exp(f1);
  //out << "P(th^{2} | E) " <<
  Elow = pow(10.0,Elow+3.0); // units of GeV;
  Ehigh = pow(10.0,Ehigh+3.0); // units of GeV;

  out << Elow << " " <<  Ehigh << " ";
  double* p = f1->GetParameters();
  double* dp = f1->GetParErrors();
  
  // double N = p[3]*(p[0] + p[1]*p[3] + 2*p[2]*(p[3]**2.0));
  N = f1->Integral(0,10.0); // integrate to effectivly inf
  for(int i=0; i<f1->GetNpar(); i++)
    {
      if(i<3 && N != 0.0 )
	{
	  p[i] = p[i]/N;
	  dp[i] = dp[i]/N;
	}  
      out << p[i];
      out << "+/-" << dp[i] << " ";
    }
  out << f1->GetChisquare() << " " << f1->GetNDF() << " " << endl; 
}
void normalizePol2Exp(TF1* &f)
{
  // Do not use. 
  if( f == NULL ){ return(NULL); }
  double p[] = f->GetParameters();
  double dp[] = f->GetParErrors();
  double N = p[3]*(p[0] + p[1]*p[3] + 2*p[2]*(p[3]**2.0));
  for(int i=0; i<3; i++)
    {
      f->SetParameter(i,p[i]/N);
      f->SetParError(i,dp[i]/N);
    }

}
void getRatioONRegion(vector<double> Etr,vector<double> Th2,double radius = 0.17)
{
  ofstream fout("OnRatio.txt"); 
  TH1D* hEnergyOn  = new TH1D("hEnergyOn","hEnergyOn",100,-2,2);
  TH1D* hEnergyAll = new TH1D("hEnergyAll","hEnergyAll",100,-2,2);
  if(Etr.size() != Th2.size())
    {
      cerr << "Problem! Vector size mis-match!" << endl;
      return;
    }
  cout << "Filling 1D Histograms...";
  for(int i=0; i<Th2.size(); i++)
    {
      hEnergyAll->Fill(Etr.at(i));
      if(Th2.at(i) < radius**2.0)
	hEnergyOn->Fill(Etr.at(i));

      if(Th2.at(i) > 10000)
	cout << "Theta^2 = " << Th2.at(i) << " Energy= " << Etr.at(i) << endl;
    }
  cout << " done!" << endl;
  for(int i=0; i<=hEnergyAll->GetNbinsX(); i++)
    {
      fout << hEnergyAll->GetBinLowEdge(i) << " ";
      fout << hEnergyAll->GetBinLowEdge(i) + hEnergyAll->GetBinWidth(i) << " ";
      fout << hEnergyOn->GetBinContent(i) << " ";
      fout << hEnergyAll->GetBinContent(i) << " ";
      // fout << hEnergyAll->GetBinContent(i)/hEnergyOn->GetBinContent(i) << " ";
      fout << endl;
    }
}

TH2D* getSimDataCombTree(string inFile,vector<double> &Etr,vector<double> &Th2)
{
  double Th2Max = 0.1;
  const int nBinsTh2 = 100;
  double Elog10first = -1.2; // units log10(E_TeV)
  double Elog10last = 1.8;
  const int nBinsE = 10;
  double Elog10Tmp,Th2Tmp;

  Etr.clear();
  Th2.clear();
  VASimulationData* pSim = new VASimulationData(); 
  VAShowerData* pShower = new VAShowerData();
  TH2D* hTh2E = new TH2D("hTh2E","Theta^2 E",nBinsTh2,0,Th2Max,nBinsE,Elog10first,Elog10last);

  VARootIO io(inFile,true);
  io.loadTheRootFile();
  if(!io.IsOpen())
    {
      cout << "No Root file here, bub!" << endl;
      return(NULL);
    }
 
  TTree* t = (TTree*)io.loadAnObject("CombinedEventsTree","SelectedEvents",true);
  if( t == NULL )
    {
      cout << "No Combined Tree! What gives!?" << endl;
      return(NULL);
    }

  t->SetBranchAddress("S",&pShower);
  t->SetBranchAddress("Sim",&pSim);
  int numCutEvents = 0;
  int numLargeTh2Events = 0;
  double theta2Max = 0;
  for(int i=0; i<t->GetEntries(); i++)
    {
      t->GetEntry(i);
      if(i%100000 == 0)
	cout << "On Event: " << i << " of " << t->GetEntries() << endl;
      
      if(pSim->fEnergyGeV == 0.0)
	{
	  numCutEvents++;
	  continue;
	}
      /*      
      if( pShower->fMSL >= 1.3 || pShower->fMSW >= 1.1 || pShower->fMSL <=0.05 || pShower->fMSW <= 0.05 || pShower->fShowerMaxHeight_KM <= 7.0)
	{ 
	  numCutEvents++;
	  continue;
	}
      */
      //cout << "MSW: " << pShower->fMSW << " MSL: " << pShower->fMSL << endl;
      Elog10Tmp = TMath::Log10(pSim->fEnergyGeV/1000);
      //Elog10Tmp = TMath::Log10(pShower->fEnergy_GeV/1000);
      Th2Tmp = pShower->fTheta2_Deg2;
      if( Th2Tmp > theta2Max ){ theta2Max = Th2Tmp; }
      if( Th2Tmp > 4.0 ){ numLargeTh2Events++; }
      //if( Elog10Tmp == 0.0 || Th2Tmp != 0.0 )
      //	{
	  //if( Th2Tmp > theta2Max ){ theta2Max = Th2Tmp; }
	  //  if( Th2Tmp > 0.1 ){ numLargeTh2Events++; }    
      Etr.push_back(Elog10Tmp);
      Th2.push_back(Th2Tmp);
      hTh2E->Fill(Th2Tmp,Elog10Tmp);
      //	}
      //else
      //	{
      //	  numCutEvents++;
	  //	  cout << "True Energy: " << pSim->fEnergyGeV << endl;
	  //cout << "Rec Energy: " << pShower->fEnergy_GeV << endl;
	  // cout << "Th^2: " << pShower->fTheta2_Deg2 << endl;

      //}
    }
  cout << "Number of passing events" << Etr.size() << endl;
  cout << "Number of cut events: " << numCutEvents << endl;
  cout << "Number of large Th2 Events: " << numLargeTh2Events << endl;
  cout << "Largest Th2 value: " << theta2Max << endl;
  return(hTh2E);
}
