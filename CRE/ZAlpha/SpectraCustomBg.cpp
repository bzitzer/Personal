 void SpectraExternBg(string inFile,string externBgFileList,string outFile = "spectraExternBgOut.root",string configPath = "~/vegas/configfiles/vegas2_5/Soft_hFit")
{
  gROOT->ProcessLine(".L dumpEATextInput.C");
  gROOT->ProcessLine(".L lima.C");
  
  const int rebin = 4;
  const int nBins = 100;
  const double emin = 0.30;
  const double emax = 20.0; 
  bool useExternAlpha = false;
  
  TH1D* hNon = new TH1D("hNon","hNon",nBins/rebin,-2,2);
  TH1D* hNoff = new TH1D("hNoff","hNoff",nBins/rebin,-2,2);
  TH1D* hAlpha = new TH1D("hAlpha","hAlpha",nBins/rebin,-2,2);
  TH1D* hDiff = new TH1D("hExcess","hExcess",nBins/rebin,-2,2);
  TH1D* hSig = new TH1D("hSig","hSig",nBins/rebin,-2,2);

  TH1D* hDiff_dEA = new TH1D("hDiff_dEA","hDiff_dEA",nBins,-2,2);
  TH1D* hNoff_dEA = new TH1D("hNoff_dEA","hNoff_dEA",nBins,-2,2);
  TH1D* hTmp = new TH1D();
  TH1D* hFlux = new TH1D("hFlux","hFlux",nBins,-2,2);
  
  hNon->Sumw2();
  hNoff->Sumw2();
  hFlux->Sumw2();
  
  string str;
  ifstream in(externBgFileList.c_str());
  if(!in.is_open())
    {
      cout << "Missing Bg File list!" << endl;
      return;
    }
  
  TH1D* hNon_dEA = getFluxHist(inFile,configPath,hNon,hNoff,hAlpha,true,useExternAlpha,nBins);
  while( in >> str )
    {
      hTmp = getFluxHist(str,configPath,hNon,hNoff,hAlpha,false,useExternAlpha,nBins);
      hNoff_dEA->Add(hTmp);
    }

  if(hNon_dEA == NULL)
    {
      cerr << "Missing On Flux hist! " << endl;
      return;
    }
  if(hNoff_dEA == NULL)
    {
      cerr << "Missing Off Flux hist! " << endl;
      return;
    }
  
  TCanvas* c1 = new TCanvas();
  hNon_dEA->SetLineColor(kRed);
  hNoff_dEA->SetLineColor(kBlue);
  hNon_dEA->Draw();
  hNoff_dEA->Draw("same");
   // return;

  double Non,Noff,alpha;
   // Excess and significance 
  for(int i=0; i<nBins/rebin; i++)
    {
       Non = hNon->GetBinContent(i);
       Noff = hNoff->GetBinContent(i);
       if(Noff < 0.5){ continue; }
       alpha = hAlpha->GetBinContent(i)/Noff;
       hAlpha->SetBinContent(i,alpha);
       hDiff->SetBinContent(i,Non - alpha*Noff);
       hSig->SetBinContent(i,lima(Non,Noff,alpha));
       cout << "Energy: " << pow(10.0,hNon->GetBinCenter(i))  <<" Non: " << Non << " Noff: " << Noff << " alpha: " << alpha << " Excess: " << hDiff->GetBinContent(i) << " Significance: " << hSig->GetBinContent(i) << endl;
    }
  //return;

  double E,dE,dE_log;
  double Fon,Foff;
  
  for(int i=0; i<nBins; i++)
     {
       E = hNon_dEA->GetBinCenter(i);
       dE_log = hNon_dEA->GetBinWidth(i);
       dE = pow(10.0,E+dE_log/2) - pow(10.0,E-dE_log/2);
       
       Fon  = hNon_dEA->GetBinContent(i)/dE;
       Foff = hNoff_dEA->GetBinContent(i)/dE;
       hFlux->SetBinContent(i,Fon - Foff);
     }
   TCanvas* c2 = new TCanvas("c2","c2",50,50,600,500);
   c2->SetLogy(1);
   c2->SetLogx(1);
   hFlux->Scale(1.0/rebin);
   hFlux->Rebin(rebin);
   //hFlux->Draw();

   TGraphErrors* gFlux = new TGraphErrors();
   const double Nex_min = 5.0;
   const double Nsig_min = 2.0;
   double Nex,sig,err;
   int nFluxPoint = 0;
   for(int i=0; i<hFlux->GetNbinsX(); i++)
     {
       sig = hSig->GetBinContent(i);
       Nex = hDiff->GetBinContent(i);
       if( sig < Nsig_min || Nex < Nex_min ){ continue; }
       
       Non   = hNon->GetBinContent(i);
       Noff  = hNoff->GetBinContent(i);
       alpha = hAlpha->GetBinContent(i);

       err = sqrt(Non + pow(alpha,2)*Noff);
       //err = hFlux->GetBinContent(i)/err;
       err = hFlux->GetBinContent(i)*err/(Nex+DBL_EPSILON);
       E = pow(10,hFlux->GetBinCenter(i));
       gFlux->SetPoint(nFluxPoint,E,hFlux->GetBinContent(i));
       gFlux->SetPointError(nFluxPoint,0,err);
       nFluxPoint++;
       cout << "Energy: " << E << " Flux: " << hFlux->GetBinContent(i) << " +/- " << err << endl;
     }
   TF1* fSpecFit = new TF1("fSpecFit","[0]*pow(x,-[1])",emin,emax);
   //TF1* fSpecFit = new TF1("fSpecFit","[0]*pow(x,-([1]+[2]*TMath::Log10(x)))",emin,emax);
   fSpecFit->SetParameter(0,3.2e-7);
   fSpecFit->SetParameter(1,2.5);
   fSpecFit->SetParameter(2,0.0);
   //gFlux->Draw("A*");
   gFlux->GetXaxis()->SetTitle("Energy [TeV]");
   gFlux->GetYaxis()->SetTitle("Flux [m^{-2}s^{-1}sr^{-1}TeV^{-1}]");
   gFlux->Fit(fSpecFit,"","",emin,emax);
   cout << "Chi^2: " << fSpecFit->GetChisquare() << " NDF: " << fSpecFit->GetNDF() << " Prob: " << fSpecFit->GetProb() << endl;
   
   LaffertyWyattEnergies(gFlux,fSpecFit,hFlux->GetBinWidth(1),0,emin,emax);
   cout << "Chi^2: " << fSpecFit->GetChisquare() << " NDF: " << fSpecFit->GetNDF() << " Prob: " << fSpecFit->GetProb() << endl;
   double intFlux = fSpecFit->Integral(emin,emax);
   double intFluxErr = fSpecFit->IntegralError(emin,emax);
   cout << "Integral flux between " << emin << " and " << emax << " TeV: " << intFlux << " +/- " << intFluxErr << endl;

   gFlux->Draw("A*");
   //return;
   
   TGraphErrors* gResid = new TGraphErrors();
   double resid;
   for(int i=0; i<gFlux->GetN(); i++)
     {
       E = gFlux->GetX()[i];
       resid = (gFlux->GetY()[i] - fSpecFit->Eval(E))/fSpecFit->Eval(E);
       err = gFlux->GetErrorY(i)/fSpecFit->Eval(E);
       gResid->SetPoint(i,E,resid);
       gResid->SetPointError(i,0,err);
     }
   TCanvas* c3 = new TCanvas("c3","c3",50,50,700,500);
   c3->SetLogx(1);
   gResid->Draw("A*");
   gResid->GetXaxis()->SetTitle("Energy [TeV]");
   gResid->GetYaxis()->SetTitle("(Data - Fit)/Fit");
   TLine* l = new TLine(emin,0.0,emax,0);
   l->SetLineStyle(2);
   l->SetLineColor(kRed);
   l->Draw("same");

   TFile* fOut = new TFile(outFile.c_str(),"RECREATE");
   hNon->Write();
   hNoff->Write();
   hNon_dEA->Write();
   hNoff_dEA->Write();
   hFlux->Write();
   gFlux->Write("gSpectrumGraph");
   gResid->Write("gResiduals");
   //for(int i=0; i<nRuns; i++)
   //  gEA[i]->Write(os[i].str().c_str());
   hSig->Write();
   hDiff->Write();
   fOut->Close();

}

double getExternAlpha()
{
  // Code here for Glenn and Mary to use for z-alpha tables...
  return(0.1); // placeholder for now... 
}

TH1D* getFluxHist(string inFile,string configPath,TH1D* &hNon,TH1D* &hNoff,TH1D* &hAlpha,bool isOnRun = true,bool useExternAlpha = false,const int nBins = 100)
{

  TFile* f = new TFile(inFile.c_str(),"READ");
  if(!f->IsOpen() )
    {
      cerr << "Problem opening ROOT file!" << endl;
      return(NULL);
    }
  TTree* EventTree = (TTree*)gDirectory->Get("EventStatsTree");
  if( EventTree == NULL )
    {
      cout << "No Event Tree!" << endl;
      return(NULL);
    }
  TTree* RunTree = (TTree*)gDirectory->Get("RunStatsTree");
  if( RunTree == NULL )
    {
      cout << "No Run Tree!" << endl;
      return(NULL);
    }
  
  bool IsOn,IsOff;
  float Az,El,EnergyGeV;
  double liveTime,W,Noise;
  UInt_t RunID = 0;
  int RunNumRunTree = 0;
  double Omega = TMath::TwoPi()*(1.0 - TMath::Cos(0.17*TMath::DegToRad()));
  
  EventTree->SetBranchAddress("OnEvent",&IsOn);  
  EventTree->SetBranchAddress("OffEvent",&IsOff);    
  EventTree->SetBranchAddress("Weight",&W);

  EventTree->SetBranchAddress("Azimuth",&Az);
  EventTree->SetBranchAddress("Elevation",&El);
  EventTree->SetBranchAddress("Noise",&Noise);
  EventTree->SetBranchAddress("EnergyGeV",&EnergyGeV);

  RunTree->SetBranchAddress("faLiveTime",&liveTime);
  RunTree->SetBranchAddress("faRunNumber",&RunNumRunTree);
  
  EventTree->SetBranchAddress("RunNum",&RunID);

  const int nEvents = EventTree->GetEntries();
  const int nRuns   = RunTree->GetEntries();
  
  vector<TGraphAsymmErrors*> gEA;
  int j = 0;
  int numCountsInRun = 0;
  double Az_avg,El_avg,No_avg;
  double Lt_tot = 0;
  double EA;
  
  TH1D* hNon_dEA = new TH1D("hNon_dEA","hNon_dEA",nBins,-2,2);
  TH1D* hNoff_dEA = new TH1D("hNoff_dEA","hNoff_dEA",nBins,-2,2);
  
  hNon_dEA->Sumw2();
  hNoff_dEA->Sumw2();
  
  vector<double> onEnergy;
  vector<double> offEnergy;
  ostringstream os[nRuns];
    
  RunTree->GetEntry(0);
  
  for(int i=0; i<nEvents; i++)
    {
      EventTree->GetEntry(i);
      if( i%10000 == 0 )
	{ cout << "On Event: " << i << " of: " << nEvents << endl; }

      if( RunNumRunTree != RunID && RunNumRunTree != 0 )
	{
	  Az_avg = Az_avg/(numCountsInRun + DBL_EPSILON);
	  El_avg = El_avg/(numCountsInRun + DBL_EPSILON);
	  No_avg = No_avg/(numCountsInRun + DBL_EPSILON);

	  gEA.push_back(NULL);
	  cout << "Az_avg: " << Az_avg << " El_avg: " << El_avg << endl;
	  gEA[j] = getEAInterpolatedCurve(RunNumRunTree,El_avg,Az_avg,No_avg,
					  0,configPath);
	  os[j] << "gEA_" << RunNumRunTree;
	  gEA[j]->SetTitle(os[j].str().c_str());
	  
	  numCountsInRun = 0;
	  j++;
	  RunTree->GetEntry(j);
	  Lt_tot += liveTime;
	}

      
      if(useExternAlpha && !isOnRun)
	W = getExternAlpha();
      
      if(IsOn && isOnRun)
	{
	  hNon->Fill(TMath::Log10(EnergyGeV/1000.0),1.0);
	}
      if(IsOff && !isOnRun)
	{
	  hNoff->Fill(TMath::Log10(EnergyGeV/1000.0),1.0);
	  hAlpha->Fill(TMath::Log10(EnergyGeV/1000.0),W);
	}
      Az_avg += Az*TMath::RadToDeg();
      El_avg += El*TMath::RadToDeg();
      No_avg += Noise;
      
      numCountsInRun++;
    }
   Az_avg = Az_avg/(numCountsInRun + DBL_EPSILON);
   El_avg = El_avg/(numCountsInRun + DBL_EPSILON);
   No_avg = No_avg/(numCountsInRun + DBL_EPSILON);
   gEA.push_back(NULL);
   gEA[j] = getEAInterpolatedCurve(RunNumRunTree,El_avg,Az_avg,No_avg,
				   0,configPath);
   Lt_tot += liveTime;
   cout << "Number of runs: " << nRuns << " Number of EA curves: " << gEA.size() << endl;
   cout << "Starting second loop: " << endl;

   j = 0;
   // gEA[j]->Draw("A*");
   //return;
   RunNumRunTree = 0;
   
   for(int i=0; i<nEvents; i++)
     {
       EventTree->GetEntry(i);
       if( i%10000 == 0 )
	 { cout << "On Event: " << i << " of: " << nEvents << endl; }

       if( RunNumRunTree != RunID )
	 {
	   cout << RunNumRunTree << " " << RunID << endl;
	   RunTree->GetEntry(j);
	   j++;
	 }
       if(IsOn || IsOff)
	 {
	   // cout << j << endl;
	   EA = gEA.at(j-1)->Eval(TMath::Log10(EnergyGeV/1000.0));
	   //cout << EA << endl;
	   //if(EA < 0.0)
	   //  EA = 0.0;
	 }

      if(useExternAlpha)
	W = getExternAlpha();  
       
       if(IsOn && EA > 0.0)
	 hNon_dEA->Fill(TMath::Log10(EnergyGeV/1000.0),1.0/EA);
       if(IsOff && EA > 0.0)
	 hNoff_dEA->Fill(TMath::Log10(EnergyGeV/1000.0),W/EA);
     }

   cout << "live time " << Lt_tot << endl;
   hNon_dEA->Scale(1.0/Lt_tot);
   hNoff_dEA->Scale(1.0/Lt_tot);
   
   if(isOnRun)
     return(hNon_dEA);
   else
     return(hNoff_dEA);

}


void LaffertyWyattEnergies(TGraphErrors* &g,TF1* &f,double dx_log,bool drawThis = true,double emin,double emax)
{
  cout << "Finding Lafferty-Wyatt Energies..." << endl;
  if(g == NULL)
    {
      cout << "No Spectra! Exiting..." << endl;
      return;
    }
  if(f == NULL)
    {
      cout << "No Spectral fit! Exiting... " << endl;
      return;
    }
  if(g->GetN() <= 2)
    {
      cout << "Not enough data points! Exiting... " << endl;
      return;
    }
  double dx;
  double x_lw;
  double x1,x2,x,x_log;
  double fx_lw,fx_mean,dfx,dfx_min;
  const int nSteps = 1e4;
  double x_tmp;
  TGraphErrors* g_lw = (TGraphErrors*)g->Clone("g_lw");
  
  for(int i=0; i<g->GetN(); i++)
    {
      dfx_min = 1e3;
      x = g->GetX()[i];
      x_log = TMath::Log10(x); // log_TeV
      x1 = pow(10.0,x_log - dx_log/2.0);
      x2 = pow(10.0,x_log + dx_log/2.0);
      dx = x2 - x1;

      fx_mean = f->Integral(x1,x2)/dx;
      for(int j=0; j<nSteps; j++)
	{
	  x_tmp = x1 + (x2 - x1)*j/nSteps;
	  dfx = TMath::Abs(fx_mean - f->Eval(x_tmp));
	  if(dfx < dfx_min)
	    {
	      dfx_min = dfx;
	      x_lw = x_tmp;
	    }
	}
      cout << "Bin center: " << x << " Lafferty-Wyatt Energy: " << x_lw << endl; 
      g_lw->SetPoint(i,x_lw,g->GetY()[i]);
      g_lw->SetPointError(i,0,g->GetErrorY(i));
      
    }
  if(drawThis)
    {
      TCanvas* c1 = new TCanvas();
      c1->SetLogy(1);
      c1->SetLogx(1);
      g->Draw("A*");
      g_lw->SetMarkerColor(kBlue);
      g_lw->SetMarkerStyle(kCircle);
      g_lw->SetLineColor(kBlue);
      f->SetLineColor(kBlue);
      g_lw->Draw("p");
      TLegend* l = new TLegend(0.6,0.6,0.9,0.9);
      l->AddEntry(g,"Bin Center Energy","p");
      l->AddEntry(g_lw,"Lafferty-Wyatt Energy","p");
      l->Draw("same");
    }
  g = g_lw;
 
  g->Fit(f,"","",emin,emax);
} 
