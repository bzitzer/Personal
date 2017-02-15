
void dumpStg6Results_wZnCorr(string inFile,string outFile, int nEvents = 0,bool useRBM = 1,double theta2Cut = 0.03, bool NoFromStg5 = true,string Stg5Path = "/data/veritas/bzitzer/bootes_1/data/root/")
{
  gROOT->SetBatch(true);
  TH1D* hSigRF = new TH1D("SigDistRF","Significance Distrbution RF",50,-5,5);
  TH1D* hSigCBG = new TH1D("SigDistCBG","Significance Distrbution CBG",50,-5,5);
  
  // Opening files and getting VEGAS objects:
  TFile* f = new TFile(inFile.c_str(),"READ");
  if(!f->IsOpen() )
    {
      cerr << "Problem opening ROOT file!" << endl;
      return;
    }
  TTree* EventTree = (TTree*)gDirectory->Get("EventStatsTree");
  if( EventTree == NULL )
    {
      cout << "No Event Tree!" << endl;
      return;
    }
  TTree* RunTree = (TTree*)gDirectory->Get("RunStatsTree");
  if( RunTree == NULL )
    {
      cout << "No Run Tree!" << endl;
      return;
    }
  TH1F* accPlot = gDirectory->Get("RingBackgroundModelAnalysis/AcceptancePlotSmoothed");
  if(accPlot == NULL)
    {
      cout << "NO acceptance plot!" << endl;
      return;
    }
  //EventTree->SetDirectory(0);
  //RunTree->SetDirectory(0);
  VASkyMap* vaOnMap  = (VASkyMap*)gDirectory->Get("RingBackgroundModelAnalysis/fIntOnMap");
  VASkyMap* vaAccMap = (VASkyMap*)gDirectory->Get("RingBackgroundModelAnalysis/fIntAccMap");
  if(vaOnMap == NULL || vaAccMap == NULL )
    {
      cout << "Problem with skymaps!" << endl;
      return;
    }
  VACoordinatePair onCenter = vaOnMap->GetCenter();
  VACoordinatePair pt_ev;
  VACoordinatePair pt_tr;  
  VACoordinatePair pt_so;

  double TelLatRad = 5.52828386357865242e-01;
  double TelLongRad = -1.93649167430676461e+00;
  Float_t EnergyGeV,El,Az;
  double RA,Dec;
  double RA_tr,Dec_tr;
  double RA_so,Dec_so;
  double DayNS;
  UInt_t MJD;
  UInt_t RunID,RunID_old;
  Float_t El_tr,Az_tr;
  //Float_t El_check,Az_check;
  double MJDDbl;
  Double_t W;
  Double_t liveTime;
  Double_t PsiEventTree;
  Int_t RunNumRunTree;
  int NumRuns = RunTree->GetEntries(); 
  Bool_t IsOn,IsOff;
  double Noise;

  vector<VASkyMapExclusionRegion*> vaExcl;
  getExclRegions(f,vaExcl);
  cout << "Got here" << endl;

  EventTree->SetBranchAddress("RunNum",&RunID);
  EventTree->SetBranchAddress("Azimuth",&Az);
  EventTree->SetBranchAddress("Elevation",&El);
  EventTree->SetBranchAddress("Noise",&Noise);
  EventTree->SetBranchAddress("EnergyGeV",&EnergyGeV);
  
  EventTree->SetBranchAddress("TrackingAzimuth",&Az_tr);
  EventTree->SetBranchAddress("TrackingElevation",&El_tr);
  EventTree->SetBranchAddress("OnEvent",&IsOn);  
  EventTree->SetBranchAddress("OffEvent",&IsOff);    
  EventTree->SetBranchAddress("Weight",&W);
  EventTree->SetBranchAddress("Psi",&PsiEventTree);

  EventTree->SetBranchAddress("MJDDbl",&MJDDbl);
  EventTree->SetBranchAddress("DayNSDbl",&DayNS);

  RunTree->SetBranchAddress("faLiveTime",&liveTime);
  RunTree->SetBranchAddress("fSourceDecDeg",&Dec_so);
  RunTree->SetBranchAddress("fSourceRADeg",&RA_so);
  RunTree->SetBranchAddress("faRunNumber",&RunNumRunTree);

  VAAzElRADecXY coord(TelLongRad,TelLatRad);
  VATime time;
  TGraph* map = new TGraph();
  TGraph* trackError = new TGraph();
  TH2D* map2 = new TH2D("skymap","raw counts map",100,65,115,100,10,30);
  double X,Y;
  double XRot,YRot;
  double theta;
  
  // vectors to store information through two loops. _ev is for vectors storing event-by-event information, _run is for run information
  //vector<Bool_t> IsOn_ev;
  vector<double> EnergyGeV_ev,El_ev,Az_ev,El_tr_ev,No_ev;
  vector<double> RA_ev,Dec_ev,RA_tr_run,Dec_tr_run,RA_tr_tmp,Dec_tr_tmp;
  vector<double> RA_tr_ev,Dec_tr_ev,Az_tr_ev;
  vector<double> RA_so_run,Dec_so_run;
  vector<double> MJDDbl_ev;
  vector<UInt_t> RunID_ev;
  vector<double> No_run,LT_run;
  vector<double> psi_ev;
  vector<int> IsOn_ev;
  // ---------------------------------------------------------------------------------
  //  First loop: Getting events from trees, coordinate transforms, calculating noise
  // ---------------------------------------------------------------------------------

  cout << "Got here" << endl;
  if(nEvents == 0 )
    nEvents = EventTree->GetEntries();
  
  RunTree->GetEntry(0);
  int j=0;
  // RunNumRunTree = 0;
  for(int i=0; i<nEvents; i++)
    {
      EventTree->GetEntry(i);
      if( i%10000 == 0 ){ cout << "On Event: " << i << " of: " << nEvents << endl; }
      if( RunNumRunTree != RunID && RunNumRunTree != 0 )
	{
	  RA_tr_run.push_back(meanVector(RA_tr_tmp));
	  Dec_tr_run.push_back(meanVector(Dec_tr_tmp));
	  
	  RA_so_run.push_back(RA_so);
	  Dec_so_run.push_back(Dec_so);

	  RA_tr_tmp.clear();
	  Dec_tr_tmp.clear();

	  LT_run.push_back(liveTime);
	  if(NoFromStg5)
	    No_run.push_back(getAvgPedVar(RunNumRunTree,Stg5Path));
	  //else
	  //  No_run.push_back(Noise);
	  j++;
	  RunTree->GetEntry(j);

	  if(RunNumRunTree != RunID)
	    cout << "Warning! Tree mis-match!" << endl;
	}
      // cout << El << " " << El_tr << endl;
      // Energy Cut:
      // if( EnergyGeV < 300.0 ){ continue; }

      // coodinate transforms:     
      time.setFromMJDDbl(MJDDbl);
      coord.AzEl2RADec2000(Az,El,time,RA,Dec); // RA,Dec in radians
      // Az_track, El_track in degrees
      coord.AzEl2RADec2000(Az_tr*TMath::DegToRad(),El_tr*TMath::DegToRad(),time,RA_tr,Dec_tr); // RATrack,DecTrack in radians
      Az *= TMath::RadToDeg();
      El *= TMath::RadToDeg();
      RA_tr *= TMath::RadToDeg();
      Dec_tr *= TMath::RadToDeg();
      RA *= TMath::RadToDeg();
      Dec *= TMath::RadToDeg();
     
      RA_ev.push_back(RA);
      Dec_ev.push_back(Dec);

      RA_tr_tmp.push_back(RA_tr);
      Dec_tr_tmp.push_back(Dec_tr);

      RA_tr_ev.push_back(RA_tr);
      Dec_tr_ev.push_back(Dec_tr);

      El_ev.push_back(El);
      Az_ev.push_back(Az);
      El_tr_ev.push_back(El_tr);
      Az_tr_ev.push_back(Az_tr);
      if(!NoFromStg5)
	No_ev.push_back(Noise);

      RunID_ev.push_back(RunID);
      psi_ev.push_back(PsiEventTree);
      EnergyGeV_ev.push_back(EnergyGeV);
      MJDDbl_ev.push_back(MJDDbl);
      IsOn_ev.push_back(IsOn);
    } //Endoth the 1st loop.
  cout << "Size of run ID: " << RunID_ev.size() << endl;
  cout << "Size of RA: " << RA_ev.size() << endl;
  cout << "Size of RA tracking: " << RA_tr_ev.size() << endl;
  cout << "Size of El tracking: " << El_tr_ev.size() << endl;

  cout << "Size of El: " << El_ev.size() << endl;
  cout << "Size of Energy: " << EnergyGeV_ev.size() << endl;
  cout << "Size of MJD: " << MJDDbl_ev.size() << endl;
  cout << "Size of IsOn: " << IsOn_ev.size() << endl;
  cout << "Size of psi: " << psi_ev.size() << endl;
  cout << "Size of Noise: " << No_ev.size() << endl;
  //return;
  //cout << "Size of RA tracking: " << RA_tr_ev.size() << endl;
  

  RA_tr_run.push_back(meanVector(RA_tr_tmp));
  Dec_tr_run.push_back(meanVector(Dec_tr_tmp));
	  
  RA_so_run.push_back(RA_so);
  Dec_so_run.push_back(Dec_so);
  LT_run.push_back(liveTime);
  //No_run.push_back(getAvgPedVar(RunNumRunTree,Stg5Path));
  if(NoFromStg5)
    No_run.push_back(getAvgPedVar(RunNumRunTree,Stg5Path));
  //  else
  //  No_run.push_back(Noise);
  // -----------------------------------------
  //  Between loops: Making Zn acceptance map
  // -----------------------------------------
 
  string strFitFun = "pol4";
  double binSize = 0.01;
  cout << "Making zenith map... ";
  VASkyMap* znMap = makeZnMap(El_ev,El_tr_ev,RA_ev,Dec_ev,vaAccMap,binSize);
  TCanvas* c1 = new TCanvas("c1","c1",50,50,700,700);
  znMap->Draw();
  cout << "done!" << endl;


  cout << "Making ON/Acceptance map..";
  VASkyMap* ratioMap = makeRatioMap(vaOnMap,vaAccMap,binSize);
  TCanvas* c2 = new TCanvas("c2","c2",60,60,700,700);
  ratioMap->Draw();
  cout << "done!" << endl;


  cout << "Correlating maps... ";
  TGraph* gCorr = correlateZnMap(znMap,ratioMap,strFitFun);
  TF1* fFitFun = gCorr->GetFunction(strFitFun.c_str());
  cout << "done!" << endl;


  cout << "Making Zenith-corrected acceptance map... ";
  VASkyMap* znAccMap0 = makeZnAccMap(fFitFun,znMap,vaAccMap,binSize);
  cout << "done!" << endl;
  //TCanvas* c3 = new TCanvas("c3","c3",70,70,700,700);
  //znAccMap0->Draw();
  VASkyMap* znAccMap = new VASkyMap();
  // return;
  // ---------------------------------------------------------------------------------------------------------
  //  Second loop: Determining ON/OFF flags, calculating alpha (weight), putting information into ascii files
  // ---------------------------------------------------------------------------------------------------------
  int RunIDOld = 0;
  double W,offset,theta,psi,W_old;
  double cres_upper,cres_lower;
  int Non,Noff;
  double sig,ex,alpha;
  int j=0;
  bool IsInExcl;
  vector<double> Non_run,Noff_run,alpha_run;
  vector<double> sig_run,ex_run,runID_run;
  filebuf fb;
  fb.open(outFile.c_str(),ios::out);
  ostream os(&fb);
  os << "RunID LiveTime(sec) Time(MJD) RA       Dec       RA_track        Dec_track    Energy(GeV)  IsOn Weight Elevation Azimuth Noise Offset" << endl;
  os << "-------------------------------------------------------------------------------------------------------------------------------------" << endl;
  cout.precision(12);
  os.precision(7);
  nEvents = RunID_ev.size();
  for(int i=0; i<nEvents; i++)
    {
      if( i%1000 == 0 ){ cout << "On event: " << i << " of: " << nEvents << endl; }
      IsOn = false;
      IsOff = false;
      IsInExcl = false;
      //cout << RunID_ev.at(i) << endl;
      // Getting Weight:
      if( RunIDOld != RunID_ev.at(i) )
	{
	  if( j != 0 )
	    {
	      // Do stuff for the end of run accounting:
	      sig = lima(Non,Noff,W);
	      ex = Non - W*Noff;
	      sig_run.push_back(sig);
	      ex_run.push_back(ex);
	      Non_run.push_back(Non);
	      Noff_run.push_back(Noff);
	      alpha_run.push_back(W);
	      runID_run.push_back(RunIDOld);
	      cout << "Run#: " << RunIDOld;
	      cout << " Non: " << Non;
	      cout << " Noff: " << Noff;
	      cout << " alpha: " << W;
	      cout << " Significance: " << sig;
	      cout << " Excess: " << ex << endl;
	    }
	  Non = 0;
	  Noff = 0;
	  pt_tr.setCoordinates_Deg(RA_tr_run.at(j),Dec_tr_run.at(j),2);
	  pt_so.setCoordinates_Deg(RA_so_run.at(j),Dec_so_run.at(j),2);
	  offset = pt_tr.angularSeparation_Deg(pt_so);
	  cres_upper = offset + sqrt(theta2Cut);
	  cres_lower = offset - sqrt(theta2Cut);
	  znAccMap = makeZnAccMapRun(pt_tr,accPlot,znAccMap0);
	  znAccMap->Draw();
	  //return;
	  W = calcAlpha(znAccMap,pt_tr,pt_so,vaExcl,cres_upper,cres_lower,sqrt(theta2Cut));
	  cout << "Alpha for run: " << RunID_ev.at(i) << " " << W << endl;
	  j++;
	}
      // Checking ON flag:
      pt_ev.setCoordinates_Deg(RA_ev.at(i),Dec_ev.at(i),2);
      theta = pt_ev.angularSeparation_Deg(pt_so);
      
      //cout << theta << " " << theta2Cut << endl;
      if( theta < sqrt(theta2Cut) )
	IsOn = true;

      if( IsOn != IsOn_ev.at(i) )
	cout << "Warning Mis-match with the ON Flag!" << endl;

      // checking OFF Flag
      psi = pt_ev.angularSeparation_Deg(pt_tr);
      if( cres_lower < psi_ev.at(i) && psi_ev.at(i) < cres_upper  && !isInExclRegion(pt_ev,vaExcl) )
	IsOff = true;

      if(IsOn)
	{
	  Non++;
	  alpha = 0.0;
	}
      if(IsOff)
	{
	  Noff++;
	  alpha = W;
	}

      // Information from ON/OFF flags dumped into ASCII file...
      if(IsOn || IsOff)
	{
	  //cout << "Got here: " << j << endl;
	  os << RunID_ev.at(i) << " ";
	  os << LT_run.at(j-1) << " ";

	  os.precision(12);
	  os << MJDDbl_ev.at(i) << " ";
	  os.precision(9);
	  os << RA_ev.at(i) << " ";
	  os << Dec_ev.at(i) << " ";
	  os << RA_tr_ev.at(i) << " ";
	  os << Dec_tr_ev.at(i) << " ";
	  os.precision(7);
	  os << EnergyGeV_ev.at(i) << "        ";
	  os << IsOn << "    ";
	  os << alpha << "       ";
	
	  os << El_tr_ev.at(i) << " ";
	  os << Az_tr_ev.at(i) << " ";
	  //os << No_run.at(j-1) << " ";
	  if(!NoFromStg5)
	    os << No_ev.at(i) << " ";
	  else
	    os << No_run.at(j-1) << " ";
	  os << psi_ev.at(i) << " ";
	  os << endl;
	}
	  
      RunIDOld = RunID_ev.at(i);
    }
  cout << "Got here" << endl;
  sig = lima(Non,Noff,W);
  ex = Non - W*Noff;
  sig_run.push_back(sig);
  ex_run.push_back(ex);
  Non_run.push_back(Non);
  Noff_run.push_back(Noff);
  alpha_run.push_back(W);
  runID_run.push_back(RunIDOld);
  cout << "Run#: " << RunIDOld;
  cout << " Non: " << Non;
  cout << " Noff: " << Noff;
  cout << " alpha: " << W;
  cout << " Significance: " << sig;
  cout << " Excess: " << ex << endl;

  fb.close();
  /*
  cout << Non_run.size() << endl;
  cout << Noff_run.size() << endl;
  cout << LT_run.size() << endl;
  cout << alpha_run.size() << endl;
  return;
  */
  // --------------------------------
  //  Dumping results to Results.txt
  // --------------------------------
  VAStatisticsUtilitiesAnl* StatAnl = new VAStatisticsUtilitiesAnl(Non_run,Noff_run,LT_run,LT_run,alpha_run);
  
  fb.open("Results.txt",ios::out);
  ostream os(&fb);
  for(int i=0; i<Non_run.size(); i++)
    {
      os << "Results for run#: " << runID_run.at(i) << endl;
      os << "  Number of ON events: " << Non_run.at(i) << endl;
      os << "  Number of OFF events: " << Noff_run.at(i) << endl;
      os << "  Alpha: " << alpha_run.at(i) << endl;
      os << "  Exp Time: " << LT_run.at(i) << endl;
      os << "  Significance: " << sig_run.at(i) << endl;
      os << "  " << endl;
      hSigCBG->Fill(sig_run.at(i));
    }
  
  os << "---------------------------" << endl;
  os << "Final Results for all runs:" << endl;
  os << "  Number of ON events: " << sumVector(Non_run) << endl;
  os << "  Number of OFF events: " << sumVector(Noff_run) << endl;
  os << "  Mean Alpha: " << calcWeightAvgVector(alpha_run,LT_run) << endl;
  os << "  Total Exp Time: " << sumVector(LT_run) << endl;
  os << "  Excess: " << StatAnl->ExcessRate()*60.0 << " +/- " << StatAnl->ExcessRateError()*60.0 << endl;
  os << "  Generalized LiMa Significance: " << StatAnl->GeneralisedLiMa() << endl;

  // ---------------------------------------------
  //  putting useful plots, etc into Results.root
  // ---------------------------------------------
  TFile* fOut = new TFile("Results.root","RECREATE");
  StatAnl->Write();
  znAccMap0->Write();
  znMap->Write();
  ratioMap->Write();
  hSigCBG->Write();
  gCorr->Write("gCorr");
  accPlot->Write("accPlot");
  fOut->Close();
}

double calcWeightAvgVector(vector<double> x,vector<double> w)
{

  if( x.size() != w.size() )
    cout << "Warning! Vector size mis-match!" << endl;

  double x_sum,w_sum = 0;
  for(int i=0; i<x.size(); i++)
    {
      x_sum += x.at(i)*w.at(i);
      w_sum += w.at(i);
    }
  if( w_sum != 0.0 )
    return(x_sum/w_sum);
  else
    return(0.0);
}
double calcAlpha(VASkyMap* accAreaMap,VACoordinatePair pt_tr,VACoordinatePair pt_so,vector<VASkyMapExclusionRegion*> exclList,const double cres_upper,const double cres_lower,const double thetaCut = 0.17)
{
  double areaOn = 0;
  double psi;  
  double areaOff = DBL_EPSILON;
  VACoordinatePair pt;
  double theta;
  for(int i=1; i<=accAreaMap->GetNbinsX(); i++)
    for(int j=1; j<=accAreaMap->GetNbinsY(); j++)
      {
	//if(!accAreaMap->isBinFilled(i,j)){ continue; }
	pt = accAreaMap->GetBinCoordinates(i,j);
	psi = pt.angularSeparation_Deg(pt_tr);
	theta = pt.angularSeparation_Deg(pt_so);
	if( theta < thetaCut )
	  areaOn += accAreaMap->GetBinArea(i,j);
	//areaOn += accAreaMap->GetBinArea(i,j)*accAreaMap->GetBinContent(i,j);
	if( psi < cres_lower || cres_upper < psi ){ continue; }
	//if( theta <= ptExclRadius ){ continue; }
	if( isInExclRegion(pt,exclList) ){ continue; }
	areaOff += accAreaMap->GetBinArea(i,j);
	//accAreaMap->GetBinContent(i,j);
      }
  return(areaOn/areaOff);

}

VASkyMap* makeZnAccMapRun(VACoordinatePair ptTrack,TH1F* hRadAccept,VASkyMap* accMap0)
{
  VASkyMap* vaZnAccMap = (VASkyMap*)accMap0->Clone();
  double acc0,psi,acc;
  VACoordinatePair pt;
  int psiBin;
  for(int i=1; i<=accMap0->GetNbinsX(); i++)
    for(int j=1; j<=accMap0->GetNbinsY(); j++)
      {
	if(!accMap0->isBinFilled(i,j)){ continue; }
	pt = accMap0->GetBinCoordinates(i,j);
	psi = pt.angularSeparation_Deg(ptTrack);
	acc = accMap0->GetBinContent(pt);
	psiBin = hRadAccept->FindBin(psi);
	acc *= hRadAccept->GetBinContent(psiBin);
	vaZnAccMap->SetBinContent(pt,acc);
      }
  return(vaZnAccMap);
}
VASkyMap* makeZnAccMap(TF1* fFit,VASkyMap* znMap,VASkyMap* accMap,double binSize = 0.005)
{

  VASkyMap* vaZnAccMap = new VASkyMap("vaZnAccMap","vaZnAccMap",accMap->GetCenter(),accMap->GetFoVX(),binSize);
  vaZnAccMap->MakeBinAreaMap();
  double* p = fFit->GetParameters();
  const int nPar = fFit->GetNpar();
  double zn;
  VACoordinatePair pt;
  double acc;
  if(p[0] == 0.0){ return(vaZnAccMap); }
  
  for(int i=1; i<=znMap->GetNbinsX(); i++)
    for(int j=1; j<=znMap->GetNbinsY(); j++)
      {
	
        if( !znMap->isBinFilled(i,j) ){ continue; }
	pt = accMap->GetBinCoordinates(i,j);
	zn = znMap->GetBinContent(pt);
	if( zn == 0.0 ){ continue; }
	acc = 1.0;
	for(int k=1; k<nPar; k++)
	  acc += (p[k]/p[0])*pow(zn,k);

	//acc *= accMap->GetBinContent(pt);
	vaZnAccMap->SetBinContent(pt,acc);
      }
  return(vaZnAccMap);
}

TGraph* correlateZnMap(VASkyMap* znMap,VASkyMap* ratioMap,string strFitFun)
{

  TGraph* gCorr = new TGraph();
  VACoordinatePair pt;
  double zn,ratio;
  int k=0;
  for(int i=1; i<=znMap->GetNbinsX(); i++)
    for(int j=1; j<=znMap->GetNbinsY(); j++)
      {
	if(!znMap->isBinFilled(i,j)){ continue; }
	pt = znMap->GetBinCoordinates(i,j);
	zn = znMap->GetBinContent(pt);
	ratio = ratioMap->GetBinContent(pt);
	if(ratio == 0.0){ continue; }
	gCorr->SetPoint(k,zn,ratio);
	k++;
      }
  cout << "Number of graph points: " << k << endl;
  gCorr->Fit(strFitFun.c_str());
  return(gCorr);
}

VASkyMap* makeRatioMap(VASkyMap* onMap,VASkyMap* accMap,double binSize = 0.005)
{
  double on,acc;
  VASkyMap* vaOnAccMap = new VASkyMap("vaOnAccMap","vaOnAccMap",onMap->GetCenter(),onMap->GetFoVX(),binSize);
  VACoordinatePair pt;
  for(int i=1; i<=vaOnAccMap->GetNbinsX(); i++)
    for(int j=1; j<=vaOnAccMap->GetNbinsY(); j++)
      {
	//if(!onMap->isBinFilled(i,j)){ continue; }
	pt = vaOnAccMap->GetBinCoordinates(i,j);
	on  = onMap->GetBinContent(pt);
	acc = accMap->GetBinContent(pt);
	if( acc == 0.0 ){ continue; }
	vaOnAccMap->SetBinContent(pt,on/acc);
      }
  return(vaOnAccMap);


}
VASkyMap* makeZnMap(vector<double> el,vector<double> el_tr,vector<double> RA,vector<double> Dec,VASkyMap* vaMap0,double binSize = 0.005)
{
  VACoordinatePair ptCenter = vaMap0->GetCenter();
  VASkyMap* vaRawZnMap = new VASkyMap("vaRawZnMap","vaRawZnMap",vaMap0->GetCenter(),vaMap0->GetFoVX(),binSize);
  VASkyMap* vaRawOnMap = new VASkyMap("vaRawOnMap","vaRawOnMap",vaMap0->GetCenter(),vaMap0->GetFoVX(),binSize);
  VASkyMap* vaZnMap = new VASkyMap("vaZnMap","vaZnMap",vaMap0->GetCenter(),vaMap0->GetFoVX(),binSize);
  VACoordinatePair pt;
  //cout << el.size() << " " << el_tr.size() << " "  << RA.size() << " " << Dec.size() << endl; 
  if( el.size() != el_tr.size())
    cout << "Warning! elevation and tracking elevation vector size mis-match!" << endl;
  if( el.size() != RA.size())
    cout << "Warning! elevation and tracking elevation vector size mis-match!" << endl;
  if( el.size() != Dec.size())
    cout << "Warning! elevation and tracking elevation vector size mis-match!" << endl;
  
  double zn,zn_tr;
  VACoordinatePair pt2;
  const double smoothRadius = 0.17;
  for(int i=0; i<el.size(); i++)
    {
      //if(i%100 == 0){ cout << "On event: " << i << " of: " << el.size() << endl; }
      pt.setCoordinates_Deg(RA.at(i),Dec.at(i),2);
      //if( ptCenter.angularSeparation_Deg(pt) > 1.7 ){ continue; }
      zn = 90.0 - el.at(i);
      zn_tr = 90.0 - el_tr.at(i);
      vaRawZnMap->FillCoordinates(pt,zn - zn_tr);
      vaRawOnMap->FillCoordinates(pt,1.0);
	
    }
    
  VASkyMap* vaSmZnMap = smoothSkyMap(vaRawZnMap,smoothRadius);
  VASkyMap* vaSmOnMap = smoothSkyMap(vaRawOnMap,smoothRadius);

  double znRaw,onRaw;
  for(int i=1; i<=vaZnMap->GetNbinsX(); i++)
    for(int j=1; j<=vaZnMap->GetNbinsY(); j++)
      {
	znRaw = vaSmZnMap->GetBinContent(i,j);
	onRaw = vaSmOnMap->GetBinContent(i,j);
	//cout << znRaw << " " << onRaw << endl;
	if(onRaw <= 0.0){ continue; }
	vaZnMap->SetBinContent(i,j,znRaw/onRaw);
      }
  return(vaZnMap);
}

VASkyMap* smoothSkyMap(VASkyMap* map,double radius)
{
  //map->MakeBinAreaMap();
  int nBinsX,nBinsY;
  double binWidthX,binWidthY;
  int nBinsInRegion = 0;
  double sumInRegion = 0;
  VACoordinatePair pt,pt0;
  VASkyMap* map_smooth = map->Clone("map_smooth");
  /*
  for(int i=1; i<=map_smooth->GetNbinsX(); i++)
    for(int j=1; j<=map_smooth->GetNbinsY(); j++)
      map_smooth->SetBinContent(i,j,0.0);
  */
  for(int x=1; x<=map->GetNbinsX(); x++)
    for(int y=1; y<=map->GetNbinsX(); y++)
      {
	//if( !map->isBinFilled(x,y) ){ continue; }

	pt0 = map->GetBinCoordinates(x,y); 
	binWidthX = map->GetXaxis()->GetBinWidth(x);
	binWidthY = map->GetYaxis()->GetBinWidth(y);
	nBinsX = 1 + (int)radius/binWidthX;
	nBinsY = 1 + (int)radius/binWidthY;
	for(int i=-nBinsX; i<=nBinsX; i++)
	  for(int j=-nBinsY; j<=nBinsY; j++)
	    {
	      //if( map->GetBinContent(x+i,y+j) == 0.0 ){ continue; }

	      pt = map->GetBinCoordinates(x+i,y+j);
	      if( pt.angularSeparation_Deg(pt0) <= radius )
		{
		  nBinsInRegion++;
		  sumInRegion += map->GetBinContent(pt);
		}
	    }
	if(nBinsInRegion == 0)
	  map_smooth->SetBinContent(x,y,map->GetBinContent(pt0));
	else
	  map_smooth->SetBinContent(x,y,sumInRegion/nBinsInRegion);

	nBinsInRegion = 0;
	sumInRegion = 0;
      }
    
  return(map_smooth);

}


double meanVector(vector<double> vec)
{
  if(vec.size() == 0)
    return(0);

  double sum = 0;
  for(int i=0; i<vec.size(); i++)
    sum += vec.at(i);

  return(sum/vec.size());
}
double sumVector(vector<double> vec)
{
  if(vec.size() == 0)
    return(0);

  double sum = 0;
  for(int i=0; i<vec.size(); i++)
    sum += vec.at(i);

  return(sum);
}

double lima(double on, double off, double alpha)
{


  double diff = on-alpha*off;

  if (on>0 && off>0){
    double lima17 = TMath::Sqrt(2) * 
      TMath::Power(
                   on * TMath::Log( ((1.0+alpha)/alpha) * (on/(on+off))  )
                   +
                   off * TMath::Log( (1.0+alpha)*(off/(on+off)))
                   , 0.5);

  }
  else {
    cout<<"warning "<<on<<" "<<off<<endl;
    return 0;
  }
  
  if(diff<0)
    return -lima17;
  else
    return lima17;
}

double lima(vector<double> &on, vector<double> &off, vector<double> &alpha){
 vector<double> onePlusAlpha;

 double fVAOn=0;
 double fVAOff=0;

 for (unsigned i=0; i<on.size(); i++){
   fVAOn+=on.at(i);
   fVAOff+=off.at(i);
 }
  double fFirstDenominator=0;
  double fSecondDenominator=0;
  double fFirstTerm=0;
  double fSecondTerm=0;
  double fG17;
  for (unsigned int i=0 ;i <alpha.size() ;i++) {
    onePlusAlpha.push_back(1+alpha[i]);
  }
  for (unsigned int i=0 ;i <alpha.size() ;i++) {
    fFirstDenominator += (alpha[i]/onePlusAlpha[i])*(on[i]+off[i]);
    fSecondDenominator += (1/onePlusAlpha[i])*(on[i]+off[i]);
    // printf("The denominator : %f %f\n",
    //fFirstDenominator,fSecondDenominator);
  }
  double fDiff = 0;
  double fSum = 0;
  if ( fVAOn > 0 && fVAOff > 0 )
    {
      for (unsigned int i=0 ;i <alpha.size() ;i++)
        {
          fFirstTerm += on[i]*log(fVAOn/fFirstDenominator);
          fSecondTerm += off[i]*log(fVAOff/fSecondDenominator);
        } 
      fSum = fFirstTerm + fSecondTerm;
      for (unsigned int i=0 ;i <alpha.size() ;i++)
        {
          fDiff += on[i] - alpha[i]*off[i];
        }
    }
  else
    {
      fDiff = 0;
      fSum = 0;
    }
  if (fDiff >= 0 )
    {
      fG17=sqrt(2.)*sqrt(fSum);
      return fG17;
    }
  else
else
    {
      fG17=sqrt(2.)*sqrt(fSum);
      return -fG17;
    }


}

void getExclRegions(TFile* f,vector<VASkyMapExclusionRegion*> &vaSourceExcl)
{
  /*
  TFile* f = new TFile(inFile.c_str(),"READ");
  if(!f->IsOpen())
    {
      cerr << "Problem with input file! " << endl;
      return;
    }
  */
  TDirectory* RBMExclusion = (TDirectory*)gDirectory->Get("RingBackgroundModelAnalysis/ExclusionRegions");

  if( RBMExclusion == NULL )
    {
      cerr << "Problem loading the RBM exclusion directory!" << endl;
      return;
    }
  int nRegions = RBMExclusion->GetNkeys();
  VASkyMapExclusionRegion* hSourceExclusion;
  const int tmp = nRegions;
  VASkyMapExclusionRegion* exclList[tmp];
  //vector<VASkyMapExclusionRegion*> vaSourceExcl;

  TIter next(RBMExclusion->GetListOfKeys());
  TKey *key;
  int i=0;
  while(key=(TKey*)next())
    {
      hSourceExclusion = (VASkyMapExclusionRegion*)RBMExclusion->FindObjectAny(key->GetName())->Clone();
      if( hSourceExclusion != NULL)
	{
	  if( hSourceExclusion->wasUsed() )
	    {
	      cout << i << endl;
	      exclList[i] = hSourceExclusion;
	      vaSourceExcl.push_back(hSourceExclusion);
	      cout << hSourceExclusion->GetName() << endl;
	      //cout << "Exclusion Center RA: " << hSourceExclusion->center().getRA_J2000_Deg() << endl;
	      cout << "Exclusion Center RA: " << exclList[i]->center().getRA_J2000_Deg() << endl;
	      
	      cout << "Exclusion Center Dec: " << hSourceExclusion->center().getDec_J2000_Deg() << endl;	  
	      cout << "Exclusion Radius: " << hSourceExclusion->radius_Deg() << endl;
	      i++;
	    }
	}
    }
  nRegions = i;
  cout << "Number of exclusion regions: " << vaSourceExcl.size() << endl;
  // f->Close();

}

double integrateOnRegion(VACoordinatePair ptBinCenter,VASkyMap* accAreaMap,double thetaCut)
{

  double area = DBL_EPSILON;
  double theta;
  VACoordinatePair pt;
  for(int i=1; i<=accAreaMap->GetNbinsX(); i++)
    for(int j=1; j<=accAreaMap->GetNbinsY(); j++)
      {
	pt = accAreaMap->GetBinCoordinates(i,j);
	theta = pt.angularSeparation_Deg(ptBinCenter);
	if( theta < thetaCut)
	  area+= accAreaMap->GetBinArea(i,j)*accAreaMap->GetBinContent(i,j);
      }
  return(area);
}

double integrateBgRegion(VACoordinatePair ptBinCenter,VACoordinatePair ptTrack,VASkyMap* accAreaMap,double lowerRadius,double upperRadius,vector<VASkyMapExclusionRegion*> exclList,double ptExclRadius)
{
  double psi,theta;
  double area = DBL_EPSILON;
  VACoordinatePair pt;
  for(int i=1; i<=accAreaMap->GetNbinsX(); i++)
    for(int j=1; j<=accAreaMap->GetNbinsY(); j++)
      {
	pt = accAreaMap->GetBinCoordinates(i,j);
	psi = pt.angularSeparation_Deg(ptTrack);
	theta = pt.angularSeparation_Deg(ptBinCenter);
	if( psi <= lowerRadius || upperRadius <= psi ){ continue; }
	if( theta <= ptExclRadius ){ continue; }
	if( isInExclRegion(pt,exclList) ){ continue; }
	area += accAreaMap->GetBinArea(i,j)*accAreaMap->GetBinContent(i,j);
      }
  return(area);
}



bool isInExclRegion(VACoordinatePair pt,vector<VASkyMapExclusionRegion*> excl)
{

  for(int i=0; i<excl.size(); i++)
    if(excl.at(i)->isWithinRegion(pt))
      return(true);

  return(false);
}

void getPositions(string inFile,vector<double> &eventRA,vector<double> &eventDec,vector<double> &eventPsi,vector<double> &eventZn,vector<double> &eventTrZn,double &meanTrRA,double &meanTrDec)
{
  cout << "Getting positions for every event..." << endl;
  TFile* f = new TFile(inFile.c_str(),"READ");
  if(!f->IsOpen() )
    {
      cerr << "Problem opening ROOT file!" << endl;
      return;
    }
  TTree* EventTree = (TTree*)gDirectory->Get("EventStatsTree");
  if( EventTree == NULL )
    {
      cout << "No Event Tree!" << endl;
      return;
    }

  Float_t az,el,azTr,elTr;
  double psi,mjddbl;
  double ra,dec,raTr,decTr;
  vector<double> eventTrDec,eventTrRA;

  double TelLatRad = 5.52828386357865242e-01;
  double TelLongRad = -1.93649167430676461e+00;

  VAAzElRADecXY coord(TelLongRad,TelLatRad);
  VATime time;

  EventTree->SetBranchAddress("Azimuth",&az);
  EventTree->SetBranchAddress("Elevation",&el);
  EventTree->SetBranchAddress("TrackingAzimuth",&azTr);
  EventTree->SetBranchAddress("TrackingElevation",&elTr);
  EventTree->SetBranchAddress("Psi",&psi);
  EventTree->SetBranchAddress("MJDDbl",&mjddbl);

  for(int i=0; i<EventTree->GetEntries(); i++)
    {
      EventTree->GetEntry(i);
      if(i%5000 == 0){ cout << "On Event " << i << " of: " << EventTree->GetEntries() << endl; }
      time.setFromMJDDbl(mjddbl);
      azTr *= TMath::DegToRad();
      elTr *= TMath::DegToRad();
      coord.AzEl2RADec2000(az,el,time,ra,dec); // RA,Dec in radians
      // Az_track, El_track in degrees
      coord.AzEl2RADec2000(azTr,elTr,time,raTr,decTr); // RATrack,DecTrack in radians

      // converting to degrees:
      raTr *= TMath::RadToDeg();
      decTr *= TMath::RadToDeg();
      ra *= TMath::RadToDeg();
      dec *= TMath::RadToDeg();

      eventRA.push_back(ra);
      eventDec.push_back(dec);
      eventPsi.push_back(psi);
      eventTrRA.push_back(raTr);
      eventTrDec.push_back(decTr);
      eventZn.push_back(90.0 - (el*TMath::RadToDeg()) );
      eventTrZn.push_back(90.0 - (elTr*TMath::RadToDeg()) );
    }

  meanTrRA = meanVector(eventTrRA);
  meanTrDec = meanVector(eventTrDec);
  cout << "done! " << endl;
  cout << meanVector(eventRA) << endl;
  cout << meanVector(eventDec) << endl;
  f->Close();
  //f->Delete();
}


double getAvgPedVar(int RunNum,string inRootPath,int winSize = 7, int numTels = 4)
{
  ostringstream os;
  os << inRootPath << "/" << RunNum << ".stg5.root";
  cout << "Getting PedVar from: " << os.str() << endl;
  return(getAvgPedVar(os.str(),winSize,numTels));

}
double getAvgPedVar(string inRootFile,int winSize = 7, int numTels = 4)
{
  //cout << inRootFile << endl;
  if(numTels==0)
    return(0);

  VARootIO io(inRootFile.c_str(),true);
  io.loadTheRootFile();
  if(!io.IsOpen())
    {
      cout << "Warning! No root file! " << endl;
      return(0);
    }
  const VAPixelStatusData* pd = io.loadThePixelStatusData();
  const VAArrayInfo* ai = io.loadTheArrayInfo();
  const VAQStatsData* qd = io.loadTheQStatsData();
  
  double PedVar = 0;
  for(int i=0; i<numTels; i++)
    PedVar += qd->getCameraAverageTraceVarTimeIndpt(i,7,pd,ai);
  io.closeTheRootFile();
  return(PedVar/numTels);
}
