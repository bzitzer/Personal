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
#include "TGraphAsymmErrors.h"

#include "aclicPreProcCommands.h"
#include <Riostream.h>

#include "VATime.h"
#include "VACoordinatePair.h"
#include "VAAzElRADecXY.h"
#include "VASpectrumAnl.h"
#include "VAUpperLimit.h"


using namespace std;

TH1F* getCombinedAcceptance(TTree* runTree);

void vaThreeMLPlugin(string inFile,string outFile)
{

  TFile* f = new TFile(inFile.c_str(),"READ");
  if(!f->IsOpen())
    {
      cout << "Missing file!" << endl;
      return;
    }

  TH1F*  acc1D   = (TH1F*)f->Get("RingBackgroundModelAnalysis/AcceptancePlot");
  TTree* runTree = (TTree*)f->Get("RunStatsTree");
  TTree* evTree  = (TTree*)f->Get("EventStatsTree");
  TH2F*  eDisp   = (TH2F*)f->Get("Spectrum/hMigrationMatrix");
  VASpectrumAnl* spec = (VASpectrumAnl*)f->Get("Spectrum/VASpectrumAnl");
  
  if(runTree == NULL)
    {
      cout << "Missing Run Tree!" << endl;
      return;
    }
  if(evTree == NULL)
    {
      cout << "Missing Event Tree!" << endl;
      return;
    }
  if(acc1D == NULL)
    {
      cout << "Missing acceptance curve!" << endl;
      return;
    }
  if(eDisp == NULL)
    {
      cout << "Missing energy migration!" << endl;
      return;
    }
  if(spec == NULL)
    {
      cout << "Missing spectrum!" << endl;
      return;
    }

  if(acc1D->GetBinContent(1) < DBL_EPSILON)
    {
      cout << "Missing combined acceptance, making one... " << endl;
      acc1D = getCombinedAcceptance(runTree);
      cout << "done!" << endl;
    }
  
  acc1D->Smooth(4);
  acc1D->Scale(1.0/acc1D->GetBinContent(1));
  acc1D->Fit("pol4","0");
  TF1* fAccZe_0 = acc1D->GetFunction("pol4");

  TGraphAsymmErrors* gEffArea_cmSqr = (TGraphAsymmErrors*)spec->GetUpperLimit()->GetEffectiveArea();
  if(gEffArea_cmSqr == NULL)
    {
      cout << "Missing Effective Area!" << endl;
      return;
    }
  TGraphAsymmErrors* gEffArea_mSqr = new TGraphAsymmErrors();
  for(int i=0; i<gEffArea_cmSqr->GetN(); i++)
    {
      gEffArea_mSqr->SetPoint(i,gEffArea_cmSqr->GetX()[i],gEffArea_cmSqr->GetY()[i]*1e-4);
    }
  const int nRuns = runTree->GetEntries();
  ostringstream osDir[nRuns];

  double TelLatRad = 5.52828386357865242e-01;
  double TelLongRad = -1.93649167430676461e+00;
  Float_t EffArea,EnergyGeV,logEnergyTeV,El,Az;
  double RA,Dec;
  double RATrack,DecTrack;
  double DayNS;
  UInt_t MJD;
  UInt_t runNum;
  Float_t El_track,Az_track;
  Float_t El_check,Az_check;
  double MJDDbl;
  Double_t W;
  Double_t liveTime,durTime;
  Double_t PsiEventTree;

  Double_t RASource,DecSource,RAOffset,DecOffset;
  double RAError,DecError;
  double fSigRF,fSigCBG;
  int runNumRunTree;
  //bool isOn,isOff;
  Bool_t isOn,isOff;
  //  int isOn,isOff;
  Float_t thSqr;
  Double_t X,Y,X_off,Y_off,X_tr,Y_tr;
  Double_t alpha;

  VATime time;
  VAAzElRADecXY coord(TelLongRad,TelLatRad);
  VACoordinatePair pt;
  VACoordinatePair pt_tr;

  evTree->SetBranchAddress("RunNum",&runNum);
  evTree->SetBranchAddress("Azimuth",&Az);
  evTree->SetBranchAddress("Elevation",&El);
  evTree->SetBranchAddress("EnergyGeV",&EnergyGeV);

  evTree->SetBranchAddress("TrackingAzimuth",&Az_track);
  evTree->SetBranchAddress("TrackingElevation",&El_track);
  evTree->SetBranchAddress("OnEvent",&isOn);  
  evTree->SetBranchAddress("OffEvent",&isOff);    
  evTree->SetBranchAddress("Weight",&W);
  
  evTree->SetBranchAddress("ThetaSq",&thSqr);

  // evTree->SetBranchAddress("RA",&RA_fRoot);    
  // evTree->SetBranchAddress("Dec",&Dec_fRoot);
  
  evTree->SetBranchAddress("MJDDbl",&MJDDbl);
  evTree->SetBranchAddress("DayNSDbl",&DayNS);
  
  runTree->SetBranchAddress("faLiveTime",&liveTime);
  runTree->SetBranchAddress("fSourceDecDeg",&DecSource);
  runTree->SetBranchAddress("fSourceRADeg",&RASource);
  runTree->SetBranchAddress("fOffsetDecDeg",&DecOffset);
  runTree->SetBranchAddress("fOffsetRADeg",&RAOffset);
  runTree->SetBranchAddress("fSignificance",&fSigRF);
  runTree->SetBranchAddress("faRunNumber",&runNumRunTree);
  runTree->SetBranchAddress("faDuration",&durTime);
  runTree->SetBranchAddress("faLiveTime",&liveTime);
  runTree->SetBranchAddress("fAlpha",&alpha);

  TFile* fOut = new TFile(outFile.c_str(),"RECREATE");
  int j = 0;
  runTree->GetEntry(j);
  //osDir[0] << "run_" << runNumRunTree;

  
  runNumRunTree = 0;
  int runNum_old = 0;
  Float_t EnergyOn,EnergyOff;
  double dt_frac;

  TTree* runSummary[nRuns];
  runSummary[0] = new TTree("tRunSummary","tRunSummary");
  TTree* dataOn[nRuns];
  dataOn[0] = new TTree("data_on","data_on");
  TTree* dataOff[nRuns];
  dataOff[0]    = new TTree("data_off","data_off");
  
  dataOn[0]->Branch("Time",&MJDDbl,"Time/D");
  dataOn[0]->Branch("Erec",&logEnergyTeV,"Erec/F");
  dataOn[0]->Branch("theta2",&thSqr,"theta2/F");
  dataOn[j]->Branch("RA",&RA,"RA/D");
  dataOn[j]->Branch("Dec",&Dec,"Dec/D");
  dataOn[j]->Branch("Xoff",&X_off,"Xoff/D");
  dataOn[j]->Branch("Yoff",&Y_off,"Yoff/D");

  dataOff[j]->Branch("Time",&MJDDbl,"Time/D");
  dataOff[j]->Branch("Erec",&logEnergyTeV,"Erec/F");
  dataOff[j]->Branch("theta2",&thSqr,"theta2/F");
  dataOff[j]->Branch("RA",&RA,"RA/D");
  dataOff[j]->Branch("Dec",&Dec,"Dec/D");
  dataOff[j]->Branch("Xoff",&X_off,"Xoff/D");
  dataOff[j]->Branch("Yoff",&Y_off,"Yoff/D");

  runSummary[0]->Branch("DeadtimeFracOn",&dt_frac,"DeadtimeFracOn/D");
  runSummary[0]->Branch("OffNorm",&alpha,"OffNorm/D");
  runSummary[0]->Branch("TargetRAJ2000",&RASource,"TargetRAJ2000/D");
  runSummary[0]->Branch("TargetDecJ2000",&DecSource,"TargetDecJ2000/D");
  runSummary[0]->Branch("tOn",&durTime,"tOn/D");

  for(int i=0; i<evTree->GetEntries(); i++)
    {
      evTree->GetEntry(i);
      if(i%10000 == 0)
	cout << "On event: " << i << " of: " << evTree->GetEntries() << endl;

      //cout << EnergyGeV << endl;
      //cout << runNum << endl;
      if( runNum_old != runNum && runNum_old != 0 )
	{
	  dt_frac = (durTime - liveTime)/liveTime;
	  runSummary[j]->Fill();
	  
	  cout << "on run: " << j << " of: " << nRuns << endl;
	  //osDir[j] << "run_" << runNumRunTree;
	  osDir[j] << "run_" << runNum_old;
	  // TDirectory* dir = new TDirectory(osDir[j].str().c_str(),osDir[j].str().c_str());
	  fOut->mkdir(osDir[j].str().c_str());
	  fOut->cd(osDir[j].str().c_str());

	  logEnergyTeV = TMath::Log10(EnergyGeV/1.0e3);
	  runSummary[j]->Write();
	  dataOn[j]->Write();
	  dataOff[j]->Write();
	  fAccZe_0->Write("fAccZe_0");
	  eDisp->Write("hMigration");
	  gEffArea_mSqr->Write("gMeanEffectiveArea");

	  cout << "Finished Run: " << runNum_old << endl;
	  cout << " Data On entries: " << dataOn[j]->GetEntries() << endl;
	  cout << " Data Off entries: " << dataOff[j]->GetEntries() << endl;
	  cout << " Run Summary entries: " << runSummary[j]->GetEntries() << endl;

	  fOut->cd("..");
	  
	  j++;
	  runTree->GetEntry(j);
	  if( runNumRunTree != runNum )
	    {
	      cout << "Warning! Num numbers do not match! " << endl;

	    }
	 

	  runSummary[j] = new TTree("tRunSummary","tRunSummary");
	  dataOn[j]     = new TTree("data_on","data_on");
	  dataOff[j]    = new TTree("data_off","data_off");
	  
	  dataOn[j]->Branch("Time",&MJDDbl,"Time/D");
	  dataOn[j]->Branch("Erec",&EnergyGeV,"Erec/F");
	  dataOn[j]->Branch("theta2",&thSqr,"theta2/F");
	  dataOn[j]->Branch("RA",&RA,"RA/D");
	  dataOn[j]->Branch("Dec",&Dec,"Dec/D");
	  dataOn[j]->Branch("Xoff",&X_off,"Xoff/D");
	  dataOn[j]->Branch("Yoff",&Y_off,"Yoff/D");
  
	  dataOff[j]->Branch("Time",&MJDDbl,"Time/D");
	  dataOff[j]->Branch("Erec",&EnergyGeV,"Erec/F");
	  dataOff[j]->Branch("theta2",&thSqr,"theta2/F");
	  dataOff[j]->Branch("RA",&RA,"RA/D");
	  dataOff[j]->Branch("Dec",&Dec,"Dec/D");
	  dataOff[j]->Branch("Xoff",&X_off,"Xoff/D");
	  dataOff[j]->Branch("Yoff",&Y_off,"Yoff/D");

	  runSummary[j]->Branch("DeadtimeFrac",&dt_frac,"DeadtimeFrac/D");
	  runSummary[j]->Branch("OffNorm",&alpha,"OffNorm/D");
	  runSummary[j]->Branch("TargetRAJ2000",&RASource,"TargetRAJ2000/D");
	  runSummary[j]->Branch("TargetDecJ2000",&DecSource,"TargetDecJ2000/D");
	  runSummary[j]->Branch("tOn",&durTime,"tOn/D");

  
	}

      //      if(isOn)
      //	cout << EnergyGeV << endl;
      if(isOn || isOff)
	{
	  time.setFromMJDDbl(MJDDbl);
	  coord.setRASource2000(RASource*TMath::DegToRad());
	  coord.setDecSource2000(DecSource*TMath::DegToRad());

	  coord.AzEl2RADec2000(Az,El,time,RA,Dec);
	  RA  *= TMath::RadToDeg();
	  Dec *= TMath::RadToDeg();
	  //pt.
	  coord.AzElToXY(Az,El,time,X,Y);

	  Az_track *= TMath::DegToRad();	  
	  El_track *= TMath::DegToRad();

	  coord.AzElToXY(Az_track,El_track,time,X_tr,Y_tr);
	  //coord.AzEl2RADec2000(Az_track,El_track,time,RA_track,Dec_track);
	  X_off = X - X_tr;
	  Y_off = Y - Y_tr;
	}


      if(isOn)
	{
	  EnergyOn = EnergyGeV;
	  dataOn[j]->Fill();
	}
      if(isOff)
	{
	  EnergyOff = EnergyGeV;
	  dataOff[j]->Fill();
	}
      runNum_old = runNum;
    }

  dt_frac = (durTime - liveTime)/liveTime;
  runSummary[j]->Fill();
  
  cout << "on run: " << j << " of: " << nRuns << endl;
  if( runNumRunTree != runNum )
    {
      cout << "Warning! Num numbers do not match! " << endl;
      
    }
  osDir[j] << "run_" << runNumRunTree;
  // TDirectory* dir = new TDirectory(osDir[j].str().c_str(),osDir[j].str().c_str());
  fOut->mkdir(osDir[j].str().c_str());
  fOut->cd(osDir[j].str().c_str());
  runSummary[j]->Write();
  dataOn[j]->Write();
  dataOff[j]->Write();
  fAccZe_0->Write("fAccZe_0");
  eDisp->Write("hMigration");
  gEffArea_mSqr->Write("gMeanEffectiveArea");

  cout << "Finished Run: " << runNum_old << endl;
  cout << " Data On entries: " << dataOn[j]->GetEntries() << endl;
  cout << " Data Off entries: " << dataOff[j]->GetEntries() << endl;
  cout << " Run Summary entries: " << runSummary[j]->GetEntries() << endl;


}

TH1F* getCombinedAcceptance(TTree* runTree)
{
  const int nRuns = runTree->GetEntries();
  ostringstream os[nRuns];
  TH1F* acc[nRuns];
  //TH1F* accTot = new TH1F("Combined Acceptance Curve","Combined Accpetance",29,0,1.7*1.7);
  double liveTime,totLiveTime;
  vector<double> lt_vec;
  int runNum;
  runTree->SetBranchAddress("faLiveTime",&liveTime);  
  runTree->SetBranchAddress("faRunNumber",&runNum);
  int nBins;
  double xmin,xmax;
  for(int i=0; i<nRuns; i++)
    {
      runTree->GetEntry(i);
      os[i] << "RingBackgroundModelAnalysis/AcceptanceCurves/AcceptanceForRun"
	    << runNum;
      acc[i] = (TH1F*)gDirectory->Get(os[i].str().c_str());
      if(acc[i] == NULL)
	{
	  cout << "Warning! Missing Acceptance for run: " << runNum << endl; 
	  continue;
	}
      totLiveTime += liveTime;
      lt_vec.push_back(liveTime);
      nBins = acc[i]->GetNbinsX();
      xmax  = acc[i]->GetXaxis()->GetXmax();
      xmin  = acc[i]->GetXaxis()->GetXmin();
    }
  TH1F* accTot = new TH1F("Combined Acceptance Curve","Combined Accpetance",nBins,xmin,xmax);
  if(totLiveTime == 0)
    {
      cout << "Warning! Total live time is zero!" << endl;
      return(NULL);
    }
  
  for(int i=0; i<nRuns; i++)
    {
      if(acc[i] == NULL){ continue; }
      accTot->Add(acc[i],lt_vec.at(i)/totLiveTime);
    }
  return(accTot);
}
