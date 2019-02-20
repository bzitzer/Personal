
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

#include "aclicPreProcCommands.h"
#include <Riostream.h>

#include "VASkyMapExclusionRegion.h"
#include "VASkyMap.h"
#include "VATime.h"
#include "VAAzElRADecXY.h"
#include "VACoordinatePair.h"
#include "VAStatisticsUtilitiesAnl.h"

using namespace std;

VATime getAvgTime(TTree* tRun);

void getPositionInfo(TTree* tRun,double &source_ra,double &source_dec,double &offset_ra,double &offset_dec,double &exp);

void rotateSkyMap(string inOnFile,string inOffFile,VASkyMap* &alphaMapRot,VASkyMap* &offMapRot,VASkyMap* &onMap);

void getExclRegions(string inFile,vector<VASkyMapExclusionRegion> &vaSourceExcl);

void combineMatchedRuns(string s6RootFileList,bool useModLiMa = false)
{

  vector<string> onFile,offFile;

  string onFileTmp,offFileTmp;
  
  vector<VASkyMap*> onMap;
  vector<VASkyMap*> offMap;
  vector<VASkyMap*> alphaMap;
  vector<VASkyMapExclusionRegion> vaExclReg;
  
  ifstream in(s6RootFileList.c_str());
  int numPairs = 0;
  while(in >> onFileTmp >> offFileTmp)
    {
      cout << onFileTmp << endl;
      onFile.push_back(onFileTmp);
      offFile.push_back(offFileTmp);
      cout << "On File: " << numPairs+1 << " " << onFile.at(numPairs) << endl;
      cout << "Off File: " << numPairs+1 << " " << offFile.at(numPairs) << endl; 
      numPairs++;
    }
  cout << "Number of pairs: " << numPairs << endl;
  for(int i=0; i<numPairs; i++)
    {

      onMap.push_back(new VASkyMap());
      offMap.push_back(new VASkyMap());
      alphaMap.push_back(new VASkyMap());
      rotateSkyMap(onFile.at(i),offFile.at(i),
		   alphaMap.at(i),offMap.at(i),onMap.at(i));
      if(i==0)
	getExclRegions(onFile.at(i),vaExclReg);
    }
  /*
  TCanvas* c1 = new TCanvas("c1","c1",50,50,500,500);
  onMap.at(0)->Draw();
  TCanvas* c2 = new TCanvas("c2","c2",60,60,500,500);
  offMap.at(0)->Draw();
  TCanvas* c3 = new TCanvas("c3","c3",70,70,500,500);
  alphaMap.at(0)->Draw();
  
  return;
  */
  VAStatisticsUtilitiesAnl* vaStat = new VAStatisticsUtilitiesAnl();
  vector<double> on;
  vector<double> off;
  vector<double> alpha;
  vector<double> exp;
  vector<double> sigma_alpha;
  double ex = 0;
  double sig = 0;
  int numBins = 0;
  VASkyMap* sigMap = (VASkyMap*)onMap.at(0)->Clone("fSigMap");
  VASkyMap* exMap  = (VASkyMap*)onMap.at(0)->Clone("fExMap");
  sigMap->SetTitle("Significance Map");
  exMap->SetTitle("Excess Map");

  on.clear();
  off.clear();
  alpha.clear();
  exp.clear();

  //return;
  TH1D* hSig = new TH1D("hSig","hSigDist",160,-6,10);
  TH1D* hSigNExcl = new TH1D("hSigNExcl","hSigDistNExcl",160,-6,10);
  
  bool isWithinExclReg = false;
  VACoordinatePair pt;
  VACoordinatePair ptMaxSig;
  const int nExclRegion = vaExclReg.size();

  double maxSig = 0.0;
  double maxExc = 0.0;
  double on_sum = 0;
  double off_sum = 0;
  double alpha_eff = 0.0;
  double on_sum_max = 0;
  double off_sum_max = 0;
  double alpha_eff_max = 0;
  double sigma_alpha_max = 0;
  double sigma_alpha_eff = 0;
  for(int i=1; i<=onMap.at(0)->GetNbinsX(); i++)
    for(int j=1; j<=onMap.at(0)->GetNbinsY(); j++)
      {
	
	if(numBins%10000 == 0)
	  cout << "On bin: " << numBins << " of: " << onMap.at(0)->GetNbinsX()*onMap.at(0)->GetNbinsY() << endl;
	
	numBins++;
	ex = 0.0;
	on_sum = 0;
	off_sum = 0;
	alpha_eff = 0;
	for(int k=0; k<numPairs; k++)
	  {
	    if( !onMap.at(k)->isBinFilled(i,j) ){ continue; }
	    if( TMath::Abs(alphaMap.at(k)->GetBinContent(i,j)) < DBL_EPSILON ){ continue; }
	    on.push_back(onMap.at(k)->GetBinContent(i,j));
	    off.push_back(offMap.at(k)->GetBinContent(i,j));
	    alpha.push_back(alphaMap.at(k)->GetBinContent(i,j));
	    if(useModLiMa)
	      {
		//sigma_alpha.push_back(alphaMap.at(k)->GetBinError(i,j));
		sigma_alpha.push_back(
				      alphaMap.at(k)->GetHistogram()->GetBinError(i,j));
		//cout << "Alpha: " << alpha.back() << " sigma alpha: " << sigma_alpha.back() <<
		//  " sqrt(alpha): << " << sqrt(alpha.back()) << endl;
	      }
	    exp.push_back(1.0);
	    on_sum  += onMap.at(k)->GetBinContent(i,j);
	    off_sum += offMap.at(k)->GetBinContent(i,j);
	    alpha_eff +=  offMap.at(k)->GetBinContent(i,j)*alphaMap.at(k)->GetBinContent(i,j);
	  }
	if(on.size() == 0 ){ continue; }

	alpha_eff = alpha_eff/(off_sum + DBL_EPSILON);
	
	vaStat = new VAStatisticsUtilitiesAnl(on_sum,off_sum,1.0,1.0,alpha_eff);
	if(!useModLiMa)
	  {
	    sig = vaStat->ExcessSignificance();
	  }
	else
	  {
	    sigma_alpha_eff = 0;
	    for(int k=0; k<on.size(); k++)
	      {
		sigma_alpha_eff += pow(off.at(k)*sigma_alpha.at(k),2.0);
		sigma_alpha_eff += off.at(k)*pow(alpha.at(k) - alpha_eff,2.0);
	      }
	    sigma_alpha_eff = sqrt(sigma_alpha_eff)/(off_sum + DBL_EPSILON);	    
	    sig = vaStat->ExcessSignificanceModified(sigma_alpha_eff);
	  }
	for(int k=0; k<on.size(); k++)
	  ex += on.at(k) - alpha.at(k)*off.at(k);

	sigMap->SetBinContent(i,j,sig);
	exMap->SetBinContent(i,j,vaStat->ExcessRate());

	if(useModLiMa)
	  hSig->Fill(sig);
	else
	  hSig->Fill(sig);

	isWithinExclReg = false;
	pt = onMap.at(0)->GetBinCoordinates(i,j);
	
	for(int k=0; k<nExclRegion; k++)
	  {
	    if( vaExclReg.at(k).isWithinRegion(pt) )
	      isWithinExclReg = true;
	  }
	if(!isWithinExclReg)
	  hSigNExcl->Fill(sig);

	if( sig > maxSig )
	  {
	    maxSig = sig;
	    ptMaxSig = pt;
	    on_sum_max = on_sum;
	    off_sum_max = off_sum;
	    alpha_eff_max = alpha_eff;
	    if(useModLiMa)
	      sigma_alpha_max = sigma_alpha_eff;
	  }
	on.clear();
	off.clear();
	alpha.clear();
	exp.clear();
	sigma_alpha.clear();
	ex = 0.0;

      }
  cout << "Max. Significance : " << maxSig << endl;
  cout << "Sum ON counts: " << on_sum_max << " Sum OFF counts: " << off_sum_max << " eff. alpha: " << alpha_eff_max << endl;
  if(useModLiMa)
    cout << "Sigma alpha: " << sigma_alpha_max << endl;

  if(!useModLiMa)
    for(int k = 0; k<onMap.size(); k++)
      cout << "On counts: " << onMap.at(k)->GetBinContent(ptMaxSig) << " Off counts: " << offMap.at(k)->GetBinContent(ptMaxSig) << " alpha: " << alphaMap.at(k)->GetBinContent(ptMaxSig) << " " << endl;
  else
    {
      int i,j;
      for(int k=0; k<onMap.size(); k++)
	{
	  alphaMap.at(k)->FindBinPosition(ptMaxSig,i,j);
	  cout << "On counts: " << onMap.at(k)->GetBinContent(ptMaxSig) << " Off counts: " << offMap.at(k)->GetBinContent(ptMaxSig) << " alpha: " << alphaMap.at(k)->GetBinContent(ptMaxSig) << " sigma_alpha: " << alphaMap.at(k)->GetBinError(i,j) << endl;
	}

    }

  
  hSigNExcl->SetLineColor(kGreen+2);
  
  TCanvas* c1 = new TCanvas("c1","c1",50,50,600,600);
  sigMap->setScale(-6,6);
  sigMap->Draw();
  for(int i=0; i<nExclRegion; i++)
    {
      vaExclReg.at(i).SetDrawColor(kBlack);
      vaExclReg.at(i).Draw(sigMap);
    }
  TCanvas* c2 = new TCanvas("c2","c2",60,60,600,600);
  exMap->Draw();
   for(int i=0; i<nExclRegion; i++)
    {
      vaExclReg.at(i).SetDrawColor(kBlack);
      vaExclReg.at(i).Draw(exMap);
    }
  TCanvas* c3 = new TCanvas("c3","c3",70,70,600,600);
  c3->SetLogy(1);
  hSig->Draw();
  hSigNExcl->Draw("same");
  hSigNExcl->Fit("gaus");
  
}

void rotateSkyMap(string inOnFile,string inOffFile,VASkyMap* &alphaMapRot,VASkyMap* &offMapRot,VASkyMap* &onMap)
{

  TFile* f1 = new TFile(inOnFile.c_str());
  
  if(!f1->IsOpen())
    {
      cout << "Blah Blah file not open blah blah..." << endl;
      return;
    }

  onMap = (VASkyMap*)gDirectory->Get("RingBackgroundModelAnalysis/fIntOnMap");
  TTree* tRunOn = (TTree*)gDirectory->Get("RunStatsTree")->Clone("tRunOn");
  VATime tOn = getAvgTime(tRunOn);
    
  if(tRunOn == NULL)
    {
      cout << "Can not find run tree for ON run! " << endl;
      return;
    }

  //  tRunOn->SetDirectory(0);
  //f1->Close();
  
  double x_coord,y_coord;
  double ra, dec;
  double noff,alpha;
  double source_ra,source_dec;
  double offset_ra,offset_dec;
  double off_source_ra,off_source_dec;
  double off_offset_ra,off_offset_dec;
  double exp_on,exp_off;
  cout << "Got here!" << endl;
  /*
  tRunOn->SetBranchAddress("fSourceRADeg" ,&source_ra);
  tRunOn->SetBranchAddress("fSourceDecDeg",&source_dec);
  tRunOn->SetBranchAddress("fOffsetRADeg" ,&offset_ra);
  tRunOn->SetBranchAddress("fOffsetDecDeg",&offset_dec);
  tRunOn->GetEntry(0);
  */
  //tRunOn->Delete();
  //f1->Close();
  //  tRunOn->SetDirectory(0);
  cout << "On Source RA: " << source_ra << " On Source Dec: " << source_dec << endl;
  cout << "On Offset RA: " << offset_ra << " On Offset Dec: " << offset_dec << endl;
  
  TFile* f2 = new TFile(inOffFile.c_str());
  if(!f2->IsOpen())
    {
      cout << "Blah Blah file not open blah blah..." << endl;
      return;
    }
  
  VASkyMap* offMap = (VASkyMap*)gDirectory->Get("RingBackgroundModelAnalysis/fRingBackground");
  //VASkyMap* offMap = (VASkyMap*)gDirectory->Get("RingBackgroundModelAnalysis/fIntOnMap");
  
  VASkyMap* alphaMap = (VASkyMap*)gDirectory->Get("RingBackgroundModelAnalysis/fAlphaMap");

  offMapRot = (VASkyMap*)onMap->Clone("offMapRot");
  alphaMapRot = (VASkyMap*)onMap->Clone("alphaMapRot");

  offMapRot->SetTitle("Off Map");
  alphaMapRot->SetTitle("Alpha Map");
  
  VACoordinatePair pt;
  VACoordinatePair ptCenterOn  = onMap->GetCenter();
  VACoordinatePair ptCenterOff = offMap->GetCenter();
  VACoordinatePair pt2;
  
  const double TelLatRad = 5.52828386357865242e-01;
  const double TelLongRad = -1.93649167430676461e+00;

  VAAzElRADecXY coord(TelLongRad,TelLatRad);

  TTree* tRunOff = (TTree*)gDirectory->Get("RunStatsTree")->Clone("tRunOff");

  if(tRunOff == NULL)
    {
      cout << "Can not find run tree for OFF run! " << endl;
      return;
    }

  //tRunOff->SetDirectory(0);
  
  VATime t   = getAvgTime(tRunOff);
 
  for(int i=1; i<=onMap->GetNbinsX(); i++)
    for(int j=1; j<=onMap->GetNbinsY(); j++)
      {
	if(!onMap->isBinFilled(i,j)){ continue; }
	offMapRot->SetBinContent(i,j,0.0);
	alphaMapRot->SetBinContent(i,j,0.0);

      }
  
  getPositionInfo(tRunOn,source_ra,source_dec,offset_ra,offset_dec,exp_on);  
  getPositionInfo(tRunOff,off_source_ra,off_source_dec,off_offset_ra,off_offset_dec,exp_off);

  double s = exp_on/(exp_off + DBL_EPSILON);
  cout << "Got here 2" << endl;
  
  
  for(int i=1; i<=onMap->GetNbinsX(); i++)
    for(int j=1; j<=onMap->GetNbinsY(); j++)
      {
	if(!onMap->isBinFilled(i,j)){ continue; }
	pt = onMap->GetBinCoordinates(i,j);
	//x_coord = pt.getRA_J2000_Deg()  - ( offset_ra  + source_ra  );
	//y_coord = pt.getDec_J2000_Deg() - ( offset_dec + source_dec );
	//	cout << "x_coord: " << x_coord << " y_coord: " << y_coord << endl;
	//coord.XY2RADec2000(x_coord, y_coord, t,
	//		   ptCenterOff.getRA_J2000_Deg() + off_offset_ra,
	//		   ptCenterOff.getDec_J2000_Deg() + off_offset_dec,
	//		   ra, dec);
	//	coord.RADec2000ToXY(pt.getRA_J2000_Deg(),pt.getDec_J2000_Deg(),tOn,
	//		    ptCenterOn.getRA_J2000_Deg()+offset_ra,ptCenterOn.getDec_J2000_Deg()+offset_dec,
	//		    x_coord,y_coord);
	coord.RADec2000ToXY(pt.getRA_J2000_Rad(),pt.getDec_J2000_Rad(),tOn,
			    ptCenterOn.getRA_J2000_Rad()+offset_ra*TMath::DegToRad(),ptCenterOn.getDec_J2000_Rad()+offset_dec*TMath::DegToRad(),
			    x_coord,y_coord);
	//ra  = x_coord + (off_offset_ra  + off_source_ra);
	//dec = y_coord + (off_offset_dec + off_source_dec);
	//ra  += +1.0*off_offset_ra;
	//dec += +1.0*off_offset_dec;
	//cout << "x coord: " << x_coord << " y coord: " << y_coord << endl;
	coord.XY2RADec2000(x_coord,y_coord,t,
			   ptCenterOff.getRA_J2000_Rad()+off_offset_ra*TMath::DegToRad(),ptCenterOff.getDec_J2000_Rad()+off_offset_dec*TMath::DegToRad(),
			   ra,dec);
	//	cout << "ra: " << ra*TMath::RadToDeg() << " dec: " << dec*TMath::RadToDeg() << endl;
	pt2.setCoordinates_Deg(ra*TMath::RadToDeg() ,dec*TMath::RadToDeg(), VACoordinates::J2000);
	//pt2.setCoordinates_Deg(ptCenterOff.getRA_J2000_Deg() + x_coord,
	//		       ptCenterOff.getDec_J2000_Deg() + y_coord,
	//		       VACoordinates::J2000);
	noff  = offMap->GetBinContent(pt2);
	alpha = s*alphaMap->GetBinContent(pt2);
	//if(noff != 0)
	//cout << "alpha: " << alpha << " noff: " << noff << endl;
	offMapRot->SetBinContent(i,j,noff);
	alphaMapRot->SetBinContent(i,j,alpha);
	alphaMapRot->SetBinError(i,j,s*alphaMap->GetBinError(i,j));
	//alphaMapRot->SetBinContent(i,j,1.0);
	//cout << off_offset_ra << " " << off_offset_dec << endl;
      }
  cout << "On Source RA: " << source_ra << " On Source Dec: " << source_dec << endl;
  cout << "On Offset RA: " << offset_ra << " On Offset Dec: " << offset_dec << endl;
  cout << "Off Source RA: " << off_source_ra << " Off Source Dec: " << off_source_dec << endl;
  cout << "Off Offset RA: " << off_offset_ra << " Off Offset Dec: " << off_offset_dec << endl;    

}

void getPositionInfo(TTree* tRun,double &source_ra,double &source_dec,double &offset_ra,double &offset_dec,double &exp)
{
  tRun->SetBranchAddress("faDuration"   ,&exp);
  tRun->SetBranchAddress("fSourceRADeg" ,&source_ra);
  tRun->SetBranchAddress("fSourceDecDeg",&source_dec);
  tRun->SetBranchAddress("fOffsetRADeg" ,&offset_ra);
  tRun->SetBranchAddress("fOffsetDecDeg",&offset_dec);
  tRun->GetEntry(0);

}

VATime getAvgTime(TTree* tRun)
{
  double runStart,runEnd;
  tRun->SetBranchAddress("faRunStartMJD",&runStart);
  tRun->SetBranchAddress("faRunEndMJD",&runEnd);
  tRun->GetEntry(0);
  double tMid = (runStart + runEnd)/2.0;

  VATime t;
  
  t.setFromMJDDbl(tMid);

  return(t);
}

void getExclRegions(string inFile,vector<VASkyMapExclusionRegion> &vaSourceExcl)
{
  
  TFile* f = new TFile(inFile.c_str(),"READ");
  if(!f->IsOpen())
    {
      cerr << "Problem with input file! " << endl;
      return;
    }
  
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
	      vaSourceExcl.push_back(*hSourceExclusion);
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
  f->Close();

}
