
#include <vector>
#include <sstream>
#include <iostream>
#include <fstream>

#include <TH1D.h>
#include <TH2D.h>
#include <TGraph.h>
#include <TGraph2D.h>
#include <TGraphAsymmErrors.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TAxis.h>
#include <TF1.h>
#include <TF2.h>
#include <TLegend.h>
#include <TLine.h>
#include <TFile.h>
#include <TTree.h>
//#include "VASlalib.h"
//#include "VASlalib_abr.cpp"
//#include "ZitzEventWeighting.cpp"


enum dwarf_t {segue_1,draco,ursa_minor,bootes,crab};
enum decay_t {WW, ZZ, bbar, tt, ee, uu};

TH2D* renormalizeDisp(TH2D** h,int numFilled);

TTree* readEventFile(string inFile = "../Conventional/Pass5f/segue1_eventList_pass5f_wZnCorr.txt",dwarf_t dwarf = segue_1);

//TTree* addDetResponse(TTree* t,string inEAFile,string dispPath = "../Conventional/Pass5f/segue1_bias/");

//TH2D* readDispFile(string dispFile,string mcBinFile,string recBinFile);
TGraphAsymmErrors* getEffectiveArea(string inFile,int runNum,double &exp);

TH2D* getEnergyDispersion(string path,int runNum);

//double* readBinFile(string inFile,const int nBins = 100);

//double* readBinFile2(string inFile,int &N);
double* readBinFile2(string inFile,int &N,int nHeaderLines = 2);

TH1D* convolveEnergyDisp(TH1D* h,TH2D* hDisp);

TH1D* convolveEDisp_varBinWidth(TH1D* h,TH2D* hDisp);

TH1D* calcPLOnDist(TGraphAsymmErrors* g,double exp,double k,double index = 2.5);

double min_b(double Noff,double Non,double g,double alpha,vector<double> h_E,vector<double> g_E,TGraph* &g1);

double LLBg(double b,double g,double Non,double Noff,double alpha,vector<double> h,vector<double> p);
