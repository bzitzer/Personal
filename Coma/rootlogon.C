{

cout<<"--------------------------"<<endl;
cout<<"-                        -"<<endl;
cout<<"-     VEGAS MACROS       -"<<endl;
cout<<"-                        -"<<endl;
cout<<"--------------------------"<<endl;
cout<<endl;

#include <vector>
#include "TTree.h"
#include "TH1.h"
cout<<"- Loading Shared Library -"<<endl;
gSystem->Load("libTreePlayer.so");
gSystem->Load("libPhysics.so");
gSystem->Load("/homes/gretel/bzitzer/vegas/vegas_liMaDev/common/lib/libSP24sharedLite.so");
gSystem->Load("/homes/gretel/bzitzer/vegas/vegas_liMaDev/resultsExtractor/lib/libStage6shared.so");
gSystem->Load("/homes/gretel/bzitzer/vegas/vegas_liMaDev/showerReconstruction2/lib/libStage4.so");
cout<<"-     Loading Macros     -"<<endl;

 gSystem->AddIncludePath("-Wno-unused -Wno-shadow -Wno-unused-parameter");


gROOT->ProcessLine(".L /homes/gretel/bzitzer/vegas/vegas_liMaDev/common/include/VACommon.h");
gROOT->ProcessLine(".include /homes/gretel/bzitzer/vegas/vegas_liMaDev/common/include/");
gROOT->ProcessLine(".include /homes/gretel/bzitzer/vegas/vegas_liMaDev/resultsExtractor/include/");
gROOT->ProcessLine(".include /homes/gretel/bzitzer/vegas/vegas_liMaDev/cfitsio/include/");
gROOT->ProcessLine(".include /homes/gretel/bzitzer/vegas/vegas_liMaDev/resultsExtractor/src/");

gEnv->SetValue("TGraph.BrowseOption", "ALP*");
gStyle->SetFrameBorderMode(0); 
gStyle->SetCanvasDefW(1000);
gStyle->SetCanvasColor(0);
gStyle->SetTitleX(0.5);
gStyle->SetTitleAlign(23);
gStyle->SetTitleFillColor(0);
gStyle->SetHistLineWidth(2);
gStyle->SetStatColor(0);
cout<<endl;
cout<<"--------------------------"<<endl;
cout<<"-   Loading Complete     -"<<endl; 
cout<<"--------------------------"<<endl;
cout<<endl;
cout<<"You can query what each of the above macros does by doing eg"<<endl;
cout<<"dumpConfigFile()"<<endl;
cout<<"This will tell you what this macro does, and show you how to use it"<<endl;
cout<<"Feel free to report any errors/unexpected behaviour."<<endl;
cout<<"Errors will not be fixed unless you report them!"<<endl;

}


