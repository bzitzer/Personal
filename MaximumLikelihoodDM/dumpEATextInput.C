

void dumpEATextInput(string InFile,string configPath,bool useMCE = 0, bool testMode = 0,bool useOlder = 0)
{
  ifstream ifs(InFile.c_str());
  if(!ifs.is_open())
    {
      cerr << "Problem with input file! " << endl;
      return;
    }

  ostringstream OutFileNm;
  if(useMCE)
    OutFileNm << InFile << ".mce";
  else
    OutFileNm << InFile << ".rece";

  cout << OutFileNm.str() << endl;

  ofstream out(OutFileNm.str().c_str());
  out << "Run#   Live Time(min)   Energy(GeV)   EA(m^2)   -DeltaEA(m^2)   +DeltaEA(m^2)" << endl;
  out << "-----------------------------------------------------------------------------" << endl;

  char tmpString[250];
  ifs.getline(tmpString,250);
  cout << tmpString << endl;
  ifs.getline(tmpString,250);
  cout << tmpString << endl;

  int RunID;
  float LT,EventTime,RA,Dec,Energy,W,El,Az,N,Offset;
  float RA_tr,Dec_tr;
  bool IsOn;
  vector<float> El_vec,Az_vec,N_vec;
  vector<float> El_vec_bg,Az_vec_bg,N_vec_bg;
  int RunID_old = 0 ;
  float LT_old;
  float El_m,Az_m,N_m;
  TGraphAsymmErrors* g = new TGraphAsymmErrors();
  TGraphAsymmErrors* gAvg = new TGraphAsymmErrors();
  //TGraphAsymmErrors* gAvg_tmp = new TGraphAsymmErrors();
  
  double TotLT = 0;
  int nLines = 0;
  int nRuns = 0;
  while(!ifs.eof())
    {
      
      nLines++;
      if( testMode && nLines > 5e3){ break; }
      if(!useOlder)
	ifs >> RunID >> LT >> EventTime >> RA >> Dec >> RA_tr >> Dec_tr >> Energy >> IsOn >> W >> El >> Az >> N >> Offset;
      else
	ifs >> RunID >> LT >> EventTime >> RA >> Dec >> Energy >> IsOn >> W >> El >> Az >> N >> Offset;
      //if(singleRun !=0 && singleRun != RunID){ continue; }
      //cout << RunID << endl;
      if(RunID != RunID_old && RunID_old !=0)
	{	
	  nRuns++;
	  //TGraphAsymmErrors* gAvg_tmp = new TGraphAsymmErrors();
	  cout << "Finished reading run number " << RunID_old << endl;
	 
	  if( El_vec.size() == 0 || Az_vec.size() == 0 || N_vec.size() == 0 )
	    {
	      El_m = AvgVec(El_vec_bg);
	      Az_m = AvgVec(Az_vec_bg);
	      N_m = AvgVec(N_vec_bg);
	      cout << "Run with no ON events! " << endl;
	    }
	  else
	    {
	      El_m = AvgVec(El_vec);
	      Az_m = AvgVec(Az_vec);
	      N_m = AvgVec(N_vec);	  
	      cout << El_m << " " << Az_m << " " << N_m << endl;
	    }
	  g = getEAInterpolatedCurve(RunID_old,El_m,Az_m,N_m,useMCE,configPath);
	  El_vec.clear();
	  Az_vec.clear();
	  N_vec.clear();
	  El_vec_bg.clear();
	  Az_vec_bg.clear();
	  N_vec_bg.clear();
	  
	  WriteTGraphToText(OutFileNm.str(),g,RunID_old,LT_old);
	  cout << "Run Number: " <<  RunID_old << " Number of points: " << g->GetN() <<" Number of average points: " << gAvg->GetN() << endl;
	  //cout << TotLT << endl;
	  //cout << LT_old << endl;
	  cout << nRuns << endl;
	  if( nRuns == 1 )
	    gAvg = g;
	  else
	    {
	      cout << "Got here 5" << endl;
	      AvgEACurve3(gAvg,g,TotLT,LT_old);	       
	      //AvgEACurve3(g,gAvg,LT_old,TotLT);
	  
	      //   gAvg = gAvg_tmp;
	      // delete gAvg_tmp;
	    }
	  cout << "Sum LT: " << TotLT << " Run LT: " << LT_old << " Avg. EA @ 1TeV " << gAvg->Eval(0) << " Run EA @ 1TeV: " << g->Eval(0) << endl;
	  TotLT += LT_old;
	}
       
      if(IsOn)
	{
	  El_vec.push_back(El);
	  Az_vec.push_back(Az);
	  N_vec.push_back(N);
	}
      else
	{
	  El_vec_bg.push_back(El);
	  Az_vec_bg.push_back(Az);
	  N_vec_bg.push_back(N);
	}
      RunID_old = RunID;
      LT_old = LT;
    }

  cout << "reading file complete" << endl;
  if(El_vec.size() != 0  && Az_vec.size() != 0 && N_vec.size() != 0 )
    {
      El_m = AvgVec(El_vec);
      Az_m = AvgVec(Az_vec);
      N_m = AvgVec(N_vec);
      cout << El_m << " " << Az_m << " " << N_m << endl;
      g = getEAInterpolatedCurve(RunID_old,El_m,Az_m,N_m,useMCE,configPath);
      if(nRuns == 0)
	gAvg = g;
      else
	AvgEACurve3(gAvg,g,TotLT,LT_old);
      TotLT += LT;
      bool didWrite = !WriteTGraphToText(OutFileNm.str(),g,RunID,LT);
      if(!didWrite)
	{
	  cerr << "Problem Writting output file!" << endl;
	}
    }
  
  bool didWrite = !WriteTGraphToText(OutFileNm.str(),gAvg,0,TotLT);
  if(!didWrite)
    {
      cerr << "Problem Writting output file!" << endl;
    }
  gAvg->Draw("AE*");
  //c1->SetLogy(1);
  out.close();
}

void AvgEACurve3(TGraphAsymmErrors* &g1, TGraphAsymmErrors* g2,double w1,double w2)
{
cout << "Got Here..." << endl;
  cout << w1 << " " << w2 << endl;
  cout << g1->GetN() << " " << g2->GetN() << endl;
  

  TGraph* g2Low = new TGraph();
  TGraph* g2High = new TGraph();
  TGraph* g1Low = new TGraph();
  TGraph* g1High = new TGraph();
  
  if( g1 == NULL && g2 == NULL) { return(NULL); }
  if( g1 == NULL) { return(g2); }
  if( g2 == NULL) { return(g1); }
  // if( w1 == w2 ){ return(g1); }
  
  int n1 = g1->GetN();
  int n2 = g2->GetN();
  
  if( w1 == 0 && w2 ==0 ){ return(NULL); }
  if( w1 == 0 ){ return(g2); }
  if( w2 == 0 ){ return(g1); }
  // if( n1 == 0 && n2 == 0 ){ return(g); }
  //if( n1 == 0 ) { return(g2); }
  //if( n2 == 0 ) { return(g1); }

  // looping to get TGraphs of Y +/- errorbars:
  
  double* X1 = g1->GetX();
  double* Y1 = g1->GetY();
  double* X2 = g2->GetX();
  double* Y2 = g2->GetY();
  cout << n1 << " " << n2 << endl;

  for(int i=0; i<n1; i++)
    {
      g1Low->SetPoint(i,X1[i],Y1[i] - g1->GetErrorYlow(i));
      g1High->SetPoint(i,X1[i],Y1[i] + g1->GetErrorYhigh(i));
      //cout << "Got here 2" << endl;
    }
  for(int i=0; i<n2; i++)
    {
      g2Low->SetPoint(i,X2[i],Y2[i] - g2->GetErrorYlow(i));
      g2High->SetPoint(i,X2[i],Y2[i] + g2->GetErrorYhigh(i));
      //cout << "Got here 3" << endl;
    }

  // getting errorbar arrays:
  double x,y,y1,y2;
  double y1Low,y2Low,y1High,y2High;
  double yL,yH;
  
  for(int i=0; i<n1; i++)
    {
      //    cout << "Got here 4" << endl;
      if( g1->GetX()[i] != g2->GetX()[i] )
	{
	  cout << "X positions not equal!" << endl;
	  cout << g1->GetX()[i] << " " <<  g2->GetX()[i] << endl;
	}
     
      x = g1->GetX()[i];
      y1 = g1->GetY()[i];
      y2 = g2->Eval(x);
      y = (y1*w1 + y2*w2)/(w1 + w2);
      //cout << y1 << " " << y2 << endl;
      g1->SetPoint(i,x,y);

      cout << i << " " << x << " " << y << endl;

      y1Low = g1Low->GetY()[i];
      y2Low = g2Low->Eval(x);
      yL = (y1Low*w1 + y2Low*w2)/(w1 + w2);

      y1High = g1High->GetY()[i];
      y2High = g2High->Eval(x);
      yH = (y1High*w1 + y2High*w2)/(w1 + w2);
      
      g1->SetPointEYlow(i,y - yL);
      g1->SetPointEYhigh(i,yH - y);
    }
  cout << "---------------------" << endl;
  cout << "Number of Eff Area Points (Avg): " << g1->GetN() << endl;
  cout << "---------------------" << endl;
  //  return(g);

}

TGraphAsymmErrors* AvgEACurve2(TGraphAsymmErrors* g1, TGraphAsymmErrors* g2,double w1,double w2)
{
  cout << "Got Here..." << endl;
  cout << w1 << " " << w2 << endl;
  cout << g1->GetN() << " " << g2->GetN() << endl;
  

  TGraph* g2Low = new TGraph();
  TGraph* g2High = new TGraph();
  TGraph* g1Low = new TGraph();
  TGraph* g1High = new TGraph();
  
  if( g1 == NULL && g2 == NULL) { return(g); }
  if( g1 == NULL) { return(g2); }
  if( g2 == NULL) { return(g1); }
  // if( w1 == w2 ){ return(g1); }
  
  const int n1 = g1->GetN();
  const int n2 = g2->GetN();
  TGraphAsymmErrors* g = new TGraphAsymmErrors(n1);
  if( w1 == 0 && w2 ==0 ){ return(g); }
  if( w1 == 0 ){ return(g2); }
  if( w2 == 0 ){ return(g1); }
  // if( n1 == 0 && n2 == 0 ){ return(g); }
  //if( n1 == 0 ) { return(g2); }
  //if( n2 == 0 ) { return(g1); }

  // looping to get TGraphs of Y +/- errorbars:
  
  double* X1 = g1->GetX();
  double* Y1 = g1->GetY();
  double* X2 = g2->GetX();
  double* Y2 = g2->GetY();
  
  for(int i=0; i<n1; i++)
    {
      g1Low->SetPoint(i,X1[i],Y1[i] - g1->GetErrorYlow(i));
      g1High->SetPoint(i,X1[i],Y1[i] + g1->GetErrorYhigh(i));
    }
  for(int i=0; i<n2; i++)
    {
      g2Low->SetPoint(i,X2[i],Y2[i] - g2->GetErrorYlow(i));
      g2High->SetPoint(i,X2[i],Y2[i] + g2->GetErrorYhigh(i));
    }

  // getting errorbar arrays:
  double x,y,y1,y2;
  double y1Low,y2Low,y1High,y2High;
  double yL,yH;
  for(int i=0; i<n1; i++)
    {
      /*
      if( g1->GetX()[i] != g2->GetX()[i] )
	{
	  cout << "X positions not equal!" << endl;
	  cout << g1->GetX()[i] << " " <<  g2->GetX()[i] << endl;
	}
     */
      x = g1->GetX()[i];
      y1 = g1->GetY()[i];
      y2 = g2->Eval(x);
      y = (y1*w1 + y2*w2)/(w1 + w2);
      g->SetPoint(i,x,y);
      cout << i << " " << x << " " << y << endl;
      y1Low = g1Low->GetY()[i];
      y2Low = g2Low->Eval(x);
      yL = (y1Low*w1 + y2Low*w2)/(w1 + w2);

      y1High = g1High->GetY()[i];
      y2High = g2High->Eval(x);
      yH = (y1High*w1 + y2High*w2)/(w1 + w2);
      
      g->SetPointEYlow(i,y - yL);
      g->SetPointEYhigh(i,yH - y);
    }
  cout << "---------------------" << endl;
  cout << "Number of Eff Area Points (Avg): " << g->GetN() << endl;
  cout << "---------------------" << endl;
  //  return(g);

}

TGraphAsymmErrors* AvgEACurve(TGraphAsymmErrors* g1,TGraphAsymmErrors* g2,double w1,double w2)
{
  TGraphAsymmErrors* g = new TGraphAsymmErrors();

  if( g1 == NULL && g2 == NULL) { return(NULL); }
  if( g1 == NULL) { return(g2); }
  if( g2 == NULL) { return(g1); }
  if( w1 == w2 ){ return(g1); }
  
  const int n1 = g1->GetN();
  const int n2 = g2->GetN();
  int n;
  if( n1 == 0 && n2 == 0 ){ return(g); }
  if( n1 == 0 ) { return(g2); }
  if( n2 == 0 ) { return(g1); }
  if( n1 != n2 ){ cout << "Warning uneven number of graph points..." << endl; }
  if( n1 <= n2 ){ n = n1; }
  
  if( n2 < n1 )
    { 
      n = n2;
      // Swap TGraph Objects, number in g1 should always be smaller     
      TGraphAsymmErrors* gTmp = g1;
      g1 = g2;
      g2 = gTmp;     
    }
  
  int diff = n2 - n;
  double* x1 = g1->GetX();
  double* x2 = g2->GetX();	
  double* y1 = g1->GetY();
  double* y2 = g2->GetY();
 
  double* y = new double[n];
  double* x = new double[n];
  
  double* yH = new double[n];
  double* yL = new double[n];
  double* y1H = new double[n];
  double* y1L = new double[n];
  double* y2H = new double[n];
  double* y2L = new double[n];

  int j = 0; // j is the number of times the graphs have been shifted to match each other
  
  for(int i=0; i<n; i++)
    {
      // cout << x1[i] << " " << x2[i] << endl;
      while(x1[i] != x2[i+j] )
	{ 
	  cout << "Warning! x positions not equal!" << endl;
	  cout << x1[i] << " " << x2[i+j] << endl;
	  j++;
	  if( j > diff ){ break; }
	} 
      x[i] = (x1[i] + x2[i])/2;
      y[i] = (y2[i+j]*w2 + y1[i]*w1)/(w2 + w1);
      g->SetPoint(i,x[i],y[i]);
      // doing errorbars:
      y1H[i] = y1[i] + g1->GetErrorYhigh(i); 
      y1L[i] = y1[i] - g1->GetErrorYlow(i); 
      y2H[i] = y2[i+j] + g1->GetErrorYhigh(i+j); 
      y2L[i] = y2[i+j] - g1->GetErrorYlow(i+j); 

      yH[i] = (y2H[i]*w2 + y1H[i]*w1)/(w2 + w1);
      yL[i] = (y2L[i]*w2 + y1L[i]*w1)/(w2 + w1);

      g->SetPointEYlow(i,0);
      g->SetPointEYhigh(i,0);
    }
	  
  /*
  for(int i=0; i<n; i++)
    {
      x[i] = x1[i];
      y2[i] = g2->Eval(x[i]);

      y1L[i] = y1[i] + g1->GetErrorYhigh(i); 
      y1L[i] = y1[i] - g1->GetErrorYlow(i); 
    }
  */
  cout << "Average done!" << endl;
  cout << n << " " << j << endl;
  return(g);

}

bool WriteTGraphToText(string outFileNm,TGraphAsymmErrors* g,int RunID = 0,float LT = 0)
{

  ofstream out("tmp.txt",ios_base::out);
  ofstream out2("tmp2.txt",ios_base::out);
 
  double E_logTeV,EA;
  double E_GeV;
  double dEA_low,dEA_high;

  string outFileNm2 = outFileNm;
  outFileNm2.append(".ch");
  
  for(int i=0; i<g->GetN(); i++)
    {
     
      g->GetPoint(i,E_logTeV,EA);
      E_GeV = pow(10.0,E_logTeV + 3.0);

      //if(EA == 0.0){ continue; }
      dEA_low = g->GetErrorYlow(i);
      dEA_high = g->GetErrorYhigh(i);
      out.precision(6);
      out << RunID << "   ";
      out << LT << "       ";
      out.precision(8);
      out << E_GeV << "   ";
      out << EA << "   ";
      out << -1.0*dEA_low << "   ";
      out << dEA_high << "   ";
      out << endl;
    }
 
  ostringstream os;
  os << "cat tmp.txt >>" << outFileNm;
  cout << RunID << " EA Curve Written to file" << endl;
  ostringstream os2;
  out2 << RunID << " " << g->GetN() << endl;
  os2 << "cat tmp2.txt >> " << outFileNm2;
  system(os2.str().c_str());
  return(system(os.str().c_str()));
  
}

float AvgVec(vector<float> vec)
{
  if( vec.size() == 0 ){ return(0); }

  float sum = 0;
  float avg;
  for(int i=0; i<vec.size(); i++)
    sum += vec.at(i);

  avg = sum/vec.size();
  return(avg);
}

//TGraphAsymmErrors* getEAInterpolatedCurve(int RunID, float El, float Az, float No, bool useMCE = 0, bool useKascade = 0,bool useHFit = 1)
TGraphAsymmErrors* getEAInterpolatedCurve(int RunID, float El, float Az, float No, bool useMCE = 0, string configPath = "/veritas/configfiles/vegas2_5/Soft")
{
  TGraphAsymmErrors* g = new TGraphAsymmErrors();
  ostringstream os;
  /*  
if(!useKascade && !useHFit)
    os << "bash ~/vegas/scripts/PickEATable.bash " << RunID << " /veritas/configfiles/vegas2_5/Soft > EATmp.txt";
  else if(useKascade && !useHFit)
    os << "bash PickEATableKas.bash " << RunID << " > EATmp.txt";
  else if(useHFit && !useKascade)
    os << "bash /veritas/scripts/PickEATable.bash " << RunID << " /veritas/configfiles/vegas2_5/Soft_hFit > EATmp.txt";
  */
  os << "bash ~/vegas/scripts/PickEATable.bash " << RunID << " " << configPath << " > EATmp.txt ";
  cout << os.str() << endl;
  if( system(os.str().c_str()) )
    {
      cerr << "Problem picking correct EA Table!" << endl;
      return;
    }
  ifstream ifs("EATmp.txt");
  string EAFileNm;
  ifs >> EAFileNm;
  cout << EAFileNm << endl;
 
  VARootIO* EAio = new VARootIO(EAFileNm,true);
  EAio->loadTheRootFile();
  
  if( !EAio->IsOpen() )
    {
      cerr << "Problem opening EA tables! " << endl;
      return;
    }

  vector<float> AzList,ZnList,OffsetList,NoList;
  getListOfEAParams(EAio,AzList,ZnList,OffsetList,NoList);
  
  cout << "Range of Azimuth: ";
  for(int i=0; i<AzList.size(); i++ ){ cout << AzList.at(i) << " "; }
  cout << endl;

  cout << "Range of Zenith: ";
  for(int i=0; i<ZnList.size(); i++ ){ cout << ZnList.at(i) << " "; }
  cout << endl;

  cout << "Range of Noise: ";
  for(int i=0; i<NoList.size(); i++ ){ cout << NoList.at(i) << " "; }
  cout << endl;

  if(OffsetList.size() != 0)
    {
      cout << "Range of Offsets: ";
      for(int i=0; i<OffsetList.size(); i++ ){ cout << OffsetList.at(i) << " "; }
      cout << endl;
    }
  
  float Zn = 90.0 - El;
  float AzL,AzH,ZnL,ZnH,NoL,NoH;
  getEAParamLimits(AzList,Az,AzL,AzH,1);
  getEAParamLimits(ZnList,Zn,ZnL,ZnH);
  getEAParamLimits(NoList,No,NoL,NoH);
  /*
  cout << "Azimuth " << Az << " is between " << AzL << " and " << AzH << endl;
  cout << "Zenith " << Zn << " is between " << ZnL << " and " << ZnH << endl;
  cout << "Noise " << No << " is between " << NoL << " and " << NoH << endl;
  */
  //if(!useKascade)
  //  {
  TGraphAsymmErrors* gEACurve_AzLZnLNoL = getEACurve(EAio,AzL,ZnL,NoL,useMCE);
  TGraphAsymmErrors* gEACurve_AzLZnLNoH = getEACurve(EAio,AzL,ZnL,NoH,useMCE);
  TGraphAsymmErrors* gEACurve_AzLZnHNoL = getEACurve(EAio,AzL,ZnH,NoL,useMCE);
  TGraphAsymmErrors* gEACurve_AzLZnHNoH = getEACurve(EAio,AzL,ZnH,NoH,useMCE);
  TGraphAsymmErrors* gEACurve_AzHZnLNoL = getEACurve(EAio,AzH,ZnL,NoL,useMCE);
  TGraphAsymmErrors* gEACurve_AzHZnLNoH = getEACurve(EAio,AzH,ZnL,NoH,useMCE);
  TGraphAsymmErrors* gEACurve_AzHZnHNoL = getEACurve(EAio,AzH,ZnH,NoL,useMCE);
  TGraphAsymmErrors* gEACurve_AzHZnHNoH = getEACurve(EAio,AzH,ZnH,NoH,useMCE);
      /*   }
      //else
      // {
      // Kascade EA tables use all offsets
      TGraphAsymmErrors* gEACurve_AzLZnLNoL = getOffsetEACurve(EAio,AzL,ZnL,NoL);
      TGraphAsymmErrors* gEACurve_AzLZnLNoH = getOffsetEACurve(EAio,AzL,ZnL,NoH);
      TGraphAsymmErrors* gEACurve_AzLZnHNoL = getOffsetEACurve(EAio,AzL,ZnH,NoL);
      TGraphAsymmErrors* gEACurve_AzLZnHNoH = getOffsetEACurve(EAio,AzL,ZnH,NoH);
      TGraphAsymmErrors* gEACurve_AzHZnLNoL = getOffsetEACurve(EAio,AzH,ZnL,NoL);
      TGraphAsymmErrors* gEACurve_AzHZnLNoH = getOffsetEACurve(EAio,AzH,ZnL,NoH);
      TGraphAsymmErrors* gEACurve_AzHZnHNoL = getOffsetEACurve(EAio,AzH,ZnH,NoL);
      TGraphAsymmErrors* gEACurve_AzHZnHNoH = getOffsetEACurve(EAio,AzH,ZnH,NoH);

    }
      */
  // Interp. using Cos(Zn) instead of Zn:
  
  ZnL = TMath::Cos(TMath::DegToRad()*ZnL);
  Zn  = TMath::Cos(TMath::DegToRad()*Zn);
  ZnH = TMath::Cos(TMath::DegToRad()*ZnH);
  
  // Using Guilliaume's nomenclature for interpolations: 
 
  TGraphAsymmErrors* gEACurve_a = InterpolateGraphs2(gEACurve_AzLZnLNoL,gEACurve_AzHZnLNoL,AzL,AzH,Az); // interp. Az w/ low Zn and low noise
  TGraphAsymmErrors* gEACurve_b = InterpolateGraphs2(gEACurve_AzLZnHNoL,gEACurve_AzHZnHNoL,AzL,AzH,Az); // interp. Az w/ high Zn and low noise
  // cout << "Got Here 1" << endl;
  TGraphAsymmErrors* gEACurve_c = InterpolateGraphs2(gEACurve_a,gEACurve_b,ZnL,ZnH,Zn); // interp. Zn w/ low noise
 
  TGraphAsymmErrors* gEACurve_d = InterpolateGraphs2(gEACurve_AzLZnLNoH,gEACurve_AzHZnLNoH,AzL,AzH,Az); // interp. Az w/ low Zn and high noise
  TGraphAsymmErrors* gEACurve_e = InterpolateGraphs2(gEACurve_AzLZnHNoH,gEACurve_AzHZnHNoH,AzL,AzH,Az); // interp. Az w/ high Zn and high noise
  //cout << "Got Here 2" << endl;
  TGraphAsymmErrors* gEACurve_f = InterpolateGraphs2(gEACurve_d,gEACurve_e,ZnL,ZnH,Zn); // interp. Zn w/ high noise
 
  TGraphAsymmErrors* gEACurve_final = InterpolateGraphs2(gEACurve_c,gEACurve_f,NoL,NoH,No); // interp. Noise
  //cout << "Got Final Interpolation" << endl;
  //gEACurve_AzLZnLNoL->Draw("A*");
  gEACurve_final->GetXaxis()->SetTitle("Energy ( log_{10}TeV )");
  gEACurve_final->GetYaxis()->SetTitle("Effective area ( m^{2} )");

  //gEACurve_final->Draw("AE*");
  /*
  gEACurve_AzLZnLNoL->Draw("*");
  gEACurve_AzHZnHNoH->Draw("*");
  gEACurve_AzLZnLNoL->SetMarkerColor(kBlue);
  gEACurve_AzHZnHNoH->SetMarkerColor(kRed);

  gEACurve_AzLZnLNoL->SetLineColor(kBlue);
  gEACurve_AzHZnHNoH->SetLineColor(kRed);
  */

  //c1->SetLogy(1);
  EAio->closeTheRootFile();
  return(gEACurve_final);

}


TGraphAsymmErrors* InterpolateGraphs2(TGraphAsymmErrors* g1, TGraphAsymmErrors* g2,double w1,double w2,double w)
{
  TGraphAsymmErrors* g = new TGraphAsymmErrors();

  TGraph* g2Low = new TGraph();
  TGraph* g2High = new TGraph();
  TGraph* g1Low = new TGraph();
  TGraph* g1High = new TGraph();

  if( g1 == NULL && g2 == NULL) { return(NULL); }
  if( g1 == NULL) { return(g2); }
  if( g2 == NULL) { return(g1); }
  if( w1 == w2 ){ return(g1); }
  if( w1 == w )
    { 
      cout << "No Interp. required: " << w1 << " = " << w << endl;
      return(g1); 
    }
  if( w2 == w )
    { 
      cout << "No Interp. required: " << w2 << " = " << w << endl;
      return(g2); 
    }

  cout << "g1 number of points: " << g1->GetN() << " g2 number of points: " << g2->GetN() << endl;

  const int n1 = g1->GetN();
  const int n2 = g2->GetN();
  if( n1 == 0 && n2 == 0 ){ return(g); }
  if( n1 == 0 ) { return(g2); }
  if( n2 == 0 ) { return(g1); }
  // looping to get TGraphs of errorbars:
  double* X1 = g1->GetX();
  double* Y1 = g1->GetY();
  double* X2 = g2->GetX();
  double* Y2 = g2->GetY();
  
  for(int i=0; i<n1; i++)
    {
      g1Low->SetPoint(i,X1[i],Y1[i] - g1->GetErrorYlow(i));
      g1High->SetPoint(i,X1[i],Y1[i] + g1->GetErrorYhigh(i));
    }
  for(int i=0; i<n2; i++)
    {
      g2Low->SetPoint(i,X2[i],Y2[i] - g2->GetErrorYlow(i));
      g2High->SetPoint(i,X2[i],Y2[i] + g2->GetErrorYhigh(i));
    }

  // getting errorbar arrays:
  double x,y,y1,y2;
  double y1Low,y2Low,y1High,y2High;
  double yL,yH;
  for(int i=0; i<n1; i++)
    {/*
      if( g1->GetX()[i] != g2->GetX()[i] )
	{
	  cout << "X positions not equal!" << endl;
	  cout << g1->GetX()[i] << " " <<  g2->GetX()[i] << endl;
	}
     */
      x = g1->GetX()[i];
      y1 = g1->GetY()[i];
      y2 = g2->Eval(x);
      y = ((y2 - y1)/(w2 - w1))*(w -  w1) + y1;
      g->SetPoint(i,x,y);

      y1Low = g1Low->GetY()[i];
      y2Low = g2Low->Eval(x);
      yL = ((y2Low - y1Low)/(w2 - w1))*(w -  w1) + y1Low;

      y1High = g1High->GetY()[i];
      y2High = g2High->Eval(x);
      yH = ((y2High - y1High)/(w2 - w1))*(w -  w1) + y1High;
      
      g->SetPointEYlow(i,y - yL);
      g->SetPointEYhigh(i,yH - y);
    }
  return(g);

}

TGraphAsymmErrors* InterpolateGraphs(TGraphAsymmErrors* g1, TGraphAsymmErrors* g2,double w1,double w2,double w)
{

  TGraphAsymmErrors* g = new TGraphAsymmErrors();

  if( g1 == NULL && g2 == NULL) { return(NULL); }
  if( g1 == NULL) { return(g2); }
  if( g2 == NULL) { return(g1); }
  if( w1 == w2 ){ return(g1); }
  if( w1 == w )
    { 
      cout << "No Interp. required: " << w1 << " = " << w << endl;
      return(g1); 
    }
  if( w2 == w )
    { 
      cout << "No Interp. required: " << w2 << " = " << w << endl;
      return(g2); 
    }
  
  const int n1 = g1->GetN();
  const int n2 = g2->GetN();
  int n;
  if( n1 == 0 && n2 == 0 ){ return(g); }
  if( n1 == 0 ) { return(g2); }
  if( n2 == 0 ) { return(g1); }
  if( n1 != n2 ){ cout << "Warning uneven number of graph points..." << endl; }
  if( n1 <= n2 ){ n = n1; }
  if( n2 < n1 )
    { 
      n = n2;
      // Swap TGraph Objects, number in g1 should always be smaller
      
      TGraphAsymmErrors* gTmp = g1;
      g1 = g2;
      g2 = gTmp;
      double wTmp = w1;
      w1 = w2;
      w2 = wTmp;
      
    }
  int diff = n2 - n;
  double* x1 = g1->GetX();
  double* x2 = g2->GetX();	
  double* y1 = g1->GetY();
  double* y2 = g2->GetY();
 
  double* y = new double[n];
  double* x = new double[n];
  
  double* yH = new double[n];
  double* yL = new double[n];
  double* y1H = new double[n];
  double* y1L = new double[n];
  double* y2H = new double[n];
  double* y2L = new double[n];

  int j = 0; // j is the number of times the graphs have been shifted to match each other
  for(int i=0; i<n; i++)
    {
      // cout << x1[i] << " " << x2[i] << endl;
      while(x1[i] != x2[i+j] )
	{ 
	  cout << "Warning! x positions not equal!" << endl;
	  cout << x1[i] << " " << x2[i+j] << endl;
	  j++;
	  if( j > diff ){ break; }
	} 
      x[i] = (x1[i] + x2[i])/2;
      y[i] = ((y2[i+j] - y1[i])/(w2 - w1))*(w - w1) + y1[i];
      g->SetPoint(i,x[i],y[i]);
      // doing errorbars:
      y1H[i] = y1[i] + g1->GetErrorYhigh(i); 
      y1L[i] = y1[i] - g1->GetErrorYlow(i); 
      y2H[i] = y2[i+j] + g1->GetErrorYhigh(i+j); 
      y2L[i] = y2[i+j] - g1->GetErrorYlow(i+j); 

      yH[i] = ((y2H[i] - y1H[i])/(w2 - w1))*(w - w1) + y1H[i];
      yL[i] = ((y2L[i] - y1L[i])/(w2 - w1))*(w - w1) + y1L[i];

      g->SetPointEYlow(i,y[i] - yL[i]);
      g->SetPointEYhigh(i,yH[i] - y[i]);
    }	  

  cout << "Interpolation done!" << endl;
  cout << n << " " << j << endl;
  return(g);
}

TGraphAsymmErrors* getEACurve(VARootIO* EAio,float Az,float Zn,float No,bool useMCE = 0)
{
  if(useMCE)
    cout << "Using EA's generated using MC energies..." << endl;

  ostringstream os;
  TGraphAsymmErrors* g = new TGraphAsymmErrors();
  VAEffectiveArea* EffArea;
  // for Nov2010 sims it seems that for Zn = 0, simulations fixed at Az = 0. BJZ hack 1/15/2013.
  if( Zn == 0.0 ){ Az = 0.0; }

  os << "EffectiveArea_Azimuth_" << (int)Az << "_Zenith_" << (int)Zn << "_Noise_" << No;
  cout << os.str() << endl;
  EffArea = (VAEffectiveArea*)EAio->loadAnObject(os.str(),"effective_areas",true);
  //if( EffArea == NULL && Zn = 0.0 ){ return(getEACurve(EAio,0,0,No)); }
  if( EffArea == NULL )
    {
      cout << "Problem getting effective area, checking aboslute offset..." << endl;
      //os << "_AbsoluteOffset_0.5";
      // BJZ hack for Galactic Center:
      os << "_AbsoluteOffset_0.75";
      cout << os.str() << endl;
      EffArea = (VAEffectiveArea*)EAio->loadAnObject(os.str(),"effective_areas",true);
      if( EffArea == NULL )
	{
	  cout << "Problem getting effective area for this event!" << endl;
	  ofstream out("missingEA.txt",std::ofstream::app);
	  out << os.str() << " " << endl;
	  return(NULL);
	}
    }
  if(useMCE)
    g = EffArea->pfEffArea_MC; // MC energies
  else
    g = EffArea->pfEffArea_Rec; // reconstructed energies
 
  if( g == NULL )
    {
      cout << "Problem getting TGraph for this event!" << endl;
      ifstream out("missingEA.txt",std::ifstream::app);
      out << os.str() << " " << endl;
      return(NULL);
    }
  cout << "Number of Eff Area TGraph points: " << g->GetN() << endl;
  //return(g);
  // BJZ adding code to get number of thrown/cut events:
  /*
  TH1F* hEnergyThrown = EffArea->pfEnergyMC_thrown;
  TH1F* hEnergyMCCut  = EffArea->pfEnergyMC_afterCut;
  for(int i=0; i<hEnergyMCCut->GetNbinsX(); i++)
    {
      //      cout << "Bin Center Thrown: " << hEnergyThrown->GetBinCenter(i) << endl;
      cout << "Bin Center Cut: " << 10.0**(hEnergyMCCut->GetBinCenter(i)) << endl;
      cout << "Number Thrown: " << hEnergyThrown->GetBinContent(i) << endl;
      cout << "Number after cuts: " << hEnergyMCCut->GetBinContent(i) << endl;
    }
  */
  return(g);
}

TGraphAsymmErrors* getOffsetEACurve(VARootIO* EAio,float Az,float Zn,float No,float Offset = 0.5)
{
  ostringstream os;
  TGraphAsymmErrors* g = new TGraphAsymmErrors();
  VAEffectiveArea* EffArea;
  os << "EffectiveArea_Azimuth_" << (int)Az << "_Zenith_" << (int)Zn << "_Noise_" << No << "_AbsoluteOffset_"<< Offset;
  cout << os.str() << endl;
  EffArea = (VAEffectiveArea*)EAio->loadAnObject(os.str(),"effective_areas",true);
  if( EffArea == NULL )
    {
      cout << "Problem getting effective area for this event!" << endl;
      return(NULL);
    }
  g = EffArea->pfEffArea_Rec;
  if( g == NULL )
    {
      cout << "Problem getting TGraph for this event!" << endl;
      return(NULL);
    }
  cout << "Number of Eff Area TGraph points: " << g->GetN() << endl;
  return(g);
}


void getEAParamLimits(vector<float> vec,float realVal,float &lowerVal,float &upperVal,bool wrap = 0)
{
  float diff;
  float min_diff1 = 1000;
  float min_diff2 = 1000;
  float diff1,diff2;

  // getting lower limit
  for(int i=0; i<vec.size(); i++)
    {
      diff = realVal - vec.at(i);
      if( diff < min_diff1 && diff > 0 )
	{
	  lowerVal = vec.at(i);
	  min_diff1 = diff;
	}
    }
  // getting upper limit
  for(int i=0; i<vec.size(); i++)
    {
      diff = vec.at(i) - realVal;
      if( diff <= min_diff2 && diff >= 0 )
	{
	  upperVal = vec.at(i);
	  min_diff2 = diff;
	}
    }
  /*
  for(int i=0; i<vec.size(); i++)
    {
      for(int j=0; j<i; j++)
	{
	  //diff = TMath::Abs(realVal - vec.at(i)) + TMath::Abs(realVal - vec.at(j)); 
	  diff2 = realVal - vec.at(j);
	  diff1 =  vec.at(i) - realVal;
	  if( diff1 < min_diff1 && diff2 < min_diff2 && diff1 > 0.0 & diff2 > 0.0)
	    {
	      min_diff1 = diff1;
	      min_diff2 = diff2;
	      
	      lowerVal = vec.at(j);
	      upperVal = vec.at(i);
	    }
	}
    }
  */
  //if( lowerVal <= realVal && upperVal <= realVal && wrap==false )
  //   lowerVal = upperVal;
  // if( realVal <= lowerVal && realVal <= upperVal && wrap==false )
  //   upperVal = lowerVal;
  if( lowerVal < realVal && upperVal < realVal && wrap==true )
    upperVal = vec.front();
  if( realVal < lowerVal && realVal < upperVal && wrap==true )
    lowerVal = vec.back();

}

void getListOfEAParams(VARootIO* effEAFile,vector<float> &AzList,vector<float> &ZnList,vector<float> &OffsetList,vector<float> &NoiseList)
{
 
  ZnList.clear();
  AzList.clear();
  OffsetList.clear();
  NoiseList.clear();
  float Value;
  char Name[50];
  
  //gDirectory->cd("effective_areas");
  TTree* ParTree = (TTree*)gDirectory->Get("effective_areas/paramTree");
  if( ParTree == NULL )
    {
      cerr << "Problem with parameter tree! " << endl;
      return;
    }

  ParTree->SetBranchAddress("value",&Value);
  ParTree->SetBranchAddress("name",&Name);
  cout << "Number of entires in parameter tree: " << ParTree->GetEntries() << endl;

  for(int i=0; i< ParTree->GetEntries(); i++)
    {
      ParTree->GetEntry(i);
      //cout << Name << endl;
      if( !strcmp(Name,"Azimuth") ){ AzList.push_back(Value); }
      if( !strcmp(Name,"Zenith") ){ ZnList.push_back(Value); }
      if( !strcmp(Name,"Noise") ){ NoiseList.push_back(Value); }
      if( !strcmp(Name,"AbsoluteOffset") ){ OffsetList.push_back(Value); }
    }
  cout << "Number of Azimuths: " << AzList.size() << endl;
  cout << "Number of Zeniths: " << ZnList.size() << endl;
  cout << "Number of Noise levels: " << NoiseList.size() << endl;
  cout << "Number of Absolute Offsets: " << OffsetList.size() << endl;
}

