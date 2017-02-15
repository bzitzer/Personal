
void dumpBiasTextInput(string InFile,string dir,string configDir,bool testMode = 0,bool useOlder = 0)
{
  ifstream ifs(InFile.c_str());
  if(!ifs.is_open())
    {
      cerr << "Problem with input file! " << endl;
      return;
    }

  ostringstream OutFileNm;
  OutFileNm << InFile << ".mat";
  cout << OutFileNm.str() << endl;

  ofstream out(OutFileNm.str().c_str());

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

  double TotLT = 0;
  int nLines = 0;
  int nRuns = 0;
  TH2F* hMat = new TH2F();
  while(!ifs.eof())
    {
      
      nLines++;
      if( testMode && nLines > 2e3){ break; }
      if(!useOlder)
	ifs >> RunID >> LT >> EventTime >> RA >> Dec >> RA_tr >> Dec_tr >> Energy >> IsOn >> W >> El >> Az >> N >> Offset;
      else
	ifs >> RunID >> LT >> EventTime >> RA >> Dec >> Energy >> IsOn >> W >> El >> Az >> N >> Offset;
    
      if(RunID != RunID_old && RunID_old !=0 )
	{	
	  nRuns++;
	
	  cout << "Finished reading run number " << RunID_old << endl;
	 
	  if( El_vec.size() == 0 || Az_vec.size() == 0 || N_vec.size() == 0 )
	    {
	      cout << "Run with no ON events! " << endl;
	      El_m = AvgVec(El_vec_bg);
              Az_m = AvgVec(Az_vec_bg);
              N_m = AvgVec(N_vec_bg);

	    }
	  else
	    {
	      El_m = AvgVec(El_vec);
	      Az_m = AvgVec(Az_vec);
	      N_m = AvgVec(N_vec);	  
	    }
	  cout << El_m << " " << Az_m << " " << N_m << endl;
	  El_vec.clear();
	  Az_vec.clear();
	  N_vec.clear();
	  El_vec_bg.clear();
	  Az_vec_bg.clear();
	  N_vec_bg.clear();
          
	  hMat = getInterpolatedBiasMatrix(RunID_old,El_m,Az_m,N_m,configDir);
	  writeBiasMatrix2(hMat,RunID_old,dir);
	    
	  if(testMode)
	    //	      makeEnergyBiasCurve(hMat);
	  //  cout << "Run Number: " <<  RunID_old << " Number of points: " << g->GetN() <<" Number of average points: " << gAvg->GetN() << endl;
	  cout << nRuns << endl;
	
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

      TotLT += LT;
      hMat = getInterpolatedBiasMatrix(RunID,El_m,Az_m,N_m,configDir);
      writeBiasMatrix2(hMat,RunID,dir);
      /*
      if(!didWrite)
	{
	  cerr << "Problem Writting output file!" << endl;
	}
      */
    }  
}

float AvgVec(vector<float> vec)
{
  float sum = 0;
  for(int i=0; i<vec.size(); i++)
    sum += vec.at(i);

  return(sum/vec.size());
}


TGraphErrors* makeEnergyResPlot(TH2F* h,double EBiasThresh=0.1)
{
  // Sanity check for me...
 
  double E_tr,E_rec;
  double w;
  TFitResultPtr fitPtr;
  double* fitPar;
  double* fitErr;
  double mean,sigma;
  TGraphErrors* g = new TGraphErrors();
  int k=0;
  for(int i=0; i<h->GetNbinsX(); i++)
    {
      E_tr = 10**(h->GetXaxis()->GetBinCenter(i));
      TH1F* hBiasEtr = new TH1F("h","",200,-1,1);
      for(int j=0; j<h->GetNbinsY(); j++)
	{
	  w = h->GetBinContent(i,j);
	  E_rec = 10**(h->GetYaxis()->GetBinCenter(j));
	  hBiasEtr->Fill((E_rec - E_tr)/E_tr,w);
	}
      if(hBiasEtr->GetSumOfWeights() <= 0){ continue; }
      /*
      fitPtr = hBiasEtr->Fit("gaus","SN E LL 0");
      fitPar = fitPtr.Get()->GetParams();
      fitErr = fitPtr.Get()->GetErrors();
      mean = fitPar[1];
      sigma = fitPar[2];
      */
      g->SetPoint(k,TMath::Log10(E_tr),hBiasEtr->GetRMS());
      g->SetPointError(k,0,hBiasEtr->GetRMSError());
      //cout << sigma << " " << fitErr[2] << endl;
      hBiasEtr->Delete();
      k++;
    }
  
  /*
    TCanvas* c1 = new TCanvas("EBias","Energy Bias",50,50,600,600);
    g->GetYaxis()->SetRangeUser(-1.8,1.8);
    g->GetXaxis()->SetRangeUser(-1.2,2.0);
  
    g->Draw("AE*");
    TLine* l = new TLine(-1.5,0,2,0);
    l->SetLineColor(kRed);
    l->SetLineStyle(5);
    l->Draw("same");
  
    g->GetYaxis()->SetTitle("(E_{rec} - E_{tr})/E_{tr}");
    g->GetXaxis()->SetTitle("True Energy [log_{10} TeV]");
    TF1* f = new TF1("","[0]*x + [1]",-1.2,-0.8);
    f->SetLineColor(kBlue);
    g->Fit(f,"R");
    for(int i=0; i<g->GetN(); i++)
    {
    if(g->GetY()[i] >= EBiasThresh )
    cout << 10.0**(g->GetX()[i] + 3.0) << " " << g->GetY()[i] << endl;
    }
  */
  g->GetXaxis()->SetTitle("True Energy [log_{10} TeV]");
  g->GetYaxis()->SetTitle("Energy Resolution");
  return(g);
}


TGraphErrors* makeEnergyBiasPlot(TH2F* h,double EBiasThresh=0.1)
{
  // Sanity check for me...
 
  double E_tr,E_rec;
  double w;
  TFitResultPtr fitPtr;
  double* fitPar;
  double* fitErr;
  double mean,sigma;
  TGraphErrors* g = new TGraphErrors();

  int k = 0;
  for(int i=0; i<h->GetNbinsX(); i++)
    {
      E_tr = 10**(h->GetXaxis()->GetBinCenter(i));
      TH1F* hBiasEtr = new TH1F("h","",20,-1,1);
      for(int j=0; j<h->GetNbinsY(); j++)
	{
	  w = h->GetBinContent(i,j);
	  E_rec = 10**(h->GetYaxis()->GetBinCenter(j));
	  hBiasEtr->Fill((E_rec - E_tr)/E_tr,w);
	}
      if(hBiasEtr->GetSumOfWeights() <= 0){ continue; }
      /*
      fitPtr = hBiasEtr->Fit("gaus","E N");
      fitPar = fitPtr.Get()->GetParams();
      fitErr = fitPtr.Get()->GetErrors();
      mean = fitPar[1];
      sigma = fitPar[2];
      */
      g->SetPoint(k,TMath::Log10(E_tr),hBiasEtr->GetMean());
      g->SetPointError(k,0,hBiasEtr->GetMeanError());
      hBiasEtr->Delete();
      k++;
    }
  g->GetXaxis()->SetTitle("E_{tr} [log_{10} TeV]");
  g->GetYaxis()->SetTitle("(E_{rec} - E_{tr})/E_{tr}");

  /*
  TCanvas* c1 = new TCanvas("EBias","Energy Bias",50,50,600,600);
   g->GetYaxis()->SetRangeUser(-1.8,1.8);
  g->GetXaxis()->SetRangeUser(-1.2,2.0);
  
  g->Draw("AE*");
  TLine* l = new TLine(-1.5,0,2,0);
  l->SetLineColor(kRed);
  l->SetLineStyle(5);
  l->Draw("same");
  
  g->GetYaxis()->SetTitle("(E_{rec} - E_{tr})/E_{tr}");
  g->GetXaxis()->SetTitle("True Energy [log_{10} TeV]");
  TF1* f = new TF1("","[0]*x + [1]",-1.2,-0.8);
  f->SetLineColor(kBlue);
  g->Fit(f,"R");
  for(int i=0; i<g->GetN(); i++)
    {
      if(g->GetY()[i] >= EBiasThresh )
	cout << 10.0**(g->GetX()[i] + 3.0) << " " << g->GetY()[i] << endl;
    }
  */
  return(g);
}


TGraphErrors* makeEnergyBiasCurve(TH2F* h,double EBiasThresh=0.1)
{
  // Sanity check for me...
 
  double E_tr,E_rec;
  double w;
  TFitResultPtr fitPtr;
  double* fitPar;
  double mean,sigma;
  TGraphErrors* g = new TGraphErrors();
  for(int i=0; i<h->GetNbinsX(); i++)
    {
      E_tr = 10**(h->GetXaxis()->GetBinCenter(i));
      TH1F* hBiasEtr = new TH1F("h","",100,-1,1);
      for(int j=0; j<h->GetNbinsY(); j++)
	{
	  w = h->GetBinContent(i,j);
	  E_rec = 10**(h->GetYaxis()->GetBinCenter(j));
	  hBiasEtr->Fill((E_rec - E_tr)/E_tr,w);
	}
      if(hBiasEtr->GetSumOfWeights() <= 0){ continue; }
      fitPtr = hBiasEtr->Fit("gaus","SN LL Q");
      fitPar = fitPtr.Get()->GetParams();
      mean = fitPar[1];
      sigma = fitPar[2];
      g->SetPoint(i,TMath::Log10(E_tr),mean);
      g->SetPointError(i,0,sigma);
      hBiasEtr->Delete();
    }
  TCanvas* c1 = new TCanvas("EBias","Energy Bias",50,50,600,600);
   g->GetYaxis()->SetRangeUser(-1.8,1.8);
  g->GetXaxis()->SetRangeUser(-1.2,2.0);
  /*
  g->Draw("AE*");
  TLine* l = new TLine(-1.5,0,2,0);
  l->SetLineColor(kRed);
  l->SetLineStyle(5);
  l->Draw("same");
  */
  g->GetYaxis()->SetTitle("(E_{rec} - E_{tr})/E_{tr}");
  g->GetXaxis()->SetTitle("True Energy [log_{10} TeV]");
  TF1* f = new TF1("","[0]*x + [1]",-1.2,-0.8);
  f->SetLineColor(kBlue);
  g->Fit(f,"R");
  for(int i=0; i<g->GetN(); i++)
    {
      if(g->GetY()[i] >= EBiasThresh )
	cout << 10.0**(g->GetX()[i] + 3.0) << " " << g->GetY()[i] << endl;
    }
  return(g);
}




void writeBiasMatrix2(TH2F* h, int runID,string dir = "Pass4d/")
{
  if(h==NULL)
    {
      cout << "No Histogram, no output! Those are the rules!" << endl;
      return;
    }

  ostringstream os1;
  ostringstream os2;
  ostringstream os3;

  os1 << dir << runID << ".EtrBin.txt";
  os2 << dir << runID << ".ErecBin.txt";
  os3 << dir << runID << ".mat.txt";

  ofstream out1(os1.str().c_str());
  ofstream out2(os2.str().c_str());
  ofstream out(os3.str().c_str());
  
  out.width(8);
  out.precision(8);
  out << scientific;
  
  writeBinEdges(h,out1,1);
  writeBinEdges(h,out2,0);

  out << "P(E_rec|E_tr) For run number: " << runID;
  out << " (rows are true energy bins, cols are recon. energy bins)" << endl;
  out << "-----------------------------------------------------------------" << endl;

  float norm = 0;
  float P;
  
  
  for(int i=0; i<h->GetNbinsX(); i++)
    {
      // cout << "got here!" << endl;
      for(int j=0; j<h->GetNbinsY(); j++)
	norm += h->GetBinContent(i,j);

      for(int j=0; j<h->GetNbinsY(); j++)
	if( norm == 0 )
	  out << 0.0 << "   ";
	else
	  {
	    P = h->GetBinContent(i,j)/norm;
	    out << P << "   ";
	  }
      out << endl;
      norm = 0;
    }
  /*
  ostringstream com;
  com << "cat " << os3.str() << endl;
  system(com.str().c_str());
  */
}
void writeBinEdges(TH2F* h, const ostream& out = cout, bool IsTrue)
{
  out.width(6);
  out.precision(6);
  out.fixed;
  out << "Bin Edges for ";
  if(IsTrue)
    out << "True";
  else
    out << "Reconstructed";

  out << " Energies" << endl;
  out << "----------------------------------------" << endl;
 
  TAxis* tAx = new TAxis();
  if(IsTrue)
    tAx = h->GetXaxis();
  else
    tAx = h->GetYaxis();

  for(int i=0; i<tAx->GetNbins(); i++)
    out << pow(10.0,tAx->GetBinLowEdge(i) + 3.0) << endl;

  out << pow(10.0,tAx->GetBinUpEdge(i)+3) << endl;


}


void writeBiasMatrix(TH2F* h, const ostream& out = cout,int runID = 56666)
{
  out << "got here" << endl;
  double ETrueLow,ETrueMid,ETrueHigh;
  double ERecLow,ERecMid,ERecHigh;

  double norm = 0;
  double P;
  out.precision(8);
  for(int i=0; i<h->GetNbinsX(); i++)
    {
      if(i==0)
	{
	  out << "Run# | Row# | E_tr_low | E_tr_mid | E_tr_high | E_rec_low | E_rec_mid | E_rec_high | P(E_rec|E_tr)";
	  // for(int j=0; j<h->GetNbinsY(); j++)
	  //  out << 10.0**(h->GetYaxis()->GetBinCenter(j)+3.0) << "   ";
	  
	  out << endl;
	  out << "--------------------------------------------------------------------------------------------------" << endl;
	}

      for(int j=0; j<h->GetNbinsY(); j++)
	norm += h->GetBinContent(i,j);

      if(norm == 0 ){ continue; }

      ETrueLow = 10.0**(h->GetXaxis()->GetBinLowEdge(i) + 3.0);
      ETrueHigh = 10.0**(h->GetXaxis()->GetBinUpEdge(i) + 3.0);
      ETrueMid = 10.0**(h->GetXaxis()->GetBinCenter(i) + 3.0);     
         
      for(int j=0; j<h->GetNbinsY(); j++)
	{
	  P = h->GetBinContent(i,j)/norm;
	  if(P==0){ continue; }
	  out << runID << "   ";
	  out << i << "    "; 
	  out << ETrueLow << "    ";
	  out << ETrueMid << "    ";
	  out << ETrueHigh << "    ";

	  ERecLow = 10.0**(h->GetXaxis()->GetBinLowEdge(j) + 3.0);
	  ERecHigh = 10.0**(h->GetXaxis()->GetBinUpEdge(j) + 3.0);
	  ERecMid = 10.0**(h->GetXaxis()->GetBinCenter(j) + 3.0); 

	  out << ERecLow << "   ";
	  out << ERecMid << "   ";
	  out << ERecHigh << "   ";
	  
	  out << P << " ";
	  out << endl;
	}
     
      out << endl;
      norm = 0;
    }

}



TH2F* getInterpolatedBiasMatrix(int RunID,float El,float Az,float No,string configDir,bool useOffsets = false)
{
 ostringstream os;
  if(!useOffsets)
    os << "bash ~/vegas/scripts/PickEATable.bash " << RunID << " " << configDir << " > EATmp.txt";
  else
    os << "bash PickEATableKas.bash " << RunID << " > EATmp.txt";

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

  float Zn = 90.0 - El;
  float AzL,AzH,ZnL,ZnH,NoL,NoH;
  getEAParamLimits(AzList,Az,AzL,AzH,1);
  getEAParamLimits(ZnList,Zn,ZnL,ZnH);
  getEAParamLimits(NoList,No,NoL,NoH);

  cout << "Azimuth " << Az << " is between " << AzL << " and " << AzH << endl;
  cout << "Zenith " << Zn << " is between " << ZnL << " and " << ZnH << endl;
  cout << "Noise " << No << " is between " << NoL << " and " << NoH << endl;
  //  return(NULL);
  if(!useOffsets)
    {
      TH2F* hBiasMat_AzLZnLNoL = getBiasMatrix(EAio,AzL,ZnL,NoL);
      TH2F* hBiasMat_AzLZnLNoH = getBiasMatrix(EAio,AzL,ZnL,NoH);
      TH2F* hBiasMat_AzLZnHNoL = getBiasMatrix(EAio,AzL,ZnH,NoL);
      TH2F* hBiasMat_AzLZnHNoH = getBiasMatrix(EAio,AzL,ZnH,NoH);

      TH2F* hBiasMat_AzHZnLNoL = getBiasMatrix(EAio,AzH,ZnL,NoL);
      TH2F* hBiasMat_AzHZnLNoH = getBiasMatrix(EAio,AzH,ZnL,NoH);
      TH2F* hBiasMat_AzHZnHNoL = getBiasMatrix(EAio,AzH,ZnH,NoL);
      TH2F* hBiasMat_AzHZnHNoH = getBiasMatrix(EAio,AzH,ZnH,NoH);
    }
  else
    return(NULL); // not yet implimented. Maybe it never will...


  //TCanvas* c1 = new TCanvas("c1","",50,50,600,600);
 
  TH2F* h = interpBiasMatrix(hBiasMat_AzLZnLNoL,hBiasMat_AzLZnLNoH,1,2,5);

    // Using Guilliaume's nomenclature for interpolations: 
  TH2F* hBiasMat_a = interpBiasMatrix(hBiasMat_AzLZnLNoL,hBiasMat_AzHZnLNoL,AzL,AzH,Az); // interp. Az w/ low Zn and low noise
  TH2F* hBiasMat_b = interpBiasMatrix(hBiasMat_AzLZnHNoL,hBiasMat_AzHZnHNoL,AzL,AzH,Az); // interp. Az w/ high Zn and low noise
  // cout << "Got Here 1" << endl;
  TH2F* hBiasMat_c = interpBiasMatrix(hBiasMat_a,hBiasMat_b,ZnL,ZnH,Zn); // interp. Zn w/ low noise
 
  TH2F* hBiasMat_d = interpBiasMatrix(hBiasMat_AzLZnLNoH,hBiasMat_AzHZnLNoH,AzL,AzH,Az); // interp. Az w/ low Zn and high noise
  TH2F* hBiasMat_e = interpBiasMatrix(hBiasMat_AzLZnHNoH,hBiasMat_AzHZnHNoH,AzL,AzH,Az); // interp. Az w/ high Zn and high noise
  //cout << "Got Here 2" << endl;
  TH2F* hBiasMat_f = interpBiasMatrix(hBiasMat_d,hBiasMat_e,ZnL,ZnH,Zn); // interp. Zn w/ high noise
 
  TH2F* hBiasMat_final = interpBiasMatrix(hBiasMat_c,hBiasMat_f,NoL,NoH,No); // interp. Noise
  //  hBiasMat_final->Draw("colz");
  return(hBiasMat_final);

}
TH2F* interpBiasMatrix(TH2F* h1,TH2F* h2,float w1,float w2,float w)
{
  if(h1==NULL){return(h2);}
  if(h2==NULL){return(h1);}
  if(TMath::IsNaN(h1->GetSumOfWeights())==1)
    {
      cout << "Warning! Histogram possibly not filled!" << endl;
      return(h2);
    }
  if(TMath::IsNaN(h2->GetSumOfWeights())==1)
    {
      cout << "Warning! Histogram possibly not filled!" << endl;
      return(h1);
    }

  //cout << TMath::IsNaN(h1->GetSumOfWeights()) << endl;
  //cout << TMath::IsNaN(h2->GetSumOfWeights()) << endl;
  if(w == w1 )
    {
      cout << "No Interp requried: " << w << " = " << w1 << endl;
      return(h1);
    }
  if(w == w2 )
    {
      cout << "No Interp requried: " << w << " = " << w2 << endl;
      return(h2);
    }
  if( w1 == w2 )
    {
      cout << "Weights are equal, perhaps at bounds for that value?" << endl;
      return(h1);
    }
  int BinNum;

  TH2F* h = (TH2F*)h1->Clone();
  double P_tmp;
  double P1,P2;
  for(int i=0; i<h1->GetNbinsX(); i++)
    for(int j=0; j<h1->GetNbinsY(); j++)
      {
	//cout << h1->GetXaxis()->GetBinCenter(i) << " " <<  h2->GetXaxis()->GetBinCenter(i) << endl;
	//cout << h1->GetYaxis()->GetBinCenter(j) << " " <<  h2->GetYaxis()->GetBinCenter(j) << endl;
	if(h1->GetXaxis()->GetBinCenter(i) != h2->GetXaxis()->GetBinCenter(i))
	  cout << "Warning! Bins don't match!" << endl;
	if(h1->GetYaxis()->GetBinCenter(j) != h2->GetYaxis()->GetBinCenter(j))
	  cout << "Warning! Bins don't match!" << endl;
	P1 = h1->GetBinContent(i,j);
	P2 = h2->GetBinContent(i,j);
	P_tmp = ((P2 - P1)/(w2 - w1))*(w - w1) + P1;
	h->SetBinContent(i,j,P_tmp);
      }

  if( h->GetSumOfWeights() != 1.0 )
    cout << "Warning! Normalization problem! N=" << h->GetSumOfWeights() << endl;

  return(h);
}


TH2F* getBiasMatrix(string inFile,float Az,float Zn,float No)
{
  VARootIO* EAio = new VARootIO(inFile.c_str(),true);
  EAio->loadTheRootFile();

  if( !EAio->IsOpen() )
    {
      cerr << "Problem opening EA tables! " << endl;
      return(NULL);
    }

  return(getBiasMatrix(EAio,Az,Zn,No));
  EAio->closeTheRootFile();
}

TH2F* getBiasMatrix(VARootIO* EAio,float Az,float Zn,float No)
{
  ostringstream os;
  TH2F* hMatrix;
  VAEffectiveArea* EffArea;
  // for Nov2010 sims it seems that for Zn = 0, simulations fixed at Az = 0. BJZ hack 1/15/2013.
  if( Zn == 0.0 ){ Az = 0.0; }

  os << "EffectiveArea_Azimuth_" << (int)Az << "_Zenith_" << (int)Zn << "_Noise_" << No;
  cout << os.str() << endl;
  EffArea = (VAEffectiveArea*)EAio->loadAnObject(os.str(),"effective_areas",true);
  if( EffArea == NULL )
    {
      cout << "Checking Absolute Offset..." << endl;
      os << "_AbsoluteOffset_0.5";
      EffArea = (VAEffectiveArea*)EAio->loadAnObject(os.str(),"effective_areas",true);
      if( EffArea == NULL )
	{
	  cout << "Problem getting effective area object for this event!" << endl;
	  return(NULL);
	}
    }
  hMatrix = EffArea->pfEnergy_Rec_VS_MC_2D;
  //hMatrix = EffArea->pfEnergyBias_2D; // If I want energy bias information instead...
  if( hMatrix == NULL )
    {
      cout << "Problem getting TH2F for this event!" << endl;
      return(NULL);
    }
  //cout << "Number of Eff Area TGraph points: " << g->GetN() << endl;
  int N = hMatrix->GetEntries();
  if( N == 0 )
    {
      cout << "Warning! Histogram not filled!" << endl;

    }
  hMatrix->Scale(1.0/N);
  //cout << hMatrix->GetNbinsX() << " " << hMatrix->GetNbinsY() << endl;
  return(hMatrix);
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

  if( lowerVal < realVal && upperVal < realVal && wrap==true )
    upperVal = vec.front();
  if( realVal < lowerVal && realVal < upperVal && wrap==true )
    lowerVal = vec.back();
 
  if( lowerVal < realVal && upperVal < realVal && wrap==false )
    upperVal = vec.back();
  if( realVal < lowerVal && realVal < upperVal && wrap==false )
    lowerVal = vec.front();

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
