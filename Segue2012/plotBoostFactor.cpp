void plotBoostFactor(string inFile)
{
  // Recreate figure 6 of the Segue 1 paper...
  TFile* f = new TFile(inFile.c_str(),"READ");
  TGraph* g_uu = gDirectory->Get("g_uu");
  if( g_uu == NULL)
    {
      cout << "Issue with TGraph!" << endl;
      return;
    }
  TGraph* g_bf = new TGraph();
  double bf,m,sigv;
  for(int i=0; i<g_uu->GetN(); i++)
    {
      g_uu->GetPoint(i,m,sigv);
      bf = sigv/3.0e-26;
      g_bf->SetPoint(i,m/1000,bf);
    }
  TCanvas* c1 = new TCanvas("c1","c1",50,50,700,500);
  c1->SetLogy(1);
  g_bf->GetXaxis()->SetLimits(0,4);
  g_bf->GetXaxis()->SetTitle("m_{#chi} [TeV]");
  g_bf->GetYaxis()->SetLimits(100,4e4);
  g_bf->GetYaxis()->SetTitle("B_{F}");
  g_bf->Draw("AL");

}
