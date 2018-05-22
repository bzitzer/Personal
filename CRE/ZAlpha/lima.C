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
