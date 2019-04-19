double lima_with_syst(double on, double off, double alpha0,double sigma_alpha)
{
  
  double a = 1;
  double b = 1.0 - alpha0;
  double c = sigma_alpha*sigma_alpha*off - alpha0;
  double d = -1.0*on*sigma_alpha*sigma_alpha;

  vector<double> r = solve_cubic(a,b,c,d);

  double alpha = 0;
  double l_max = 0;
  double l;
  for(int i=0; i<r.size(); i++)
    {
      if( r.at(i) < 0.0 )
	{
	  continue;
	}
      l = on*TMath::Log(r.at(i)) + (on + off)*TMath::Log((on+off)/(r.at(i)+1))
	- 0.5*TMath::Power((r.at(i) - alpha0)/sigma_alpha,2.0);
      //cout << l << " " << r.at(i) << endl;
      if( l > l_max )
	{
	  alpha = r.at(i);
	}
    }
  if( alpha <= 0.0 )
    {
      cout << "Warning! Negative alpha! on="<< on << " off= " << " alpha= "<< alpha <<endl;
      return(0);
    }
  double diff = on - alpha0 * off;
	
  if(on > 0 && off > 0)
    {
      double lima17 = TMath::Sqrt(2) *
	TMath::Power(
		     on * TMath::Log(((1.0 + alpha) / alpha) * (on / (on + off)))
		     +off * TMath::Log((1.0 + alpha) * (off / (on + off))), 0.5);
      lima17 = sqrt(lima17*lima17 + TMath::Power((alpha - alpha0)/sigma_alpha,2.0));
    }
  else
    {
      cout << "warning " << on << " " << off << endl;
      return 0;
    }
  
  if(diff < 0)
    {
      return -lima17;
    }
  else
    {
      return lima17;
    }
}


vector<double> solve_cubic(double a, double b, double c, double d)
{
  vector<double> r;
  if( a == 0)
    {
      return(r);
    }
  double p = (3*a*c - pow(b,2.0))/(3*pow(a,2.0));
  double q = (2*pow(b,3.0) - 9*a*b*c + 27*a*a*d)/(27*pow(a,3.0));
  
  if( q*q - 4*pow(p,3.0)/27 < 0.0 )
    {
      cout << "Warning! q^2 +4p^2/27 < 0" << endl;
    }
  if( p >= 0.0 )
    {
      cout << "Warning! p >= 0 " << endl;
    }
  
  double t0 = 2*sqrt(-p/3.0)*TMath::Cos((1.0/3)*TMath::ACos(3*q*sqrt(-3.0/p)/(2*p)));

  double t1 = 2*sqrt(-p/3.0)*TMath::Cos((1.0/3)*TMath::ACos(3*q*sqrt(-3.0/p)/(2*p)) - TMath::TwoPi()/3.0);

  double t2 = 2*sqrt(-p/3.0)*TMath::Cos((1.0/3)*TMath::ACos(3*q*sqrt(-3.0/p)/(2*p)) - 2.0*TMath::TwoPi()/3.0);

  //  cout << "t0="<<t0 << " t1="<<t1<< " t2="<<t2<< endl;
  if(!TMath::IsNaN(t0))
    {
      double x0 = t0 - b/(3*a);
      r.push_back(x0);
    }
  if(!TMath::IsNaN(t1))
    {
      double x1 = t1 - b/(3*a);
      r.push_back(x1);
    }
  if(!TMath::IsNaN(t2))
    {
      double x2 = t2 - b/(3*a);
      r.push_back(x2);
    }
  return(r);
}
