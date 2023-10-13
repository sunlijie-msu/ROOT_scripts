double R = 100;			//implantation rate in Hz
double Thalf = 27.3;		//half life in seconds
double T = 100*Thalf;		//total measurement time in seconds
double lambda = log(2)/Thalf;	//decay constant in 1/s
double t_dead = 2;		//period of dead time between implantation and measurement periods.

double calculateTotalCounts(double ti, double td, bool simplified = false)
//Calculates the total counts over pertiod T (global), with implantation time ti and decay period td
{
  if (simplified)
    return T/(ti+td)*R/lambda*(1-exp(-lambda*ti))*(1-exp(-lambda*td));    

  double N = 0;				//current number of atmos in volume
  double counts = 0;			//total decays during decay period
  double t = 0;				//current time
  
  //start calculation
  while ((t+ti+td)<T)
  {
    N += R/lambda*(1-exp(-ti*lambda));	//add atoms for implantation period
    t += ti;				//add implantation time
    counts += N*exp(-t_dead*lambda)*(1-exp(-(td-t_dead)*lambda));	//add number of decays during decay period
    N -= N*(1-exp(-td*lambda));		//subtract the decays from the atoms number
    t += td;				//add decay time
  }

  //deal with last period
  if ((t+ti)<T)
  {
    N += R/lambda*(1-exp(-ti*lambda));	//add atoms for implantation period
    t += ti;				//add implantation time
    counts += N*(1-exp(-(T-t)*lambda));	//add number of decays during remaining time
  }

  return counts;  
}

double findBestRatio(double ti = Thalf, double rmin = 0.1, double rmax =2, bool simplified = false)
//Find the best ti/td ratio between rmin and rmax
{
  double bestCounts = -1;		//best counts so far;
  double bestRatio = -1;		//best ratio so far;
  double step = (rmax-rmin)/1e2;	//search step

  double r;				//current ratio
  double counts;			//current counts

  //start search
  for (r=rmin;r<rmax;r+=step)
  {
    counts = calculateTotalCounts(ti,r*ti,simplified);	//calculate counts with current ratio
    if (counts>bestCounts)			//update if better than current best ratio
    {
      bestCounts = counts;
      bestRatio = r;
    }
  }
  
  return bestRatio; 			//return best ratio
}

double findBestTi(double tmin = Thalf/10.0, double tmax = Thalf*10, double R = 1.0, bool simplified = false)
//Find best ti period with given ratio R. Rturns best ti/Thalf
{
  double bestCounts = -1;		//best counts so far
  double bestTi = -1;			//best Ti so far
  double step = (tmax-tmin)/1e2;	//search step
  printf("step = %f T_1/2\n",step/Thalf);
  double ti;				//current ti
  double counts;			//current counts

  //start search
  for (ti=tmin;ti<tmax;ti+=step)
  {
    printf("\rcurrent ti %.4f T_1/2",ti/Thalf); 
    counts = calculateTotalCounts(ti,R*ti,simplified);	//calculate counts with current ti
    if (counts>bestCounts)
    {
      bestCounts = counts;
      bestTi = ti;
    }
  }
  printf("\n");
  return bestTi/Thalf;			//return best ti in units of Thalf
}

double findBestCombination(double &bestRatio, double &bestTi, double tmin = 0.1*Thalf, double tmax = 3*Thalf, double rmin = 0.1, double rmax = 2, bool simplified = false)
//search the phase spcae of ti and r to find the best combination
{
  double bestCounts = -1;		//best counts so far
  double step = (tmax-tmin)/1e3;	//search step
  printf("step = %f T_1/2\n",step/Thalf);
  double ti,r;				//current ti
  double counts;			//current counts

  //start search
  for (ti=tmin;ti<tmax;ti+=step)
  {
    printf("\rcurrent ti %.4f T_1/2",ti/Thalf); 
    r = findBestRatio(ti,rmin,rmax,simplified);
    counts = calculateTotalCounts(ti,r*ti,simplified);	//calculate counts with current ti
    if (counts>bestCounts)
    {
      bestCounts = counts;
      bestTi = ti;
      bestRatio = r;

    }
  }
  printf("efficiency = %.2e\n",bestCounts/(T*R));
  return bestCounts/(T*R);  
}

void calculateImplantationDecayPeriods()
{
  printf("implantation/decay periods effect on counting rate\n");
  printf("implantation (ms)   decay (ms)   counts/beam particle\n");
  double ti,td,best,counts;
  t_dead = 0;
  best = findBestCombination(ti,td);
    
  for (t_dead=10e-3;t_dead<60e-3;t_dead+=10e-3)
  {
    printf("waiting time between implantation and measurement: %.0f ms\n",t_dead/1e-3);    
    for (ti=0.15;ti<0.51;ti+=0.05)
      for (td=0.15;td<0.51;td+=0.05)
      {
        counts = calculateTotalCounts(ti,td)/(T*R);
        if (counts/best>0.3)
          printf("%.0f %.0f %.2f\n",ti/1e-3,td/1e-3,counts);    
      }
  }
}
