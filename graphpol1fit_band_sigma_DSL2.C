#include <iostream>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <time.h>
#include <map>
#include <TROOT.h>
#include <TSystem.h>
#include <TCutG.h>
#include "TChain.h"
#include "TStyle.h"
#include "TLine.h"
#include "TCanvas.h"
#include "TEllipse.h"
#include "TString.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TSpectrum.h"
#include "TVirtualFitter.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "TMatrixDSymfwd.h"
#include "TRandom.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH1D.h"
#include "TH2.h"
#include "TF1.h"
#include "TMath.h"
#include "TPaveText.h"
#include "TPad.h"
#include "TAxis.h"
#include "TClass.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "string.h"
#include "stdlib.h"
#include "stdio.h"
#include "TLegend.h"
using namespace std;
void graphpol1fit_band_sigma_DSL2()// the best band explanation
	//graph-pol1 fit sigma (with confidence band) of GRIFFIN obtained from gausnerfcpol1_peakfit_response_DSL.C for 6 unshifited peaks
{//also get the sigma at five 31S gamma energies of interest (1248,2234,3076,4971,5156) for the next simulation and uncertainty analysis.
	time_t tim;
	struct tm *at;
	char now[80];
	int irawroot;
	char rawrootname[300];
	double slope[300],slopeerr[300];
	double intercept[300],intercepterr[300];
	double p7[300],p6[300],p5[300],p4[300],p3[300],p2[300],p1[300],p0[300];//for output
	double p7err[300],p6err[300],p5err[300],p4err[300],p3err[300],p2err[300],p1err[300],p0err[300];
	double p01err[300];
	const int ID1=0;// no need to change idetector=ID1//which detector
	const int ID2=0;// no need to change idetector<=ID2//which detector modify

	int binwidth=1;
	//	int binnum=(upperlimit-lowerlimit)/binwidth;
	int minrange=0,maxrange=0,minbin,maxbin;
	double sigma=1,thresh=0.8,sigmab=20;//adjust
	double gaplow=70.,gaphigh=70.;//fitting range随分辨不同调整//
	unsigned long idetector;
	int jj,iipeak,ibin,k,ii;
	char paraprint[100],b_name[200],histo_name[200],h_name[200],hfit_name[200],hcali_name[200],tauflag[200];
	char command1[100],command2[100];
	//string paraprint[30],h_name[80],histo_name[80],b_name[80],hfit_name[80],hcali_name[80];
	TCanvas *canvaspeak[300];
	TCanvas *canvascali[300];
	TH1F *histo[ID2+1];//TH1F peak search+gauss fit,creat histograms
	const int nummax=10; //modify
	int peaknum=nummax;
	TF1 *total[nummax];//creat function
	TF1 *p[nummax], *g[nummax], *b[nummax];
	TGraph *graph1[300], *graph1bandhigh[300], *graph1bandlow[300];//TGraph
	TGraph *graph2[300], *graph2bandhigh[300], *graph2bandlow[300];//TGraph
	TGraph *graph3[300], *graph3bandhigh[300], *graph3bandlow[300];//TGraph
	TGraph *graph4[300], *graph4bandhigh[300], *graph4bandlow[300];//TGraph
	double energy[ID2+1][nummax],energyerr[ID2+1][nummax];
	double peaky[nummax],peakyerr[nummax];
	double sig[ID2+1][nummax],sigerr[ID2+1][nummax];
	double tau[ID2+1][nummax],tauerr[ID2+1][nummax];
	double residual[ID2+1][nummax],residualerr[ID2+1][nummax];
	double range_sig1[ID2+1], range_sig2[ID2+1], range_tau1[ID2+1], range_tau2[ID2+1];
	double parchi[ID2+1],parNDF[ID2+1],p_value[ID2+1];
	// 	double *energy;//if you don't know how many peaks will be found, use this
	// 	double *peaky;//if you don't know how many peaks will be found, use this
//	double energylit[6]={451.7, 493.3, 944.9, 1460.820, 1612.4, 2614.511};//25Si
//	double energyliterr[6]={0.5, 0.7, 0.5, 0.005, 0.5, 0.010};//25Si
// 	double energylit[6]={450.7, 1460.820, 1599, 2614.511, 2908, 7801};//23Al
// 	double energyliterr[6]={0.15, 0.005, 2, 0.010, 3, 2};//23Al

	const int Npoint=17;
	double energy_point[Npoint] = {500,1000,1500,2000,2500,3000,3500,4000,4500,5000,5500,6000,6500,7000,7500,8000,8500};// all NuDat values, set location of point for single value modify
	double err_point[Npoint];  // error on the function at point x0 for single value
	
//  	const int Npoint=500;
//  	double energy_point[Npoint];// set location of point for drawing
//  	double err_point[Npoint];  // error on the function at point x0 for drawing
//  	for(jj=0;jj<Npoint;jj++){	energy_point[jj]=jj*10;}

	double peakmin;
	double chtemp;
	double par[nummax][7],parerr[nummax][7]; //for input
	
// 	par[0][6]= 0.00048726;	parerr[0][6]=0.00000025;
// 	par[0][5]=-0.01184952;	parerr[0][5]=0.00000197;
// 	par[0][4]= 0.11877336;	parerr[0][4]=0.00001408;
// 	par[0][3]=-0.61313810;	parerr[0][3]=0.00009402;
// 	par[0][2]= 1.33576400;	parerr[0][2]=0.00060392;
// 	par[0][1]= 0.95942945;	parerr[0][1]=0.00366288;
// 	par[0][0]=-8.24241330;	parerr[0][0]=0.01691285; //useless for now
//	par[0][7]=1;	parerr[0][7]=0.1;//C useless for now

	int lowerlimit[ID2+1][nummax],upperlimit[ID2+1][nummax];
	for(idetector=ID1;idetector<=ID2;idetector++)//which detector
	{
		peaknum=nummax;
//		if(idetector==8||idetector==11||idetector==14|idetector==15) continue;
// 		if(idetector==5)peaknum=nummax-1;
// 		if(idetector==0)peaknum=nummax-4;
 		sprintf(b_name,"%s","D:/X/out/DSL2/DSL2_sigma_tau_addback.dat");
 		ifstream infile(b_name,ios::in);//The input data that need to be fitted
		for(iipeak=0;iipeak<peaknum;iipeak++)//which peak in a detector modify
		{
// 			energy[0][0]=278.997;	energyerr[0][0]=0.5264; sig[0][0]=1.1572; sigerr[0][0]=0.0396;
// 			energy[0][1]=547.655;	energyerr[0][1]=0.4874; sig[0][1]=1.3567;	sigerr[0][1]=0.0715;
// 			energy[0][2]=2813.990;	energyerr[0][2]=0.1093;	sig[0][2]=1.3881;	sigerr[0][2]=0.0347;
// 			energy[0][3]=3597.020;	energyerr[0][3]=0.2993;	sig[0][3]=1.8728;	sigerr[0][3]=0.0948;
// 			energy[0][4]=3735.890;	energyerr[0][4]=0.2103; sig[0][4]=1.6225;	sigerr[0][4]=0.0700;
// 			energy[0][5]=6128.490;	energyerr[0][5]=0.2495; sig[0][5]=1.8517;	sigerr[0][5]=0.0778;
			infile>>energy[idetector][iipeak]>>energyerr[idetector][iipeak]>>sig[idetector][iipeak]>>sigerr[idetector][iipeak]>>tau[idetector][iipeak]>>tauerr[idetector][iipeak];
			//energyerr[idetector][iipeak] = energyerr[idetector][iipeak] * 100; // for fun
			//sigerr[idetector][iipeak] = sigerr[idetector][iipeak] * 5; // for fun
			cout<<"GRIFFIN_"<<idetector<<'	'<<energy[idetector][iipeak]<<'	'<<energyerr[idetector][iipeak]<<'	'<<sig[idetector][iipeak]<<'	'<<sigerr[idetector][iipeak]<<'	'<<tau[idetector][iipeak]<<'	'<<tauerr[idetector][iipeak]<<endl;
		}
		if(idetector==0) {range_sig1[idetector]=0; range_sig2[idetector]=5.4; range_tau1[idetector]=0; range_tau2[idetector]=8.0;}
		//if(idetector==1) {range_sig1[idetector]=0; range_sig2[idetector]=4; range_tau1[idetector]=0; range_tau2[idetector]=3.5;}
	}

	ofstream outfile("D:/X/out/DSL2/peakcali.dat",ios::out);
	ofstream outfile2("D:/X/out/DSL2/peakcalipara.dat",ios::out);


	for(idetector=0;idetector<=0;idetector++)//don't change if there is only one detector
	{
		peaknum=nummax-0;// nummax-1 for fitting, nummax-0 for figure
// 		if(idetector==8||idetector==11||idetector==14|idetector==15) continue;
// 		if(idetector==5)peaknum=nummax-1;
// 		if(idetector==0)peaknum=nummax-4;


		//************ sigma fit by linear func *****************************************

		//sprintf(hcali_name,"%s%d","Efficiency of SeGA",idetector);
		sprintf(hcali_name,"%s","GRIFFIN0deg_Addback");
		canvascali[idetector]=new TCanvas(hcali_name,hcali_name,1300,700);//建立画布
		canvascali[idetector]->cd();//进入画布
		TPad* pad1 = new TPad("pad1", "The pad 70% of the height", 0.0, 0.42, 1.0, 1.0);
		// xlow, ylow, xup, yup
		TPad* pad2 = new TPad("pad2", "The pad 30% of the height", 0.0, 0.0, 1.0, 0.42);

		// Set margins for pad1
		pad1->SetTopMargin(0.04);  // relative to pad1 height
		pad1->SetBottomMargin(0.05); // relative to pad1 height
		pad1->SetLeftMargin(0.08);  // relative to pad1 width
		pad1->SetRightMargin(0.04); // relative to pad1 width

		// Set margins for pad2
		pad2->SetTopMargin(0.04);    // relative to pad2 height
		pad2->SetBottomMargin(0.4); // relative to pad2 height
		pad2->SetLeftMargin(0.08);    // relative to pad2 width
		pad2->SetRightMargin(0.04);   // relative to pad2 width

		pad1->Draw();
		pad1->cd();
		pad1->SetFrameLineWidth(2);

		//		graph1[idetector]=new TGraph(peaknum,energy,sig);//TGraph *gr1=new TGraph(n,x,y);
		graph1[idetector]= new TGraphErrors(peaknum,energy[idetector],sig[idetector],energyerr[idetector],sigerr[idetector]);//画error bars TGraph(n,x,y,ex,ey);
		graph1[idetector]->SetTitle("");
		graph1[idetector]->GetXaxis()->SetTitle("E_{#gamma} (keV)");
		graph1[idetector]->GetYaxis()->SetTitle("#sigma");
		graph1[idetector]->GetXaxis()->CenterTitle();//居中
		graph1[idetector]->GetYaxis()->CenterTitle();//居中
		graph1[idetector]->GetXaxis()->SetLabelFont(132);//坐标字体
		graph1[idetector]->GetYaxis()->SetLabelFont(132);//坐标字体
		graph1[idetector]->GetXaxis()->SetLabelSize(0.135);
		graph1[idetector]->GetYaxis()->SetLabelSize(0.135);
		graph1[idetector]->GetXaxis()->SetLabelOffset(0.07);//坐标偏移
		graph1[idetector]->GetXaxis()->SetTitleOffset(1.0);//轴名偏移
		graph1[idetector]->GetXaxis()->SetTitleFont(132);//轴名字体
		graph1[idetector]->GetYaxis()->SetTitleFont(132);//轴名字体
		graph1[idetector]->GetYaxis()->SetTitleOffset(0.295);//轴名偏移
		graph1[idetector]->GetXaxis()->SetTitleSize(0.135);
		graph1[idetector]->GetYaxis()->SetTitleSize(0.135);
		graph1[idetector]->GetYaxis()->SetNdivisions(505);
		graph1[idetector]->GetYaxis()->SetTickLength(0.010);
		graph1[idetector]->GetXaxis()->SetLimits(0, 9000);
		graph1[idetector]->GetYaxis()->SetRangeUser(range_sig1[idetector], range_sig2[idetector]);
		graph1[idetector]->SetLineWidth(2);
		graph1[idetector]->SetStats(0);
		graph1[idetector]->SetMarkerStyle(7);
		graph1[idetector]->SetMarkerColor(1);



		TF1 *pol1 = new TF1("pol1","pol1",0, 8998);
		pol1->SetLineWidth(0);
		pol1->SetNpx(89980);
		pol1->SetParNames("p0","p1");//y=p1*x+p0
		//graph1[idetector]->Fit("pol1");//pol1 can be used directly without TF1 constructor in CINT
		TFitResultPtr r_sig = graph1[idetector]->Fit("pol1","MES");//"S" means the result of the fit is returned in the TFitResultPtr
		//If the option "S" is instead used, TFitResultPtr contains the TFitResult and behaves as a smart pointer to it.
		//The fit parameters, error and chi2 (but not covariance matrix) can be retrieved also from the fitted function.
		TH1D *hint_sig = new TH1D("hint_sig", "Fitted exponential with conf.band", 8998, 0, 8998);//Create a histogram to hold the confidence intervals
		TVirtualFitter * fitter = TVirtualFitter::GetFitter();//The method TVirtualFitter::GetFitter())->Get the parameters of your fitting function after having it fitted to an histogram.
		fitter->GetConfidenceIntervals(hint_sig, 0.683);//By default the intervals are inflated using the chi2/ndf value of the fit if a chi2 fit is performed
		//confidence interval for the colored band: 1σ confidence interval: P=0.683, 1σ confidence interval: P=0.95, 3σ confidence interval: P=0.997
		//hint_sig will contain the CL result that you can draw on top of your fitted graph.
		//where hint_sig will hold the errors and could superimpose it on the same canvas where you plot central values.
		//The method TVirtualFitter::GetFitter())->GetConfidenceIntervals(hint_sig, 0.683 or 0.95 or 0.997) computes the confidence level of your model (fitting) function after having it fitted to an histogram.
		hint_sig->SetStats(kFALSE);
		hint_sig->SetFillColor(kRed-9);//confidence band [colored] shows 1 or 2 or 3 standard deviation uncertainty
		//TMatrixDSym cov = r_sig->GetCovarianceMatrix();//useless for now
		TMatrixD cov = r_sig->GetCovarianceMatrix();//error matrix
		TMatrixD cor = r_sig->GetCorrelationMatrix();//parameter correlation coefficients
		cov.Print();
		cor.Print();

		graph1[idetector]->Draw("AP");//"A": Axis are drawn around the graph, "P": The current marker is plotted at each point
		hint_sig->Draw("e3 same");//draw the confidence band [colored] showing 1 or 2 or 3 standard deviation uncertainty
		graph1[idetector]->Draw("P same");//draw the points again above error band
		pad1->RedrawAxis();

		p0[idetector]=pol1->GetParameter(0);//intercept
		p1[idetector]=pol1->GetParameter(1);//slope
		p0err[idetector]=pol1->GetParError(0);
		//cout<<p0err[idetector]<<endl;
		p0err[idetector]=sqrt(fitter->GetCovarianceMatrixElement(0,0));//the same as pol1->GetParError(0);
		//cout<<p0err[idetector]<<endl;
		p1err[idetector]=pol1->GetParError(1);
		//cout<<p1err[idetector]<<endl;
		p1err[idetector]=sqrt(fitter->GetCovarianceMatrixElement(1,1));//the same pol1->GetParError(1);
		//cout<<p1err[idetector]<<endl;
		p01err[idetector]=fitter->GetCovarianceMatrixElement(0,1);//this one cannot be obtained from pol1->GetParError
		//cout<<p01err[idetector]<<endl;
		parchi[idetector]=pol1->GetChisquare();
		parNDF[idetector]=pol1->GetNDF();
		p_value[idetector]=pol1->GetProb();//This probability is not the “probability that your fit is good.” If you did many fake experiments (draw many random samples of data points from the assumed distribution (your fit function)), this is the percentage of experiments that would give χ2 values ≥ to the one you got in this experiment.

		r_sig->GetConfidenceIntervals(Npoint, 1, 1, energy_point, err_point, 0.683, true);//get the error of sigma/tau at energies of interest, this command is independent of the colored band
		//(Number of x points, 1, 1, x, err, confidence level, false); norm is a flag to control if the intervals need to be inflated by the chi2/ndf value. true is inflated, false is not inflated.

		for(iipeak=0;iipeak<peaknum;iipeak++)//get sigma/tau at energies for plotting the residuals
		{
			residual[idetector][iipeak]= sig[idetector][iipeak]-pol1->Eval(energy[idetector][iipeak]);
			residualerr[idetector][iipeak]=sigerr[idetector][iipeak];
			outfile<<"Eg=	"<<energy[idetector][iipeak]<<"	residual=	"<<residual[idetector][iipeak]<<endl;
		}

		for(iipeak=0;iipeak<Npoint;iipeak++)//get sigma/tau at energies of interest for output file
		{
			outfile<<"GRIFFIN_"<<idetector<<"	Eg=	"<<energy_point[iipeak]<<"	sig=	"<<pol1->Eval(energy_point[iipeak])<<"	err_sig=	"<<err_point[iipeak]<<endl;
// 			outfile<<"err_sig from band=	"<<err_point[iipeak]<<endl;
// 			outfile<<"err_sig from cov=	"<<sqrt(energy_point[iipeak]*energy_point[iipeak]*p1err[idetector]*p1err[idetector]+p0err[idetector]*p0err[idetector]+2*energy_point[iipeak]*p01err[idetector])<<endl;
// 			outfile<<"err_sig from cov inflated=	"<<sqrt(energy_point[iipeak]*energy_point[iipeak]*p1err[idetector]*p1err[idetector]+p0err[idetector]*p0err[idetector]+2*energy_point[iipeak]*p01err[idetector])*sqrt(parchi[idetector]/parNDF[idetector])<<endl;
				//GRIFFIN_0       Eg=     4971    sig=    1.71437
				//err_sig from band=      0.0403738
				//err_sig from cov=       0.0403221
				//err_sig from band inflated=	0.106545
				//err_sig from cov inflated=      0.0932589
		}
		
// 		TVirtualFitter * fitter = TVirtualFitter::GetFitter();//for peak fit, useless here
// 		p0[idetector]=fitter->GetParameter(0);
// 		p1[idetector]=fitter->GetParameter(1);
// 		double E_mu_mu = fitter->GetCovarianceMatrixElement(0,0);
// 		cout<<sqrt(E_mu_mu)<<endl;
// 		double E_mu_mu = fitter->GetCovarianceMatrixElement(1,1);
// 		cout<<sqrt(E_mu_mu)<<endl;
// 		double E_mu_mu = fitter->GetCovarianceMatrixElement(0,1);
// 		cout<<sqrt(abs(E_mu_mu))<<endl;
// 		double errylow[Npoint],erryhigh[Npoint];
// 		for(iipeak=0;iipeak<Npoint;iipeak++)
// 		{
// 			//x[iipeak]=energy_point[iipeak];
// 			//y[iipeak]=pol1->Eval(energy_point[iipeak]);//is equal to y[iipeak]=p1[idetector]*x[iipeak]+p0[idetector];
// 			//dy[iipeak]=sqrt(x[iipeak]*x[iipeak]*p1err[idetector]*p1err[idetector]+p0err[idetector]*p0err[idetector]+p1[idetector]*p1[idetector]*dx[iipeak]*dx[iipeak]);
// 			//dy[iipeak]=err_point[iipeak];
// 			//cout<<"	"<<energy_point[iipeak]<<"	"<<err_point[iipeak]<<endl;
// 			errylow[iipeak]=pol1->Eval(energy_point[iipeak])-err_point[iipeak];
// 			erryhigh[iipeak]=pol1->Eval(energy_point[iipeak])+err_point[iipeak];
// 		}
// 		graph1bandlow[idetector]= new TGraph(Npoint,energy_point,errylow);//画error bars TGraph(n,x,y,ex,ey);
// 		graph1bandlow[idetector]->SetLineColor(2);
// 		graph1bandlow[idetector]->SetLineWidth(1);
// 		graph1bandlow[idetector]->Draw("C");//"C" A smooth Curve is drawn
// 		graph1bandhigh[idetector]= new TGraph(Npoint,energy_point,erryhigh);//画error bars TGraph(n,x,y,ex,ey);
// 		graph1bandhigh[idetector]->SetLineColor(2);
// 		graph1bandhigh[idetector]->SetLineWidth(1);
// 		graph1bandhigh[idetector]->Draw("C");//"C" A smooth Curve is drawn

		TPaveText *textpol1 = new TPaveText(0.69,0.10,0.95,0.45,"brNDC");//left, down, right, up
		textpol1->SetBorderSize(1);
		textpol1->SetFillColor(0);
		textpol1->SetTextAlign(12);//align = 10*HorizontalAlign + VerticalAlign, 12 means水平左对齐、垂直居中对齐
		textpol1->SetTextFont(132);//font = 10 * fontID + precision, 12+2,12 means Symbol; 13+2, 13 means News Times Roman
		sprintf(paraprint,"y=p1*x+p0");
		textpol1->AddText(paraprint);
		sprintf(paraprint,"p1=%.9f+/-%.9f",p1[idetector],p1err[idetector]);
		textpol1->AddText(paraprint);
		sprintf(paraprint,"p0=%.5f+/-%.5f",p0[idetector],p0err[idetector]);
		textpol1->AddText(paraprint);
		sprintf(paraprint,"Chi2=%.2f",parchi[idetector]);
		textpol1->AddText(paraprint);
		sprintf(paraprint,"NDF=%.0f",parNDF[idetector]);
		textpol1->AddText(paraprint);
		sprintf(paraprint,"p-val=%e",p_value[idetector]);
		textpol1->AddText(paraprint);
		textpol1->Draw();

		outfile2<<"GRIFFIN_"<<idetector<<"	y=p1*x+p0	p1=	"<<setprecision(8)<<p1[idetector]<<"	+/-	"<<p1err[idetector]<<"	p0=	"<<p0[idetector]<<"	+/-	"<<p0err[idetector]<<"	Chi2=	"<<parchi[idetector]<<"	ndf=	"<<parNDF[idetector]<<"	p-value=	"<<p_value[idetector];
// 		outfile2<<"SeGA_"<<idetector<<"	y=pol1*x	"<<endl;
// 		outfile2<<"p7=	"<<setprecision(8)<<p0[idetector]<<"	+/-	"<<p0err[idetector]<<endl;
		//		outfile2<<"Chi2=	"<<parchi[idetector]<<"	ndf=	"<<parNDF[idetector]<<"	p-value=	"<<p_value[idetector];

		
		graph2[idetector]= new TGraphErrors(peaknum,energy[idetector],residual[idetector],energyerr[idetector],residualerr[idetector]);//画error bars TGraph(n,x,y,ex,ey);
		graph2[idetector]->SetTitle("");
		graph2[idetector]->GetXaxis()->SetTitle("E_{#gamma} (keV)");//轴名
		graph2[idetector]->GetYaxis()->SetTitle("Data #minus Fit");//轴名
		graph2[idetector]->GetXaxis()->CenterTitle();//居中
		//graph2[idetector]->GetYaxis()->CenterTitle();//居中
		graph2[idetector]->GetXaxis()->SetLabelFont(132);//坐标字体
		graph2[idetector]->GetYaxis()->SetLabelFont(132);//坐标字体
		graph2[idetector]->GetXaxis()->SetLabelSize(0.18);
		graph2[idetector]->GetYaxis()->SetLabelSize(0.18);
		graph2[idetector]->GetXaxis()->SetTitleFont(132);//轴名字体
		graph2[idetector]->GetYaxis()->SetTitleFont(132);//轴名字体
		graph2[idetector]->GetXaxis()->SetTitleOffset(1.10);//轴名偏移
		graph2[idetector]->GetYaxis()->SetTitleOffset(0.22);//轴名偏移
		graph2[idetector]->GetXaxis()->SetTitleSize(0.18);
		graph2[idetector]->GetYaxis()->SetTitleSize(0.18);
		graph2[idetector]->GetYaxis()->SetNdivisions(105);
		graph2[idetector]->GetYaxis()->SetTickLength(0.010);
		graph2[idetector]->SetStats(0);
		graph2[idetector]->GetXaxis()->SetLimits(0, 9000);
		graph2[idetector]->GetXaxis()->SetRangeUser(0, 9000);
		graph2[idetector]->GetYaxis()->SetRangeUser(-3, 3); 
		graph2[idetector]->SetLineWidth(2);
		graph2[idetector]->SetLineColor(kBlack);
		graph2[idetector]->SetMarkerStyle(7);
		graph2[idetector]->SetMarkerColor(kBlack);
		canvascali[idetector]->cd();//进入画布
		pad2->SetFrameLineWidth(2);
		pad2->Draw();
		pad2->RedrawAxis();
		pad2->cd();
		graph2[idetector]->Draw("AP");//"A": Axis are drawn around the graph, "P": The current marker is plotted at each point
		TLine* T1 = new TLine(0, 0, 9000, 0); //x1,y1,x2,y2
		T1->Draw("R");

		//sprintf(hcali_name,"%s.png",hcali_name);
		sprintf(hcali_name, "%s", "D:/X/out/DSL2/sigma_band.png");
		canvascali[idetector]->SaveAs(hcali_name);//存图
	}//for (idetector=0;idetector<ID;idetector++)
}//peakcali main