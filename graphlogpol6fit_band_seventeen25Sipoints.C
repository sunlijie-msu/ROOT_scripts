#include <iostream>
#include <fstream>
#include <iomanip.h>
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
void graphlogpol6fit_band_seventeen25Sipoints()// graph-logpol6 fit the simulated SeGA efficiency points to get the curve (with colored bands)
{
	time_t tim;
	struct tm *at;
	char now[80];
	int irawroot;
	char rawrootname[300];
	double slope[300],slopeerr[300];
	double intercept[300],intercepterr[300];
	double p7[300],p6[300],p5[300],p4[300],p3[300],p2[300],p1[300],p0[300];//for output
	double p7err[300],p6err[300],p5err[300],p4err[300],p3err[300],p2err[300],p1err[300],p0err[300];
	const int ID1=0;// no need to change i=ID1//which detector
	const int ID2=15;// no need to change i<=ID2//which detector

	int binwidth=1;
	//	int binnum=(upperlimit-lowerlimit)/binwidth;
	int minrange=0,maxrange=0,minbin,maxbin;
	double sigma=1,thresh=0.8,sigmab=20;//adjust
	double gaplow=70.,gaphigh=70.;//fitting range随分辨不同调整//
	unsigned long i;
	int jj,ii,ibin,k,ipeak;
	char paraprint[100],b_name[200],histo_name[200],h_name[200],hfit_name[200],hcali_name[200],tauflag[200];
	char command1[100],command2[100];
	//string paraprint[30],h_name[80],histo_name[80],b_name[80],hfit_name[80],hcali_name[80];
	TCanvas *canvaspeak[300];
	TCanvas *canvascali[300];
	TH1F *histo[ID2+1];//TH1F peak search+gauss fit,creat histograms
	const int nummax=17;//search peaknumber can be many, fitnumber can be limited to 3//no need to modify
	int peaknum=nummax;
	TF1 *total[nummax];//creat function
	TF1 *p[nummax], *g[nummax], *b[nummax];
	TGraph *graph1[300], *graph1bandhigh[300], *graph1bandlow[300];//TGraph
	TGraph *graph2[300], *graph2bandhigh[300], *graph2bandlow[300];//TGraph
	TGraph *graph3[300], *graph3bandhigh[300], *graph3bandlow[300];//TGraph
	TGraph *graph4[300], *graph4bandhigh[300], *graph4bandlow[300];//TGraph
	double energy[ID2+1][nummax],energyerr[ID2+1][nummax];
	double peaky[nummax],peakyerr[nummax];
	double eff[ID2+1][nummax],efferr[ID2+1][nummax];
	double residual[ID2+1][nummax],residualerr[ID2+1][nummax];
	double range_sig1[ID2+1], range_sig2[ID2+1], range_tau1[ID2+1], range_tau2[ID2+1];
	// 	double *energy;//if you don't know how many peaks will be found, use this
	// 	double *peaky;//if you don't know how many peaks will be found, use this
//	double energylit[6]={451.7, 493.3, 944.9, 1460.820, 1612.4, 2614.511};//25Si
//	double energyliterr[6]={0.5, 0.7, 0.5, 0.005, 0.5, 0.010};//25Si
// 	double energylit[6]={450.7, 1460.820, 1599, 2614.511, 2908, 7801};//23Al
// 	double energyliterr[6]={0.15, 0.005, 2, 0.010, 3, 2};//23Al

	const int Npoint=17;
	double energy_point[Npoint] = {452,493,845,884,945,1338,1369,1612,1789,2222,2673,2754,2870,4238,6288,6955,7900};// all NuDat values, set location of point for single value modify
	double err_point[Npoint];  // error on the function at point x0 for single value
	
//  	const int Npoint=500;
//  	double energy_point[Npoint];// set location of point for drawing
//  	double err_point[Npoint];  // error on the function at point x0 for drawing
//  	for(jj=0;jj<Npoint;jj++){	energy_point[jj]=jj*10;}

	double peakmin;
	double chtemp;
	double par[nummax][7],parerr[nummax][7]; //for input
	double parchi[ID2+1],parNDF[ID2+1],p_value[ID2+1];

	//source at 25Si SeGA simulated points 12SeGA
	energy[0][0]=452;	energyerr[0][0]=0.5;	eff[0][0]=0.0574;	efferr[0][0]=0.0004;
	energy[0][1]=493;	energyerr[0][1]=0.7;	eff[0][1]=0.0548165;	efferr[0][1]=0.0003476;
	energy[0][2]=845;	energyerr[0][2]=0.7;	eff[0][2]=0.0421328;	efferr[0][2]=0.0002907;
	energy[0][3]=884;	energyerr[0][3]=0.8;	eff[0][3]=0.0416398;	efferr[0][3]=0.0002871;
	energy[0][4]=945;	energyerr[0][4]=0.5;	eff[0][4]=0.0398055;	efferr[0][4]=0.0002811;
	energy[0][5]=1338;	energyerr[0][5]=0.7;	eff[0][5]=0.0338215;	efferr[0][5]=0.0002570;
	energy[0][6]=1369;	energyerr[0][6]=0.005;	eff[0][6]=0.0330812;	efferr[0][6]=0.0002551;
	energy[0][7]=1612;	energyerr[0][7]=0.5;	eff[0][7]=0.0295800;	efferr[0][7]=0.0002370;
	energy[0][8]=1789;	energyerr[0][8]=0.5;	eff[0][8]=0.0282285;	efferr[0][8]=0.0002304;
	energy[0][9]=2222;	energyerr[0][9]=0.8;	eff[0][9]=0.0239624;	efferr[0][9]=0.0002148;
	energy[0][10]=2673;	energyerr[0][10]=0.6;	eff[0][10]=0.0211531;	efferr[0][10]=0.0002014;
	energy[0][11]=2754;	energyerr[0][11]=0.011;	eff[0][11]=0.0207001;	efferr[0][11]=0.0001943;
	energy[0][12]=2870;	energyerr[0][12]=0.06;	eff[0][12]=0.0198263;	efferr[0][12]=0.0001887;
	energy[0][13]=4238;	energyerr[0][13]=0.06;	eff[0][13]=0.0140199;	efferr[0][13]=0.0001609;
	energy[0][14]=6288;	energyerr[0][14]=2;	eff[0][14]=0.0090627;	efferr[0][14]=0.0001397;
	energy[0][15]=6955;	energyerr[0][15]=2;	eff[0][15]=0.0079365;	efferr[0][15]=0.0001276;
	energy[0][16]=7900;	energyerr[0][16]=2;	eff[0][16]=0.0067176;	efferr[0][16]=0.0001071;
	//source at 25Si SeGA simulated points 11SeGA
// 	energy[0][0]=452;	energyerr[0][0]=0.5;	eff[0][0]=0.0529;	efferr[0][0]=0.0003;
// 	energy[0][1]=493;	energyerr[0][1]=0.7;	eff[0][1]=0.0503727;	efferr[0][1]=0.0003334;
// 	energy[0][2]=845;	energyerr[0][2]=0.7;	eff[0][2]=0.0387974;	efferr[0][2]=0.0002789;
// 	energy[0][3]=884;	energyerr[0][3]=0.8;	eff[0][3]=0.0383053;	efferr[0][3]=0.0002753;
// 	energy[0][4]=945;	energyerr[0][4]=0.5;	eff[0][4]=0.0365968;	efferr[0][4]=0.0002695;
// 	energy[0][5]=1338;	energyerr[0][5]=0.7;	eff[0][5]=0.0310930;	efferr[0][5]=0.0002464;
// 	energy[0][6]=1369;	energyerr[0][6]=0.005;	eff[0][6]=0.0304623;	efferr[0][6]=0.0002447;
// 	energy[0][7]=1612;	energyerr[0][7]=0.5;	eff[0][7]=0.0273165;	efferr[0][7]=0.0002277;
// 	energy[0][8]=1789;	energyerr[0][8]=0.5;	eff[0][8]=0.0261392;	efferr[0][8]=0.0002215;
// 	energy[0][9]=2222;	energyerr[0][9]=0.8;	eff[0][9]=0.0220312;	efferr[0][9]=0.0002060;
// 	energy[0][10]=2673;	energyerr[0][10]=0.6;	eff[0][10]=0.0194582;	efferr[0][10]=0.0001932;
// 	energy[0][11]=2754;	energyerr[0][11]=0.011;	eff[0][11]=0.0190426;	efferr[0][11]=0.0001865;
// 	energy[0][12]=2870;	energyerr[0][12]=0.06;	eff[0][12]=0.0182325;	efferr[0][12]=0.0001809;
// 	energy[0][13]=4238;	energyerr[0][13]=0.06;	eff[0][13]=0.0128265;	efferr[0][13]=0.0001539;
// 	energy[0][14]=6288;	energyerr[0][14]=2;	eff[0][14]=0.0082518;	efferr[0][14]=0.0001335;
// 	energy[0][15]=6955;	energyerr[0][15]=2;	eff[0][15]=0.0073202;	efferr[0][15]=0.0001223;
// 	energy[0][16]=7900;	energyerr[0][16]=2;	eff[0][16]=0.0061710;	efferr[0][16]=0.0001026;



//	ifstream infile("C:/Si24/Si22peakcali/Clover.dat",ios::in);//The data that need to be fitted
//	int lowerlimit[ID2+1][nummax],upperlimit[ID2+1][nummax];
// 	for(i=ID1;i<=ID2;i++)//which detector
// 	{
// 		peaknum=nummax;
// 		if(i==8||i==11||i==14|i==15) continue;
// 		if(i==5)peaknum=nummax-1;//modify
// 		if(i==0)peaknum=nummax-4;
//		sprintf(b_name,"%s%d%s","D:/X/out/Si25/SeGA_",i,"_sigma_freetau.dat");//modify 12 peaks
//		sprintf(b_name,"%s%d%s","D:/X/out/Si25/SeGA_",i,"_sigma_tau_onlybg.dat");//only 8 beta-delayed gamma peaks
//		ifstream infile(b_name,ios::in);//The data that need to be fitted
// 		for(ii=0;ii<peaknum;ii++)//which peak in one SeGA detector modify =0<12
// 		{
// 			infile>>energy[i][ii]>>energyerr[i][ii]>>eff[i][ii]>>efferr[i][ii]>>residual[i][ii]>>residualerr[i][ii];
// 			//cout<<"SeGA_"<<i<<'	'<<energy[i][ii]<<'	'<<energyerr[i][ii]<<'	'<<eff[i][ii]<<'	'<<efferr[i][ii]<<'	'<<residual[i][ii]<<'	'<<residualerr[i][ii]<<endl;
// 		}
// 		if(i==0) {range_sig1[i]=1; range_sig2[i]=2; range_tau1[i]=0; range_tau2[i]=1;}
// 		if(i==1) {range_sig1[i]=0; range_sig2[i]=4; range_tau1[i]=0; range_tau2[i]=3.5;}
// 		if(i==2) {range_sig1[i]=0; range_sig2[i]=4; range_tau1[i]=0; range_tau2[i]=4;}
// 		if(i==3) {range_sig1[i]=0; range_sig2[i]=4; range_tau1[i]=0; range_tau2[i]=3;}
// 		if(i==4) {range_sig1[i]=0; range_sig2[i]=4; range_tau1[i]=0; range_tau2[i]=3;}
// 		if(i==5) {range_sig1[i]=0; range_sig2[i]=2; range_tau1[i]=0; range_tau2[i]=1.6;}
// 		if(i==6) {range_sig1[i]=0; range_sig2[i]=4; range_tau1[i]=0; range_tau2[i]=3;}
// 		if(i==7) {range_sig1[i]=0; range_sig2[i]=4; range_tau1[i]=0; range_tau2[i]=2.5;}
// 		if(i==9) {range_sig1[i]=0; range_sig2[i]=4; range_tau1[i]=-2.5; range_tau2[i]=4;}
// 		if(i==10) {range_sig1[i]=0; range_sig2[i]=4.5; range_tau1[i]=0; range_tau2[i]=4;}
// 		if(i==12) {range_sig1[i]=-6; range_sig2[i]=10; range_tau1[i]=-3; range_tau2[i]=5;}
// 		if(i==13) {range_sig1[i]=-4; range_sig2[i]=7; range_tau1[i]=-2; range_tau2[i]=6;}
// 	}
	sprintf(tauflag,"%s","taufree");
	ofstream outfile("D:/X/out/Si25/outfile/peakcali.dat",ios::out);
	ofstream outfile2("D:/X/out/Si25/outfile/peakcalipara.dat",ios::out);




	for(i=0;i<=0;i++)//don't change
	{
		peaknum=nummax-0;//modify nummax-1 for fitting, nummax-0 for figure
// 		if(i==8||i==11||i==14|i==15) continue;
// 		if(i==5)peaknum=nummax-1;//modify
// 		if(i==0)peaknum=nummax-4;


		//************ sigma fit by linear func *****************************************

		//sprintf(hcali_name,"%s%d","Efficiency of SeGA",i);
		sprintf(hcali_name,"%s","Efficiency of SeGA");
		canvascali[i]=new TCanvas(hcali_name,hcali_name,1000,700);//建立画布
		canvascali[i]->cd();//进入画布
		TPad *pad1 = new TPad("pad1", "The pad 70% of the height",0.0,0.3,1.0,1.0);// Double_t xlow, Double_t ylow, Double_t xup, Double_t yup,
		TPad *pad2 = new TPad("pad2", "The pad 30% of the height",0.0,0.0,1.0,0.3);
		pad1->SetTopMargin(0.03);
		pad1->SetRightMargin(0.02);
		pad1->SetLeftMargin(0.08);
		pad1->SetBottomMargin(0.08);
		//pad1->SetBorderMode(0);

		pad2->SetTopMargin(0.03);
		pad2->SetRightMargin(0.02);
		pad2->SetLeftMargin(0.08);
		pad2->SetBottomMargin(0.21);
		//pad2->SetBorderMode(0);

		pad1->Draw();
		pad2->Draw();
		pad1->cd();
		gStyle->SetOptTitle(0);
		//canvascali[i]->SetGrid();//显示网格
		//		graph1[i]=new TGraph(peaknum,energy,eff);//TGraph *gr1=new TGraph(n,x,y);
		graph1[i]= new TGraphErrors(peaknum,energy[i],eff[i],energyerr[i],efferr[i]);//画error bars TGraph(n,x,y,ex,ey);
		graph1[i]->SetTitle(hcali_name);
		graph1[i]->GetXaxis()->SetTitle("E_{#gamma} (keV)");//轴名
		graph1[i]->GetYaxis()->SetTitle("Relative Efficiency  (arbitrary units)");//轴名
		graph1[i]->GetXaxis()->CenterTitle();//居中
		graph1[i]->GetYaxis()->CenterTitle();//居中
		graph1[i]->GetXaxis()->SetLabelFont(132);//坐标字体
		graph1[i]->GetYaxis()->SetLabelFont(132);//坐标字体
		graph1[i]->GetXaxis()->SetTitleFont(132);//轴名字体
		graph1[i]->GetYaxis()->SetTitleFont(132);//轴名字体
		//graph1[i]->GetYaxis()->SetLabelSize(0.05);//坐标字号
		//graph1[i]->GetYaxis()->SetTitleSize(0.05);//轴名字号
		graph1[i]->GetXaxis()->SetTitleOffset(1.2);//轴名偏移
		graph1[i]->GetYaxis()->SetTitleOffset(1.0);//轴名偏移
		graph1[i]->GetXaxis()->SetRangeUser(120,8500);
		//graph1[i]->GetYaxis()->SetRangeUser(0.00,0.06);
		graph1[i]->SetMarkerStyle(21);
		graph1[i]->SetMarkerColor(1);

		TF1 *logpol6 = new TF1("logpol6","exp([0]+[1]*log(x)+[2]*pow(log(x),2)+[3]*pow(log(x),3)+[4]*pow(log(x),4)+[5]*pow(log(x),5)+[6]*pow(log(x),6))",10,10000);
		logpol6->SetNpx(19800);
		logpol6->SetParNames("p0","p1","p2","p3","p4","p5","p6");//y=exp(p0+p1*lnx+p2*(lnx)^2+...
		logpol6->SetParameter(0, 2.1558592);
		//logpol6->SetParError(0,0.005598814);
		//logpol6->SetParLimits(0,1.777853,1.777853);
		logpol6->SetParameter(1,-1.0714089);
		//logpol6->SetParError(1,0.0009436133);
		//logpol6->SetParLimits(1,-1.175923,-1.175923);
		logpol6->SetParameter(2,1.1115575);
		//logpol6->SetParError(2,0.0001292948);
		//logpol6->SetParLimits(2,1.106058,1.106058);
		logpol6->SetParameter(3,-0.60187579);
		//logpol6->SetParError(3,1.688012e-05);
		//logpol6->SetParLimits(3,-0.6004948,-0.6004948);
		logpol6->SetParameter(4,0.12472777);
		//logpol6->SetParError(4,2.128892e-06);
		//logpol6->SetParLimits(4,0.1250527,0.1250527);
		logpol6->SetParameter(5,-0.011253297);
		//logpol6->SetParError(5,2.523441e-07);
		//logpol6->SetParLimits(5,-0.0112302,-0.0112302);
		logpol6->SetParameter(6,0.00037136058);
		//logpol6->SetParError(6,2.741992e-08);
		//logpol6->SetParLimits(6,0.0003661849,0.0003661849);
		//graph1[i]->Fit("logpol6");//logpol6 can be used directly without TF1 constructor in CINT

		TFitResultPtr r_sig = graph1[i]->Fit("logpol6","MS");//"S" means the result of the fit is returned in the TFitResultPtr
		r_sig = graph1[i]->Fit("logpol6","MS");//"S" means the result of the fit is returned in the TFitResultPtr
		r_sig = graph1[i]->Fit("logpol6","MS");//"S" means the result of the fit is returned in the TFitResultPtr
		//If the option "S" is instead used, TFitResultPtr contains the TFitResult and behaves as a smart pointer to it.
		//The fit parameters, error and chi2 (but not covariance matrix) can be retrieved also from the fitted function.
		TH1D *hint_sig = new TH1D("hint_sig", "Fitted exponential with conf.band", 9000, 0, 9000);//Create a histogram to hold the confidence intervals
		(TVirtualFitter::GetFitter())->GetConfidenceIntervals(hint_sig, 0.683);//By default the intervals are corrected using the chi2/ndf value of the fit if a chi2 fit is performed
		//hint_sig will contain the CL result that you can draw on top of your fitted graph.
		//where hint_sig will hold the errors and could superimpose it on the same canvas where you plot central values.
		//The method TVirtualFitter::GetFitter())->GetConfidenceIntervals(hint_sig, 0.95) computes the confidence level of your model (fitting) function after having it fitted to an histogram.
		hint_sig->SetStats(kFALSE);
		hint_sig->SetFillColor(kRed-9);
		//TMatrixDSym cov = r_sig->GetCovarianceMatrix();//useless for now
		TMatrixD cov = r_sig->GetCorrelationMatrix();
		TMatrixD cor = r_sig->GetCovarianceMatrix();
		cov.Print();
		cor.Print();
		r_sig->GetConfidenceIntervals(Npoint, 1, 1, energy_point, err_point, 0.683, true);//(Number of x points, 1, 1, x, err, confidence level, false); norm is a flag to control if the intervals need to be normalized to the chi2/ndf value. modify
		for(ii=0;ii<Npoint;ii++)
		{
			residual[i][ii]=(logpol6->Eval(energy_point[ii])-eff[i][ii])/logpol6->Eval(energy_point[ii]);
			residualerr[i][ii]=efferr[i][ii]/logpol6->Eval(energy_point[ii]);
			outfile<<"Eg=	"<<energy_point[ii]<<"	eff=	"<<logpol6->Eval(energy_point[ii])<<"	relative residual=	"<<residual[i][ii]<<"	err=	"<<err_point[ii]<<endl;
		}
		graph1[i]->Draw("AP");//"A": Axis are drawn around the graph, "P": The current marker is plotted at each point
//		hint_sig->Draw("e3 same");
		graph1[i]->Draw("P same");//draw the points again above error band

		p0[i]=logpol6->GetParameter(0);
		p1[i]=logpol6->GetParameter(1);
		p2[i]=logpol6->GetParameter(2);
		p3[i]=logpol6->GetParameter(3);
		p4[i]=logpol6->GetParameter(4);
		p5[i]=logpol6->GetParameter(5);
		p6[i]=logpol6->GetParameter(6);

		p0err[i]=logpol6->GetParError(0);
		p1err[i]=logpol6->GetParError(1);
		p2err[i]=logpol6->GetParError(2);
		p3err[i]=logpol6->GetParError(3);
		p4err[i]=logpol6->GetParError(4);
		p5err[i]=logpol6->GetParError(5);
		p6err[i]=logpol6->GetParError(6);

		parchi[i]=logpol6->GetChisquare();
		parNDF[i]=logpol6->GetNDF();
		p_value[i]=logpol6->GetProb();//This probability is not the “probability that your fit is good.” If you did many fake experiments (draw many random samples of data points from the assumed distribution (your fit function)), this is the percentage of experiments that would give χ2 values ≥ to the one you got in this experiment.
		
// 		TVirtualFitter * fitter = TVirtualFitter::GetFitter();
// 		p0[i]=fitter->GetParameter(0);
// 		p1[i]=fitter->GetParameter(1);
// 		double E_mu_mu = fitter->GetCovarianceMatrixElement(0,0);
// 		cout<<sqrt(E_mu_mu)<<endl;
// 		double E_mu_mu = fitter->GetCovarianceMatrixElement(1,1);
// 		cout<<sqrt(E_mu_mu)<<endl;
// 		double E_mu_mu = fitter->GetCovarianceMatrixElement(0,1);
// 		cout<<sqrt(abs(E_mu_mu))<<endl;
// 		double errylow[Npoint],erryhigh[Npoint];
// 		for(ii=0;ii<Npoint;ii++)
// 		{
// 			//x[ii]=energy_point[ii];
// 			//y[ii]=logpol6->Eval(energy_point[ii]);//is equal to y[ii]=p1[i]*x[ii]+p0[i];
// 			//dy[ii]=sqrt(x[ii]*x[ii]*p1err[i]*p1err[i]+p0err[i]*p0err[i]+p1[i]*p1[i]*dx[ii]*dx[ii]);
// 			//dy[ii]=err_point[ii];
// 			//cout<<"	"<<energy_point[ii]<<"	"<<err_point[ii]<<endl;
// 			errylow[ii]=logpol6->Eval(energy_point[ii])-err_point[ii];
// 			erryhigh[ii]=logpol6->Eval(energy_point[ii])+err_point[ii];
// 		}
// 		graph1bandlow[i]= new TGraph(Npoint,energy_point,errylow);//画error bars TGraph(n,x,y,ex,ey);
// 		graph1bandlow[i]->SetLineColor(2);
// 		graph1bandlow[i]->SetLineWidth(1);
// 		graph1bandlow[i]->Draw("C");//"C" A smooth Curve is drawn
// 		graph1bandhigh[i]= new TGraph(Npoint,energy_point,erryhigh);//画error bars TGraph(n,x,y,ex,ey);
// 		graph1bandhigh[i]->SetLineColor(2);
// 		graph1bandhigh[i]->SetLineWidth(1);
// 		graph1bandhigh[i]->Draw("C");//"C" A smooth Curve is drawn

		TPaveText *textpol6 = new TPaveText(0.64,0.44,0.96,0.95,"brNDC");//left, down, right, up
		textpol6->SetBorderSize(1);
		textpol6->SetFillColor(0);
		textpol6->SetTextAlign(12);//align = 10*HorizontalAlign + VerticalAlign, 12 means水平左对齐、垂直居中对齐
		textpol6->SetTextFont(132);//font = 10 * fontID + precision, 12+2,12 means Symbol; 13+2, 13 means News Times Roman
		sprintf(paraprint,"y=exp[sum(p_i ln(E)^i)]");
		textpol6->AddText(paraprint);
		sprintf(paraprint,"p6=%e+/-%e",p6[i],p6err[i]);
		textpol6->AddText(paraprint);
		sprintf(paraprint,"p5=%e+/-%e",p5[i],p5err[i]);
		textpol6->AddText(paraprint);
		sprintf(paraprint,"p4=%e+/-%e",p4[i],p4err[i]);
		textpol6->AddText(paraprint);
		sprintf(paraprint,"p3=%e+/-%e",p3[i],p3err[i]);
		textpol6->AddText(paraprint);
		sprintf(paraprint,"p2=%e+/-%e",p2[i],p2err[i]);
		textpol6->AddText(paraprint);
		sprintf(paraprint,"p1=%e+/-%e",p1[i],p1err[i]);
		textpol6->AddText(paraprint);
		sprintf(paraprint,"p0=%e+/-%e",p0[i],p0err[i]);
		textpol6->AddText(paraprint);
		sprintf(paraprint,"Chi2=%.2f",parchi[i]);
		textpol6->AddText(paraprint);
		sprintf(paraprint,"NDF=%.0f",parNDF[i]);
		textpol6->AddText(paraprint);
		sprintf(paraprint,"p-val=%e",p_value[i]);
		textpol6->AddText(paraprint);
		textpol6->Draw();
		//sprintf(hcali_name,"%s.png",hcali_name);
		sprintf(hcali_name,"%s%s%s%s%s","D:/X/out/responseforDoppler/",hcali_name,"_",tauflag,"_band.png");
		canvascali[i]->SaveAs(hcali_name);//存图
		outfile2<<"SeGA_"<<i<<"	y=logpol6	"<<endl;
		outfile2<<"p6=	"<<setprecision(8)<<p6[i]<<"	+/-	"<<p6err[i]<<endl;
		outfile2<<"p5=	"<<setprecision(8)<<p5[i]<<"	+/-	"<<p5err[i]<<endl;
		outfile2<<"p4=	"<<setprecision(8)<<p4[i]<<"	+/-	"<<p4err[i]<<endl;
		outfile2<<"p3=	"<<setprecision(8)<<p3[i]<<"	+/-	"<<p3err[i]<<endl;
		outfile2<<"p2=	"<<setprecision(8)<<p2[i]<<"	+/-	"<<p2err[i]<<endl;
		outfile2<<"p1=	"<<setprecision(8)<<p1[i]<<"	+/-	"<<p1err[i]<<endl;
		outfile2<<"p0=	"<<setprecision(8)<<p0[i]<<"	+/-	"<<p0err[i]<<endl;
		outfile2<<"Chi2=	"<<parchi[i]<<"	ndf=	"<<parNDF[i]<<"	p-value=	"<<p_value[i];

		pad2->cd();
		graph2[i]= new TGraphErrors(peaknum,energy[i],residual[i],energyerr[i],residualerr[i]);//画error bars TGraph(n,x,y,ex,ey);
		graph2[i]->GetXaxis()->SetTitle("E_{#gamma} (keV)");//轴名
		graph2[i]->GetYaxis()->SetTitle("Relative residual");//轴名
		graph2[i]->GetXaxis()->CenterTitle();//居中
		graph2[i]->GetYaxis()->CenterTitle();//居中
		graph2[i]->GetXaxis()->SetLabelFont(132);//坐标字体
		graph2[i]->GetYaxis()->SetLabelFont(132);//坐标字体
		graph2[i]->GetXaxis()->SetTitleFont(132);//轴名字体
		graph2[i]->GetYaxis()->SetTitleFont(132);//轴名字体
		//graph1[i]->GetYaxis()->SetLabelSize(0.05);//坐标字号
		//graph1[i]->GetYaxis()->SetTitleSize(0.05);//轴名字号
		graph2[i]->GetXaxis()->SetRangeUser(120,8500);
		graph2[i]->GetYaxis()->SetRangeUser(-0.2,0.2);
		graph2[i]->GetXaxis()->SetTitleOffset(1.3);
		graph2[i]->GetYaxis()->SetTitleOffset(0.4);
		graph2[i]->GetXaxis()->SetTitleSize(0.08);
		graph2[i]->GetXaxis()->SetLabelOffset(0.015);
		graph2[i]->GetXaxis()->SetLabelSize(0.08);
		graph2[i]->GetYaxis()->SetLabelSize(0.08);
		graph2[i]->GetYaxis()->SetTitleSize(0.08);
		graph2[i]->GetXaxis()->SetNdivisions(520);//n = n1 + 100*n2 + 10000*n3
		graph2[i]->GetXaxis()->SetNdivisions(10,10,1);
		graph2[i]->GetYaxis()->SetNdivisions(505);
		graph2[i]->SetMarkerStyle(21);
		graph2[i]->SetMarkerColor(1);
		graph2[i]->Draw("AP");//"A": Axis are drawn around the graph, "P": The current marker is plotted at each point
		TLine *T1=new TLine(120,0,8500,0);
		T1->Draw("R");
/*
		//************ residual fit by sqrt + const func *****************************************

		sprintf(hcali_name,"%s%d","Tau of SeGA",i);
		canvascali[i]=new TCanvas(hcali_name,hcali_name,900,600);//建立画布
		canvascali[i]->cd();//进入画布
		//canvascali[i]->SetGrid();//显示网格
		//		graph1[i]=new TGraph(peaknum,energy,residual);//TGraph *gr1=new TGraph(n,x,y);
		graph3[i]= new TGraphErrors(peaknum,energy[i],residual[i],energyerr[i],residualerr[i]);//画error bars TGraph(n,x,y,ex,ey);
		graph3[i]->SetTitle(hcali_name);
		graph3[i]->GetXaxis()->SetTitle("Energy");//轴名
		graph3[i]->GetYaxis()->SetTitle("Tau");//轴名
		graph3[i]->GetXaxis()->CenterTitle();//居中
		graph3[i]->GetYaxis()->CenterTitle();//居中
		graph3[i]->GetXaxis()->SetLabelFont(132);//坐标字体
		graph3[i]->GetYaxis()->SetLabelFont(132);//坐标字体
		graph3[i]->GetXaxis()->SetTitleFont(132);//轴名字体
		graph3[i]->GetYaxis()->SetTitleFont(132);//轴名字体
		//graph1[i]->GetYaxis()->SetLabelSize(0.05);//坐标字号
		//graph1[i]->GetYaxis()->SetTitleSize(0.05);//轴名字号
		graph3[i]->GetXaxis()->SetTitleOffset(1.2);//轴名偏移
		graph3[i]->GetYaxis()->SetTitleOffset(1.3);//轴名偏移
		graph3[i]->GetYaxis()->SetRangeUser(range_tau1[i],range_tau2[i]);
		graph3[i]->SetMarkerStyle(21);
		graph3[i]->SetMarkerColor(1);

		TF1 *func_sqrt_const=new TF1("func_sqrt_const","[0]+[1]*sqrt(x)",0,60000);//多项式拟合，调节合适的拟合range
		func_sqrt_const->SetParNames("p0","p1");//y=p1*sqrt(x)+p0
//		graph3[i]->Fit("func_sqrt_const");//logpol6 can be used directly without TF1 constructor in CINT

		TFitResultPtr r_tau_sqrt_const = graph3[i]->Fit("func_sqrt_const","MS");//"S" means the result of the fit is returned in the TFitResultPtr
		//If the option "S" is instead used, TFitResultPtr contains the TFitResult and behaves as a smart pointer to it.
		//The fit parameters, error and chi2 (but not covariance matrix) can be retrieved also from the fitted function.
		TH1D *hint_tau_sqrt_const = new TH1D("hint_tau_sqrt_const", "Fitted exponential with conf.band", 9000, 0, 9000);//Create a histogram to hold the confidence intervals
		(TVirtualFitter::GetFitter())->GetConfidenceIntervals(hint_tau_sqrt_const, 0.683);//By default the intervals are corrected using the chi2/ndf value of the fit if a chi2 fit is performed
		//hint_tau_sqrt_const will contain the CL result that you can draw on top of your fitted graph.
		//where hint_sig will hold the errors and could superimpose it on the same canvas where you plot central values.
		//The method TVirtualFitter::GetFitter())->GetConfidenceIntervals(hint_sig, 0.95) computes the confidence level of your model (fitting) function after having it fitted to an histogram.
		hint_tau_sqrt_const->SetStats(kFALSE);
		hint_tau_sqrt_const->SetFillColor(kRed-9);
		//TMatrixDSym cov = r_tau_sqrt_const->GetCovarianceMatrix();//useless for now
		TMatrixD cov = r_tau_sqrt_const->GetCorrelationMatrix();
		TMatrixD cor = r_tau_sqrt_const->GetCovarianceMatrix();
		cov.Print();
		cor.Print();
		r_tau_sqrt_const->GetConfidenceIntervals(Npoint, 1, 1, energy_point, err_point, 0.683, false);//(Number of x points, 1, 1, x, err, confidence level, norm); norm is a flag to control if the intervals need to be normalized to the chi2/ndf value.
		for(ii=0;ii<Npoint;ii++)
		{
			outfile<<"SeGA_"<<i<<"	Eg=	"<<energy_point[ii]<<"	residual=	"<<func_sqrt_const->Eval(energy_point[ii])<<"	err_tau=	"<<err_point[ii]<<endl;
		}
		graph3[i]->Draw("AP");//"A": Axis are drawn around the graph, "P": The current marker is plotted at each point
		hint_tau_sqrt_const->Draw("e3 same");
		graph3[i]->Draw("P same");//draw the points again above error band

		p0[i]=func_sqrt_const->GetParameter(0);
		p0err[i]=func_sqrt_const->GetParError(0);
		p1[i]=func_sqrt_const->GetParameter(1);
		p1err[i]=func_sqrt_const->GetParError(1);
		parchi[i]=func_sqrt_const->GetChisquare();
		parNDF[i]=func_sqrt_const->GetNDF();
		p_value[i]=func_sqrt_const->GetProb();//This probability is not the “probability that your fit is good.” If you did many fake experiments (draw many random samples of data points from the assumed distribution (your fit function)), this is the percentage of experiments that would give χ2 values ≥ to the one you got in this experiment.
		double errylow[Npoint],erryhigh[Npoint];
		for(ii=0;ii<Npoint;ii++)
		{
			//x[ii]=energy_point[ii];
			//y[ii]=func_sqrt_const->Eval(energy_point[ii]);//is equal to y[ii]=p1[i]*sqrt(x[ii])+p0[i];
			//dy[ii]=sqrt(sqrt(x[ii])*sqrt(x[ii])*p1err[i]*p1err[i]+p0err[i]*p0err[i]);
			//dy[ii]=err_point[ii];
			//cout<<"	"<<energy_point[ii]<<"	"<<err_point[ii]<<endl;
			errylow[ii]=func_sqrt_const->Eval(energy_point[ii])-err_point[ii];
			erryhigh[ii]=func_sqrt_const->Eval(energy_point[ii])+err_point[ii];
		}
		graph3bandlow[i]= new TGraph(Npoint,energy_point,errylow);//画error bars TGraph(n,x,y,ex,ey);
		graph3bandlow[i]->SetLineColor(4);
		graph3bandlow[i]->SetLineWidth(1);
		graph3bandlow[i]->Draw("C");
		graph3bandhigh[i]= new TGraph(Npoint,energy_point,erryhigh);//画error bars TGraph(n,x,y,ex,ey);
		graph3bandhigh[i]->SetLineColor(4);
		graph3bandhigh[i]->SetLineWidth(1);
		graph3bandhigh[i]->Draw("C");

		TPaveText *textpol3 = new TPaveText(0.10,0.64,0.34,0.89,"brNDC");//left, down, right, up
		textpol3->SetBorderSize(1);
		textpol3->SetFillColor(0);
		textpol3->SetTextAlign(12);//align = 10*HorizontalAlign + VerticalAlign, 12 means水平左对齐、垂直居中对齐
		textpol3->SetTextFont(132);//font = 10 * fontID + precision, 12+2,12 means Symbol; 13+2, 13 means News Times Roman
		sprintf(paraprint,"y=p1*sqrt(x)+p0");
		textpol3->AddText(paraprint);
		sprintf(paraprint,"p1=%.5f+/-%.5f",p1[i],p1err[i]);
		textpol3->AddText(paraprint);
		sprintf(paraprint,"p0=%.5f+/-%.5f",p0[i],p0err[i]);
		textpol3->AddText(paraprint);
		sprintf(paraprint,"Chisquare=%.2f",parchi[i]);
		textpol3->AddText(paraprint);
		sprintf(paraprint,"NDF=%.0f",parNDF[i]);
		textpol3->AddText(paraprint);
		sprintf(paraprint,"p-value=%e",p_value[i]);
		textpol3->AddText(paraprint);
		textpol3->Draw();
		//sprintf(hcali_name,"%s.png",hcali_name);
		sprintf(hcali_name,"%s%s%s%s%s","D:/X/out/responseforDoppler/",hcali_name,"_",tauflag,"_band.png");
		canvascali[i]->SaveAs(hcali_name);//存图
		outfile2<<"	SeGA_"<<i<<"	y=p1*sqrt(x)+p0	taup1=	"<<setprecision(8)<<p1[i]<<"	+/-	"<<p1err[i]<<"	taup0=	"<<p0[i]<<"	+/-	"<<p0err[i]<<"	Chi2=	"<<parchi[i]<<"	ndf=	"<<parNDF[i]<<"	p-value=	"<<p_value[i]<<endl;





		//************ residual fit by quadratic func *****************************************

		sprintf(hcali_name,"%s%d","AlSi_tau_pol2_SeGA_",i);
		canvascali[i]=new TCanvas(hcali_name,hcali_name,900,600);//建立画布
		canvascali[i]->cd();//进入画布
		//canvascali[i]->SetGrid();//显示网格
		//		graph1[i]=new TGraph(peaknum,energy,residual);//TGraph *gr1=new TGraph(n,x,y);
		graph2[i]= new TGraphErrors(peaknum,energy[i],residual[i],energyerr[i],residualerr[i]);//画error bars TGraph(n,x,y,ex,ey);
		graph2[i]->SetTitle(hcali_name);
		graph2[i]->GetXaxis()->SetTitle("Energy");//轴名
		graph2[i]->GetYaxis()->SetTitle("Tau");//轴名
		graph2[i]->GetXaxis()->CenterTitle();//居中
		graph2[i]->GetYaxis()->CenterTitle();//居中
		graph2[i]->GetXaxis()->SetLabelFont(132);//坐标字体
		graph2[i]->GetYaxis()->SetLabelFont(132);//坐标字体
		graph2[i]->GetXaxis()->SetTitleFont(132);//轴名字体
		graph2[i]->GetYaxis()->SetTitleFont(132);//轴名字体
		//graph1[i]->GetYaxis()->SetLabelSize(0.05);//坐标字号
		//graph1[i]->GetYaxis()->SetTitleSize(0.05);//轴名字号
		graph2[i]->GetXaxis()->SetTitleOffset(1.2);//轴名偏移
		graph2[i]->GetYaxis()->SetTitleOffset(1.3);//轴名偏移
		graph2[i]->GetYaxis()->SetRangeUser(range_tau1[i],range_tau2[i]);
		graph2[i]->SetMarkerStyle(21);
		graph2[i]->SetMarkerColor(1);
		TF1 *pol2=new TF1("pol2","pol2",0,60000);//多项式拟合，调节合适的拟合range
		pol2->SetParNames("p0","p1","p2");//y=p2*x^2+p1*x+p0
//		graph2[i]->Fit("pol2");//logpol6 can be used directly without TF1 constructor in CINT

		TFitResultPtr r_tau_pol2 = graph2[i]->Fit("pol2","MS");//"S" means the result of the fit is returned in the TFitResultPtr
		//If the option "S" is instead used, TFitResultPtr contains the TFitResult and behaves as a smart pointer to it.
		//The fit parameters, error and chi2 (but not covariance matrix) can be retrieved also from the fitted function.
		TH1D *hint_tau_pol2 = new TH1D("hint_tau_pol2", "Fitted exponential with conf.band", 9000, 0, 9000);//Create a histogram to hold the confidence intervals
		(TVirtualFitter::GetFitter())->GetConfidenceIntervals(hint_tau_pol2, 0.683);//By default the intervals are corrected using the chi2/ndf value of the fit if a chi2 fit is performed
		//hint_tau_pol2 will contain the CL result that you can draw on top of your fitted graph.
		//where hint_sig will hold the errors and could superimpose it on the same canvas where you plot central values.
		//The method TVirtualFitter::GetFitter())->GetConfidenceIntervals(hint_sig, 0.95) computes the confidence level of your model (fitting) function after having it fitted to an histogram.
		hint_tau_pol2->SetStats(kFALSE);
		hint_tau_pol2->SetFillColor(kRed-9);
		TMatrixDSym cov = r_tau_pol2->GetCovarianceMatrix();//useless for now
		r_tau_pol2->GetConfidenceIntervals(Npoint, 1, 1, energy_point, err_point, 0.683, false);//(Number of x points, 1, 1, x, err, confidence level, false); norm is a flag to control if the intervals need to be normalized to the chi2/ndf value.
		outfile<<"Eg=	"<<energy_point[3]<<"	residual=	"<<pol2->Eval(energy_point[3])<<"	err_tau=	"<<err_point[3]<<endl;
		graph2[i]->Draw("AP");//"A": Axis are drawn around the graph, "P": The current marker is plotted at each point
		hint_tau_pol2->Draw("e3 same");
		graph2[i]->Draw("P same");//draw the points again above error band

		p0[i]=pol2->GetParameter(0);
		p1[i]=pol2->GetParameter(1);
		p2[i]=pol2->GetParameter(2);
		p0err[i]=pol2->GetParError(0);
		p1err[i]=pol2->GetParError(1);
		p2err[i]=pol2->GetParError(2);
		parchi[i]=pol2->GetChisquare();
		parNDF[i]=pol2->GetNDF();
		p_value[i]=pol2->GetProb();//This probability is not the “probability that your fit is good.” If you did many fake experiments (draw many random samples of data points from the assumed distribution (your fit function)), this is the percentage of experiments that would give χ2 values ≥ to the one you got in this experiment.
		double errylow[Npoint],erryhigh[Npoint];
		for(ii=0;ii<Npoint;ii++)
		{
			//x[ii]=energy_point[ii];
			//y[ii]=pol2->Eval(energy_point[ii]);//is equal to y[ii]=p2[i]*x[ii]*x[ii]+p1[i]*x[ii]+p0[i];
			//dy[ii]=sqrt(x[ii]*x[ii]*x[ii]*x[ii]*p2err[i]*p2err[i]+x[ii]*x[ii]*p1err[i]*p1err[i]+p0err[i]*p0err[i]+(2*p2[i]+p1[i])*(2*p2[i]+p1[i])*dx[ii]);
			//dy[ii]=err_point[ii];
			//cout<<"	"<<energy_point[ii]<<"	"<<err_point[ii]<<endl;
			errylow[ii]=pol2->Eval(energy_point[ii])-err_point[ii];
			erryhigh[ii]=pol2->Eval(energy_point[ii])+err_point[ii];
		}
		graph2bandlow[i]= new TGraph(Npoint,energy_point,errylow);//画error bars TGraph(n,x,y,ex,ey);
		graph2bandlow[i]->SetLineColor(2);
		graph2bandlow[i]->SetLineWidth(1);
		graph2bandlow[i]->Draw("C");
		graph2bandhigh[i]= new TGraph(Npoint,energy_point,erryhigh);//画error bars TGraph(n,x,y,ex,ey);
		graph2bandhigh[i]->SetLineColor(2);
		graph2bandhigh[i]->SetLineWidth(1);
		graph2bandhigh[i]->Draw("C");

		TPaveText *textpol2 = new TPaveText(0.10,0.64,0.34,0.89,"brNDC");//left, down, right, up
		textpol2->SetBorderSize(1);
		textpol2->SetFillColor(0);
		textpol2->SetTextAlign(12);//align = 10*HorizontalAlign + VerticalAlign, 12 means水平左对齐、垂直居中对齐
		textpol2->SetTextFont(132);//font = 10 * fontID + precision, 12+2,12 means Symbol; 13+2, 13 means News Times Roman
		sprintf(paraprint,"y=p2*x^2+p1*x+p0");
		textpol2->AddText(paraprint);
		sprintf(paraprint,"p2=%.9f+/-%.9f",p2[i],p2err[i]);
		textpol2->AddText(paraprint);
		sprintf(paraprint,"p1=%.5f+/-%.5f",p1[i],p1err[i]);
		textpol2->AddText(paraprint);
		sprintf(paraprint,"p0=%.5f+/-%.5f",p0[i],p0err[i]);
		textpol2->AddText(paraprint);
		sprintf(paraprint,"Chi2=%.2f",parchi[i]);
		textpol2->AddText(paraprint);
		sprintf(paraprint,"NDF=%.0f",parNDF[i]);
		textpol2->AddText(paraprint);
		sprintf(paraprint,"p-val=%e",p_value[i]);
		textpol2->AddText(paraprint);
		textpol2->Draw();
		//sprintf(hcali_name,"%s.png",hcali_name);
		sprintf(hcali_name,"%s%s%s%s%s","D:/X/out/responseforDoppler/",hcali_name,"_",tauflag,"_band.png");
		canvascali[i]->SaveAs(hcali_name);//存图
		outfile2<<"	SeGA_"<<i<<"	y=p2*x^2+p1*x+p0	taup2=	"<<setprecision(11)<<p2[i]<<"	+/-	"<<p2err[i]<<"	taup1=	"<<setprecision(8)<<p1[i]<<"	+/-	"<<p1err[i]<<"	taup0=	"<<p0[i]<<"	+/-	"<<p0err[i]<<"	Chi2=	"<<parchi[i]<<"	ndf=	"<<parNDF[i]<<"	p-value=	"<<p_value[i];





		//************ residual fit by sqrt func *****************************************

		sprintf(hcali_name,"%s%d","AlSi_tau_sqrt_SeGA_",i);
		canvascali[i]=new TCanvas(hcali_name,hcali_name,900,600);//建立画布
		canvascali[i]->cd();//进入画布
		//canvascali[i]->SetGrid();//显示网格
		//		graph1[i]=new TGraph(peaknum,energy,residual);//TGraph *gr1=new TGraph(n,x,y);
		graph4[i]= new TGraphErrors(peaknum,energy[i],residual[i],energyerr[i],residualerr[i]);//画error bars TGraph(n,x,y,ex,ey);
		graph4[i]->SetTitle(hcali_name);
		graph4[i]->GetXaxis()->SetTitle("Energy");//轴名
		graph4[i]->GetYaxis()->SetTitle("Tau");//轴名
		graph4[i]->GetXaxis()->CenterTitle();//居中
		graph4[i]->GetYaxis()->CenterTitle();//居中
		graph4[i]->GetXaxis()->SetLabelFont(132);//坐标字体
		graph4[i]->GetYaxis()->SetLabelFont(132);//坐标字体
		graph4[i]->GetXaxis()->SetTitleFont(132);//轴名字体
		graph4[i]->GetYaxis()->SetTitleFont(132);//轴名字体
		//graph1[i]->GetYaxis()->SetLabelSize(0.05);//坐标字号
		//graph1[i]->GetYaxis()->SetTitleSize(0.05);//轴名字号
		graph4[i]->GetXaxis()->SetTitleOffset(1.2);//轴名偏移
		graph4[i]->GetYaxis()->SetTitleOffset(1.3);//轴名偏移
		graph4[i]->GetYaxis()->SetRangeUser(range_tau1[i],range_tau2[i]);
		graph4[i]->SetMarkerStyle(21);
		graph4[i]->SetMarkerColor(1);

		TF1 *func_sqrt=new TF1("func_sqrt","[0]*sqrt(x)",0,60000);//多项式拟合，调节合适的拟合range
		func_sqrt->SetParNames("p0");//y=p0*sqrt(x)
//		graph4[i]->Fit("func_sqrt");//logpol6 can be used directly without TF1 constructor in CINT

		TFitResultPtr r_tau_sqrt = graph4[i]->Fit("func_sqrt","MS");//"S" means the result of the fit is returned in the TFitResultPtr
		//If the option "S" is instead used, TFitResultPtr contains the TFitResult and behaves as a smart pointer to it.
		//The fit parameters, error and chi2 (but not covariance matrix) can be retrieved also from the fitted function.
		TH1D *hint_tau_sqrt = new TH1D("hint_tau_sqrt", "Fitted exponential with conf.band", 9000, 0, 9000);//Create a histogram to hold the confidence intervals
		(TVirtualFitter::GetFitter())->GetConfidenceIntervals(hint_tau_sqrt, 0.683);//By default the intervals are corrected using the chi2/ndf value of the fit if a chi2 fit is performed
		//hint_tau_sqrt will contain the CL result that you can draw on top of your fitted graph.
		//where hint_tau_sqrt will hold the errors and could superimpose it on the same canvas where you plot central values.
		//The method TVirtualFitter::GetFitter())->GetConfidenceIntervals(hint_sig, 0.95) computes the confidence level of your model (fitting) function after having it fitted to an histogram.
		hint_tau_sqrt->SetStats(kFALSE);
		hint_tau_sqrt->SetFillColor(kRed-9);
		TMatrixDSym cov = r_tau_sqrt->GetCovarianceMatrix();
		r_tau_sqrt->GetConfidenceIntervals(Npoint, 1, 1, energy_point, err_point, 0.683, false);//(Number of x points, 1, 1, x, err, confidence level, norm); norm is a flag to control if the intervals need to be normalized to the chi2/ndf value.
		outfile<<"Eg=	"<<energy_point[3]<<"	residual=	"<<func_sqrt->Eval(energy_point[3])<<"	err_tau=	"<<err_point[3]<<endl;
		graph4[i]->Draw("AP");//"A": Axis are drawn around the graph, "P": The current marker is plotted at each point
		hint_tau_sqrt->Draw("e3 same");
		graph4[i]->Draw("P same");//draw the points again above error band

		p0[i]=func_sqrt->GetParameter(0);
		p0err[i]=func_sqrt->GetParError(0);
		parchi[i]=func_sqrt->GetChisquare();
		parNDF[i]=func_sqrt->GetNDF();
		p_value[i]=func_sqrt->GetProb();//This probability is not the “probability that your fit is good.” If you did many fake experiments (draw many random samples of data points from the assumed distribution (your fit function)), this is the percentage of experiments that would give χ2 values ≥ to the one you got in this experiment.
		double errylow[Npoint],erryhigh[Npoint];
		for(ii=0;ii<Npoint;ii++)
		{
			//x[ii]=energy_point[ii];
			//y[ii]=func_sqrt->Eval(energy_point[ii]);//is equal to y[ii]=p1[i]*sqrt(x[ii])+p0[i];
			//dy[ii]=sqrt(sqrt(x[ii])*sqrt(x[ii])*p1err[i]*p1err[i]+p0err[i]*p0err[i]);
			//dy[ii]=err_point[ii];
			//cout<<"	"<<energy_point[ii]<<"	"<<err_point[ii]<<endl;
			errylow[ii]=func_sqrt->Eval(energy_point[ii])-err_point[ii];
			erryhigh[ii]=func_sqrt->Eval(energy_point[ii])+err_point[ii];
		}
		graph4bandlow[i]= new TGraph(Npoint,energy_point,errylow);//画error bars TGraph(n,x,y,ex,ey);
		graph4bandlow[i]->SetLineColor(2);
		graph4bandlow[i]->SetLineWidth(1);
		graph4bandlow[i]->Draw("C");
		graph4bandhigh[i]= new TGraph(Npoint,energy_point,erryhigh);//画error bars TGraph(n,x,y,ex,ey);
		graph4bandhigh[i]->SetLineColor(2);
		graph4bandhigh[i]->SetLineWidth(1);
		graph4bandhigh[i]->Draw("C");

		TPaveText *textpol4 = new TPaveText(0.10,0.64,0.34,0.89,"brNDC");//left, down, right, up
		textpol4->SetBorderSize(1);
		textpol4->SetFillColor(0);
		textpol4->SetTextAlign(12);//align = 10*HorizontalAlign + VerticalAlign, 12 means水平左对齐、垂直居中对齐
		textpol4->SetTextFont(132);//font = 10 * fontID + precision, 12+2,12 means Symbol; 13+2, 13 means News Times Roman
		sprintf(paraprint,"y=p0*sqrt(x)");
		textpol4->AddText(paraprint);
		sprintf(paraprint,"p0=%.5f+/-%.5f",p0[i],p0err[i]);
		textpol4->AddText(paraprint);
		sprintf(paraprint,"Chisquare=%.2f",parchi[i]);
		textpol4->AddText(paraprint);
		sprintf(paraprint,"NDF=%.0f",parNDF[i]);
		textpol4->AddText(paraprint);
		sprintf(paraprint,"p-value=%e",p_value[i]);
		textpol4->AddText(paraprint);
		textpol4->Draw();
		//sprintf(hcali_name,"%s.png",hcali_name);
		sprintf(hcali_name,"%s%s%s%s%s","D:/X/out/responseforDoppler/",hcali_name,"_",tauflag,"_band.png");
		canvascali[i]->SaveAs(hcali_name);//存图
		outfile2<<"	SeGA_"<<i<<"	y=p0*sqrt(x)	taup0=	"<<p0[i]<<"	+/-	"<<p0err[i]<<"	Chi2=	"<<parchi[i]<<"	ndf=	"<<parNDF[i]<<"	p-value=	"<<p_value[i]<<endl;
*/

	}//for (i=0;i<ID;i++)
	//}//for(irawroot=runstart; irawroot<=runstop; irawroot++)
}//peakcali main