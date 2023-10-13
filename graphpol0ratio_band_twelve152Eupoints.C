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
void graphpol0ratio_band_twelve152Eupoints()// graph-pol0 fit find the scaling factor by fitting the ratios of SeGA efficiency data points to simulation points.
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
	const int nummax=12;//search peaknumber can be many, fitnumber can be limited to 3//no need to modify
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

	const int Npoint=12;
	double energy_point[Npoint] = {122,245,344,411,444,779,867,964,1112,1213,1299,1408};// all NuDat values, set location of point for single value modify
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
//	
//	par[0][7]=1;	parerr[0][7]=0.1;//C useless for now

	double parchi[ID2+1],parNDF[ID2+1],p_value[ID2+1];

	// data points / simulated points upstream
// 	energy[0][0]=245;	energyerr[0][0]=0.0008;	eff[0][0]=1.1334917;	efferr[0][0]=0.0135023;
// 	energy[0][1]=344;	energyerr[0][1]=0.0012;	eff[0][1]=1.2056473;	efferr[0][1]=0.0100275;
// 	energy[0][2]=411;	energyerr[0][2]=0.0012;	eff[0][2]=1.1980085;	efferr[0][2]=0.0287613;
// 	energy[0][3]=444;	energyerr[0][3]=0.003;	eff[0][3]=1.2499955;	efferr[0][3]=0.0252096;
// 	energy[0][4]=779;	energyerr[0][4]=0.0024;	eff[0][4]=1.2898608;	efferr[0][4]=0.0155404;
// 	energy[0][5]=867;	energyerr[0][5]=0.003;	eff[0][5]=1.2460410;	efferr[0][5]=0.0258933;
// 	energy[0][6]=964;	energyerr[0][6]=0.018;	eff[0][6]=1.2867404;	efferr[0][6]=0.0150528;
// 	energy[0][7]=1112;	energyerr[0][7]=0.003;	eff[0][7]=1.3181519;	efferr[0][7]=0.0166318;
// 	energy[0][8]=1213;	energyerr[0][8]=0.0011;	eff[0][8]=1.3519620;	efferr[0][8]=0.0561612;
// 	energy[0][9]=1299;	energyerr[0][9]=0.008;	eff[0][9]=1.2377750;	efferr[0][9]=0.0430786;
// 	energy[0][10]=1408;	energyerr[0][10]=0.003;	eff[0][10]=1.3153758;	efferr[0][10]=0.0146292;


	// data points / simulated points center
// 	energy[0][0]=245;	energyerr[0][0]=0.0008;	eff[0][0]=1.1424415;	efferr[0][0]=0.0099944;
// 	energy[0][1]=344;	energyerr[0][1]=0.0012;	eff[0][1]=1.1537502;	efferr[0][1]=0.0081442;
// 	energy[0][2]=411;	energyerr[0][2]=0.0012;	eff[0][2]=1.1289797;	efferr[0][2]=0.0157035;
// 	energy[0][3]=444;	energyerr[0][3]=0.003;	eff[0][3]=1.1507738;	efferr[0][3]=0.0138569;
// 	energy[0][4]=779;	energyerr[0][4]=0.0024;	eff[0][4]=1.1573512;	efferr[0][4]=0.0099531;
// 	energy[0][5]=867;	energyerr[0][5]=0.003;	eff[0][5]=1.1364312;	efferr[0][5]=0.0148093;
// 	energy[0][6]=964;	energyerr[0][6]=0.018;	eff[0][6]=1.1556650;	efferr[0][6]=0.0095907;
// 	energy[0][7]=1112;	energyerr[0][7]=0.003;	eff[0][7]=1.2079218;	efferr[0][7]=0.0107188;
// 	energy[0][8]=1213;	energyerr[0][8]=0.0011;	eff[0][8]=1.1120834;	efferr[0][8]=0.0259362;
// 	energy[0][9]=1299;	energyerr[0][9]=0.008;	eff[0][9]=1.0842734;	efferr[0][9]=0.0219037;
// 	energy[0][10]=1408;	energyerr[0][10]=0.003;	eff[0][10]=1.1809978;	efferr[0][10]=0.0096613;


	// data points / simulated points downstream
// 	energy[0][0]=245;	energyerr[0][0]=0.0008;	eff[0][0]=0.9708478;	efferr[0][0]=0.0113462;
// 	energy[0][1]=344;	energyerr[0][1]=0.0012;	eff[0][1]=1.0358245;	efferr[0][1]=0.0085228;
// 	energy[0][2]=411;	energyerr[0][2]=0.0012;	eff[0][2]=1.0232191;	efferr[0][2]=0.0235027;
// 	energy[0][3]=444;	energyerr[0][3]=0.003;	eff[0][3]=1.0401022;	efferr[0][3]=0.0202101;
// 	energy[0][4]=779;	energyerr[0][4]=0.0024;	eff[0][4]=1.0831699;	efferr[0][4]=0.0126985;
// 	energy[0][5]=867;	energyerr[0][5]=0.003;	eff[0][5]=1.0980590;	efferr[0][5]=0.0224501;
// 	energy[0][6]=964;	energyerr[0][6]=0.018;	eff[0][6]=1.0985774;	efferr[0][6]=0.0125971;
// 	energy[0][7]=1112;	energyerr[0][7]=0.003;	eff[0][7]=1.1439810;	efferr[0][7]=0.0141348;
// 	energy[0][8]=1213;	energyerr[0][8]=0.0011;	eff[0][8]=1.1369743;	efferr[0][8]=0.0455561;
// 	energy[0][9]=1299;	energyerr[0][9]=0.008;	eff[0][9]=1.0271111;	efferr[0][9]=0.0338235;
// 	energy[0][10]=1408;	energyerr[0][10]=0.003;	eff[0][10]=1.1402474;	efferr[0][10]=0.0124771;

	// data points / simulated points center Timi
// 	energy[0][0]=122;	energyerr[0][0]=0.0003;	eff[0][0]=1.2044;	efferr[0][0]=0.0083;
// 	energy[0][1]=245;	energyerr[0][1]=0.0008;	eff[0][1]=1.1293030;	efferr[0][1]=0.0098630;
// 	energy[0][2]=344;	energyerr[0][2]=0.0012;	eff[0][2]=1.1455815;	efferr[0][2]=0.0080807;
// 	energy[0][3]=411;	energyerr[0][3]=0.0012;	eff[0][3]=1.0907080;	efferr[0][3]=0.0149309;
// 	energy[0][4]=444;	energyerr[0][4]=0.003;	eff[0][4]=1.1369190;	efferr[0][4]=0.0136259;
// 	energy[0][5]=779;	energyerr[0][5]=0.0024;	eff[0][5]=1.1494941;	efferr[0][5]=0.0098670;
// 	energy[0][6]=867;	energyerr[0][6]=0.003;	eff[0][6]=1.1127502;	efferr[0][6]=0.0143527;
// 	energy[0][7]=964;	energyerr[0][7]=0.018;	eff[0][7]=1.1571845;	efferr[0][7]=0.0095998;
// 	energy[0][8]=1112;	energyerr[0][8]=0.003;	eff[0][8]=1.2004595;	efferr[0][8]=0.0106312;
// 	energy[0][9]=1213;	energyerr[0][9]=0.0011;	eff[0][9]=1.0902919;	efferr[0][9]=0.0251024;
// 	energy[0][10]=1299;	energyerr[0][10]=0.008;	eff[0][10]=1.1259356;	efferr[0][10]=0.0231225;
// 	energy[0][11]=1408;	energyerr[0][11]=0.003;	eff[0][11]=1.1587510;	efferr[0][11]=0.0094370;

	// data points / simulated points center Mike
	energy[0][0]=122;	energyerr[0][0]=0.0003;	eff[0][0]=1.3433;	efferr[0][0]=0.0093;
	energy[0][1]=245;	energyerr[0][1]=0.0008;	eff[0][1]=1.1457611;	efferr[0][1]=0.0100234;
	energy[0][2]=344;	energyerr[0][2]=0.0012;	eff[0][2]=1.1571026;	efferr[0][2]=0.0081679;
	energy[0][3]=411;	energyerr[0][3]=0.0012;	eff[0][3]=1.1322601;	efferr[0][3]=0.0157492;
	energy[0][4]=444;	energyerr[0][4]=0.003;	eff[0][4]=1.1541176;	efferr[0][4]=0.0138972;
	energy[0][5]=779;	energyerr[0][5]=0.0024;	eff[0][5]=1.1607141;	efferr[0][5]=0.0099820;
	energy[0][6]=867;	energyerr[0][6]=0.003;	eff[0][6]=1.1397333;	efferr[0][6]=0.0148523;
	energy[0][7]=964;	energyerr[0][7]=0.018;	eff[0][7]=1.1590229;	efferr[0][7]=0.0096186;
	energy[0][8]=1112;	energyerr[0][8]=0.003;	eff[0][8]=1.2114315;	efferr[0][8]=0.0107500;
	energy[0][9]=1213;	energyerr[0][9]=0.0011;	eff[0][9]=1.1153148;	efferr[0][9]=0.0260116;
	energy[0][10]=1299;	energyerr[0][10]=0.008;	eff[0][10]=1.0874239;	efferr[0][10]=0.0219674;
	energy[0][11]=1408;	energyerr[0][11]=0.003;	eff[0][11]=1.1844294;	efferr[0][11]=0.0096893;

	// Timi simulated points / Mike simulated points center
// 	energy[0][0]=122;	energyerr[0][0]=0.0003;	eff[0][0]=1.1170;	efferr[0][0]=0.0080;
// 	energy[0][1]=245;	energyerr[0][1]=0.0008;	eff[0][1]=1.0176832;	efferr[0][1]=0.0094732;
// 	energy[0][2]=344;	energyerr[0][2]=0.0012;	eff[0][2]=1.0091053;	efferr[0][2]=0.0075579;
// 	energy[0][3]=411;	energyerr[0][3]=0.0012;	eff[0][3]=1.0362483;	efferr[0][3]=0.0175527;
// 	energy[0][4]=444;	energyerr[0][4]=0.003;	eff[0][4]=1.0133515;	efferr[0][4]=0.0147300;
// 	energy[0][5]=779;	energyerr[0][5]=0.0024;	eff[0][5]=1.0124184;	efferr[0][5]=0.0098175;
// 	energy[0][6]=867;	energyerr[0][6]=0.003;	eff[0][6]=1.0254849;	efferr[0][6]=0.0156168;
// 	energy[0][7]=964;	energyerr[0][7]=0.018;	eff[0][7]=1.0047404;	efferr[0][7]=0.0097323;
// 	energy[0][8]=1112;	energyerr[0][8]=0.003;	eff[0][8]=1.0125280;	efferr[0][8]=0.0103641;
// 	energy[0][9]=1213;	energyerr[0][9]=0.0011;	eff[0][9]=1.0227239;	efferr[0][9]=0.0290305;
// 	energy[0][10]=1299;	energyerr[0][10]=0.008;	eff[0][10]=0.9770162;	efferr[0][10]=0.0245980;
// 	energy[0][11]=1408;	energyerr[0][11]=0.003;	eff[0][11]=1.0246378;	efferr[0][11]=0.0095773;

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
		pad1->SetTopMargin(0.04);
		pad1->SetRightMargin(0.03);
		pad1->SetLeftMargin(0.08);
		pad1->SetBottomMargin(0.10);
		//pad1->SetBorderMode(0);

		pad2->SetTopMargin(0.01);
		pad2->SetRightMargin(0.03);
		pad2->SetLeftMargin(0.08);
		pad2->SetBottomMargin(0.23);
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
		graph1[i]->GetYaxis()->SetTitle("Ratio of data/simulation");//轴名
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
		graph1[i]->GetXaxis()->SetRangeUser(100,1530);
//		graph1[i]->GetYaxis()->SetRangeUser(0.9,1.5);
		graph1[i]->SetMarkerStyle(21);
		graph1[i]->SetMarkerColor(1);
		TF1 *pol0 = new TF1("pol0","pol0",10,10000);
		pol0->SetParNames("p0");
		//graph1[i]->Fit("pol0");//pol0 can be used directly without TF1 constructor in CINT

		TFitResultPtr r_sig = graph1[i]->Fit("pol0","MS");//"S" means the result of the fit is returned in the TFitResultPtr
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
			residual[i][ii]=(pol0->Eval(energy_point[ii])-eff[i][ii])/pol0->Eval(energy_point[ii]);
			residualerr[i][ii]=efferr[i][ii]/pol0->Eval(energy_point[ii]);
			outfile<<"Eg=	"<<energy_point[ii]<<"	relative residual=	"<<residual[i][ii]<<endl;
		}
		graph1[i]->Draw("AP");//"A": Axis are drawn around the graph, "P": The current marker is plotted at each point
		hint_sig->Draw("e3 same");
		graph1[i]->Draw("P same");//draw the points again above error band

		p0[i]=pol0->GetParameter(0);
		
		p0err[i]=pol0->GetParError(0);
//		p1err[i]=pol0->GetParError(1);
		parchi[i]=pol0->GetChisquare();
		parNDF[i]=pol0->GetNDF();
		p_value[i]=pol0->GetProb();//This probability is not the “probability that your fit is good.” If you did many fake experiments (draw many random samples of data points from the assumed distribution (your fit function)), this is the percentage of experiments that would give χ2 values ≥ to the one you got in this experiment.
// 		TVirtualFitter * fitter = TVirtualFitter::GetFitter();
// 		p7[i]=fitter->GetParameter(0);
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
// 			//y[ii]=pol0->Eval(energy_point[ii]);//is equal to y[ii]=p1[i]*x[ii]+p7[i];
// 			//dy[ii]=sqrt(x[ii]*x[ii]*p1err[i]*p1err[i]+p7err[i]*p7err[i]+p1[i]*p1[i]*dx[ii]*dx[ii]);
// 			//dy[ii]=err_point[ii];
// 			//cout<<"	"<<energy_point[ii]<<"	"<<err_point[ii]<<endl;
// 			errylow[ii]=pol0->Eval(energy_point[ii])-err_point[ii];
// 			erryhigh[ii]=pol0->Eval(energy_point[ii])+err_point[ii];
// 		}
// 		graph1bandlow[i]= new TGraph(Npoint,energy_point,errylow);//画error bars TGraph(n,x,y,ex,ey);
// 		graph1bandlow[i]->SetLineColor(2);
// 		graph1bandlow[i]->SetLineWidth(1);
// 		graph1bandlow[i]->Draw("C");//"C" A smooth Curve is drawn
// 		graph1bandhigh[i]= new TGraph(Npoint,energy_point,erryhigh);//画error bars TGraph(n,x,y,ex,ey);
// 		graph1bandhigh[i]->SetLineColor(2);
// 		graph1bandhigh[i]->SetLineWidth(1);
// 		graph1bandhigh[i]->Draw("C");//"C" A smooth Curve is drawn

		TPaveText *textpol1 = new TPaveText(0.69,0.64,0.96,0.95,"brNDC");//left, down, right, up
		textpol1->SetBorderSize(1);
		textpol1->SetFillColor(0);
		textpol1->SetTextAlign(12);//align = 10*HorizontalAlign + VerticalAlign, 12 means水平左对齐、垂直居中对齐
		textpol1->SetTextFont(132);//font = 10 * fontID + precision, 12+2,12 means Symbol; 13+2, 13 means News Times Roman
		sprintf(paraprint,"y=C*logpol6");
		textpol1->AddText(paraprint);
		sprintf(paraprint,"C=%.5f+/-%.5f",p0[i],p0err[i]);
		textpol1->AddText(paraprint);
		sprintf(paraprint,"Chi2=%.2f",parchi[i]);
		textpol1->AddText(paraprint);
		sprintf(paraprint,"NDF=%.0f",parNDF[i]);
		textpol1->AddText(paraprint);
		sprintf(paraprint,"p-val=%e",p_value[i]);
		textpol1->AddText(paraprint);
		textpol1->Draw();
		//sprintf(hcali_name,"%s.png",hcali_name);
		sprintf(hcali_name,"%s%s%s%s%s","D:/X/out/responseforDoppler/",hcali_name,"_",tauflag,"_band.png");
		canvascali[i]->SaveAs(hcali_name);//存图
		outfile2<<"SeGA_"<<i<<"	y=pol0*x	"<<endl;
		outfile2<<"p7=	"<<setprecision(8)<<p0[i]<<"	+/-	"<<p0err[i]<<endl;
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
		graph2[i]->GetXaxis()->SetRangeUser(100,1530);
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
		TLine *T1=new TLine(100,0,1530,0);
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
		func_sqrt_const->SetParNames("p7","p1");//y=p1*sqrt(x)+p7
//		graph3[i]->Fit("func_sqrt_const");//pol0 can be used directly without TF1 constructor in CINT

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

		p7[i]=func_sqrt_const->GetParameter(0);
		p7err[i]=func_sqrt_const->GetParError(0);
		p1[i]=func_sqrt_const->GetParameter(1);
		p1err[i]=func_sqrt_const->GetParError(1);
		parchi[i]=func_sqrt_const->GetChisquare();
		parNDF[i]=func_sqrt_const->GetNDF();
		p_value[i]=func_sqrt_const->GetProb();//This probability is not the “probability that your fit is good.” If you did many fake experiments (draw many random samples of data points from the assumed distribution (your fit function)), this is the percentage of experiments that would give χ2 values ≥ to the one you got in this experiment.
		double errylow[Npoint],erryhigh[Npoint];
		for(ii=0;ii<Npoint;ii++)
		{
			//x[ii]=energy_point[ii];
			//y[ii]=func_sqrt_const->Eval(energy_point[ii]);//is equal to y[ii]=p1[i]*sqrt(x[ii])+p7[i];
			//dy[ii]=sqrt(sqrt(x[ii])*sqrt(x[ii])*p1err[i]*p1err[i]+p7err[i]*p7err[i]);
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
		sprintf(paraprint,"y=p1*sqrt(x)+p7");
		textpol3->AddText(paraprint);
		sprintf(paraprint,"p1=%.5f+/-%.5f",p1[i],p1err[i]);
		textpol3->AddText(paraprint);
		sprintf(paraprint,"p7=%.5f+/-%.5f",p7[i],p7err[i]);
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
		outfile2<<"	SeGA_"<<i<<"	y=p1*sqrt(x)+p7	taup1=	"<<setprecision(8)<<p1[i]<<"	+/-	"<<p1err[i]<<"	taup0=	"<<p7[i]<<"	+/-	"<<p7err[i]<<"	Chi2=	"<<parchi[i]<<"	ndf=	"<<parNDF[i]<<"	p-value=	"<<p_value[i]<<endl;





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
		pol2->SetParNames("p7","p1","p2");//y=p2*x^2+p1*x+p7
//		graph2[i]->Fit("pol2");//pol0 can be used directly without TF1 constructor in CINT

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

		p7[i]=pol2->GetParameter(0);
		p1[i]=pol2->GetParameter(1);
		p2[i]=pol2->GetParameter(2);
		p7err[i]=pol2->GetParError(0);
		p1err[i]=pol2->GetParError(1);
		p2err[i]=pol2->GetParError(2);
		parchi[i]=pol2->GetChisquare();
		parNDF[i]=pol2->GetNDF();
		p_value[i]=pol2->GetProb();//This probability is not the “probability that your fit is good.” If you did many fake experiments (draw many random samples of data points from the assumed distribution (your fit function)), this is the percentage of experiments that would give χ2 values ≥ to the one you got in this experiment.
		double errylow[Npoint],erryhigh[Npoint];
		for(ii=0;ii<Npoint;ii++)
		{
			//x[ii]=energy_point[ii];
			//y[ii]=pol2->Eval(energy_point[ii]);//is equal to y[ii]=p2[i]*x[ii]*x[ii]+p1[i]*x[ii]+p7[i];
			//dy[ii]=sqrt(x[ii]*x[ii]*x[ii]*x[ii]*p2err[i]*p2err[i]+x[ii]*x[ii]*p1err[i]*p1err[i]+p7err[i]*p7err[i]+(2*p2[i]+p1[i])*(2*p2[i]+p1[i])*dx[ii]);
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
		sprintf(paraprint,"y=p2*x^2+p1*x+p7");
		textpol2->AddText(paraprint);
		sprintf(paraprint,"p2=%.9f+/-%.9f",p2[i],p2err[i]);
		textpol2->AddText(paraprint);
		sprintf(paraprint,"p1=%.5f+/-%.5f",p1[i],p1err[i]);
		textpol2->AddText(paraprint);
		sprintf(paraprint,"p7=%.5f+/-%.5f",p7[i],p7err[i]);
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
		outfile2<<"	SeGA_"<<i<<"	y=p2*x^2+p1*x+p7	taup2=	"<<setprecision(11)<<p2[i]<<"	+/-	"<<p2err[i]<<"	taup1=	"<<setprecision(8)<<p1[i]<<"	+/-	"<<p1err[i]<<"	taup0=	"<<p7[i]<<"	+/-	"<<p7err[i]<<"	Chi2=	"<<parchi[i]<<"	ndf=	"<<parNDF[i]<<"	p-value=	"<<p_value[i];





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
		func_sqrt->SetParNames("p7");//y=p7*sqrt(x)
//		graph4[i]->Fit("func_sqrt");//pol0 can be used directly without TF1 constructor in CINT

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

		p7[i]=func_sqrt->GetParameter(0);
		p7err[i]=func_sqrt->GetParError(0);
		parchi[i]=func_sqrt->GetChisquare();
		parNDF[i]=func_sqrt->GetNDF();
		p_value[i]=func_sqrt->GetProb();//This probability is not the “probability that your fit is good.” If you did many fake experiments (draw many random samples of data points from the assumed distribution (your fit function)), this is the percentage of experiments that would give χ2 values ≥ to the one you got in this experiment.
		double errylow[Npoint],erryhigh[Npoint];
		for(ii=0;ii<Npoint;ii++)
		{
			//x[ii]=energy_point[ii];
			//y[ii]=func_sqrt->Eval(energy_point[ii]);//is equal to y[ii]=p1[i]*sqrt(x[ii])+p7[i];
			//dy[ii]=sqrt(sqrt(x[ii])*sqrt(x[ii])*p1err[i]*p1err[i]+p7err[i]*p7err[i]);
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
		sprintf(paraprint,"y=p7*sqrt(x)");
		textpol4->AddText(paraprint);
		sprintf(paraprint,"p7=%.5f+/-%.5f",p7[i],p7err[i]);
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
		outfile2<<"	SeGA_"<<i<<"	y=p7*sqrt(x)	taup0=	"<<p7[i]<<"	+/-	"<<p7err[i]<<"	Chi2=	"<<parchi[i]<<"	ndf=	"<<parNDF[i]<<"	p-value=	"<<p_value[i]<<endl;
*/

	}//for (i=0;i<ID;i++)
	//}//for(irawroot=runstart; irawroot<=runstop; irawroot++)
}//peakcali main