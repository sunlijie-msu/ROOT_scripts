#include <iostream>
#include <fstream>
#include <iomanip.h>
#include <math.h>
#include <time.h>
#include <map>
#include <TROOT.h>
#include <TSystem.h>
#include <TStyle.h>
#include <TCutG.h>
#include "TChain.h"
#include "TStyle.h"
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
void GausUpperLimit()// construct a Gauss(mean,sigma) function and obtain the 90% confidence level upper limit
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
	const int ID2=0;// no need to change i<=ID2//which detector

	int binwidth=1;
	//	int binnum=(upperlimit-lowerlimit)/binwidth;
	int minrange=0,maxrange=0,minbin,maxbin;
	double sigma=1,thresh=0.8,sigmab=20;//adjust
	double gaplow=70.,gaphigh=70.;//fitting range随分辨不同调整//no need to change
	unsigned long i;
	int jj,ii,k,ipeak;
	char paraprint[100],o_name[200],b_name[200],histo_name[200],h_name[200],hfit_name[200],hcali_name[200],tauflag[200];
	char command1[100],command2[100];
	//string paraprint[30],h_name[80],histo_name[80],b_name[80],hfit_name[80],hcali_name[80];
	TCanvas *canvaspeak[300];
	TCanvas *canvascali[300];
	TH1F *histo[ID2+1];//TH1F peak search+gauss fit,creat histograms
	const int nummax=100;//search peaknumber can be many, fitnumber can be limited to 3//no need to change
	int peaknum=nummax;
	TF1 *total[nummax];//creat function
	TF1 *p[nummax], *g[nummax], *b[nummax];
	TH1F *h1[300];
	TGraph *graph1[300], *graph1bandhigh[300], *graph1bandlow[300];//TGraph
	TGraph *graph2[300], *graph2bandhigh[300], *graph2bandlow[300];//TGraph
	TGraph *graph3[300], *graph3bandhigh[300], *graph3bandlow[300];//TGraph
	TGraph *graph4[300], *graph4bandhigh[300], *graph4bandlow[300];//TGraph
	double Mean[nummax]={0},meanerr[ID2+1][nummax];
	double SigmaLow[nummax]={0},peakyerr[nummax];
	double SigmaUp[nummax]={0},Chierr[ID2+1][nummax];
	double tau[ID2+1][nummax],tauerr[ID2+1][nummax];
	double range_sig1[ID2+1], range_sig2[ID2+1], range_tau1[ID2+1], range_tau2[ID2+1];
	// 	double *Mean;//if you don't know how many peaks will be found, use this
	// 	double *peaky;//if you don't know how many peaks will be found, use this
//	double energylit[6]={451.7, 493.3, 944.9, 1460.820, 1612.4, 2614.511};//25Si
//	double energyliterr[6]={0.5, 0.7, 0.5, 0.005, 0.5, 0.010};//25Si
// 	double energylit[6]={450.7, 1460.820, 1599, 2614.511, 2908, 7801};//23Al
// 	double energyliterr[6]={0.15, 0.005, 2, 0.010, 3, 2};//23Al

// 	const int Npoint=4;
// 	double energy_point[Npoint] = {1369,2754,2870,4238};// set location of point for single value
// 	double err_point[Npoint];  // error on the function at point x0 for single value
	
// 	const int Npoint=230;
// 	double energy_point[Npoint];// set location of point for drawing
// 	double err_point[Npoint];  // error on the function at point x0 for drawing
// 	for(jj=0;jj<Npoint;jj++){	energy_point[jj]=jj*40;}

	double peakmin;
	double chtemp;
	double par[nummax][5],parerr[nummax][5];//随峰数不同改动par[peaknum][3]
	double parchi[ID2+1],parNDF[ID2+1],p_value[ID2+1];
	int Eg=1248; //DSLmodify always finite value
	//1248, 2234, 3076, 4971, 5156, 5141, 4156, 4045, 3435, 2186, 2838, 4270, 3541, 5294
	//  49,     47,      46,     42,     42,     39,      39,     39,      45,      45,     44,     41,     40,      39
	int Ea=0; // no need to change
	double Egv=0; // no need to change
	char flag[20]; // no need to change
	sprintf(flag,"base"); // no need to change

	ipeak=10; //choose a peak DSLmodify
	if (ipeak==0) { Eg=1248; Ea=49; Egv=1248; minrange=700; maxrange=5700; }
	if (ipeak==1) { Eg=2234; Ea=47; Egv=2234; minrange=120; maxrange=500; }
	if (ipeak==2) { Eg=3076; Ea=46; Egv=3076; minrange=0; maxrange=30; }
	if (ipeak==3) { Eg=4971; Ea=42; Egv=4970; minrange=0; maxrange=30; }
	if (ipeak==4) { Eg=5156; Ea=42; Egv=5156; minrange=0; maxrange=30; }
	if (ipeak==5) { Eg=5141; Ea=39; Egv=5141; minrange=0; maxrange=30; }
	if (ipeak==6) { Eg=4156; Ea=39; Egv=4156; minrange=0; maxrange=30; }
	if (ipeak==7) { Eg=3435; Ea=45; Egv=3435; minrange=1; maxrange=50; }
	if (ipeak==8) { Eg=2186; Ea=45; Egv=2186; minrange=0; maxrange=30; }
	if (ipeak==9) { Eg=2838; Ea=44; Egv=2838; minrange=0; maxrange=30; }
	if (ipeak==10) { Eg=4270; Ea=41; Egv=4270; minrange=0; maxrange=50; }
	if (ipeak==11) { Eg=3541; Ea=40; Egv=3541; minrange=0; maxrange=30; }
	if (ipeak==12) { Eg=5294; Ea=39; Egv=5293; minrange=0; maxrange=50; }

	
	sprintf(o_name,"%s%d%s%d%s%.0f%s","D:/X/out/DSL/UL/Chi2_Gamma",Eg,"_Ea",Ea,"_Egv",Egv,"slice_vTau.dat");
	cout<<o_name<<endl;
	ifstream infile(o_name,ios::in);//The data that need to be fitted
	ii=0;
	string line;
	while ( getline(infile, line) )
	{
		stringstream(line)>>Mean[ii]>>SigmaLow[ii]>>SigmaUp[ii];
		cout<<Mean[ii]<<'	'<<SigmaUp[ii]<<endl;
		ii++;
	}
	sprintf(b_name,"%s%d%s%d%s%.0f%s","D:/X/out/DSL/UL/Gaus_Gamma",Eg,"_Ea",Ea,"_Egv",Egv,"slice_vTau.dat");
	cout<<b_name<<endl;
	ofstream outfile2(b_name,ios::out);
//	ifstream infile("C:/Si24/Si22peakcali/Clover.dat",ios::in);//The data that need to be fitted
//	int lowerlimit[ID2+1][nummax],upperlimit[ID2+1][nummax];
	for(int irun=0;irun<=5;irun++)//loop each Gamma's Tau to do Chi minimization fit
	{// DSLmodify choose a flag
		if (irun==0) sprintf(flag,"base");
		if (irun==1) sprintf(flag,"bkglow");
		if (irun==2) sprintf(flag,"bkghigh");
		if (irun==3) sprintf(flag,"SP0.9");
		if (irun==4) sprintf(flag,"SP1.1");
		if (irun==5) sprintf(flag,"AC2-0.8");
		if (irun==6) sprintf(flag,"respsigma1");
		if (irun==7) sprintf(flag,"respsigma-1");
		if (irun==8) sprintf(flag,"resptau1");
		if (irun==9) sprintf(flag,"resptau-1");


		sprintf(hcali_name,"%s%d%s%d%s%.0f%s%s","Gaus_Gamma",Eg,"_Ea",Ea,"_Egv",Egv,"slice_vTau_",flag);
		canvascali[irun]=new TCanvas(hcali_name,hcali_name,900,600);//建立画布
		canvascali[irun]->cd();//进入画布
		//canvascali[i]->SetGrid();//显示网格
		TF1 *gaus=new TF1("gaus","gaus", 0,5800);
		gaus->SetParameters(100000,Mean[irun],SigmaUp[irun]);
		gaus->SetNpx(50000);
		gaus->Draw();
		gaus->GetXaxis()->SetRangeUser(0,maxrange);
		gaus->GetXaxis()->SetTitle("#tau (fs)");
		gaus->GetYaxis()->SetTitle("A.U.");
		gaus->GetXaxis()->CenterTitle();
		gaus->GetYaxis()->CenterTitle();
		gPad->SetRightMargin(0.03);

// 		for (float ibin=0;ibin<90;ibin+=0.01)
// 		{
// 			if (gaus->Integral(0,ibin)/gaus->Integral(0,maxrange)>=0.1) break;
// 		}
// 		outfile2<<ibin;
// 		for (float ibin=0;ibin<90;ibin+=0.01)
// 		{
// 			if (gaus->Integral(0,ibin)/gaus->Integral(0,maxrange)>=0.5) break;
// 		}
// 		outfile2<<"	"<<ibin;
// 		for (float ibin=0;ibin<90;ibin+=0.01)
// 		{
// 			if (gaus->Integral(0,ibin)/gaus->Integral(0,maxrange)>=0.9) break;
// 		}
// 		outfile2<<"	"<<ibin<<endl;

// 		for (float ibin=0;ibin<6000;ibin+=1)
// 		{
// 			if (gaus->Integral(0,ibin)/gaus->Integral(0,maxrange)>=0.17) break;
// 		}
// 		outfile2<<ibin;
// 		for (float ibin=0;ibin<6000;ibin+=1)
// 		{
// 			if (gaus->Integral(0,ibin)/gaus->Integral(0,maxrange)>=0.5) break;
// 		}
// 		outfile2<<"	"<<ibin;
// 		for (float ibin=0;ibin<6000;ibin+=1)
// 		{
// 			if (gaus->Integral(0,ibin)/gaus->Integral(0,maxrange)>=0.83) break;
// 		}
// 		outfile2<<"	"<<ibin<<endl;

		for (float ibin=0;ibin<600;ibin+=0.1)
		{
			if (gaus->Integral(0,ibin)/gaus->Integral(0,maxrange)>=0.17) break;
		}
		outfile2<<ibin;
		for (float ibin=0;ibin<600;ibin+=0.1)
		{
			if (gaus->Integral(0,ibin)/gaus->Integral(0,maxrange)>=0.5) break;
		}
		outfile2<<"	"<<ibin;
		for (float ibin=0;ibin<600;ibin+=0.1)
		{
			if (gaus->Integral(0,ibin)/gaus->Integral(0,maxrange)>=0.83) break;
		}
		outfile2<<"	"<<ibin<<endl;

		TPaveText *textpol6 = new TPaveText(0.30,0.70,0.95,0.90,"brNDC");//left, down, right, up
		textpol6->SetBorderSize(1);
		textpol6->SetFillColor(0);
		textpol6->SetTextAlign(12);//align = 10*HorizontalAlign + VerticalAlign, 12 means水平左对齐、垂直居中对齐
		textpol6->SetTextFont(132);//font = 10 * fontID + precision, 12+2,12 means Symbol; 13+2, 13 means News Times Roman
		sprintf(paraprint,"UL=%.2f%s",ibin," fs");
		textpol6->AddText(paraprint);
		sprintf(paraprint, "%s", hcali_name);
		textpol6->AddText(paraprint);
		textpol6->Draw();
		//sprintf(hcali_name,"%s.png",hcali_name);
		sprintf(hcali_name,"%s%s%s","D:/X/out/Chi2/",hcali_name,".png");
		canvascali[irun]->SaveAs(hcali_name);//存图

	}//irun
	//}//for(irawroot=runstart; irawroot<=runstop; irawroot++)
}//peakcali main