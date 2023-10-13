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
void response_variation()// how does response function look like using different sigma and tau
{
	time_t tim;
	struct tm *at;
	char now[80];
	int irawroot;
	char rawrootname[300];
	double slope[300],slopeerr[300];
	double intercept[300],intercepterr[300];
	const int ID1=0;// i=ID1//which detector
	const int ID2=15;// i<=ID2//which detector

	int binwidth=1;
	//	int binnum=(upperlimit-lowerlimit)/binwidth;
	int minrange=0,maxrange=0,minbin,maxbin;
	double sigma=1,thresh=0.8,sigmab=20;//adjust
	float gaplow=10.,gaphigh=5.;//fitting range随分辨不同调整
	unsigned long i;
	int jj,ii,ibin,k,ipeak;
	char paraprint[30],b_name[200],histo_name[200],h_name[200],hfit_name[200],hcali_name[200],tauflag[200];
	char command1[100],command2[100];
	//string paraprint[30],h_name[80],histo_name[80],b_name[80],hfit_name[80],hcali_name[80];
	TCanvas *canvaspeak[300];
	TCanvas *canvascali[300];
	TH1F *histo[ID2+1];//TH1F peak search+gauss fit,creat histograms
	const int nummax=100;//search peaknumber can be many, fitnumber can be limited to 3//no need to modify
	int peaknum=nummax;
	TF1 *total[nummax];//creat function
	TF1 *p[nummax], *g[nummax], *b[nummax];
	TGraph *graph1[300], *graph2[300];//TGraph刻度
	Double_t peakx[nummax],peakxerr[nummax];
	Double_t peaky[nummax],peakyerr[nummax];
	Double_t sig[nummax],sigerr[nummax];
	Double_t tau[nummax],tauerr[nummax];
	// 	float *peakx;//if you don't know how many peaks will be found, use this
	// 	float *peaky;//if you don't know how many peaks will be found, use this
//	double energylit[6]={451.7, 493.3, 944.9, 1460.820, 1612.4, 2614.511};//25Si
//	double energyliterr[6]={0.5, 0.7, 0.5, 0.005, 0.5, 0.010};//25Si
	double energylit[6]={450.7, 1460.820, 1599, 2614.511, 2908, 7801};//23Al
	double energyliterr[6]={0.15, 0.005, 2, 0.010, 3, 2};//23Al
	float peakmin;
	float chtemp;
	double par[nummax][5],parerr[nummax][5];//随峰数不同改动par[peaknum][3]//adjust
	double parchi[nummax],parNDF[nummax];

	for(i=0;i<=0;i++)//which SeGA detector modify =0<=15. 0-10 are usually pretty good to fit, 12,13 are troublesome
	{
		sprintf(hfit_name,"%s%d","Response_test_",i);
		canvaspeak[i]=new TCanvas(hfit_name,hfit_name,1200,600);//建立画布
		canvaspeak[i]->cd();//进入画布
// 		histo[i]->SetTitle(hfit_name);//图名
// 		histo[i]->GetXaxis()->SetTitle("Energy");//轴名
// 		histo[i]->GetYaxis()->SetTitle("Counts");//轴名
// 		histo[i]->GetXaxis()->CenterTitle();//居中
// 		histo[i]->GetYaxis()->CenterTitle();//居中
// 		histo[i]->GetXaxis()->SetLabelFont(132);//坐标字体
// 		histo[i]->GetYaxis()->SetLabelFont(132);//坐标字体
// 		histo[i]->GetXaxis()->SetTitleFont(132);//轴名字体
// 		histo[i]->GetYaxis()->SetTitleFont(132);//轴名字体
// 		histo[i]->GetXaxis()->SetTitleOffset(1.2);//轴名偏移
// 		histo[i]->GetYaxis()->SetTitleOffset(1.3);//轴名偏移		
// 		histo[i]->Draw();
		for(ii=0; ii<=29; ii++)//which peak in one SeGA detector modify =0<=5
		{
// 			if(ii==0)sig[ii]=0.1;
// 			else sig[ii]=sig[ii-1]+0.2;
// 			if(ii==0)tau[ii]=1;
// 			else tau[ii]=tau[ii-1];
			if(ii==0)sig[ii]=1;
			else sig[ii]=sig[ii-1];
			if(ii==0)tau[ii]=0.1;
			else tau[ii]=tau[ii-1]+0.2;
		}
		for(ii=0; ii<=29; ii++)//which peak in one SeGA detector modify =0<=5
		{
			peakx[ii]=2000;
			total[ii]=new TF1("total","[0]*x+[1]+sqrt(3.141592654/2)*[2]/[3]*[4]*exp(0.5*([4]*[4]/([3]*[3]))+(x-[5])/[3])*ROOT::Math::erfc(1/sqrt(2)*([4]/[3]+(x-[5])/[4]))",peakx[ii]-gaplow,peakx[ii]+gaphigh);//高斯拟合，调节合适的拟合range / Glassman PRC2019 low-energy tail
			g[ii]=new TF1("g","gausn",peakx[ii]-gaplow,peakx[ii]+gaphigh);// The [2]-N parameter in total is equivalent to the Constant in gausn
			p[ii]=new TF1("p","[0]*x+[1]-[0]*x-[1]+sqrt(3.141592654/2)*[2]/[3]*[4]*exp(0.5*([4]*[4]/([3]*[3]))+(x-[5])/[3])*ROOT::Math::erfc(1/sqrt(2)*([4]/[3]+(x-[5])/[4]))",peakx[ii]-gaplow,peakx[ii]+gaphigh);//pure peak
			b[ii]=new TF1("b","[0]*x+[1]",peakx[ii]-gaplow,peakx[ii]+gaphigh);//pure bkg
			total[ii]->SetNpx((gaphigh+gaplow)*100);
			g[ii]->SetNpx((gaphigh+gaplow)*100);
			p[ii]->SetNpx((gaphigh+gaplow)*100);
			b[ii]->SetNpx((gaphigh+gaplow)*100);
			
			total[ii]->SetParameters(0,0,1000,tau[ii],sig[ii],2000);//initial value [0]-A, [1]-B, [2]-N, [3]-τ, [4]-σ, [5]-μ
//			total[ii]->SetParameters(0,500,150000,0.7,1.5,3480);//initial value [0]-A, [1]-B, [2]-N, [3]-λ, [4]-σ, [5]-μ  //Glassman
// 			total[ii]->SetParLimits(0,-500,500);//Bkg A  //Glassman
// 			total[ii]->SetParLimits(1,-50000,300000);//Bkg B  //Glassman
// 			total[ii]->SetParLimits(2,3,60000);//Constant,min,max  //Glassman
// 			total[ii]->SetParLimits(3,0.0001,4);//Tau  //Glassman portal
// 			total[ii]->SetParLimits(3,0.7,0.7);//Tau  //Glassman portal
// 			total[ii]->SetParLimits(4,0.4,20);//Sigma  //Glassman
// 			total[ii]->SetParLimits(5,440,7840);//Mean  //Glassman
			total[ii]->SetParNames("BkgA","BkgB","Const*bin","Tau","Sigma","Mean");
			if(ii==0){	total[ii]->SetLineColor(1); total[ii]->Draw();}
			else total[ii]->Draw("same");
		}//for(ii=0;ii<peaknum;ii++)
	}//for (i=0;i<ID;i++)
}//peakcali main