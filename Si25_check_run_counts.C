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
void Si25_check_run_counts()// relative counts in each run.
{//also for the ratio of 401 proton/1369 gamma in peak areas
	time_t tim;
	struct tm *at;
	char now[80];
	int irawroot;
	char rawrootname[300];
	double slope[300],slopeerr[300];
	double intercept[300],intercepterr[300];
	const int ID1=0;// no need to modify i=ID1//which detector
	const int ID2=141;// no need to modify i<=ID2//which detector

	int binwidth=1;
	//	int binnum=(upperlimit-lowerlimit)/binwidth;
	float minrange=0,maxrange=0,minbin=0,maxbin=0;
	float minbinb=0,maxbinb=0;
	double sigma=1,thresh=0.8,sigmab=20;//adjust
	float gaplow=70.,gaphigh=70.;//fitting range随分辨不同调整//no need to modify
	int irawroot=0;
	int jj,ii,ibin,k,ipeak;
	char paraprint[30],b_name[200],histo_name[200],h_name[200],hfit_name[200],hcali_name[200],tauflag[200];
	char command1[100],command2[100];
	//string paraprint[30],h_name[80],histo_name[80],b_name[80],hfit_name[80],hcali_name[80];
	TCanvas *canvaspeak[300];
	TCanvas *canvascali[300];
	TH1F *histo[ID2+1];//TH1F peak search+gauss fit,creat histograms
	const int nummax=4;//search peaknumber can be many, fitnumber can be limited to 3//no need to modify
	TF1 *total[nummax];//creat function
	TF1 *p[nummax], *g[nummax], *b[nummax];
	TGraph *graph1[300], *graph2[300];//TGraph刻度
	int peaknum=nummax;
	Double_t peakx[nummax],peakxerr[nummax];
	Double_t peaky[nummax],peakyerr[nummax];
	Double_t sig[nummax],sigerr[nummax];
	Double_t tau[nummax],tauerr[nummax];
	// 	float *peakx;//if you don't know how many peaks will be found, use this
	// 	float *peaky;//if you don't know how many peaks will be found, use this
	double energylit[4]={1368.626, 2754.007, 2869.50, 4237.96};
	double energyliterr[4]={0.005, 0.011, 0.06, 0.06};
	float peakmin;
	float chtemp;
	double par[nummax][5],parerr[nummax][5];//随峰数不同改动par[peaknum][3]//adjust
	double parchi[nummax],parNDF[nummax];
	int lowerlimit[ID2+1][nummax],upperlimit[ID2+1][nummax];
	ofstream outfile("D:/X/out/outfile/peakcali.dat",ios::out);
	ofstream outfile2("D:/X/out/outfile/peakcalipara.dat",ios::out);
	//FILE *outfile=fopen ("D:/X/peakcali.dat","a");
	int runstart=123,runstop=140;
	for(irawroot=runstart; irawroot<=runstop; irawroot++)
	{
		if(irawroot==124||irawroot==126||irawroot==128||irawroot==130||irawroot==132||irawroot==133||irawroot==135||irawroot==137||irawroot==139)continue;
	sprintf(rawrootname,"%s%04d%s","D:/X/out/Moshe/secondAnalysis/calibrated/protons-",irawroot,"-all.root");
//	sprintf(rawrootname,"%s","D:/out/moshe/run-all_after_cali_25Si.root");//ungated modify
//	sprintf(rawrootname,"%s","X:/aaron/rootfiles/run-all_after_cali_23Al.root");//ungated modify
//	sprintf(rawrootname,"%s","D:/X/out/Si25/100eV-bin-individual_gated_sun_25Si_include_bad_SeGA_do_not_use_total.root");//gated modify
//	sprintf(rawrootname,"%s","D:/X/out/Si25/100eV-bin-individual_gated_sun_23Al_include_bad_SeGA_do_not_use_total.root");//gated modify
	TFile *fin = new TFile(rawrootname);//after this statement, you can use any ROOT command1 for this rootfile
	//TTree *tree = (TTree*)fin->Get("tree");
	//TTree *tree = (TTree*)fin->Get("tree");//此句可有可无
	//	unsigned long nentries=tree->GetEntries();//读事件数
	cout<<rawrootname<<endl;

	// get histogram *******************
//	for(i=ID1;i<=ID2;i++)//which detector
//	{
//		if(i==8) continue;
		sprintf(b_name,"%s","h6");//proton spectrum modify
		sprintf(b_name,"%s","gbCoincidences");//gamma-ray spectrum modify
		histo[irawroot] = (TH1F*)fin->Get(b_name);
		//histo[i]->Rebin(10);
		//histo[i]->GetYaxis()->SetRangeUser(0,4000);//zoom the axis
//	}
	// get histogram *******************
	
//	for(i=0;i<=0;i++)//which SeGA detector modify =0<=15
//	{
		sprintf(hfit_name,"%s%d","Fivepads_Rn",irawroot);
		canvaspeak[irawroot]=new TCanvas(hfit_name,hfit_name,900,600);//建立画布
//		canvaspeak[i]->Divide(2,2);//
//		canvaspeak[irawroot]->cd();
//		if(irawroot==8) continue;
		histo[irawroot]->SetTitle(hfit_name);//图名
		histo[irawroot]->GetXaxis()->SetTitle("Energy");//轴名
		histo[irawroot]->GetYaxis()->SetTitle("Counts");//轴名
		histo[irawroot]->GetXaxis()->CenterTitle();//居中
		histo[irawroot]->GetYaxis()->CenterTitle();//居中
		histo[irawroot]->GetXaxis()->SetLabelFont(132);//坐标字体
		histo[irawroot]->GetYaxis()->SetLabelFont(132);//坐标字体
		histo[irawroot]->GetXaxis()->SetTitleFont(132);//轴名字体
		histo[irawroot]->GetYaxis()->SetTitleFont(132);//轴名字体
		histo[irawroot]->GetXaxis()->SetTitleOffset(1.2);//轴名偏移
		histo[irawroot]->GetYaxis()->SetTitleOffset(1.3);//轴名偏移
		histo[irawroot]->Draw();
		
//		for(ii=0;ii<=0;ii++)//which peak in one SeGA detector modify =0<4
//		{
			float highcounts=0,lowcounts=0,peakcounts=0;
			minbinb=histo[irawroot]->FindBin(230);
			minbin=histo[irawroot]->FindBin(350);
			maxbin=histo[irawroot]->FindBin(470);
			maxbinb=histo[irawroot]->FindBin(590);

			minbinb=histo[irawroot]->FindBin(1330);
			minbin=histo[irawroot]->FindBin(1355);
			maxbin=histo[irawroot]->FindBin(1380);
			maxbinb=histo[irawroot]->FindBin(1405);

// 			minbinb=histo[i]->FindBin(431);
// 			minbin=histo[i]->FindBin(444);
// 			maxbin=histo[i]->FindBin(457);
// 			maxbinb=histo[i]->FindBin(470);

// 			minbinb=histo[i]->FindBin(2867);
// 			minbin=histo[i]->FindBin(2890);
// 			maxbin=histo[i]->FindBin(2913);
// 			maxbinb=histo[i]->FindBin(2936);

// 			minbinb=histo[i]->FindBin(7720);
// 			minbin=histo[i]->FindBin(7770);
// 			maxbin=histo[i]->FindBin(7820);
// 			maxbinb=histo[i]->FindBin(7870);

// 			minbinb=histo[i]->FindBin(918);
// 			minbin=histo[i]->FindBin(940);
// 			maxbin=histo[i]->FindBin(955);
// 			maxbinb=histo[i]->FindBin(963);

// 			minbinb=histo[i]->FindBin(1591);
// 			minbin=histo[i]->FindBin(1605);
// 			maxbin=histo[i]->FindBin(1619);
// 			maxbinb=histo[i]->FindBin(1633);

// 			minbinb=histo[i]->FindBin(2725);
// 			minbin=histo[i]->FindBin(2745);
// 			maxbin=histo[i]->FindBin(2765);
// 			maxbinb=histo[i]->FindBin(2785);

// 			minbinb=histo[i]->FindBin(4190);
// 			minbin=histo[i]->FindBin(4220);
// 			maxbin=histo[i]->FindBin(4250);
// 			maxbinb=histo[i]->FindBin(4280);
 
// 			minbinb=histo[i]->FindBin(1340);
// 			minbin=histo[i]->FindBin(1360);
// 			maxbin=histo[i]->FindBin(1380);
// 			maxbinb=histo[i]->FindBin(1400);

			for(ibin=minbin;ibin<maxbin;ibin++)
			{
				peakcounts+=histo[irawroot]->GetBinContent(ibin);
				//cout<<histo[i]->GetBinContent(ibin)<<endl;
			}
			for(ibin=minbinb;ibin<minbin;ibin++)
			{
				lowcounts+=histo[irawroot]->GetBinContent(ibin);
				//cout<<histo[i]->GetBinContent(ibin)<<endl;
			}
			for(ibin=maxbin;ibin<maxbinb;ibin++)
			{
				highcounts+=histo[irawroot]->GetBinContent(ibin);
				//cout<<histo[i]->GetBinContent(ibin)<<endl;
			}
			histo[irawroot]->GetXaxis()->SetRange(minbinb,maxbinb);
			cout<<"Five_pads"<<irawroot<<"	"<<highcounts<<"	"<<lowcounts<<"	"<<peakcounts<<endl;
			//histo[i]->Draw();
// 			if(ii==0){	gaplow=7; gaphigh=7;}//good dets 12,13's peak2 are troublesome
// 			if(ii==1){	gaplow=7; gaphigh=7;}
// 			if(ii==2){	gaplow=8; gaphigh=9;}
// 			if(ii==3){	gaplow=10; gaphigh=10;}
// 			histo[i]->GetXaxis()->SetRangeUser(peakx[ii]-gaplow,peakx[ii]+gaphigh);//zoom the axis
// 			//cout<<"************"<<peakx[ii]<<"	"<<peaky[ii]<<endl;
// 			minrange=peakx[ii]-gaplow;
// 			maxrange=peakx[ii]+gaphigh;
// 			minbin=histo[i]->FindBin(minrange);
// 			maxbin=histo[i]->FindBin(maxrange);
//		}//for(ii=0;ii<peaknum;ii++)
			//outfile<<"SeGA_"<<i<<"	Constant"<<ii<<"=	"<<par[ii][2]<<"	+/-	"<<parerr[ii][2]<<"	Mean"<<ii<<"=	"<<par[ii][5]<<"	+/-	"<<parerr[ii][5]<<"	Maximum"<<ii<<"=	"<<peakx[ii]<<"	Sigma"<<ii<<"=	"<<par[ii][4]<<"	+/-	"<<parerr[ii][4]<<"	Tau"<<ii<<"=	"<<par[ii][3]<<"	+/-	"<<parerr[ii][3]<<"	A"<<ii<<"=	"<<par[ii][0]<<"	+/-	"<<parerr[ii][0]<<"	B"<<ii<<"=	"<<par[ii][1]<<"	+/-	"<<parerr[ii][1]<<"	Chi2"<<ii<<"=	"<<parchi[ii]<<"	NDF"<<ii<<"=	"<<parNDF[ii]<<"	Area"<<ii<<"=	"<<sqrt(3.1415926536*2)*par[ii][2]*par[ii][4]<<endl;//输出文本查看
			
//	}//for (i=0;i<ID;i++)
	}//for(irawroot=runstart; irawroot<=runstop; irawroot++)
}//peakcali main