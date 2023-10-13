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
void graphpol12fit()// graph-pol1 fit sigma, pol2 fit tau of SeGA 16 channels obtained from 12 25Si,23Al peaks, graphpol12fit_band.C is the same with this code, but also consider error bars.
{
	time_t tim;
	struct tm *at;
	char now[80];
	int irawroot;
	char rawrootname[300];
	double slope[300],slopeerr[300];
	double intercept[300],intercepterr[300];
	double p2[300],p1[300],p0[300];
	double p2err[300],p1err[300],p0err[300];
	const int ID1=0;// no need to modify i=ID1//which detector
	const int ID2=15;// no need to modify i<=ID2//which detector

	int binwidth=1;
	//	int binnum=(upperlimit-lowerlimit)/binwidth;
	int minrange=0,maxrange=0,minbin,maxbin;
	double sigma=1,thresh=0.8,sigmab=20;//adjust
	float gaplow=70.,gaphigh=70.;//fitting range随分辨不同调整//no need to modify
	unsigned long i;
	int jj,ii,ibin,k,ipeak;
	char paraprint[30],b_name[200],histo_name[200],h_name[200],hfit_name[200],hcali_name[200],tauflag[200];
	char command1[100],command2[100];
	//string paraprint[30],h_name[80],histo_name[80],b_name[80],hfit_name[80],hcali_name[80];
	TCanvas *canvaspeak[300];
	TCanvas *canvascali[300];
	TH1F *histo[ID2+1];//TH1F peak search+gauss fit,creat histograms
	const int nummax=12;//search peaknumber can be many, fitnumber can be limited to 3//no need to modify
	TF1 *total[nummax];//creat function
	TF1 *p[nummax], *g[nummax], *b[nummax];
	TGraph *graph1[300], *graph2[300];//TGraph刻度
	int peaknum=nummax;
	Double_t mean[ID2+1][nummax],meanerr[ID2+1][nummax];
	Double_t peaky[nummax],peakyerr[nummax];
	Double_t sig[ID2+1][nummax],sigerr[ID2+1][nummax];
	Double_t tau[ID2+1][nummax],tauerr[ID2+1][nummax];
	// 	float *mean;//if you don't know how many peaks will be found, use this
	// 	float *peaky;//if you don't know how many peaks will be found, use this
	//	double energylit[6]={451.7, 493.3, 944.9, 1460.820, 1612.4, 2614.511};//25Si
	//	double energyliterr[6]={0.5, 0.7, 0.5, 0.005, 0.5, 0.010};//25Si
	double energylit[6]={450.7, 1460.820, 1599, 2614.511, 2908, 7801};//23Al
	double energyliterr[6]={0.15, 0.005, 2, 0.010, 3, 2};//23Al
	float peakmin;
	float chtemp;
	double par[nummax][5],parerr[nummax][5];//随峰数不同改动par[peaknum][3]//adjust
	double parchi[nummax],parNDF[nummax];

//	ifstream infile("C:/Si24/Si22peakcali/Clover.dat",ios::in);//The data that need to be fitted
//	int lowerlimit[ID2+1][nummax],upperlimit[ID2+1][nummax];
	for(i=ID1;i<=ID2;i++)//which detector
	{
		if(i==8||i==11||i==14|i==15) continue;
		sprintf(b_name,"%s%d%s","D:/X/out/Si25/SeGA_",i,"_sigma_freetau.dat");
		ifstream infile(b_name,ios::in);//The data that need to be fitted
		for(ii=0;ii<peaknum;ii++)//which peak in one SeGA detector modify =0<12
		{
			infile>>mean[i][ii]>>meanerr[i][ii]>>sig[i][ii]>>sigerr[i][ii]>>tau[i][ii]>>tauerr[i][ii];
			cout<<"SeGA_"<<i<<'	'<<mean[i][ii]<<'	'<<meanerr[i][ii]<<'	'<<sig[i][ii]<<'	'<<sigerr[i][ii]<<'	'<<tau[i][ii]<<'	'<<tauerr[i][ii]<<endl;
		}
	}
	sprintf(tauflag,"%s","taufree");
	ofstream outfile("D:/X/out/outfile/peakcalipara.dat",ios::out);

	for(i=1;i<=1;i++)//which SeGA detector modify =0<=15
	{
		if(i==8||i==11||i==14|i==15) continue;
		if(i==5)peaknum=11;//modify
		sprintf(hcali_name,"%s%d","AlSi_sigma_SeGA_",i);
		canvascali[i]=new TCanvas(hcali_name,hcali_name,900,600);//建立画布
		canvascali[i]->cd();//进入画布
		//canvascali[i]->SetGrid();//显示网格
		//		graph1[i]=new TGraph(peaknum,mean,sig);//TGraph *gr1=new TGraph(n,x,y);
		//float eenergy[2]={200,280};
		//float epeakch[2]={0,0};
		graph1[i]= new TGraphErrors(peaknum,mean[i],sig[i],meanerr[i],sigerr[i]);//画error bars TGraph(n,x,y,ex,ey);
		graph1[i]->SetTitle(hcali_name);
		graph1[i]->GetXaxis()->SetTitle("Energy");//轴名
		graph1[i]->GetYaxis()->SetTitle("Sigma");//轴名
		graph1[i]->GetXaxis()->CenterTitle();//居中
		graph1[i]->GetYaxis()->CenterTitle();//居中
		graph1[i]->GetXaxis()->SetLabelFont(132);//坐标字体
		graph1[i]->GetYaxis()->SetLabelFont(132);//坐标字体
		graph1[i]->GetXaxis()->SetTitleFont(132);//轴名字体
		graph1[i]->GetYaxis()->SetTitleFont(132);//轴名字体
		//graph1[i]->GetYaxis()->SetLabelSize(0.05);//坐标字号
		//graph1[i]->GetYaxis()->SetTitleSize(0.05);//轴名字号
		graph1[i]->GetXaxis()->SetTitleOffset(1.2);//轴名偏移
		graph1[i]->GetYaxis()->SetTitleOffset(1.3);//轴名偏移
		graph1[i]->SetMarkerStyle(21);
		graph1[i]->SetMarkerColor(1);
		TF1 *pol1=new TF1("pol1","pol1",0,60000);//多项式拟合，调节合适的拟合range
		pol1->SetParNames("p0","p1");
		graph1[i]->Fit("pol1");//pol1 can be used directly without TF1 constructor in CINT
		p0[i]=pol1->GetParameter(0);
		p1[i]=pol1->GetParameter(1);
		p0err[i]=pol1->GetParError(0);
		p1err[i]=pol1->GetParError(1);
		graph1[i]->Draw("AP"); //A-Axis around the graph1,AP is suitable
		TPaveText *textpol1 = new TPaveText(0.13,0.74,0.55,0.89,"brNDC");//left, down, right, up
		textpol1->SetBorderSize(1);
		textpol1->SetFillColor(0);
		textpol1->SetTextAlign(12);//align = 10*HorizontalAlign + VerticalAlign, 12 means水平左对齐、垂直居中对齐
		textpol1->SetTextFont(132);//font = 10 * fontID + precision, 12+2,12 means Symbol; 13+2, 13 means News Times Roman
		sprintf(paraprint,"slope(p1)=%.9f+/-%.9f",p1[i],p1err[i]);
		textpol1->AddText(paraprint);
		sprintf(paraprint,"intercept(p0)=%.5f+/-%.5f",p0[i],p0err[i]);
		textpol1->AddText(paraprint);
		textpol1->Draw();
		//sprintf(hcali_name,"%s.png",hcali_name);
		sprintf(hcali_name,"%s%s%s%s%s","D:/X/out/responseforDoppler/",hcali_name,"_",tauflag,".png");
		canvascali[i]->SaveAs(hcali_name);//存图
		outfile<<"SeGA_"<<i<<"	sigp1=	"<<p1[i]<<"	+/-	"<<p1err[i]<<"	sigp0=	"<<p0[i]<<"	+/-	"<<p0err[i];

		sprintf(hcali_name,"%s%d","AlSi_tau_SeGA_",i);
		canvascali[i]=new TCanvas(hcali_name,hcali_name,900,600);//建立画布
		canvascali[i]->cd();//进入画布
		//canvascali[i]->SetGrid();//显示网格
		//		graph1[i]=new TGraph(peaknum,mean,tau);//TGraph *gr1=new TGraph(n,x,y);
		//float eenergy[2]={200,280};
		//float epeakch[2]={0,0};
		graph2[i]= new TGraphErrors(peaknum,mean[i],tau[i],meanerr[i],tauerr[i]);//画error bars TGraph(n,x,y,ex,ey);
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
		graph2[i]->SetMarkerStyle(21);
		graph2[i]->SetMarkerColor(1);
		TF1 *pol2=new TF1("pol2","pol2",0,60000);//多项式拟合，调节合适的拟合range
		pol2->SetParNames("p0","p1","p2");
		graph2[i]->Fit("pol2");//pol1 can be used directly without TF1 constructor in CINT
		p0[i]=pol2->GetParameter(0);
		p1[i]=pol2->GetParameter(1);
		p2[i]=pol2->GetParameter(2);
		p0err[i]=pol2->GetParError(0);
		p1err[i]=pol2->GetParError(1);
		p2err[i]=pol2->GetParError(2);
		graph2[i]->Draw("AP"); //A-Axis around the graph1,AP is suitable
		TPaveText *textpol2 = new TPaveText(0.13,0.74,0.55,0.89,"brNDC");//left, down, right, up
		textpol2->SetBorderSize(1);
		textpol2->SetFillColor(0);
		textpol2->SetTextAlign(12);//align = 10*HorizontalAlign + VerticalAlign, 12 means水平左对齐、垂直居中对齐
		textpol2->SetTextFont(132);//font = 10 * fontID + precision, 12+2,12 means Symbol; 13+2, 13 means News Times Roman
		sprintf(paraprint,"p2=%.9f+/-%.9f",p2[i],p2err[i]);
		textpol2->AddText(paraprint);
		sprintf(paraprint,"p1=%.5f+/-%.5f",p1[i],p1err[i]);
		textpol2->AddText(paraprint);
		sprintf(paraprint,"p0=%.5f+/-%.5f",p0[i],p0err[i]);
		textpol2->AddText(paraprint);
		textpol2->Draw();
		//sprintf(hcali_name,"%s.png",hcali_name);
		sprintf(hcali_name,"%s%s%s%s%s","D:/X/out/responseforDoppler/",hcali_name,"_",tauflag,".png");
		canvascali[i]->SaveAs(hcali_name);//存图
		outfile<<"	SeGA_"<<i<<"	taup2=	"<<p2[i]<<"	+/-	"<<p2err[i]<<"	taup1=	"<<p1[i]<<"	+/-	"<<p1err[i]<<"	taup0=	"<<p0[i]<<"	+/-	"<<p0err[i]<<endl;

	}//for (i=0;i<ID;i++)
	//}//for(irawroot=runstart; irawroot<=runstop; irawroot++)
}//peakcali main