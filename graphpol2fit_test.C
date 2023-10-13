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
void graphpol2fit()
{
	time_t tim;
	struct tm *at;
	char now[80];
	int irawroot;
	char rawrootname[80];
	double p2[100],p2err[100];
	double p1[100],p1err[100];
	double p0[100],p0err[100];
	const int ID1=0;//adjust no use
	const int ID2=1;//adjust no use
	const int binmin=300,binmax=6000;//adjust//search for peaks in this range no use
	const int binwidth=10;// no use
	int binnum=(binmax-binmin)/binwidth;// no use
	double sigma=7,thresh=0.3,sigmab=20;//adjust
	float gaplow=100.,gaphigh=100.;//fitting range随分辨不同调整//adjust
	unsigned long i;
	int jj,ii;
	char paraprint[30],b_name[50],histo_name[50],h_name[50],hfit_name[50],hcali_name[50];
	char command[60];
	TCanvas *canvascali[300];
	TGraph *graph[300];//TGraph刻度
	int peaknum=7;//adjust

	float energy[18],channel[18];
	ifstream infile("C:/Si24/Si22peakcali/Clover.dat",ios::in);//The data that need to be fitted
	for(ii=0;ii<peaknum;ii++)
	{
		infile>>channel[ii]>>energy[ii];
		cout<<' '<<channel[ii]<<' '<<energy[ii]<<endl;
	}
	ofstream outfile("C:/Si24/Si22peakcali/Clovercali.dat",ios::out);//modify if change path
	sprintf(b_name,"%s","Clover");
	for(i=ID1;i<ID2;i++)
	{
		sprintf(hcali_name,"cali_%s%s%d%s",b_name,"[",i,"]");
		canvascali[i]=new TCanvas(hcali_name,hcali_name,900,600);//建立画布
		canvascali[i]->cd();//进入画布
		//canvascali[i]->SetGrid();//显示网格
		graph[i]=new TGraph(peaknum,channel,energy);//TGraph *gr1=new TGraph(n,x,y);
		//float eenergy[2]={200,280};
		//float epeakch[2]={0,0};
		//graph[i]= new TGraphErrors(peaknum,peakx,energy,epeakch,eenergy);//画error bars
		graph[i]->SetTitle(hcali_name);
		graph[i]->GetXaxis()->SetTitle("Channel");//轴名
		graph[i]->GetYaxis()->SetTitle("Energy");//轴名
		graph[i]->GetXaxis()->CenterTitle();//居中
		graph[i]->GetYaxis()->CenterTitle();//居中
		graph[i]->GetXaxis()->SetLabelFont(132);//坐标字体
		graph[i]->GetYaxis()->SetLabelFont(132);//坐标字体
		graph[i]->GetXaxis()->SetTitleFont(132);//轴名字体
		graph[i]->GetYaxis()->SetTitleFont(132);//轴名字体
		//graph[i]->GetYaxis()->SetLabelSize(0.05);//坐标字号
		//graph[i]->GetYaxis()->SetTitleSize(0.05);//轴名字号
		graph[i]->GetXaxis()->SetTitleOffset(1.2);//轴名偏移
		graph[i]->GetYaxis()->SetTitleOffset(1.3);//轴名偏移
		graph[i]->SetMarkerStyle(21);
		graph[i]->SetMarkerColor(1);
		TF1 *pol2=new TF1("pol2","pol2",20,5000);//多项式拟合，调节合适的拟合range
		pol2->SetParNames("p0","p1","p2");
		graph[i]->Fit("pol2");//pol2 can be used directly without TF1 constructor in CINT
		p0[i]=pol2->GetParameter(0);
		p1[i]=pol2->GetParameter(1);
		p2[i]=pol2->GetParameter(2);
		p0err[i]=pol2->GetParError(0);//Obtaining the error of the 1st parameter
		p1err[i]=pol2->GetParError(1);//Obtaining the error of the 2nd parameter
		p2err[i]=pol2->GetParError(2);//Obtaining the error of the 3rd parameter
		cout<<pol2->GetChisquare()/pol2->GetNDF();
		graph[i]->Draw("AP"); //A-Axis around the graph,AP is suitable
		TPaveText *textpol1 = new TPaveText(0.13,0.68,0.45,0.89,"brNDC");//left, down, right, up
		textpol1->SetBorderSize(1);
		textpol1->SetFillColor(0);
		textpol1->SetTextAlign(12);//align = 10*HorizontalAlign + VerticalAlign, 12 means水平左对齐、垂直居中对齐
		textpol1->SetTextFont(132);//font = 10 * fontID + precision, 12+2,12 means Symbol; 13+2, 13 means News Times Roman
		sprintf(paraprint,"p2=%.5f%s%.5f",p2[i],"+/-",p2err[i]);
		textpol1->AddText(paraprint);
		sprintf(paraprint,"p1=%.3f%s%.3f",p1[i],"+/-",p1err[i]);
		textpol1->AddText(paraprint);
		sprintf(paraprint,"p0=%.3f%s%.3f",p0[i],"+/-",p0err[i]);
		textpol1->AddText(paraprint);
		sprintf(paraprint,"y=%.5f*x^2+%.3f*x+%.3f",p2[i],p1[i],p0[i]);
		textpol1->AddText(paraprint);
		textpol1->Draw();
		//sprintf(hcali_name,"%s.png",hcali_name);
		sprintf(hcali_name,"%s%s%s","C:/Si24/Si22peakcali/",hcali_name,".png");//modify if change path
		canvascali[i]->SaveAs(hcali_name);//存图
		outfile<<b_name<<i<<" p2=	"<<p2[i]<<'	'<<p2err[i]<<endl;//输出文本查看
		outfile<<b_name<<i<<" p1=	"<<p1[i]<<'	'<<p1err[i]<<endl;//输出文本查看
		outfile<<b_name<<i<<" p0=	"<<p0[i]<<'	'<<p0err[i]<<endl;//输出文本查看
		outfile<<b_name<<i<<" y="<<p2[i]<<"*x^2+"<<p1[i]<<"*x+"<<p0[i]<<endl;//输出文本查看
	}//for (i=0;i<ID;i++)
	//}//for(irawroot=runstart; irawroot<=runstop; irawroot++)
}//peakcali main