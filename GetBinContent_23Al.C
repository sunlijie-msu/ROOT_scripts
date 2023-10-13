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
#include<TChain.h>
#include "TChain.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TEllipse.h"
#include "TString.h"
#include "TGraph.h"
#include "TSpectrum.h"
#include "TVirtualFitter.h"
#include "TRandom.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TMath.h"
#include "TPaveText.h"
#include "TPad.h"
#include "string.h"
#include "stdlib.h"
#include "stdio.h"
#include "TLegend.h"
using namespace std;
void GetBinContent_23Al()//get lots of peak integrals from lots of histograms
{
	int ii,jj,i,j,numberofbin10,numberofbin3,numberofbin4,numberofbin6,numberofbin8;
	char calrootname[80];
	char hname[80],gname[80];
	char txtname[80];
	int runstart,runstop;
// 	cout<<"input runstart:";
// 	cin>>runstart;
// 	cout<<"input runstop:";
// 	cin>>runstop;
//	for(icalroot=runstart; icalroot<=runstop; icalroot++)
	ofstream outfile("D:/X/out/Si25/outfile/peakcali.dat",ios::out);
	// 	sprintf(calrootname,"%s","D:/X/out/Moshe/source_at_0_0_20_upstream_5M.root");//up
	//	sprintf(calrootname,"%s","D:/X/out/Moshe/source_at_0_0_20_upstream_5M_noFET.root");//up
	//	sprintf(calrootname,"%s","D:/X/out/Moshe/source_at_0_0_0_center_5M.root");//center sega
	//	sprintf(calrootname,"%s","D:/X/out/Moshe/source_at_0_0_-20_downstream_5M.root");//down
	sprintf(calrootname,"%s","D:/X/out/Moshe/23Al_5M.root");//23Al
	//sprintf(calrootname,"%s","D:/X/out/Moshe/25Si_10M.root");//25Si
	TFile *fin = new TFile(calrootname);//after this statement, you can use any ROOT command for this rootfile
	TH1F *histo[17];
// 	TTree *tree = (TTree*)fin->Get("tree");
// 	long nentries=tree->GetEntries();//读事件数
	//Firstly, redefine the variables to hold the read values.
	//Float_t SCA[3];
	//tree->SetBranchAddress("SCA", SCA);
// 	ULong64_t sdc[32];
// 	tree->SetBranchAddress("sdc", sdc);//cannot be omitted
	//用SetBranchAddress函数将tree的Branch TOFC与重定义好的变量地址&TOFC联系起来
//	tree->GetEntry(nentries-1);
	//获得输入root文件的第某个entry相应的Branch变量数据
//	int Eg[12]={122,245,344,411,444,779,867,964,1112,1213,1299,1408};//152Eu
	int Eg[12]={451,1599,2904,7801,7801,7801,7801,7801,7801,7801,7801,7801};//23Al
//	int Eg[17]={452,493,845,884,945,1338,1369,1612,1789,2222,2673,2754,2870,4238,6288,6955,7900};//25Si

	i=16;
	numberofbin3=3; numberofbin4=4; numberofbin6=6; numberofbin8=8; numberofbin10=10;
	sprintf(hname,"%s","hSega");
	int Integral[17]={0};
	int bkgl[17]={0}, bkgh[17]={0};
	histo[i] = (TH1F*)fin->Get(hname);
	outfile<<hname<<endl;
	for(ii=0;ii<=11;ii++)//which peak
	{
		for (jj=-numberofbin10;jj<=numberofbin10;jj++)
		{
			Integral[ii]+=histo[i]->GetBinContent(Eg[ii]+jj);//-8+8
			//cout<<Eg[ii]+jj<<'	';
			bkgl[ii]+=histo[i]->GetBinContent(Eg[ii]-numberofbin10*2+jj-1);//-25-9
			//cout<<Eg[ii]-numberofbin10*2+jj-1<<'	';
			bkgh[ii]+=histo[i]->GetBinContent(Eg[ii]+numberofbin10*2+jj+1);//+9+25
			//cout<<Eg[ii]+numberofbin10*2+jj+1<<endl;
		}
		outfile<<"Eg=	"<<Eg[ii]<<"	Integral=	"<<Integral[ii]<<"	bkgl=	"<<bkgl[ii]<<"	bkgh=	"<<bkgh[ii]<<endl;
	}
	outfile<<endl;
	for(i=0;i<16;i++)//which SeGA detector
	{
		sprintf(hname,"%s%d","hEnergyDepositSega_",i);
		int Integral[17]={0};
		int bkgl[17]={0}, bkgh[17]={0};
		histo[i] = (TH1F*)fin->Get(hname);
		outfile<<hname<<endl;
		for(ii=0;ii<=11;ii++)//which peak
		{
			for (jj=-numberofbin10;jj<=numberofbin10;jj++)
			{
				Integral[ii]+=histo[i]->GetBinContent(Eg[ii]+jj);//-3+3
				//cout<<Eg[ii]+jj<<'	';
				bkgl[ii]+=histo[i]->GetBinContent(Eg[ii]-numberofbin10*2+jj-1);//-10-4
				//cout<<Eg[ii]-numberofbin10*2+jj-1<<'	';
				bkgh[ii]+=histo[i]->GetBinContent(Eg[ii]+numberofbin10*2+jj+1);//+4+10
				//cout<<Eg[ii]+numberofbin10*2+jj+1<<endl;
			}
			outfile<<"Eg=	"<<Eg[ii]<<"	Integral=	"<<Integral[ii]<<"	bkgl=	"<<bkgl[ii]<<"	bkgh=	"<<bkgh[ii]<<endl;
		}
// 		for(ii=9;ii<=9;ii++)//which peak
// 		{
// 			for (jj=-numberofbin4;jj<=numberofbin4;jj++)
// 			{
// 				Integral[ii]+=histo[i]->GetBinContent(Eg[ii]+jj);//-4+4
// 				//cout<<Eg[ii]+jj<<'	';
// 				bkgl[ii]+=histo[i]->GetBinContent(Eg[ii]-numberofbin4*2+jj-1);//-13-5
// 				//cout<<Eg[ii]-numberofbin10*2+jj-1<<'	';
// 				bkgh[ii]+=histo[i]->GetBinContent(Eg[ii]+numberofbin4*2+jj+1);//+5+13
// 				//cout<<Eg[ii]+numberofbin10*2+jj+1<<endl;
// 			}
// 			outfile<<"Eg=	"<<Eg[ii]<<"	Integral=	"<<Integral[ii]<<"	bkgl=	"<<bkgl[ii]<<"	bkgh=	"<<bkgh[ii]<<endl;
// 		}
// 		for(ii=10;ii<=12;ii++)//which peak
// 		{
// 			for (jj=-numberofbin6;jj<=numberofbin6;jj++)
// 			{
// 				Integral[ii]+=histo[i]->GetBinContent(Eg[ii]+jj);//-6+6
// 				//cout<<Eg[ii]+jj<<'	';
// 				bkgl[ii]+=histo[i]->GetBinContent(Eg[ii]-numberofbin6*2+jj-1);//-19-7
// 				//cout<<Eg[ii]-numberofbin10*2+jj-1<<'	';
// 				bkgh[ii]+=histo[i]->GetBinContent(Eg[ii]+numberofbin6*2+jj+1);//+7+19
// 				//cout<<Eg[ii]+numberofbin10*2+jj+1<<endl;
// 			}
// 			outfile<<"Eg=	"<<Eg[ii]<<"	Integral=	"<<Integral[ii]<<"	bkgl=	"<<bkgl[ii]<<"	bkgh=	"<<bkgh[ii]<<endl;
// 		}
// 		for(ii=13;ii<=16;ii++)//which peak
// 		{
// 			for (jj=-numberofbin10;jj<=numberofbin10;jj++)
// 			{
// 				Integral[ii]+=histo[i]->GetBinContent(Eg[ii]+jj);//-10+10
// 				//cout<<Eg[ii]+jj<<'	';
// 				bkgl[ii]+=histo[i]->GetBinContent(Eg[ii]-numberofbin10*2+jj-1);//-31-11
// 				//cout<<Eg[ii]-numberofbin10*2+jj-1<<'	';
// 				bkgh[ii]+=histo[i]->GetBinContent(Eg[ii]+numberofbin10*2+jj+1);//+11+31
// 				//cout<<Eg[ii]+numberofbin10*2+jj+1<<endl;
// 			}
// 			outfile<<"Eg=	"<<Eg[ii]<<"	Integral=	"<<Integral[ii]<<"	bkgl=	"<<bkgl[ii]<<"	bkgh=	"<<bkgh[ii]<<endl;
// 		}
		outfile<<endl;
	}

	//outfile<<"Rn"<<icalroot<<'	'<<nentries<<'	'<<sdc[0]<<'	'<<sdc[1]<<'	'<<sdc[2]<<endl;//time stamp, trigger, event
	//tree->Scan("adc[1][7]");
	fin->Close();
}