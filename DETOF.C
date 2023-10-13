#include <iostream>
#include <fstream>
#include <math.h>
#include <map>
#include <TROOT.h>
#include "TChain.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TEllipse.h"
#include "TString.h"
#include "TGraph.h"
#include "TSpectrum.h"
#include "TVirtualFitter.h"
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
#include <time.h>
#include<TChain.h>
using namespace std;

void DETOF()//draw DETOF plot
{
	time_t tim;
	struct tm *at;
	char now1[80];
	long nentries[652];
	long totalentries;
	memset(nentries,0,sizeof(nentries));
	int icalroot,ilarge;
	char calrootname[80];
	char anarootname[80];
	char resultname[80];
	int runstart,runstop;
	cout<<"input runstart:";
	cin>>runstart;
	cout<<"input runstop:";
	cin>>runstop;

	TChain* T777Chain = new TChain("T777");//root文件中的tree名
	for(icalroot=runstart; icalroot<=runstop; icalroot++)
	{
		if((icalroot==182)||(icalroot==285)||(icalroot>=252&&icalroot<=319))
		{
			nentries[icalroot]=T777Chain->GetEntries();//if wrong RN was input, nentries[wrong RN]=last nentries
			continue;
		}
		//||(icalroot>=526&&icalroot<=528)
		sprintf(calrootname,"%s%04d%s","X:/T777/S27calabcd",icalroot,".root");
		T777Chain->Add(calrootname);
		cout<<calrootname;
		nentries[icalroot]=T777Chain->GetEntries();//unsigned long是最大范围的整数，相当于ULong64_t
		cout<<" Total nentries["<<icalroot<<"]="<<nentries[icalroot]<<endl;
	}
	for(ilarge=runstop+1;ilarge<652;ilarge++)
	{
		nentries[ilarge]=500000000;//in case of no appropriate ellipse loop to get in
	}
	totalentries=nentries[icalroot-1];
	long nentriesmax;
	cout<<"input the max number of entries (not longer than 2147483647): "<<endl;
	//cin>>nentriesmax;
	nentriesmax=500000000;
	Float_t         DE1max;
	Float_t         DE2[4];
	Float_t         T1;
	Float_t         T2;
	Float_t         TDE1[4];
	Float_t         TDE2[4];
	Float_t         TOF;
	T777Chain->SetBranchAddress("DE1", &DE1max);
	T777Chain->SetBranchAddress("DE2", DE2);
	T777Chain->SetBranchAddress("TDE1", TDE1);
	T777Chain->SetBranchAddress("TDE2", TDE2);
	T777Chain->SetBranchAddress("T1", &T1);
	T777Chain->SetBranchAddress("T2", &T2);
	T777Chain->SetBranchAddress("TOF", &TOF);
	//用SetBranchAddress函数将tree的Branch TOFC与重定义好的变量地址&TOFC联系起来

	sprintf(anarootname,"%s","X:/T777/S27DE2TOF.root");
	TFile *fout = new TFile(anarootname,"RECREATE");//输出文件//It's better to define histograms and then define fout, in case of draw bugs.
	TCanvas *cDE1TOF=new TCanvas("cDE1TOF","cDE1TOF",800,600);
	//TH2F *hDE1TOF = new TH2F("hDE1TOF","hDE1TOF",400,180,230,400,100,300);//name a histogram	 (TOF,DE)
	TH2F *hDE2TOF = new TH2F("hDE2TOF","hDE2TOF",400,180,230,600,100,400);//name a histogram	 (TOF,DE)
	for(int i=0;i<totalentries;i++)
	{
		T777Chain->GetEntry(i);//获得输入root文件的第i个entry相应的Branch变量数据，then the redefined variables could be used.
		if(DE1max>100&&DE1max<300&&TOF<230&&TOF>180
			&&(DE2[0]>300||DE2[1]>300||DE2[2]>300||DE2[3]>300)
			&&((TDE1[0]-T1>-26000&&TDE1[0]-T1<-22000)||(TDE1[1]-T1>-26000&&TDE1[1]-T1<-22000)||(TDE1[2]-T1>-26000&&TDE1[2]-T1<-22000)||(TDE1[3]-T1>-26000&&TDE1[3]-T1<-22000))
			&&((TDE2[0]-T1>-26000&&TDE2[0]-T1<-22000)||(TDE2[1]-T1>-26000&&TDE2[1]-T1<-22000)||(TDE2[2]-T1>-26000&&TDE2[2]-T1<-22000)||(TDE2[3]-T1>-26000&&TDE2[3]-T1<-22000))
			)
		{
// 				//&&(pow((TOF-x1)/a1,2)+pow((DE1-y1)/b1,2)>1)
// 				//&&(pow((TOF-x2)/a2,2)+pow((DE1-y2)/b2,2)>1)
// 				)
				//hDE1TOF->Fill(TOF,DE1max);
				hDE2TOF->Fill(TOF,DE2[0]*0.12589+31.8826);
// 			else if(
// 				(pow((TOF-x1)/a1,2)+pow((DE1-y1)/b1,2)<=1&&gRandom->Uniform(0,1)>0.99)
// 				||(pow((TOF-x1)/a1/2,2)+pow((DE1-y1)/b1/2,2)<=1&&gRandom->Uniform(0,1)>0.999)
// 				||(pow((TOF-x2)/a2,2)+pow((DE1-y2)/b2,2)<=1&&gRandom->Uniform(0,1)>0.99)
// 				||(pow((TOF-x2)/a2/2,2)+pow((DE1-y2)/b2/2,2)<=1&&gRandom->Uniform(0,1)>0.999)
// 				)
// 				hDE1TOF->Fill(TOF,DE1);
// 			else if(gRandom->Uniform(0,1)>0.9)
// 				hDE1TOF->Fill(TOF,DE1);
		}
		if(i%1000000==0)
		{
			time(&tim);
			at=localtime(&tim);
			strftime(now1,79,"%Y-%m-%d %H:%M:%S",at);
			cout<<now1;
			printf(" complete %.1f%s\n",100*(float)i/(float)totalentries,"%");
			//cout<<" x1Mg="<<x1Mg<<", y1Mg="<<y1Mg<<", aMg="<<aMg<<", bMg="<<bMg<<endl;
		}
		if(i>=nentriesmax) break;
	}
	cDE1TOF->cd();
	gPad->SetLogz();
	hDE2TOF->Draw("colz");
	fout->Write();//等效于把所有的tree和一维谱都写入文件
	fout->Close();
}