#include <iostream>
#include <fstream>
#include <math.h>
#include <time.h>
#include <map>
#include <TROOT.h>
#include <TSystem.h>
#include <TStyle.h>
#include <TCutG.h>
#include <TChain.h>
#include <TMinuit.h>
#include "TChain.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TEllipse.h"
#include "TString.h"
#include "TGraph.h"
#include "TSpectrum.h"
#include "TVirtualFitter.h"
#include "TRandom.h"
#include <TRandom3.h>
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
#include "TPaletteAxis.h"
void rantest()// test random
{
	const float pi=3.141592653589793;
	time_t start,tim;
	struct tm *at;
	char now[80];
	float speed;
// 	TH1F *hxy=new TH1F("ran_Rndm","ran_Rndm",400, -2,2);
// 	TH1F *hxyz=new TH1F("costheta","costheta",400, -2,2);
// 	TH1F *h4=new TH1F("acos(costheta)","acos(costheta)",500, -1,4);
// 	TH1F *h5=new TH1F("theta","theta",500, -1,4);
// 	TH1F *h6=new TH1F("cos(theta)","cos(theta)",400, -2,2);
	TH1F *hx = new TH1F("hx","hx", 300,-1.5,1.5);
	TH1F *hy = new TH1F("hy","hy", 300,-1.5,1.5);
	TH1F *hz = new TH1F("hz","hz", 300,-1.5,1.5);
	TH2F *hxy = new TH2F("hxy","hxy", 300,-1.5,1.5,300,-1.5,1.5);
	TH3F *hxyz = new TH3F("hxyz","hxyz", 300,-1.5,1.5,300,-1.5,1.5,300,-1.5,1.5);
// 	TH1F *h4=new TH1F("Fill_Histogram_Randomfrom_Function","Fill_Histogram_Randomfrom_Function",1000, 0,1000);
// 	TH1F *h5=new TH1F("Fill_Histogram_Randomfrom_Histogram","Fill_Histogram_Randomfrom_Histogram",1000, 0,1000);
// 	TH1F *h6=new TH1F("GetRandomfrom_Histogram","GetRandomfrom_Histogram",1000, 0,1000);
// 	TH1F *h7=new TH1F("GetRandomfrom_Function","GetRandomfrom_Function",1000, 0,1000);
// 	TH1F *h8=new TH1F("GetRandomfrom_Function_Range","GetRandomfrom_Function_Range",1000, 0,1000);
//	TRandom3* ran = new TRandom3(61520154+49092137+200000);
// 	TF1 *SiDEC=new TF1("SiDEC","[0]*exp(x/(-[1]/0.693147))+[2]", 0,1000);//自定义拟合函数
// 	SiDEC->SetNpx(1000);
// 	SiDEC->SetParNames("A","T","B");
// 	SiDEC->SetParameters(500,100,0);//自定义的拟合函数必须赋初值
// 	SiDEC->SetParLimits(2,0,20);
// 	TCanvas *c1=new TCanvas("c1","c1");
// 	c1->cd();
// 	SiDEC->Draw();
// 	h4->FillRandom("SiDEC",1000000);
// 	TCanvas *c2=new TCanvas("c2","c2");
// 	c2->cd();
// 	h4->Draw();
	//TH1F *h5=new TH1F("time","time",200,0,20*T_Lifetime);

//	h5->FillRandom(h4,1000000);
//	TCanvas *c3=new TCanvas("c3","c3");
//	c3->cd();
//	h5->Draw();
	float costheta,theta;
	float phi,r;
	float x,y,z;
	long TotalCounts=100000;
	bool CosVsX = false;
	gRandom->SetSeed(0);
	TRandom3* ran = new TRandom3(61520154+49092137+200000);
	if(CosVsX)
	{
		for(long i=0;i<TotalCounts;i++)//If type "for loop" in a root window, it just execute once, instead of 1000 times. It's weird.
		{
			theta=gRandom->Uniform(0,pi);//Uniform(x_low,x_high)
			phi=gRandom->Uniform(0,2*pi);
			r=1;
			x=r*sin(theta)*cos(phi);
			y=r*sin(theta)*sin(phi);
			z=r*cos(theta);
			hx->Fill(x);hy->Fill(y);hz->Fill(z);
			hxy->Fill(x,y);//
			hxyz->Fill(x,y,z);//
			if(i==0)
			{
				time(&start);
				at=localtime(&start);
				strftime(now,79,"%Y-%m-%d %H:%M:%S",at);
				cout<<now<<" started."<<endl;
			}
			if(i%10000==0&&i!=0)
			{
				time(&tim);
				at=localtime(&tim);
				strftime(now,79,"%Y-%m-%d %H:%M:%S",at);
				cout<<now;
				speed=(float)i/(float)difftime(tim,start);
				printf(" %.1f%s%.0f%s%.0f%s\n",100*(float)i/(float)TotalCounts,"%. Speed is ",speed," e/s. Still need ",float(TotalCounts-i)/speed," s.");
			}
		}
	}
	else
		for(long i=0;i<TotalCounts;i++)//If type "for loop" in a root window, it just execute once, instead of 1000 times. It's weird.
		{
			costheta=ran->Uniform(-1,1);
			phi=ran->Uniform(0,2*pi);
			//costheta=gRandom->Uniform(-1,1);
			//phi=gRandom->Uniform(0,2*pi);
			r=1;
			x=r*sqrt(1-costheta*costheta)*cos(phi);
			y=r*sqrt(1-costheta*costheta)*sin(phi);
			z=r*costheta;
			hx->Fill(x);hy->Fill(y);hz->Fill(z);
			hxy->Fill(x,y);//
			hxyz->Fill(x,y,z);//
			if(i==0)
			{
				time(&start);
				at=localtime(&start);
				strftime(now,79,"%Y-%m-%d %H:%M:%S",at);
				cout<<now<<" started."<<endl;
			}
			if(i%10000==0&&i!=0)
			{
				time(&tim);
				at=localtime(&tim);
				strftime(now,79,"%Y-%m-%d %H:%M:%S",at);
				cout<<now;
				speed=(float)i/(float)difftime(tim,start);
				printf(" %.1f%s%.0f%s%.0f%s\n",100*(float)i/(float)TotalCounts,"%. Speed is ",speed," e/s. Still need ",float(TotalCounts-i)/speed," s.");
			}
		}

// 		rand2 = ran->Rndm();
// 		hxy->Fill(rand2);
// 		costheta=gRandom->Uniform(-1,1);
// 		hxyz->Fill(costheta);
// 		h4->Fill(acos(costheta));
// 		theta=gRandom->Uniform(0,3.141592654);
// 		h5->Fill(theta);
// 		h6->Fill(cos(theta));
// 		Double_t n=h5->GetRandom();//Return a random number distributed according the histogram bin contents.
// 		h6->Fill(n);
// 		n=SiDEC->GetRandom();//Return a random number following this function shape.
// 		h7->Fill(n);
// 		n=SiDEC->GetRandom(200,500);//Return a random number following this function shape in [xmin,xmax].
// 		h8->Fill(n);//Random numbers distributed according to a user defined function in a limited interval,
		//if(i%(counts/100)==0){cout<<j<<"%"<<" done"<<endl;j++;}
	TCanvas *cx=new TCanvas("cx","cx");
	hx->Draw();
	TCanvas *cy=new TCanvas("cy","cy");
	hy->Draw();
	TCanvas *cz=new TCanvas("cz","cz");
	hz->Draw();
 	TCanvas *cxy=new TCanvas("cxy","cxy",1000,1000);
 	hxy->Draw();
 	TCanvas *cxyz=new TCanvas("cxyz","cxyz",1000,1000);
 	hxyz->Draw();
// 	TCanvas *c6=new TCanvas("c6","c6");
// 	h8->Draw();
}