#include <iostream>
#include <iomanip.h>
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
using namespace std;
void Ran()
{
	unsigned long i,countin;
	int jj,ii,ie;
	time_t tim;
	struct tm *at;
	char now[80];
	float PrangeinSi[30],PdrangeinSi[30];
	int irawroot;
	char h_name[50];
	char rawrootname[80];
	const float pi=3.14159265;
	float x1D1=-24.7,y1D1=-21.7,z1D1=-10.149,x2D1=24.7,y2D1=27.7,z2D1=-10,x0=0,y0=0,z0=0;//DSSD1 implantation，x junction; y ohmic
	float x1D2=-24.7,y1D2=-24.7,z1D2=0,x2D2=24.7,y2D2=24.7,z2D2=0.066;//DSSD2 implantation，x junction; y ohmic
	float x1Q1=-24.7,y1Q1=-24.7,z1Q1=7,x2Q1=24.7,y2Q1=24.7,z2Q1=7.3;//QSD1 geometry
	float x1Q2=-24.7,y1Q2=-24.7,z1Q2=12,x2Q2=24.7,y2Q2=24.7,z2Q2=13.5;//QSD2 geometry
	float x1QU=-24.7,y1QU=49,z1QU=-63,x2QU=24.7,y2QU=50.5,z2QU=-13;//QSDU geometry
	float x1QD=-24.7,y1QD=-37.5,z1QD=-63,x2QD=24.7,y2QD=-36,z2QD=-13;//QSDD geometry
	float x1QR=42,y1QR=-23.7,z1QR=-50,x2QR=43.5,y2QR=25.7,z2QR=0;//QSDR geometry
	float x1QL=-43.5,y1QL=-23.7,z1QL=-50,x2QL=-42,y2QL=25.7,z2QL=0;//QSDL geometry
	float x,y,z;
	float costheta,phi,r;
	double rin,Ein,Eloss;
	sprintf(rawrootname,"V:/RIBLL2015/data24/Si24calabcd_0636_0650.root");
	TFile *fin = new TFile(rawrootname);
	TTree *T999 = (TTree*)fin->Get("T999");
	unsigned long nentries=T999->GetEntries();//读事件数
	cout<<rawrootname;
	cout<<"  Entries="<<nentries<<endl;
	TH1F *px=(TH1F*)fin->Get("hSiimpl60px");// = new TH1F("px","px",16,-23.3,23.3);//create a histogram
	TH1F *py=(TH1F*)fin->Get("hSiimpl60py");// = new TH1F("py","py",16,-23.3,23.3);//create a histogram
	TH1D *pz=(TH1D*)fin->Get("hSiimpl60pz");
// 	for(ii=0;ii<16;ii++)
// 	{
// 		sprintf(h_name,"%s%d","hSiimpl60pz",ii);
// 		pz[ii] = (TH1D*)fin->Get(h_name);
// 	}
	TCanvas *c3=new TCanvas("c3","c3",640,640);
	TH3F *h3 = new TH3F("h3","h3", 160,-23.3*2,23.3*2.2,160,-23.3*2,23.3*2,160,-23.3*3,23.3*1);
	TH2F *h2D1 = new TH2F("h2D1","h2D1", 160,-23.3*2,23.3*2,160,-23.3*2,23.3*2);
	TH1F *h1D1 = new TH1F("h1D1","h1D1", 160,-23.3*2,23.3*2);
	TH2F *h2D2 = new TH2F("h2D2","h2D2", 160,-23.3*2,23.3*2,160,-23.3*2,23.3*2);
	TH1F *h1D2 = new TH1F("h1D2","h1D2", 160,-23.3*2,23.3*2);
	TH2F *h2Q1 = new TH2F("h2Q1","h2Q1", 160,-23.3*2,23.3*2,160,-23.3*2,23.3*2);
	TH1F *h1Q1 = new TH1F("h1Q1","h1Q1", 160,-23.3*2,23.3*2);
	TH2F *h2Q2 = new TH2F("h2Q2","h2Q2", 160,-23.3*2,23.3*2,160,-23.3*2,23.3*2);
	TH1F *h1Q2 = new TH1F("h1Q2","h1Q2", 160,-23.3*2,23.3*2);
	TH2F *h2QU = new TH2F("h2QU","h2QU", 160,-23.3*2,23.3*2,160,-23.3*3,23.3*1);
	TH1F *h1QU = new TH1F("h1QU","h1QU", 160,-23.3*3,23.3*0);
	TH2F *h2QD = new TH2F("h2QD","h2QD", 160,-23.3*2,23.3*2,160,-23.3*3,23.3*1);
	TH1F *h1QD = new TH1F("h1QD","h1QD", 160,-23.3*3,23.3*0);
	TH2F *h2QR = new TH2F("h2QR","h2QR", 160,-23.3*2,23.3*2,160,-23.3*3,23.3*1);
	TH1F *h1QR = new TH1F("h1QR","h1QR", 160,-23.3*3,23.3*0.5);
	TH2F *h2QL = new TH2F("h2QL","h2QL", 160,-23.3*2,23.3*2,160,-23.3*3,23.3*1);
	TH1F *h1QL = new TH1F("h1QL","h1QL", 160,-23.3*3,23.3*0.5);
	TH1F *hEin = new TH1F("hEin","hEin", 1000,-1,4000);
	TH1F *hEloss = new TH1F("hEloss","hEloss", 1000,-1,4000);
	TH1D *hrange = new TH1D("hrange","hrange", 500,-1,20);
	TH2F *h2 = new TH2F("h2","h2", 16,-23.3,23.3,16,-23.3,23.3);
	TH1F *h1 = new TH1F("h1","h1", 66,0,0.066);
// 	h3->FillRandom("gaus",5000);
	ifstream infile("C:/Si24/Si22peakcali/EnergyRange.dat",ios::in);
	for(ii=0;ii<23;ii++)
	{
		infile>>PrangeinSi[ii]>>PdrangeinSi[ii];
	}
	for (i=0;i<10;i++)
	{
		rin=0;
		//产生[x_low,x_high]区间均匀分布的随机数
		costheta=gRandom->Uniform(-1,1);//Uniform(x_low,x_high)
		z0=pz->GetRandom()/1000;
		if(costheta>0)rin=(z2D1-z0)/costheta;
		if(costheta<0)rin=(z1D1-z0)/costheta;
		if(rin>0)hrange->Fill(rin);
	}
	for(ie=0;ie<23;ie++)
	{
		countin=0;
		for (i=0;i<10000;i++)
		{
			//产生[x_low,x_high]区间均匀分布的随机数
			costheta=gRandom->Uniform(-1,1);//Uniform(x_low,x_high)
			phi=gRandom->Uniform(0,2*pi);
			rin=hrange->GetRandom();
			Ein=58285752852*pow(rin,6) - 13979019876*pow(rin,5) + 1271120957*pow(rin,4) - 53248138*pow(rin,3) + 1084962.98*pow(rin,2) + 12271.245*rin + 13.2193;
			//Ein=0.0000000583*pow(rin,6) - 0.0000139790*pow(rin,5) + 0.0012711210*pow(rin,4) - 0.0532481391*pow(rin,3) + 1.0849630189*pow(rin,2) + 12.2712438233*rin + 13.2193162773;
			if(Ein>3000)Ein=3000;
			if(Ein=3000)hEin->Fill(Ein);
			//cout<<Ein<<' ';
			r=gRandom->Gaus(PrangeinSi[ie],PdrangeinSi[ie]);//0.1-7 MeV calculated by LISE
			//r=gRandom->Gaus(0.09205,0.00406);//3 MeV SRIM
			//x0=gRandom->Gaus(0,1);//Gaus(mu,sigma), Default mu=0, sigma=1
			// 		x0=gRandom->Uniform(-0.8,0.8);//Gaus(mu,sigma), Default mu=0, sigma=1
			// 		y0=gRandom->Uniform(-0.8,0.8);//Gaus(mu,sigma), Default mu=0, sigma=1
			// 		z0=gRandom->Uniform(-0.8,0.8);
			//☆get a random number distributed according to a function
			x0=px->GetRandom();
			y0=py->GetRandom();
 			z0=z1D1+pz->GetRandom()/1000;
 			for(;z0>-10;){z0=z1D1+pz->GetRandom()/1000;}
//			z0=gRandom->Uniform(-10.073,-10.073);
			//☆get a random number distributed according the contents of a histogram or a user defined function
			//gRandom->Rannor(x0,y0);//Generate a pair of Gaussian random numbers with mu=0 and sigma=1
			//h2->Fill(x0,y0);
			x=x0+r*sqrt(1-costheta*costheta)*cos(phi);
			y=y0+r*sqrt(1-costheta*costheta)*sin(phi);
			z=z0+r*costheta;
			//gRandom->Sphere(x,y,z,r);
			//if(x*x+y*y+z*z<=r*r)
			//rin=sqrt((x-x0)*(x-x0)+(x-x0)*(x-x0)+(x-x0)*(x-x0));
			//if(x>x2D1||x<x1D1||y>y2D1||y<y1D1||z>z2D1||z<z1D1)
			if(z<z2D1&&z>z1D1)
			{
				Eloss=58285752852*pow(r,6) - 13979019876*pow(r,5) + 1271120957*pow(r,4) - 53248138*pow(r,3) + 1084962.98*pow(r,2) + 12271.245*r + 13.2193;
				//if(Eloss>3000)Eloss=3000;
				//cout<<Eloss<<' ';
				h3->Fill(y,x,z);
				h2->Fill(y,x);
				hEloss->Fill(Eloss);
				countin++;
			}
		}
		cout<<countin/10000.<<endl;
	}

	c3->cd();
	h3->SetLineColor(2);
	h3->GetXaxis()->SetTitle("x");
	h3->GetYaxis()->SetTitle("y");
	h3->GetZaxis()->SetTitle("z");
	h3->GetXaxis()->CenterTitle();
	h3->GetYaxis()->CenterTitle();
	h3->GetZaxis()->CenterTitle();
	h3->Draw(); 
	//h2->Draw();
	//h1->Draw();
	//hEloss->Draw();
}