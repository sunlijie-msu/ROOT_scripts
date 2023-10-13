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
#include "TH3.h"
#include "TF1.h"
#include "TMath.h"
#include "TPaveText.h"
#include "TPad.h"
#include "string.h"
#include "stdlib.h"
#include "stdio.h"
#include <time.h>
using namespace std;
void RanD3effi27S()//能量相关，从D3发射一定能量质子，被各探测器探测到的概率
{
	unsigned long i,cntD1,cntD2,cntD3,cntQ1,cntQ2,cntQU,cntQD,cntQL,cntQR;
	int jj,ii,ie;
	time_t tim;
	struct tm *at;
	char now[80];
	float PrangeinSi[30],PdrangeinSi[30];
	int irawroot;
	char h_name[50];
	char anarootname[80];
	char simurootname[80];
	const float pi=3.14159265;
	float x1D1=-49.5/2.0,y1D1=-49.5/2.0,z1D1=-19-0.142/2,x2D1=49.5/2.0,y2D1=49.5/2.0,z2D1=-19+0.142/2.0;//DSSD2 implantation，x junction; x ohmic
	float x1D2=-49.5/2.0,y1D2=-49.5/2.0,z1D2=0-0.040/2,x2D2=49.5/2.0,y2D2=49.5/2.0,z2D2=0+0.040/2.0,x0=0,y0=0,z0=0;//DSSD2 implantation，x junction; y ohmic
	float x1D3=-49.5/2.0,y1D3=-49.5/2.0,z1D3=19-0.304/2.0,x2D3=49.5/2.0,y2D3=49.5/2.0,z2D3=19+0.304/2.0;//DSSD3 implantation，x junction; y ohmic
	float x1Q1=-24.7,y1Q1=-24.7,z1Q1=27,x2Q1=24.7,y2Q1=24.7,z2Q1=27.3;//QSD1 geometry
	float x1Q2=-24.7,y1Q2=-24.7,z1Q2=32,x2Q2=24.7,y2Q2=24.7,z2Q2=33.5;//QSD2 geometry
	float x1QU=-24.7,y1QU=49,z1QU=-63,x2QU=24.7,y2QU=50.5,z2QU=-13;//QSDU geometry
	float x1QD=-24.7,y1QD=-37.5,z1QD=-63,x2QD=24.7,y2QD=-36,z2QD=-13;//QSDD geometry
	float x1QR=42,y1QR=-23.7,z1QR=-50,x2QR=43.5,y2QR=25.7,z2QR=0;//QSDR geometry
	float x1QL=-43.5,y1QL=-23.7,z1QL=-50,x2QL=-42,y2QL=25.7,z2QL=0;//QSDL geometry
	//upstream x1y1z1, downstream x2y2z2
	float x,y,z,xd,yd,zd;
	float costheta,phi,r;
	double rin,Ein,Eloss;
	sprintf(anarootname,"D:/X/RIBLL2017/S27_0154_0345_p230_c250_t6T_E.root");
	TFile *fin = new TFile(anarootname);
	TTree *T999 = (TTree*)fin->Get("T999");
	unsigned long nentries=T999->GetEntries();//读事件数
	cout<<anarootname;
	cout<<"  Entries="<<nentries<<endl;
	TH1F *pxD2=(TH1F*)fin->Get("himpl40px");// = new TH1F("pxD3","pxD3",16,-23.3,23.3);//create a histogram
	TH1F *pyD2=(TH1F*)fin->Get("himpl40py");// = new TH1F("pyD3","pyD3",16,-23.3,23.3);//create a histogram
	TH1D *pzD2=(TH1D*)fin->Get("himpl40pz");
	TH1F *pxD3=(TH1F*)fin->Get("himpl304px");// = new TH1F("pxD3","pxD3",16,-23.3,23.3);//create a histogram
	TH1F *pyD3=(TH1F*)fin->Get("himpl304py");// = new TH1F("pyD3","pyD3",16,-23.3,23.3);//create a histogram
	TH1D *pzD3=(TH1D*)fin->Get("himpl304pz");
	// 	for(ii=0;ii<16;ii++)
	// 	{
	// 		sprintf(h_name,"%s%d","hSiimpl60pz",ii);
	// 		pzD3[ii] = (TH1D*)fin->Get(h_name);
	// 	}
	sprintf(simurootname,"%s","D:/X/RIBLL2017/S27simulateD3.root");
	TFile *fout = new TFile(simurootname,"RECREATE");//输出文件//It's better to define histograms and then define fout, in case of draw bugs.
	TTree *T111 = new TTree("T111","T111");//or TTree *T888 = new TTree("T888","Treetitle");
	TCanvas *c3=new TCanvas("c3","c3",640,640);
	// 	TH3F *h3 = new TH3F("h3","h3", 500,-2,2,500,-2,2,500,-2,2);
	// 	TH2F *h2Q2 = new TH2F("h2Q2","h2Q2", 500,-2,2,500,-2,2);
	// 	TH1F *hzQ2 = new TH1F("hzQ2","hzQ2", 500,-2,2);
	//TH3F *h3 = new TH3F("h3","h3", 16,-23.3,23.3,16,-23.3,23.3,66,0,0.066);
	TH3F *h3 = new TH3F("h3","h3", 100,-50,50,100,-50,50,100,-50,50);
	TH2F *hxyD1 = new TH2F("hxyD1","hxyD1",160,-24.75-1,24.75+1,160,-24.75-1,24.75+1);
	TH1F *hzD1 = new TH1F("hzD1","hzD1", 142+20,0-0.010,0.142+0.010);
	TH2F *hxyD2 = new TH2F("hxyD2","hxyD2",160,-24.75-1,24.75+1,160,-24.75-1,24.75+1);
	TH1F *hzD2 = new TH1F("hzD2","hzD2", 40+20,0-0.005,0.040+0.005);
	TH2F *hxyD3 = new TH2F("hxyD3","hxyD3",160,-24.75-1,24.75+1,160,-24.75-1,24.75+1);
	TH1F *hzD3 = new TH1F("hzD3","hzD3",304+20,0-0.010,0.304+0.010);
	TH2F *hxyQ1 = new TH2F("hxyQ1","hxyQ1", 160,-23.3*2,23.3*2,160,-23.3*2,23.3*2);
	TH1F *hzQ1 = new TH1F("hzQ1","hzQ1", 160,-23.3*2,23.3*2);
	//TH2F *h2Q2 = new TH2F("h2Q2","h2Q2", 160,-23.3*2,23.3*2,160,-23.3*2,23.3*2);
	//TH1F *h1Q2 = new TH1F("h1Q2","h1Q2", 160,-23.3*2,23.3*2);
	TH2F *hxyQU = new TH2F("hxyQU","hxyQU", 160,-23.3*2,23.3*2,160,-23.3*3,23.3*1);
	TH1F *hzQU = new TH1F("hzQU","hzQU", 160,-23.3*3,23.3*0);
	//TH2F *h2QD = new TH2F("h2QD","h2QD", 160,-23.3*2,23.3*2,160,-23.3*3,23.3*1);
	//TH1F *h1QD = new TH1F("h1QD","h1QD", 160,-23.3*3,23.3*0);
	//TH2F *h2QR = new TH2F("h2QR","h2QR", 160,-23.3*2,23.3*2,160,-23.3*3,23.3*1);
	//TH1F *h1QR = new TH1F("h1QR","h1QR", 160,-23.3*3,23.3*0.5);
	//TH2F *h2QL = new TH2F("h2QL","h2QL", 160,-23.3*2,23.3*2,160,-23.3*3,23.3*1);
	//TH1F *h1QL = new TH1F("h1QL","h1QL", 160,-23.3*3,23.3*0.5);
	TH1F *hEin = new TH1F("hEin","hEin", 1000,-1,4000);
	TH1F *hEloss = new TH1F("hEloss","hEloss", 1000,-1,4000);
	TH1D *hrange = new TH1D("hrange","hrange", 500,-1,20);
	//TH2F *h2Q2 = new TH2F("h2Q2","h2Q2", 16,-23.3,23.3,16,-23.3,23.3);
	// 	h3->FillRandom("gaus",5000);
	ofstream outfile("C:/Si24/Si22peakcali/EffiD3.dat",ios::out);
	ifstream infile("C:/Si24/Si22peakcali/EnergyRange.dat",ios::in);
	for(ii=0;ii<23;ii++)
	{
		infile>>PrangeinSi[ii]>>PdrangeinSi[ii];//质子在硅探测器中射程能损关系，射程、射程歧离0.1-10 MeV calculated by LISE
	}
	// 	for (i=0;i<1000000;i++)
	// 	{
	// 		rin=0;
	// 		//产生[x_low,x_high]区间均匀分布的随机数
	// 		costheta=gRandom->Uniform(-1,1);//Uniform(x_low,x_high)
	// 		z0=pzD3->GetRandom()/1000;
	// 		if(costheta>0)rin=(z2D3-z0)/costheta;
	// 		if(costheta<0)rin=(z1D3-z0)/costheta;
	// 		if(rin>0)hrange->Fill(rin);
	// 	}
	for(ie=0;ie<23;ie++)
	{
		cntD1=0; cntD2=0; cntD3=0; cntQ1=0; cntQ2=0; cntQU=0; cntQD=0; cntQR=0; cntQL=0;
		for (i=0;i<1000000;i++)
		{
			x=500; y=500; z=500;
			//产生[x_low,x_high]区间均匀分布的随机数
			costheta=gRandom->Uniform(-1,1);//Uniform(x_low,x_high)
			phi=gRandom->Uniform(0,2*pi);
			r=gRandom->Gaus(1*PrangeinSi[ie],PdrangeinSi[ie]);//0.1-10 MeV calculated by LISE
			//rin=hrange->GetRandom();
			//Ein=58285752852*pow(rin,6) - 13979019876*pow(rin,5) + 1271120957*pow(rin,4) - 53248138*pow(rin,3) + 1084962.98*pow(rin,2) + 12271.245*rin + 13.2193;
			//Ein=0.0000000583*pow(rin,6) - 0.0000139790*pow(rin,5) + 0.0012711210*pow(rin,4) - 0.0532481391*pow(rin,3) + 1.0849630189*pow(rin,2) + 12.2712438233*rin + 13.2193162773;
			//if(Ein>3000)Ein=3000;
			//if(Ein=3000)hEin->Fill(Ein);
			//cout<<Ein<<' ';
			//r=gRandom->Gaus(prangeinSi[ie],0.00);//3 MeV LISE
			//r=gRandom->Gaus(0.09205,0.00406);//3 MeV SRIM
			//x0=gRandom->Gaus(0,1);//Gaus(mu,sigma), Default mu=0, sigma=1
			// 		x0=gRandom->Uniform(-0.8,0.8);//Gaus(mu,sigma), Default mu=0, sigma=1
			// 		y0=gRandom->Uniform(-0.8,0.8);//Gaus(mu,sigma), Default mu=0, sigma=1
			// 		z0=gRandom->Uniform(-0.8,0.8);
			//☆get a random number distributed according to a function
			//x0=gRandom->Uniform(7,8);
			// 			x0=pxD2->GetRandom();
			// 			y0=pyD2->GetRandom();
			// 			z0=z1D2+pzD2->GetRandom()/1000;
			// 			//cout<<z0<<' ';
			// 			for(;z0>-10;){z0=z1D2+pzD2->GetRandom()/1000;}
			x0=pxD3->GetRandom();
			y0=pyD3->GetRandom();
			z0=pzD3->GetRandom()/1000+z1D3;
			//for(;z0>0.066;){z0=z1D3+pzD3->GetRandom()/1000;}
			//☆get a random number distributed according the contents of a histogram or a user defined function
			//gRandom->Rannor(x0,y0);//Generate a pair of Gaussian random numbers with mu=0 and sigma=1
			x=x0+r*sqrt(1-costheta*costheta)*cos(phi);
			y=y0+r*sqrt(1-costheta*costheta)*sin(phi);
			z=z0+r*costheta;
			if(z<z2D3&&z>z1D3)
			{
				cntD3++;//stop in D3
				hxyD3->Fill(x,y);//final xy position in D3
				hzD3->Fill(z-z1D3);//final z position in D3
				h3->Fill(x,y,z);
			}
			else if(z>=z2D3)//escape out of D3, towards z+, r should be lengthened
			{
				x=x+(z1Q1-z2D3)*tan(acos(costheta))*cos(phi);
				y=y+(z1Q1-z2D3)*tan(acos(costheta))*sin(phi);
				z=z+(z1Q1-z2D3);
				if(x<x2Q1&&x>x1Q1&&y<y2Q1&&y>y1Q1&&z<z2Q1&&z>z1Q1)//if stop in Q1
				{
					cntQ1++;
					hxyQ1->Fill(x,y);//final xy position in Q1
				}
				else if(z>=z2Q1)//if pass through Q1, r should be lengthened
				{
					x=x+(z1Q2-z2Q1)*tan(acos(costheta))*cos(phi);
					y=y+(z1Q2-z2Q1)*tan(acos(costheta))*sin(phi);
					z=z+(z1Q2-z2Q1);
					if(x<x2Q2&&x>x1Q2&&y<y2Q2&&y>y1Q2&&z<z2Q2&&z>z1Q2)//if stop in Q2
						cntQ2++;
				}
			}
			else if(z<=z1D3)//escape out of D3, towards z-, r should be lengthened
			{
				x=x-(z1D3-z2D2)*tan(acos(costheta))*cos(phi);
				y=y-(z1D3-z2D2)*tan(acos(costheta))*sin(phi);
				z=z-(z1D3-z2D2);
				if(x<x2D2&&x>x1D2&&y<y2D2&&y>y1D2&&z<z2D2&&z>z1D2)//if stop in D2
				{
					cntD2++;//stop in D2
					hxyD2->Fill(x,y);//final xy position in D2
					hzD2->Fill(z-z1D2);//final z position in D2
					h3->Fill(x,y,z);
				}
				else if(z<=z1D2)//if pass through D2, r should be lengthened
				{
					x=x-(z1D2-z2D1)*tan(acos(costheta))*cos(phi);
					y=y-(z1D2-z2D1)*tan(acos(costheta))*sin(phi);
					z=z-(z1D2-z2D1);
					if(x<x2D1&&x>x1D1&&y<y2D1&&y>y1D1&&z<z2D1&&z>z1D1)//if stop in D1
					{
						cntD1++;//stop in D1
						hxyD1->Fill(x,y);//final xy position in D1
						hzD1->Fill(z-z1D1);//final z position in D1
						h3->Fill(x,y,z);//final z position in D1
					}
				}
				else if(z<=z1D1)//if pass through D1, no need to judge whether stop or transmit QSDs
				{
					yd=y1QU-y0;//QSDU useless for 27S
					if(tan(acos(costheta))>-1000&&tan(acos(costheta))<1000)
					{
						x=x0+yd/tan(phi);//x0+r*sqrt(1-costheta*costheta)*cos(phi);
						z=z0+yd/tan(acos(costheta))/sin(phi);//y0+r*sqrt(1-costheta*costheta)*sin(phi);
					}
					if(tan(acos(costheta))>-1000&&tan(acos(costheta))<1000&&x<x2QU&&x>x1QU&&z<z2QU&&z>z1QU)
					{
						//h3->Fill(y1QU,x,z);
						hxyQU->Fill(x,z);
						hzQU->Fill(z);
						cntQU++;
					}
				}
			}
		}//for (i=0;i<10000000;i++)
		outfile<<"ie=	"<<ie<<"	cntD1=	"<<cntD1<<"	cntD2=	"<<cntD2<<"	cntD3=	"<<cntD3<<endl;
		//divided by 2, since every xd,yd,zd has a mirror value. tan(theta)-->(-1000,1000) in case of variable overflow
		//restrict phi for Q_LRUD, or restrict theta for Q2 can remove the mirror double
	}//(ie=0;ie<23;ie++)
	//T111->Fill();
	c3->cd();
	h3->SetLineColor(2);
	h3->GetXaxis()->SetTitle("x");
	h3->GetYaxis()->SetTitle("y");
	h3->GetZaxis()->SetTitle("z");
	h3->GetXaxis()->CenterTitle();
	h3->GetYaxis()->CenterTitle();
	h3->GetZaxis()->CenterTitle();
	h3->Draw(); 
	fout->Write();
	fout->Close();
}