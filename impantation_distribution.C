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
void impantation_distribution()//能量无关，全穿出，周边立体角
{
	unsigned long i,cntD1,cntD2,cntQ1,cntQ2,cntQU,cntQD,cntQL,cntQR;
	int jj,ii,ie;
	time_t tim;
	struct tm *at;
	char now[80];
	float prangeinSi[20];
	int irawroot;
	char h_name[50];
	char anarootname[80];
	char simurootname[80];
	const float pi=3.14159265;
// 	float x1D1=-24.7,y1D1=-21.7,z1D1=-10.149,x2D1=24.7,y2D1=27.7,z2D1=-10,x0=0,y0=0,z0=0;//DSSD1 implantation，x junction; y ohmic
// 	float x1D2=-24.7,y1D2=-24.7,z1D2=0,x2D2=24.7,y2D2=24.7,z2D2=0.066;//DSSD2 implantation，x junction; y ohmic
// 	float x1Q1=-24.7,y1Q1=-24.7,z1Q1=7,x2Q1=24.7,y2Q1=24.7,z2Q1=7.3;//QSD1 geometry
// 	float x1Q2=-24.7,y1Q2=-24.7,z1Q2=12,x2Q2=24.7,y2Q2=24.7,z2Q2=13.5;//QSD2 geometry
// 	float x1QU=-24.7,y1QU=49,z1QU=-63,x2QU=24.7,y2QU=50.5,z2QU=-13;//QSDU geometry
// 	float x1QD=-24.7,y1QD=-37.5,z1QD=-63,x2QD=24.7,y2QD=-36,z2QD=-13;//QSDD geometry
// 	float x1QR=42,y1QR=-23.7,z1QR=-50,x2QR=43.5,y2QR=25.7,z2QR=0;//QSDR geometry
// 	float x1QL=-43.5,y1QL=-23.7,z1QL=-50,x2QL=-42,y2QL=25.7,z2QL=0;//QSDL geometry
	//upstream x1y1z1, downstream x2y2z2
// 	float x,y,z,xd,yd,zd;
// 	float costheta,phi,r;
// 	double rin,Ein,Eloss;
// 	sprintf(anarootname,"V:/RIBLL2015/data24/Si24calabcd_0636_0650.root");
// 	TFile *fin = new TFile(anarootname);
// 	TTree *T999 = (TTree*)fin->Get("T999");
// 	unsigned long nentries=T999->GetEntries();//读事件数
// 	cout<<anarootname;
// 	cout<<"  Entries="<<nentries<<endl;
// 	TH1F *pxD1=(TH1F*)fin->Get("hSiimpl300px");// = new TH1F("pxD2","pxD2",16,-23.3,23.3);//create a histogram
// 	TH1F *pyD1=(TH1F*)fin->Get("hSiimpl300py");// = new TH1F("pyD2","pyD2",16,-23.3,23.3);//create a histogram
// 	TH1D *pzD1=(TH1D*)fin->Get("hSiimpl300pz");
// 	TH1F *pxD2=(TH1F*)fin->Get("hSiimpl60px");// = new TH1F("pxD2","pxD2",16,-23.3,23.3);//create a histogram
// 	TH1F *pyD2=(TH1F*)fin->Get("hSiimpl60py");// = new TH1F("pyD2","pyD2",16,-23.3,23.3);//create a histogram
// 	TH1D *pzD2=(TH1D*)fin->Get("hSiimpl60pz");
// // 	for(ii=0;ii<16;ii++)
// // 	{
// // 		sprintf(h_name,"%s%d","hSiimpl60pz",ii);
// // 		pzD2[ii] = (TH1D*)fin->Get(h_name);
// // 	}
	sprintf(simurootname,"%s","C:/Si24/Si22peakcali/S27simulate.root");
	TFile *fout = new TFile(simurootname,"RECREATE");//输出文件
	//TTree *T111 = new TTree("T111","T111");
	//TCanvas *c3=new TCanvas("c3","c3",640,640);
	TH2F *hx2y2 = new TH2F("hx2y2","hx2y2", 50,-25,25,50,-25,25);
	TH1F *hx2 = new TH1F("hx2","hx2", 50,-25,25);
	TH1F *hy2 = new TH1F("hy2","hy2", 50,-25,25);
	TH1F *hz1 = new TH1F("hz1","hz1", 60,9.910,9.970);
	TH1F *hz2 = new TH1F("hz2","hz2", 60,9.970,10.030);
	TH1F *hz3 = new TH1F("hz3","hz3", 300,10.030,10.330);
	for (int ie=0;ie<100000;ie++)
	{
		hx2y2->Fill(gRandom->Gaus(0,12),gRandom->Gaus(0,12));//Add one value at a time
		hx2->Fill(gRandom->Gaus(0,12));//Add one value at a time
		hy2->Fill(gRandom->Gaus(0,12));//Add one value at a time
		hz1->Fill(gRandom->Gaus(10.000,0.020));//Add one value at a time
		hz2->Fill(gRandom->Gaus(10.000,0.020));//Add one value at a time
		hz3->Fill(gRandom->Gaus(10.000,0.020));//Add one value at a time
	}
	fout->Write();
	fout->Close();
}