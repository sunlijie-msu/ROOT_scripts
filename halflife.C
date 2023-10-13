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
void halflife()//20Mg Tfit no daughter, almost useless
{
	double A,T,B; double Ae,Te,Be; double Bg;
	float ch;
	TFile *fin = new TFile("D:/X/RIBLL2017/data27/T999/S27_0154_0154first.root");//after this statement, you can use any ROOT command for this rootfile
	ofstream outfile("C:/Si24/Si22peakcali/Tdfit.dat",ios::out);
	//TH1F*hT=new TH1F(*TD40);//copy any histogram
	TH1F *T142_40 = new TH1F("T142_40","T142_40",450,1,4501);
	TH1F *T142_A = new TH1F("T142_A","T142_A",4500,1,9001); //与质子谱符合
	TH1F *T40_A = new TH1F("T40_A","T40_A",4500,1,9001); //与质子谱符合
	TH1F *T142_B = new TH1F("T142_B","T142_B",4500,1,9001); //与质子谱符合
	TH1F *T40_B = new TH1F("T40_B","T40_B",4500,1,9001); //与质子谱符合
	T999->Draw("T142>>T142_A","D142Ane+D142Bne>1000&&D142Ane+D142Bne<13600"); //与质子谱符合
	T999->Draw("T40>>T40_A","D40Ane+D40Bne>1000&&D40Ane+D40Bne<13600"); //与质子谱符合，T正式结果
	T999->Draw("T142>>T142_B","D142Ane+D142Bne>1000&&D142Ane+D142Bne<13600"); //与质子谱符合
	T999->Draw("T40>>T40_B","D40Ane+D40Bne>1000&&D40Ane+D40Bne<13600"); //与质子谱符合
	T142_40->Add(T142_40,T142_A); //与质子谱符合
	T142_40->Add(T142_40,T40_A); //与质子谱符合
	T142_40->Add(T142_40,T142_B); //与质子谱符合
	T142_40->Add(T142_40,T40_B); //与质子谱符合
	TF1 *SiDEC=new TF1("SiDEC","[0]*exp(x/(-[1]/0.693147))+[2]",1,5000);//自定义拟合函数
	SiDEC->SetParNames("A","T","B");
	SiDEC->SetParameters(300,120,1.1);//自定义的拟合函数必须赋初值
	SiDEC->SetParLimits(2,0,20);
// 	TH1F *T142_40 = new TH1F("T142_40","T142_40",4200,1,8401);
// 	T142_40->Add(T142_40,TD142);//
// 	T142_40->Add(T142_40,TD40);
// 	T142_40->Rebin(5);
// 	TF1 *SiDEC=new TF1("SiDEC","[0]*exp(x/(-[1]/0.693147))+[2]",1,9000);//自定义拟合函数
// 	SiDEC->SetParNames("A","T","B");
// 	SiDEC->SetParameters(2000,120,0.673795);//自定义的拟合函数必须赋初值
// 	//SiDEC->SetParLimits(2,0.673795,0.673795);//fix B
	TF1 *SiBG=new TF1("SiBG","[0]",1,9000);
	SiBG->SetParNames("Bg");
	for(int ii=10;ii<11;ii++)
	{
		ch=ii*5*90;
		T142_40->Fit("SiDEC","L","",1,ch);//specify a range in the Fit, recommended
		//T142_40->Fit("SiBG","L","",4200,ch);//specify a range in the Fit, recommended
		A=SiDEC->GetParameter(0);
		T=SiDEC->GetParameter(1);
		B=SiDEC->GetParameter(2);
		Ae=SiDEC->GetParError(0);
		Te=SiDEC->GetParError(1);
		Be=SiDEC->GetParError(2);
		outfile<<ii*5<<"	A=	"<<A<<"	Ae=	"<<Ae<<"	T=	"<<T<<"	Te=	"<<Te<<"	B=	"<<B<<"	Be=	"<<Be<<endl;
		//outfile<<Bg<<endl;
	}
}//Tf main