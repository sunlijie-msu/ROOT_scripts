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
#include "TMatrixDSymfwd.h"
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
void likelihood()
{
	double xmin=0; double xmax=100; double xcenter;
	double binmin=1,binmax=100;
	double myvalue;
	TH1F *h1=new TH1F("h1","h1",100,xmin,xmax);
	for (long j=0;j<100000;j++)
	{
		myvalue = gRandom->Gaus(40,7);
		h1->Fill(myvalue);//Add one x-value at a time
	}
	TH1F *h2 = new TH1F("h2","h2",100,xmin,xmax);
	TF1 *f1=new TF1("f1","gaus",xmin,xmax);
// 	for(int i=1;i<=50;i++)
// 	{
// 		h1->SetBinContent(i,counts[i-1]);
// 	}
// 	xcenter=h1->GetBinCenter(h1->GetMaximumBin());
// 	xmin=xcenter-50;
// 	xmax=xcenter+50;
// 	binmin=h1->FindBin(xmin);
// 	binmax=h1->FindBin(xmax);
// 	
	//for(int k=0;k<10;k++)
	{
	TCanvas* c1   = new TCanvas("c1","c1",1600,800);
	c1->cd();
	TPad *pad1 = new TPad("pad1", "The pad 70% of the height",0.0,0.3,1.0,1.0);// Double_t xlow, Double_t ylow, Double_t xup, Double_t yup,
	TPad *pad2 = new TPad("pad2", "The pad 30% of the height",0.0,0.0,1.0,0.3);
	pad1->SetTopMargin(0.02);
	pad1->SetRightMargin(0.02);
	pad1->SetLeftMargin(0.04);
	pad1->SetBottomMargin(0.08);
	//pad1->SetBorderMode(0);
	pad2->SetTopMargin(0.01);
	pad2->SetRightMargin(0.02);
	pad2->SetLeftMargin(0.04);
	pad2->SetBottomMargin(0.2);
	//pad2->SetBorderMode(0);
	pad1->Draw();
	pad2->Draw();
	pad1->cd();
	gStyle->SetOptTitle(0);
	gStyle->SetFitFormat("7.6g");
	h1->Sumw2(kFALSE);
	h1->SetBinErrorOption(TH1::kPoisson);//TH1::kNormal or TH1::kPoisson
//	h1->GetXaxis()->SetRange(binmin,binmax);
	h1->SetLineWidth(2);
	h1->Draw("E");
	gStyle->SetOptFit(1111);
	h1->SetStats();
	//f1->SetParLimits(1,35+k,35+k+0.00001);
	//h1->Fit(f1,"L","",xmin,xmax);
	TFitResultPtr r = h1->Fit(f1,"LS");
	//cout<<"Chi2="<<r->Chi2()<<endl;
	//cout<<"L=	"<<r->MinFcnValue()<<endl;
	//double N = h1->GetXaxis()->GetNbins()
	r->Print("V");
	pad2->cd();
	
	h2->GetXaxis()->SetLabelFont(63);
	h2->GetXaxis()->SetLabelSize(16);
//	h2->GetXaxis()->SetTitle("channel energy (keV)");
	h2->GetYaxis()->SetLabelFont(63);
	h2->GetYaxis()->SetLabelSize(16);
	h2->SetLineWidth(2);
	double chisquareNeyman=0,chisquarePearson=0,err_ydata=0,likelihood=0,residual,yfit,ydata;
	for (int i=binmin;i<=binmax;i++)
	{
		ydata=h1->GetBinContent(i);
		err_ydata=h1->GetBinError(i);
		yfit=f1->Eval(h1->GetBinCenter(i));
		residual=ydata-yfit;
		h2->SetBinContent(i,residual);
		if(ydata!=0) {chisquareNeyman += residual*residual/ydata;}//Neyman Chi2 method
		if(yfit!=0) {chisquarePearson += residual*residual/yfit;}//Pearson Chi2 method
		if(yfit!=0&&ydata!=0) {likelihood+=yfit-ydata+ydata*log(ydata/yfit);}
		if(yfit!=0&&ydata==0) {likelihood+=yfit-ydata;}
		//if(yfit!=0&&ydata!=0) 	cout<<h1->GetBinCenter(i)<<"	"<<yfit-ydata+ydata*log(ydata/yfit)<<endl;
	}
//	cout<<"chisquareNeyman=	"<<chisquareNeyman<<endl;
//	cout<<"chisquarePearson=	"<<chisquarePearson<<endl;
	cout<<"-lnL=	"<<likelihood<<endl;
//	cout<<"NDF=	"<<f1->GetNDF()<<endl;
//	h2->GetXaxis()->SetRange(binmin,binmax);
	h2->Draw("E");
	c1->cd();
	}
}
