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
void residuals1()
{
	double xmin=0; double xmax=5000; double xcenter;
	double binmin,binmax;
	double counts[50]={5, 6, 5, 6, 5, 6, 5, 6, 5, 6, 5, 6, 5, 6, 5, 6, 5, 6, 5, 6, 5, 9, 20, 30, 40, 45, 40, 30, 20, 10, 6, 5, 6, 5, 6, 5, 6, 5, 6, 5, 6, 5, 6, 5, 6, 5, 6, 5, 6, 7};
	TH1F *h1=new TH1F("h1","h1",2500,xmin,xmax);
	TH1F *h2 = new TH1F("h2","residuals",2500,xmin,xmax);
	TF1 *f1=new TF1("f1","gaus",xmin,xmax);
	for(int i=1;i<=50;i++)
	{
		h1->SetBinContent(i+1000,counts[i-1]);
	}
	xcenter=h1->GetBinCenter(h1->GetMaximumBin());
	xmin=xcenter-50;
	xmax=xcenter+50;
	binmin=h1->FindBin(xmin);
	binmax=h1->FindBin(xmax);
	
	gStyle->SetOptStat(0);
	TCanvas *c1 = new TCanvas("c1","c1",800,800);
	TPad *pad1 = new TPad("pad1","pad1",0,0.33,1,1);
	TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.33);
	pad1->SetBottomMargin(0.00001);
	pad1->SetBorderMode(0);
	pad1->SetLogy();
	pad2->SetTopMargin(0.00001);
	pad2->SetBottomMargin(0.1);
	pad2->SetBorderMode(0);
	pad1->Draw();
	pad2->Draw();
	pad1->cd();
	h1->GetXaxis()->SetRange(binmin,binmax);
	h1->Draw("E");
	gStyle->SetOptFit(1111);
	h1->SetStats();
	h1->Fit(f1,"E","",xmin,xmax);h1->Fit(f1,"E","",xmin,xmax);h1->Fit(f1,"E","",xmin,xmax);
	//double N = h1->GetXaxis()->GetNbins();

	pad2->cd();
	
	h2->GetXaxis()->SetLabelFont(63);
	h2->GetXaxis()->SetLabelSize(16);
	h2->GetXaxis()->SetTitle("channel energy (keV)");
	h2->GetYaxis()->SetLabelFont(63);
	h2->GetYaxis()->SetLabelSize(16);
	for (Int_t i=binmin;i<=binmax;i++)
	{
		Double_t diff = h1->GetBinContent(i)-f1->Eval(h1->GetBinCenter(i));
		h2->SetBinContent(i,diff);
	}
	h2->GetXaxis()->SetRange(binmin,binmax);
	h2->Draw("E");
	c1->cd();
}
