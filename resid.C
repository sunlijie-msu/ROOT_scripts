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
void getbincontent()
{
	double counts[50]={5, 6, 5, 6, 5, 6, 5, 6, 5, 6, 5, 6, 5, 6, 5, 6, 5, 6, 5, 6, 5, 9, 20, 30, 40, 45, 40, 30, 20, 10, 6, 5, 6, 5, 6, 5, 6, 5, 6, 5, 6, 5, 6, 5, 6, 5, 6, 5, 6, 7};
	double x[50];
	double diff;
	TCanvas* c1   = new TCanvas("c1","c1",1000,700);
	TH1F *h1=new TH1F("h1","h1",50,1900,2000);
	TH1F *h2=new TH1F("h2","h2",50,1900,2000);
	TF1 *f1=new TF1("f1","gaus",1900,2000);
	for(int i=1;i<=50;i++)
	{
		h1->SetBinContent(i,counts[i-1]);
	}
	c1->cd();
	TPad *pad1 = new TPad("pad1", "The pad 70% of the height",0.0,0.3,1.0,1.0);// Double_t xlow, Double_t ylow, Double_t xup, Double_t yup,
	TPad *pad2 = new TPad("pad2", "The pad 30% of the height",0.0,0.0,1.0,0.3);
	pad1->Draw();
	pad2->Draw();
	pad1->cd();
	h1->Draw("E");
	h1->Fit(f1,"ME");

	for(int i=1;i<=50;i++)
	{
		diff=h1->GetBinContent(i)-f1->Eval(h1->GetBinCenter(i));//Fit-Data
		h2->SetBinContent(i,diff);
	}
	pad2->cd();
	h2->Draw("E");
}