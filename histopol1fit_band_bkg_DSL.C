#include <iostream>
#include <fstream>
#include <iomanip.h>
#include <math.h>
#include <time.h>
#include <map>
#include <TROOT.h>
#include <TSystem.h>
#include <TCutG.h>
#include "TChain.h"
#include "TStyle.h"
#include "TLine.h"
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
// example to illustrate how to fit excluding points in a given range

// 	const double E0_gamma=1248; //Ex=1248
// 	const int startpeak=1234; //has to be even number
// 	const int endpeak=1352; //has to be even number
// 	const int up=200; //E0_gamma+up has to be even number, if E0 is odd, up should be odd
// 	const int down=100; //E0_gamma-down has to be even number, if E0 is odd, down should be odd
// 	const int Ea=49;
// 	double T0_lifetime = 0;// =0 read in simulation files   // This is Excited State Mean Lifetime (fs), Mean Lifetime τ=T/ln(2); half-life T=τ*ln(2)
// 	int Lifetimestep = 100; // read in simulation files
// 	double centroid=E0_gamma; //centroid=E0_gamma; for even energy 1246, 1248, 1250, etc; centroid=E0_gamma-1; for odd energy 1245, 1247, 1249, etc;
// 	int Yaxisrange=60;

// 	const double E0_gamma=2234; //Ex=2234
// 	const int startpeak=2222; //has to be even number
// 	const int endpeak=2420; //has to be even number
// 	const int up=320; //E0_gamma+up has to be even number, if E0 is odd, up should be odd
// 	const int down=140; //E0_gamma-down has to be even number, if E0 is odd, down should be odd
// 	const int Ea=47;
// 	double T0_lifetime = 80; //=80
// 	int Lifetimestep = 10;
// 	double centroid=E0_gamma;
// 	int Yaxisrange=20;

// const double E0_gamma=3076; //Ex=3076
// const int startpeak=3246; //has to be even number
// const int endpeak=3326; //has to be even number
// const int up=354; //E0_gamma+up has to be even number, if E0 is odd, up should be odd
// const int down=-70; //E0_gamma-down has to be even number, if E0 is odd, down should be odd
// const int Ea=46;
// double T0_lifetime = 0;
// int Lifetimestep = 1;
// double centroid=E0_gamma;
// int Yaxisrange=15;

// 	const double E0_gamma=4971; //Ex=4971
// 	const int startpeak=5268; //has to be even number
// 	const int endpeak=5370; //has to be even number
// 	const int up=531; //E0_gamma+up has to be even number, if E0 is odd, up should be odd
// 	const int down=-153; //E0_gamma-down has to be even number, if E0 is odd, down should be odd
// 	const int Ea=42;
// 	double T0_lifetime = 0;
// 	int Lifetimestep = 1;
// 	double centroid=E0_gamma-1;
// 	int Yaxisrange=15;

// 	const double E0_gamma=5156; //Ex=5156
// 	const int startpeak=5450; //has to be even number
// 	const int endpeak=5584; //has to be even number
// 	const int up=570; //E0_gamma+up has to be even number, if E0 is odd, up should be odd
// 	const int down=-198; //E0_gamma-down has to be even number, if E0 is odd, down should be odd
// 	const int Ea=42;
// 	double T0_lifetime = 0;
// 	int Lifetimestep = 1;
// 	double centroid=E0_gamma;
// 	int Yaxisrange=15;

// 	const double E0_gamma=3435; //Ex=3435
// 	const int startpeak=3636; //has to be even number
// 	const int endpeak=3716; //has to be even number
// 	const int up=381; //E0_gamma+up has to be even number, if E0 is odd, up should be odd
// 	const int down=-81; //E0_gamma-down has to be even number, if E0 is odd, down should be odd
// 	const int Ea=45;
// 	double T0_lifetime = 0;
// 	int Lifetimestep = 1;
// 	double centroid=E0_gamma;
// 	int Yaxisrange=15;

// 	const double E0_gamma=2186; //Ex=3435
// 	const int startpeak=2314; //has to be even number
// 	const int endpeak=2364; //has to be even number
// 	const int up=242; //E0_gamma+up has to be even number, if E0 is odd, up should be odd
// 	const int down=-40; //E0_gamma-down has to be even number, if E0 is odd, down should be odd
// 	const int Ea=45;
// 	double T0_lifetime = 0;
// 	int Lifetimestep = 1;
// 	double centroid=E0_gamma;
// 	int Yaxisrange=15;

// 	const double E0_gamma=2838; //Ex=4085
// 	const int startpeak=2998; //has to be even number
// 	const int endpeak=3070; //has to be even number
// 	const int up=316; //E0_gamma+up has to be even number, if E0 is odd, up should be odd
// 	const int down=-66; //E0_gamma-down has to be even number, if E0 is odd, down should be odd
// 	const int Ea=44;
// 	double T0_lifetime = 0;
// 	int Lifetimestep = 1;
// 	double centroid=E0_gamma;
// 	int Yaxisrange=15;

// 	const double E0_gamma=4270; //Ex=5518
// 	const int startpeak=4520; //has to be even number
// 	const int endpeak=4610; //has to be even number
// 	const int up=452; //E0_gamma+up has to be even number, if E0 is odd, up should be odd
// 	const int down=-126; //E0_gamma-down has to be even number, if E0 is odd, down should be odd
// 	const int Ea=41;
// 	double T0_lifetime = 0;
// 	int Lifetimestep = 1;
// 	double centroid=E0_gamma;
// 	int Yaxisrange=15;

// 	const double E0_gamma=3541; //Ex=5775
// 	const int startpeak=3742; //has to be even number
// 	const int endpeak=3834; //has to be even number
// 	const int up=415; //E0_gamma+up has to be even number, if E0 is odd, up should be odd
// 	const int down=-71; //E0_gamma-down has to be even number, if E0 is odd, down should be odd
// 	const int Ea=40;
// 	double T0_lifetime = 0;
// 	int Lifetimestep = 1;
// 	double centroid=E0_gamma;
// 	int Yaxisrange=15;

// 	const double E0_gamma=5141; //Ex=6390
// 	const int startpeak=5464; //has to be even number
// 	const int endpeak=5554; //has to be even number
// 	const int up=551; //E0_gamma+up has to be even number, if E0 is odd, up should be odd
// 	const int down=-181; //E0_gamma-down has to be even number, if E0 is odd, down should be odd
// 	const int Ea=39;
// 	double T0_lifetime = 0;
// 	int Lifetimestep = 1;
// 	double centroid=E0_gamma;
// 	int Yaxisrange=15;

// 	const double E0_gamma=4156; //Ex=6390
// 	const int startpeak=4404; //has to be even number
// 	const int endpeak=4494; //has to be even number
// 	const int up=450; //E0_gamma+up has to be even number, if E0 is odd, up should be odd
// 	const int down=-110; //E0_gamma-down has to be even number, if E0 is odd, down should be odd
// 	const int Ea=39;
// 	double T0_lifetime = 0;
// 	int Lifetimestep = 1;
// 	double centroid=E0_gamma;
// 	int Yaxisrange=15;

// 	const double E0_gamma=5294; //Ex=6542
// 	const int startpeak=5600; //has to be even number
// 	const int endpeak=5720; //has to be even number
// 	const int up=560; //E0_gamma+up has to be even number, if E0 is odd, up should be odd
// 	const int down=-150; //E0_gamma-down has to be even number, if E0 is odd, down should be odd
// 	const int Ea=39;
// 	double T0_lifetime = 0;
// 	int Lifetimestep = 1;
// 	double centroid=E0_gamma-1;
// 	int Yaxisrange=15;

double fitrangelow1, fitrangelow2, fitrangehigh1, fitrangehigh2;
const double Ebin_gamma=E0_gamma;
double minrange=Ebin_gamma-down;
double maxrange=Ebin_gamma+up;
fitrangelow1=minrange; fitrangelow2=startpeak; fitrangehigh1=endpeak; fitrangehigh2=maxrange;

Bool_t reject;
Double_t fline(Double_t *x, Double_t *par)
{ // Double_t fcn(Double_t *x, Double_t *params)
	if (reject && x[0] > startpeak && x[0] < endpeak) {
		TF1::RejectPoint();
		return 0;
	}
	return par[0] + par[1]*x[0];
}
void histopol1fit_band_bkg_DSL()// example to illustrate how to fit excluding points in a given range
{
	char rawrootname[300],hname[300];
	sprintf(rawrootname,"%s%.0f%s","D:/X/out/S31_Gamma", E0_gamma, ".root");//modify formal Eg
	TFile *fin = new TFile(rawrootname);//after this statement, you can use any ROOT command1 for this rootfile
	//TTree *tree = (TTree*)fin->Get("tree");

	TF1 *f1=new TF1("f1","[0]+x*[1]",0,9000);
	TF1 *f2=new TF1("f2","[0]+x*[1]",0,9000);
	sprintf(hname,"%s%.0f","S31_Gamma", E0_gamma);
	TCanvas *canvash=new TCanvas(hname,hname,1000,400);
	canvash->cd();
	canvash->SetRightMargin(0.02);
	hSega->Draw("e");
	hSega->SetTitle(hname);

	TF1 *fl = new TF1("fl",fline,fitrangelow1,fitrangehigh2,2); //npar=2 is the number of free parameters used by the function 
	//we want to fit only the linear background excluding the signal area
	reject = kTRUE;
	hSega->Fit(fl,"L0"); // "1" will draw a full curve, "0" will not draw the fitted curve in the excluded region
	// option "L" a likelihood fit is used instead of the default chi2 square fit. 
	reject = kFALSE;

	//store 3 separate functions for visualization (three bands)
	TH1D *hint_tau3 = new TH1D("hint_tau3", "Fitted func with conf.band", (fitrangehigh2-fitrangelow1)*10, fitrangelow1,fitrangehigh2);
	(TVirtualFitter::GetFitter())->GetConfidenceIntervals(hint_tau3, 0.683);
	hint_tau3->SetStats(kFALSE);
	hint_tau3->SetFillColor(kGreen-9);
	hint_tau3->Draw("e3 same");

	TH1D *hint_tau1 = new TH1D("hint_tau1", "Fitted func with conf.band", (fitrangelow2-fitrangelow1)*10, fitrangelow1,fitrangelow2);
	(TVirtualFitter::GetFitter())->GetConfidenceIntervals(hint_tau1, 0.683);
	hint_tau1->SetStats(kFALSE);
	hint_tau1->SetFillColor(kRed-9);
	hint_tau1->Draw("e3 same");

	TH1D *hint_tau2 = new TH1D("hint_tau2", "Fitted func with conf.band", (fitrangehigh2-fitrangehigh1)*10, fitrangehigh1,fitrangehigh2);
	(TVirtualFitter::GetFitter())->GetConfidenceIntervals(hint_tau2, 0.683);
	hint_tau2->SetStats(kFALSE);
	hint_tau2->SetFillColor(kRed-9);
	hint_tau2->Draw("e3 same");

	//store 2 separate functions for visualization (two lines)
	TF1 *fleft = new TF1("fleft",fline,fitrangelow1,fitrangelow2,2);
	fleft->SetParameters(fl->GetParameters());
	hSega->GetListOfFunctions()->Add(fleft);
	TF1 *fright = new TF1("fright",fline,fitrangehigh1,fitrangehigh2,2);
	fright->SetParameters(fl->GetParameters());
	hSega->GetListOfFunctions()->Add(fright);

	char paraprint[100];
	double parChi, parNDF, p_value, p0, p1, p0err, p1err;
	p0 = fl->GetParameter(0);
	p1 = fl->GetParameter(1);
	p0err = fl->GetParError(0);
	p1err = fl->GetParError(1);
	parChi = fl->GetChisquare();
	parNDF = fl->GetNDF();
	p_value = fl->GetProb();
	TPaveText *textchi = new TPaveText(0.12,0.60,0.38,0.98,"brNDC");//left, down, right, up
	textchi->SetBorderSize(1);
	textchi->SetFillColor(0);
	textchi->SetTextAlign(12);//align = 10*HorizontalAlign + VerticalAlign, 12 means水平左对齐、垂直居中对齐
	textchi->SetTextFont(132);//font = 10 * fontID + precision, 12+2,12 means Symbol; 13+2, 13 means News Times Roman
	sprintf(paraprint,"Chi2=%.2f",parChi);
	textchi->AddText(paraprint);
	sprintf(paraprint,"NDF=%.0f",parNDF);
	textchi->AddText(paraprint);
	sprintf(paraprint,"p-val=%f",p_value);
	textchi->AddText(paraprint);
	sprintf(paraprint,"y=p1*x+p0");
	textchi->AddText(paraprint);
	sprintf(paraprint,"p1=%.8f+/-%.8f",p1,p1err);
	textchi->AddText(paraprint);
	sprintf(paraprint,"p0=%.5f+/-%.5f",p0,p0err);
	textchi->AddText(paraprint);
	textchi->Draw();
}//peakcali main