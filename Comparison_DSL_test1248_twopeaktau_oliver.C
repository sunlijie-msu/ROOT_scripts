#include "TH1F.h"
//#include <cmath> //can't use pow() with this header
#include <stdlib.h>
#include "TMinuit.h"
#include "TFumili.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TStyle.h"
#include "TLine.h"
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
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "TRandom.h"
#include <TRandom3.h>
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
#include "TLegend.h"
#include "TPaletteAxis.h"
using namespace std;
//main_readhists; fcn(); comparehists(); main_draw_save;//search m o d i f y to change peak

TFile *fin_simu0, *fin_simu1, *fin_data;
TH1F *hFit,*hSega,*hSega_data,*bkgshort,*background;//

const double Binsperkev=0.5; //Number of bins per keV, the only place to change binwidth, all other variables is related to Binsperkev // 10 for NSCL, 2 for RIBLL, useless for DSL
const int binwidth=2; //binwidths in units of keV
const int factor_rebin=2; //simu and data Rebin factor
float bkgDown=0.99;
float bkgUp=1.01;
double ULbkg=1.1;
double LLbkg=0.9;

//Enable one of the energies and others should be commented out.

	const double E0_gamma=1248; //Ex=1248
	const int startpeak=1234; //has to be even number
	const int endpeak=1352; //has to be even number
	const int up=200; //E0_gamma+up has to be even number, if E0 is odd, up should be odd
	const int down=100; //E0_gamma-down has to be even number, if E0 is odd, down should be odd
	const int Ea=49;
	double T0_lifetime = 0;// =0 read in simulation files   // This is Excited State Mean Lifetime (fs), Mean Lifetime τ=T/ln(2); half-life T=τ*ln(2)
	int Lifetimestep = 1; // read in simulation files
	double centroid=1248.40; //centroid=E0_gamma; for even energy 1246, 1248, 1250, etc; centroid=E0_gamma-1; for odd energy 1245, 1247, 1249, etc;
	double Eg_uncertainty = 0.20;
	int Yaxisrange=200;
	ULbkg=3.8/3.2;
	LLbkg=2.7/3.2;

// 	const double E0_gamma=2234; //Ex=2234
// 	const int startpeak=2222; //has to be even number
// 	const int endpeak=2420; //has to be even number
// 	const int up=320; //E0_gamma+up has to be even number, if E0 is odd, up should be odd
// 	const int down=140; //E0_gamma-down has to be even number, if E0 is odd, down should be odd
// 	const int Ea=47;
// 	double T0_lifetime = 0; //=80
// 	int Lifetimestep = 1;
// 	double centroid=2233.97;
// 	double Eg_uncertainty = 0.20;
// 	int Yaxisrange=20;
// 	ULbkg=3.3/2.85;
// 	LLbkg=2.5/2.85;


// 	const double E0_gamma=3076; //Ex=3076
// 	const int startpeak=3246; //has to be even number
// 	const int endpeak=3326; //has to be even number
// 	const int up=354; //E0_gamma+up has to be even number, if E0 is odd, up should be odd
// 	const int down=-70; //E0_gamma-down has to be even number, if E0 is odd, down should be odd
// 	const int Ea=46;
// 	double T0_lifetime = 0;
// 	int Lifetimestep = 1;
// 	double centroid = 3076.24;
// 	double Eg_uncertainty = 0.2;
// 	int Yaxisrange=15;
// 	ULbkg=1.69/1.5;
// 	LLbkg=1.3/1.5;

// 	const double E0_gamma=4970; //Ex=4970
// 	const int startpeak=5268; //has to be even number
// 	const int endpeak=5370; //has to be even number
// 	const int up=532; //E0_gamma+up has to be even number, if E0 is odd, up should be odd
// 	const int down=-154; //E0_gamma-down has to be even number, if E0 is odd, down should be odd
// 	const int Ea=42;
// 	double T0_lifetime = 0;
// 	int Lifetimestep = 1;
// 	double centroid=4970.20;
// 	double Eg_uncertainty = 0.90;
// 	int Yaxisrange=15;
// 	ULbkg=1.8/1.6;
// 	LLbkg=1.4/1.6;

// 	const double E0_gamma=5156; //Ex=5156
// 	const int startpeak=5450; //has to be even number
// 	const int endpeak=5584; //has to be even number
// 	const int up=570; //E0_gamma+up has to be even number, if E0 is odd, up should be odd
// 	const int down=-198; //E0_gamma-down has to be even number, if E0 is odd, down should be odd
// 	const int Ea=42;
// 	double T0_lifetime = 0;
// 	int Lifetimestep = 1;
// 	double centroid=5155.70;
// 	double Eg_uncertainty = 0.60;
// 	int Yaxisrange=15;
// 	ULbkg=1.07/0.93;
// 	LLbkg=0.8/0.93;








// 	const double E0_gamma=3284; //Ex=3284
// 	const int up=336;
// 	const int down=-70;
// 	const int Ea=45;
// 	double T0_lifetime = 88;
// 	int Lifetimestep = 1;
// 	double centroid=E0_gamma;
// 	int Yaxisrange=15;

// 	const double E0_gamma=2035; //Ex=3284
// 	const int up=232;
// 	const int down=30;
// 	const int Ea=45;
// 	double T0_lifetime = 88;
// 	int Lifetimestep = 1;
// 	double centroid=E0_gamma;
// 	int Yaxisrange=15;

// 	const double E0_gamma=1050; //Ex=3284
// 	const int up=162;
// 	const int down=30;
// 	const int Ea=45;
// 	double T0_lifetime = 88;
// 	int Lifetimestep = 1;
// 	double centroid=E0_gamma;
// 	int Yaxisrange=15;

// 	const double E0_gamma=2101; //Ex=3349
// 	const int up=232;
// 	const int down=-54;
// 	const int Ea=45;
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
// 	double centroid=3434.70;
// 	double Eg_uncertainty = 0.32;
// 	int Yaxisrange=15;
// 	ULbkg=1.3/1.1;
// 	LLbkg=0.9/1.1;

// 	const double E0_gamma=2186; //Ex=3435
// 	const int startpeak=2314; //has to be even number
// 	const int endpeak=2364; //has to be even number
// 	const int up=242; //E0_gamma+up has to be even number, if E0 is odd, up should be odd
// 	const int down=-40; //E0_gamma-down has to be even number, if E0 is odd, down should be odd
// 	const int Ea=45;
// 	double T0_lifetime = 0;
// 	int Lifetimestep = 1;
// 	double centroid=2186.33;
// 	double Eg_uncertainty = 0.33;
// 	int Yaxisrange=15;
// 	ULbkg=5.38/4.6;
// 	LLbkg=3.9/4.6;

// 	const double E0_gamma = 2838; //Ex=4085
// 	const int startpeak = 2998; //has to be even number
// 	const int endpeak = 3070; //has to be even number
// 	const int up = 316; //E0_gamma+up has to be even number, if E0 is odd, up should be odd
// 	const int down = -66; //E0_gamma-down has to be even number, if E0 is odd, down should be odd
// 	const int Ea = 44;
// 	double T0_lifetime = 0;
// 	int Lifetimestep = 1;
// 	double centroid = 2837.60;
// 	double Eg_uncertainty = 0.32;
// 	int Yaxisrange = 15;
// 	ULbkg = 3.3/2.8;
// 	LLbkg = 2.3/2.8;

// 	const double E0_gamma=1852; //Ex=4085
// 	const int up=220;
// 	const int down=-34;
// 	const int Ea=44;
// 	double T0_lifetime = 0;
// 	int Lifetimestep = 1;
// 	double centroid=E0_gamma;
// 	int Yaxisrange=15;

// 	const double E0_gamma=4208; //Ex=4208
// 	const int up=420;
// 	const int down=-158;
// 	const int Ea=43;
// 	double T0_lifetime = 0;
// 	int Lifetimestep = 1;
// 	double centroid=E0_gamma;
// 	int Yaxisrange=15;

// 	const double E0_gamma=2959; //Ex=4208
// 	const int up=316;
// 	const int down=-90;
// 	const int Ea=43;
// 	double T0_lifetime = 0;
// 	int Lifetimestep = 1;
// 	double centroid=E0_gamma;
// 	int Yaxisrange=15;

// 	const double E0_gamma=1165; //Ex=4450
// 	const int up=172;
// 	const int down=90;
// 	const int Ea=43;
// 	double T0_lifetime = 415;
// 	int Lifetimestep = 1;
// 	double centroid=E0_gamma;
// 	int Yaxisrange=25;

// 	const double E0_gamma=4519; //Ex=4519
// 	const int up=440;
// 	const int down=-178;
// 	const int Ea=43;
// 	double T0_lifetime = 0;
// 	int Lifetimestep = 1;
// 	double centroid=E0_gamma;
// 	int Yaxisrange=15;

// 	const double E0_gamma=3334; //Ex=4583
// 	const int up=346;
// 	const int down=-90;
// 	const int Ea=43;
// 	double T0_lifetime = 0;
// 	int Lifetimestep = 1;
// 	double centroid=E0_gamma;
// 	int Yaxisrange=20;

// 	const double E0_gamma=1299; //Ex=4583
// 	const int up=172;
// 	const int down=-20;
// 	const int Ea=43;
// 	double T0_lifetime = 0;
// 	int Lifetimestep = 1;
// 	double centroid=E0_gamma;
// 	int Yaxisrange=15;

// 	const double E0_gamma=1232; //Ex=4583
// 	const int up=172;
// 	const int down=10;
// 	const int Ea=43;
// 	double T0_lifetime = 68;
// 	int Lifetimestep = 1;
// 	double centroid=E0_gamma;
// 	int Yaxisrange=15;

// 	const double E0_gamma=4717; //Ex=4717
// 	const int up=470;
// 	const int down=-198;
// 	const int Ea=42;
// 	double T0_lifetime = 0;
// 	int Lifetimestep = 1;
// 	double centroid=E0_gamma;
// 	int Yaxisrange=15;

// 	const double E0_gamma=2484; //Ex=4717
// 	const int up=272;
// 	const int down=-74;
// 	const int Ea=42;
// 	double T0_lifetime = 0;
// 	int Lifetimestep = 1;
// 	double centroid=E0_gamma;
// 	int Yaxisrange=15;

// 	const double E0_gamma=1434; //Ex=4717
// 	const int up=170;
// 	const int down=-14;
// 	const int Ea=42;
// 	double T0_lifetime = 0;
// 	int Lifetimestep = 1;
// 	double centroid=E0_gamma;
// 	int Yaxisrange=15;

// 	const double E0_gamma=4866; //Ex=4866
// 	const int up=468;
// 	const int down=-190;
// 	const int Ea=42;
// 	double T0_lifetime = 0;
// 	int Lifetimestep = 1;
// 	double centroid=E0_gamma;
// 	int Yaxisrange=15;

// 	const double E0_gamma=3617; //Ex=4866
// 	const int up=366;
// 	const int down=-130;
// 	const int Ea=42;
// 	double T0_lifetime = 0;
// 	int Lifetimestep = 1;
// 	double centroid=E0_gamma;
// 	int Yaxisrange=15;
	
// 	const double E0_gamma=3773; //Ex=5022
// 	const int up=384;
// 	const int down=-130;
// 	const int Ea=42;
// 	double T0_lifetime = 0;
// 	int Lifetimestep = 1;
// 	double centroid=E0_gamma;
// 	int Yaxisrange=15;

// 	const double E0_gamma=1738; //Ex=5022
// 	const int up=210;
// 	const int down=-34;
// 	const int Ea=42;
// 	double T0_lifetime = 0;
// 	int Lifetimestep = 1;
// 	double centroid=E0_gamma;
// 	int Yaxisrange=15;

// 	const double E0_gamma=1672; //Ex=5022
// 	const int up=210;
// 	const int down=-34;
// 	const int Ea=42;
// 	double T0_lifetime = 0;
// 	int Lifetimestep = 1;
// 	double centroid=E0_gamma;
// 	int Yaxisrange=15;

// 	const double E0_gamma=1950; //Ex=5302
// 	const int up=230;
// 	const int down=-44;
// 	const int Ea=41;
// 	double T0_lifetime = 0;
// 	int Lifetimestep = 1;
// 	double centroid=E0_gamma;
// 	int Yaxisrange=15;

// 	const double E0_gamma=2050; //Ex=5402
// 	const int up=230;
// 	const int down=-44;
// 	const int Ea=41;
// 	double T0_lifetime = 0;
// 	int Lifetimestep = 1;
// 	double centroid=E0_gamma;
// 	int Yaxisrange=15;

// 	const double E0_gamma=5435; //Ex=5435
// 	const int up=520;
// 	const int down=-240;
// 	const int Ea=41;
// 	double T0_lifetime = 0;
// 	int Lifetimestep = 1;
// 	double centroid=E0_gamma;
// 	int Yaxisrange=15;
	
// 	const double E0_gamma=4187; //Ex=5435
// 	const int up=440;
// 	const int down=-160;
// 	const int Ea=41;
// 	double T0_lifetime = 0;
// 	int Lifetimestep = 1;
// 	double centroid=E0_gamma;
// 	int Yaxisrange=15;
	
// 	const double E0_gamma=2167; //Ex=5518
// 	const int up=242;
// 	const int down=-84;
// 	const int Ea=41;
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
// 	double centroid = 4269.50;
// 	double Eg_uncertainty = 0.30;
// 	int Yaxisrange = 15;
// 	ULbkg = 2.9/2.55;
// 	LLbkg = 2.2/2.55;
	
// 	const double E0_gamma=2325; //Ex=5676
// 	const int up=242;
// 	const int down=-84;
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
// 	double centroid = 3541.10;
// 	double Eg_uncertainty = 0.27;
// 	int Yaxisrange = 15;
// 	ULbkg = 7.2/6.3;
// 	LLbkg = 5.6/6.3;

// 	const double E0_gamma=5824; //Ex=5824
// 	const int up=550;
// 	const int down=-280;
// 	const int Ea=40;
// 	double T0_lifetime = 0;
// 	int Lifetimestep = 1;
// 	double centroid=E0_gamma;
// 	int Yaxisrange=15;

// 	const double E0_gamma=4576; //Ex=5824
// 	const int up=442;
// 	const int down=-184;
// 	const int Ea=40;
// 	double T0_lifetime = 0;
// 	int Lifetimestep = 1;
// 	double centroid=E0_gamma;
// 	int Yaxisrange=15;

// 	const double E0_gamma=3590; //Ex=5824
// 	const int up=366;
// 	const int down=-130;
// 	const int Ea=40;
// 	double T0_lifetime = 0;
// 	int Lifetimestep = 1;
// 	double centroid=E0_gamma;
// 	int Yaxisrange=15;

// 	const double E0_gamma=2540; //Ex=5824
// 	const int up=262;
// 	const int down=-104;
// 	const int Ea=40;
// 	double T0_lifetime = 0;
// 	int Lifetimestep = 1;
// 	double centroid=E0_gamma;
// 	int Yaxisrange=15;

// 	const double E0_gamma=5890; //Ex=5890
// 	const int up=550;
// 	const int down=-280;
// 	const int Ea=40;
// 	double T0_lifetime = 0;
// 	int Lifetimestep = 1;
// 	double centroid=E0_gamma;
// 	int Yaxisrange=15;

// 	const double E0_gamma=3656; //Ex=5890
// 	const int up=440;
// 	const int down=-110;
// 	const int Ea=40;
// 	double T0_lifetime = 0;
// 	int Lifetimestep = 1;
// 	double centroid=E0_gamma;
// 	int Yaxisrange=15;

// 	const double E0_gamma=2606; //Ex=5890
// 	const int up=280;
// 	const int down=-100;
// 	const int Ea=40;
// 	double T0_lifetime = 0;
// 	int Lifetimestep = 1;
// 	double centroid=E0_gamma;
// 	int Yaxisrange=30;

// 	const double E0_gamma=1394; //Ex=5977
// 	const int up=156;
// 	const int down=-4;
// 	const int Ea=40;
// 	double T0_lifetime = 0;
// 	int Lifetimestep = 1;
// 	double centroid=E0_gamma;
// 	int Yaxisrange=15;

// 	const double E0_gamma=2780; //Ex=6129
// 	const int up=280;
// 	const int down=-90;
// 	const int Ea=40;
// 	double T0_lifetime = 0;
// 	int Lifetimestep = 1;
// 	double centroid=E0_gamma;
// 	int Yaxisrange=30;

// 	const double E0_gamma=4890; //Ex=6138
// 	const int up=500;
// 	const int down=-210;
// 	const int Ea=40;
// 	double T0_lifetime = 0;
// 	int Lifetimestep = 1;
// 	double centroid=E0_gamma;
// 	int Yaxisrange=15;

// 	const double E0_gamma=2786; //Ex=6138
// 	const int up=290;
// 	const int down=-90;
// 	const int Ea=40;
// 	double T0_lifetime = 0;
// 	int Lifetimestep = 1;
// 	double centroid=E0_gamma;
// 	int Yaxisrange=30;

// 	const double E0_gamma=3922; //Ex=6158
// 	const int up=440;
// 	const int down=-110;
// 	const int Ea=40;
// 	double T0_lifetime = 0;
// 	int Lifetimestep = 1;
// 	double centroid=E0_gamma;
// 	int Yaxisrange=15;

// 	const double E0_gamma=2874; //Ex=6158
// 	const int up=290;
// 	const int down=-90;
// 	const int Ea=40;
// 	double T0_lifetime = 0;
// 	int Lifetimestep = 1;
// 	double centroid=E0_gamma;
// 	int Yaxisrange=30;

// 	const double E0_gamma=1708; //Ex=6158
// 	const int up=210;
// 	const int down=-34;
// 	const int Ea=40;
// 	double T0_lifetime = 27;
// 	int Lifetimestep = 1;
// 	double centroid=E0_gamma;
// 	int Yaxisrange=15;

// 	const double E0_gamma=6255; //Ex=6255
// 	const int up=580;
// 	const int down=-320;
// 	const int Ea=39;
// 	double T0_lifetime = 0;
// 	int Lifetimestep = 1;
// 	double centroid=E0_gamma;
// 	int Yaxisrange=15;

// 	const double E0_gamma=2971; //Ex=6255
// 	const int up=316;
// 	const int down=-100;
// 	const int Ea=39;
// 	double T0_lifetime = 0;
// 	int Lifetimestep = 1;
// 	double centroid=E0_gamma;
// 	int Yaxisrange=15;

// 	const double E0_gamma=6268; //Ex=6268
// 	const int up=580;
// 	const int down=-320;
// 	const int Ea=39;
// 	double T0_lifetime = 0;
// 	int Lifetimestep = 1;
// 	double centroid=E0_gamma;
// 	int Yaxisrange=15;

// 	const double E0_gamma=6279; //Ex=6279
// 	const int up=580;
// 	const int down=-320;
// 	const int Ea=39;
// 	double T0_lifetime = 0;
// 	int Lifetimestep = 1;
// 	double centroid=E0_gamma;
// 	int Yaxisrange=15;

// 	const double E0_gamma=5030; //Ex=6279
// 	const int up=500;
// 	const int down=-210;
// 	const int Ea=39;
// 	double T0_lifetime = 0;
// 	int Lifetimestep = 1;
// 	double centroid=E0_gamma;
// 	int Yaxisrange=15;

// 	const double E0_gamma=4045; //Ex=6279
// 	const int up=420;
// 	const int down=-150;
// 	const int Ea=39;
// 	double T0_lifetime = 0;
// 	int Lifetimestep = 1;
// 	double centroid=E0_gamma;
// 	int Yaxisrange=15;

// 	const double E0_gamma=2995; //Ex=6279
// 	const int up=316;
// 	const int down=-100;
// 	const int Ea=39;
// 	double T0_lifetime = 0;
// 	int Lifetimestep = 1;
// 	double centroid=E0_gamma;
// 	int Yaxisrange=20;

// 	const double E0_gamma=6348; //Ex=6348
// 	const int up=580;
// 	const int down=-320;
// 	const int Ea=39;
// 	double T0_lifetime = 0;
// 	int Lifetimestep = 1;
// 	double centroid=E0_gamma;
// 	int Yaxisrange=15;

// 	const double E0_gamma=5100; //Ex=6348
// 	const int up=500;
// 	const int down=-210;
// 	const int Ea=39;
// 	double T0_lifetime = 0;
// 	int Lifetimestep = 1;
// 	double centroid=E0_gamma;
// 	int Yaxisrange=15;
	
// 	const double E0_gamma=4114; //Ex=6348
// 	const int up=420;
// 	const int down=-178;
// 	const int Ea=39;
// 	double T0_lifetime = 0;
// 	int Lifetimestep = 1;
// 	double centroid=E0_gamma;
// 	int Yaxisrange=15;

// 	const double E0_gamma=6390; //Ex=6390
// 	const int up=580;
// 	const int down=-320;
// 	const int Ea=39;
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
// 	double centroid=5141.30;
// 	double Eg_uncertainty = 0.60;
// 	int Yaxisrange=15;
// 	ULbkg=2.4/2.1;
// 	LLbkg=1.8/2.1;

// 	const double E0_gamma=4156; //Ex=6390
// 	const int startpeak=4404; //has to be even number
// 	const int endpeak=4494; //has to be even number
// 	const int up=450; //E0_gamma+up has to be even number, if E0 is odd, up should be odd
// 	const int down=-110; //E0_gamma-down has to be even number, if E0 is odd, down should be odd
// 	const int Ea=39;
// 	double T0_lifetime = 0;
// 	int Lifetimestep = 1;
// 	double centroid=4155.84;
// 	double Eg_uncertainty = 0.31;
// 	int Yaxisrange=15;
// 	ULbkg=3.9/3.5;
// 	LLbkg=3.1/3.5;

// 	const double E0_gamma=3314; //Ex=6390
// 	const int up=330;
// 	const int down=-120;
// 	const int Ea=39;
// 	double T0_lifetime = 0;
// 	int Lifetimestep = 1;
// 	double centroid=E0_gamma;
// 	int Yaxisrange=15;

// 	const double E0_gamma=3106; //Ex=6390
// 	const int up=316;
// 	const int down=-100;
// 	const int Ea=39;
// 	double T0_lifetime = 0;
// 	int Lifetimestep = 1;
// 	double centroid=E0_gamma;
// 	int Yaxisrange=15;

// 	const double E0_gamma=2183; //Ex=6390
// 	const int up=276;
// 	const int down=-40;
// 	const int Ea=39;
// 	double T0_lifetime = 0;
// 	int Lifetimestep = 1;
// 	double centroid=E0_gamma;
// 	int Yaxisrange=15;

// 	const double E0_gamma=3043; //Ex=6394
// 	const int up=320;
// 	const int down=-130;
// 	const int Ea=39;
// 	double T0_lifetime = 0;
// 	int Lifetimestep = 1;
// 	double centroid=E0_gamma;
// 	int Yaxisrange=15;

// 	const double E0_gamma=1091; //Ex=6394
// 	const int up=150;
// 	const int down=-12;
// 	const int Ea=39;
// 	double T0_lifetime = 0;
// 	int Lifetimestep = 1;
// 	double centroid=E0_gamma;
// 	int Yaxisrange=15;

// 	const double E0_gamma=6402; //Ex=6402
// 	const int up=580;
// 	const int down=-320;
// 	const int Ea=39;
// 	double T0_lifetime = 0;
// 	int Lifetimestep = 1;
// 	double centroid=E0_gamma;
// 	int Yaxisrange=15;

// 	const double E0_gamma=5154; //Ex=6402
// 	const int up=500;
// 	const int down=-210;
// 	const int Ea=39;
// 	double T0_lifetime = 0;
// 	int Lifetimestep = 1;
// 	double centroid=E0_gamma;
// 	int Yaxisrange=15;

// 	const double E0_gamma=4168; //Ex=6402
// 	const int up=420;
// 	const int down=-178;
// 	const int Ea=39;
// 	double T0_lifetime = 0;
// 	int Lifetimestep = 1;
// 	double centroid=E0_gamma;
// 	int Yaxisrange=15;

// 	const double E0_gamma=3326; //Ex=6402
// 	const int up=330;
// 	const int down=-120;
// 	const int Ea=39;
// 	double T0_lifetime = 0;
// 	int Lifetimestep = 1;
// 	double centroid=E0_gamma;
// 	int Yaxisrange=15;

// 	const double E0_gamma=3118; //Ex=6402
// 	const int up=316;
// 	const int down=-100;
// 	const int Ea=39;
// 	double T0_lifetime = 0;
// 	int Lifetimestep = 1;
// 	double centroid=E0_gamma;
// 	int Yaxisrange=15;

// 	const double E0_gamma=3052; //Ex=6402
// 	const int up=320;
// 	const int down=-120;
// 	const int Ea=39;
// 	double T0_lifetime = 0;
// 	int Lifetimestep = 1;
// 	double centroid=E0_gamma;
// 	int Yaxisrange=15;

// 	const double E0_gamma=6540; //Ex=6542
// 	const int up=580;
// 	const int down=-320;
// 	const int Ea=39;
// 	double T0_lifetime = 0;
// 	int Lifetimestep = 1;
// 	double centroid=E0_gamma;
// 	int Yaxisrange=15;

// 	const double E0_gamma=5293; //Ex=6542
// 	const int startpeak=5600; //has to be even number
// 	const int endpeak=5720; //has to be even number
// 	const int up=561; //E0_gamma+up has to be even number, if E0 is odd, up should be odd
// 	const int down=-151; //E0_gamma-down has to be even number, if E0 is odd, down should be odd
// 	const int Ea=39;
// 	double T0_lifetime = 0;
// 	int Lifetimestep = 1;
// 	double centroid=5292.50;
// 	double Eg_uncertainty = 0.40;
// 	int Yaxisrange=15;
// 	ULbkg=2.0/1.75;
// 	LLbkg=1.5/1.75;

// 	const double E0_gamma=4306; //Ex=6542
// 	const int up=420;
// 	const int down=-178;
// 	const int Ea=39;
// 	double T0_lifetime = 0;
// 	int Lifetimestep = 1;
// 	double centroid=E0_gamma;
// 	int Yaxisrange=15;

// 	const double E0_gamma=3464; //Ex=6542
// 	const int up=340;
// 	const int down=-130;
// 	const int Ea=39;
// 	double T0_lifetime = 0;
// 	int Lifetimestep = 1;
// 	double centroid=E0_gamma;
// 	int Yaxisrange=15;

// 	const double E0_gamma=3256; //Ex=6542
// 	const int up=340;
// 	const int down=-130;
// 	const int Ea=39;
// 	double T0_lifetime = 0;
// 	int Lifetimestep = 1;
// 	double centroid=E0_gamma;
// 	int Yaxisrange=15;

// 	const double E0_gamma=3190; //Ex=6542
// 	const int up=330;
// 	const int down=-120;
// 	const int Ea=39;
// 	double T0_lifetime = 0;
// 	int Lifetimestep = 1;
// 	double centroid=E0_gamma;
// 	int Yaxisrange=15;

//	const double E0_gamma=2000; //DSLmodify //E0+i*binwidth, a baseline E0 determined by the file name of the simuroot.root, 2752 4236 2868 1366, decide which simuroot to use// if you want to use for() loop and E0_gamma is related to iscanE, E0_gamma cannot be global variable out of the main(), as loop has to be in the main().
	//E0-Ebin to shift the simulated spectrum.h E0-Ebin>0 is moving simulated spectrum towards higher energy
	const double Ebin_gamma=E0_gamma; //modify center for drawing figures: 2754 4238 2870 1369, don't change to other values // if you want to use for() loop and E0_gamma is related to iscanE, E0_gamma cannot be global variable out of the main(), as loop has to be in the main().




	const double minrange=Ebin_gamma-down; //always 100-keV window
	const double maxrange=Ebin_gamma+up; //always 100-keV window
	const double minrangeb=Ebin_gamma-down; //longer range to fit the background correctly
	const double maxrangeb=Ebin_gamma+up; //longer range to fit the background correctly
	const double minrangezoom=Ebin_gamma-down; //shorter range to calculate chi2 25 30 30 40
	const double maxrangezoom=Ebin_gamma+up; //shorter range to calculate chi2 25 30 30 40


	const int Nbins=(maxrange-minrange)/binwidth;//make sure this is an integer // there is another global Nbins at the beginning because fcn() also needs an Nbins.
	const int Nbinsb=(maxrangeb-minrangeb)/binwidth;//make sure this is an integer // there is another global Nbins at the beginning, you'd better keep them consistent
	const int Nbinsz=(maxrangezoom-minrangezoom)/binwidth;//make sure this is an integer // there is another global Nbins at the beginning, you'd better keep them consistent
	const int Nbinspeak=(endpeak-startpeak)/binwidth;//make sure this is an integer // there is another global Nbins at the beginning, you'd better keep them consistent

	Bool_t reject;
//double E0_gamma=2753.0+double(iscanE)*2.0/10.0;                  // Energy of Gamma-ray (keV) useless
const int nBranches=2;                  //13,19;2,3 Number of proton Branches feeding daughter modify
TH1F *proton_branch[nBranches],*proton_copy[nBranches],*maximum_branch[nBranches],*totalmax;
//Lit 1369-keV γ 20
// double E_CoM[nBranches] = {4.556,4.303,4.2583,3.597,3.466,3.342,3.237,3.021,2.486,2.453,2.165,1.272,0.9434,0.55,1.805,1.501,1.04,1.685,1.396,0.000};         //Center of Mass energy (MeV) in descending order
// double feed[nBranches]={1.28/216.80,3.32/216.80,100/216.80,10.86/216.80,34.5/216.80,6.57/216.80,4.15/216.80,3.74/216.80,0.96/216.80,0.4/216.80,17.2/216.80,2.26/216.80,17/216.80,2.5/216.80,6.73/216.80,2.9/216.80,1.53/216.80,0.2/216.80,0.61/216.80,0.09/216.80};		//beta feeding intensity of each proton, if the sum of all the feeding intensities is less than 1, the extra events will assume a En=0.
//Lit 2754-keV γ 4
double E_CoM[nBranches] = {1};         //Center of Mass energy (MeV) in descending order
double feed[nBranches]={1};		//beta feeding intensity of each proton, if the sum of all the feeding intensities is less than 1, the extra events will assume a En=0.
//Lit 2870-keV γ 3
// double E_CoM[nBranches] = {1.685,1.396,0.000};         //Center of Mass energy (MeV) in descending order
// double feed[nBranches]={0.93/(2.89+0.93+0.034),2.89/(2.89+0.93+0.034),0.034/(2.89+0.93+0.034)};		//beta feeding intensity of each proton, if the sum of all the feeding intensities is less than 1, the extra events will assume a En=0.
//Lit 4238-keV γ 3
// double E_CoM[nBranches] = {1.685,1.396,0.000};         //Center of Mass energy (MeV) in descending order
// double feed[nBranches]={0.93/(2.89+0.93+0.039),2.89/(2.89+0.93+0.039),0.039/(2.89+0.93+0.039)};		//beta feeding intensity of each proton, if the sum of all the feeding intensities is less than 1, the extra events will assume a En=0.
//modify
//Thomas 1369-keV γ 14
// double E_CoM[nBranches] = {4.545,4.252,3.610,3.463,3.231,2.98,2.162,1.268,0.943,0.555,1.804,1.489,1.377,0.000};         //Center of Mass energy (MeV) in descending order
// double feed[nBranches]={6.6/208.3,100/208.3,5.9/208.3,28.1/208.3,5.4/208.3,1.7/208.3,18.1/208.3,6.1/208.3,17.1/208.3,7.2/208.3,6.1/208.3,5.0/208.3,0.91/208.3,0.09/208.3};		//beta feeding intensity of each proton, if the sum of all the feeding intensities is less than 1, the extra events will assume a En=0.
//Thomas 2754-keV γ 3
// double E_CoM[nBranches] = {1.804,1.489,0.000};         //Center of Mass energy (MeV) in descending order
// double feed[nBranches]={6.1/(5+6.1+0.519),5.0/(5+6.1+0.519),0.519/(5+6.1+0.519)};		//beta feeding intensity of each proton, if the sum of all the feeding intensities is less than 1, the extra events will assume a En=0.
//Thomas 2870-keV γ 2
// double E_CoM[nBranches] = {1.377,0.000};         //Center of Mass energy (MeV) in descending order
// double feed[nBranches]={4.3/(4.3+0.038),0.038/(4.3+0.038)};		//beta feeding intensity of each proton, if the sum of all the feeding intensities is less than 1, the extra events will assume a En=0.
//Thomas 4238-keV γ 2
// double E_CoM[nBranches] = {1.377,0.000};         //Center of Mass energy (MeV) in descending order
// double feed[nBranches]={4.3/(4.3+0.043),0.043/(4.3+0.043)};		//beta feeding intensity of each proton, if the sum of all the feeding intensities is less than 1, the extra events will assume a En=0.

// int Flag_SP = 1;                         //stopping power default=1
// int Flag_Det = 1;                      //Detector resolution. No f1 sampling when Flag_Det is false, will be very fast. default=1
// int Flag_Ang = 1;                       //Tests random cos vs random between -1 and 1 for angular dist (false tests cos and true tests random X(-1,1) default=1
// int Flag_Escape = 0;					//If take recoil escaping into consideration default=0
// int Flag_sig = 0;					//If take σ error into consideration default=0, -1 or 1, minus err or plus err
// int Flag_tau = 0;					//If take τ error into consideration default=0, -1 or 1, minus err or plus err
// int Flag_st = 0;					//If take stopping power error into consideration default=0, -1 or 1, minus err or plus err
// int NumberofSteps=200;	//Divide the entire lifetime into several steps
double Chi2=0;//Chi2 is global varible, use Chi2 to extract the minimized Chi2 obtained in fcn(), and output Chi2 in main().
long double LBayesian=1;
//	TCanvas *c1;

//void readHists()
//generate example histograms
//{

//hSega->GetXaxis()->SetRangeUser(E0_gamma-50,E0_gamma+50);
//background->GetXaxis()->SetRangeUser(minrange,maxrange);
// 	int minim=0; //Only change this if you need to shift the histogram some amount of bins. Wouldn't recommend for bins larger than 0.1 keV
// 	if(minim>=0)
// 	{
// 		for(int i=0;i<(p1->GetNbinsX()-minim);i++)
// 		{
// 			p1->SetBinContent(p1->GetNbinsX()-i,p1->GetBinContent(p1->GetNbinsX()-i-minim));
// 			p2->SetBinContent(p1->GetNbinsX()-i,p2->GetBinContent(p1->GetNbinsX()-i-minim));
// 			p3->SetBinContent(p1->GetNbinsX()-i,p3->GetBinContent(p1->GetNbinsX()-i-minim));
// 			//p4->SetBinContent(p1->GetNbinsX()-i,p4->GetBinContent(p1->GetNbinsX()-i-minim));
// 			//p5->SetBinContent(p1->GetNbinsX()-i,p5->GetBinContent(p1->GetNbinsX()-i-minim));
// 			//p6->SetBinContent(p1->GetNbinsX()-i,p6->GetBinContent(p1->GetNbinsX()-i-minim));
// 		}
// 	}
// 	else
// 	{
// 		for(int i=0;i<(p1->GetNbinsX()+minim);i++)
// 		{
// 			p1->SetBinContent(i,p1->GetBinContent(i-minim));
// 			p2->SetBinContent(i,p2->GetBinContent(i-minim));
// 			p3->SetBinContent(i,p3->GetBinContent(i-minim));
// 			//p4->SetBinContent(i,p4->GetBinContent(i-minim));
// 		}
// 	}

//}

void Comparison_DSL_test1248_twopeaktau()
{
	double flag_Eg, flag_bkg, flag_SP, flag_AC;
	
	char b_name[200];//3*3*3*3*16=1296
	//1 Tau * 3 Eg * 3 SP * 3 AC = 27 output files
	//1 Tau * 3 Eg * 3 bkg * 3 SP * 3 AC = 81 output files
	//0-15 fs: 1296 runs, Eg 0-431, Eg 432-863, Eg 864-1295

	//2186, 3435
	//2 fs step, 0-30 fs: 1296 runs, Eg 0-431, Eg 432-863, Eg 864-1295

//	int ifolder = 0;	for (flag_Eg = centroid - Eg_uncertainty; flag_Eg <= centroid + 0; flag_Eg += Eg_uncertainty)	// (1)-(2) first step 2186, 3435
//	int ifolder = 864; flag_Eg = centroid + Eg_uncertainty; //(3) // second step 2186, 3435

	//3076, 4970, 5156, 4156, 5141, 5293
	//1 fs step, 0-20 fs: 1701 runs, Eg- 0-566, Eg 567-1133, Eg+ 1134-1700

	//2234
	//20 fs step, 0-600 fs: 2511 runs, Eg- 0-836, Eg 837-1673, Eg+ 1674-2510
	//20 fs step, 0-500 fs: 2106 runs, Eg- 0-701, Eg 702-1403, Eg+ 1404-2105

//		supplement 30,40,50,60 fs: start from 1701
//  	int ifolder = 1701; // 1701-1781 for 30fs, 1782-1862 for 40fs, 1863-1943 for 50 fs, 1944-2024 for 60 fs supplemental step 4156

//	for (flag_Eg = centroid - Eg_uncertainty; flag_Eg <= centroid + Eg_uncertainty; flag_Eg += Eg_uncertainty)
//		supplemental step 4156
//	int ifolder = 864; flag_Eg = centroid + Eg_uncertainty; //(3) // original step 4156	

//	Eg1248, 200 fs step, 600-4000 fs: (4000-600)/200=17,+1; 4200-6000 fs: (6000-4200)/200=9,+1
//	1458 runs, Eg- 0-485, Eg 486-971, Eg+ 972-1457, supplemental Eg 10*27*3=810 runs,  1458-2267, 1458-1727, 1728-1997, 1998-2267
//  int ifolder = 0;
//  for (flag_Eg = centroid - Eg_uncertainty; flag_Eg <= centroid + 0; flag_Eg += Eg_uncertainty)
// 	int ifolder = 972;
// 	for (flag_Eg = centroid + Eg_uncertainty; flag_Eg <= centroid + Eg_uncertainty; flag_Eg += Eg_uncertainty)
//  int ifolder = 1458;
//  for (flag_Eg = centroid - Eg_uncertainty; flag_Eg <= centroid + Eg_uncertainty; flag_Eg += Eg_uncertainty)

//	choose one flag_Eg, this script can't do more than 989 files once
//	int ifolder = 0; flag_Eg = centroid - Eg_uncertainty; //(1)  // original step 3076, 4156
//	int ifolder = 567; flag_Eg = centroid; //(2) // original step 4156
//	int ifolder = 1134; flag_Eg = centroid + Eg_uncertainty; //(3) // original step 4156

//	int ifolder = 0; flag_Eg = centroid - Eg_uncertainty; //(1)  // original step 2234
//	int ifolder = 702; flag_Eg = centroid; //(2) // original step 2234
//	int ifolder = 1404; flag_Eg = centroid + Eg_uncertainty; //(3) // original step 2234
	
//	int ifolder = 0; // original step 2838, 4270, 3541. 0-647, 8*27*3=648
// 	int ifolder = 648; // supplemental step 2838, 4270, 3541. 648-972, 4*27*3=324
// 	for (flag_Eg = centroid - Eg_uncertainty; flag_Eg <= centroid + Eg_uncertainty; flag_Eg += Eg_uncertainty)
	
//	int ifolder = 2268; flag_Eg = centroid - Eg_uncertainty*5*2; //1248 supplement 1246.40
//	int ifolder = 3024; flag_Eg = centroid - Eg_uncertainty*5; //1248 supplement 1247.40
	int ifolder =0; flag_Eg = centroid;
	{
		for (flag_bkg = 1.0; flag_bkg <= 1.0; flag_bkg += 0.1)
		{
			for (flag_SP = 1.0; flag_SP <= 1.0; flag_SP += 0.1)
			{
				for (flag_AC = 0.0; flag_AC <= 0.0; flag_AC += 1.0)
				{		

//  			FILE *outfile=fopen ("D:/X/out/Compare.dat","a"); //ofstream doesn't work. output test
//  			FILE *outfile2=fopen ("D:/X/out/Compare2.dat","a"); //ofstream doesn't work. output test
// 			sprintf(b_name,"%s%d%s%d%s%.0f%s%s%s","D:/X/out/DSL/Chi2/Chi2_Gamma",E0_gamma,"_Ea",Ea,"_Egv",centroid,"slice_vTau_",flag_bkg,".dat");
// 			FILE *outfile=fopen (b_name,"w"); //ofstream doesn't work
// 			sprintf(b_name,"%s%d%s%d%s%.0f%s%s%s","D:/X/out/DSL/LBayesian/LBayesian_Gamma",E0_gamma,"_Ea",Ea,"_Egv",centroid,"slice_vTau_",flag_bkg,".dat");
// 			FILE *outfile2=fopen (b_name,"w"); //ofstream doesn't work
					//DSLportal

// 	for(int iscantau=0;iscantau<=30;iscantau+=2)// all step 2186, 3435
//	for(int iscantau=0;iscantau<=20;iscantau+=1) // original step 3076, 4156
//	for(int iscantau=30;iscantau<=60;iscantau+=10) // supplemental step 4156
//	for(int iscantau=0;iscantau<=500;iscantau+=20) // original step 2234 (no supplemental step)
//	for(int iscantau=0;iscantau<=0;iscantau+=200) // original step 1248
		int iscantau = 5000;
		int plus_tau = 0;
// 		flag_SP = 0.5;
// 		int second_SP = 2.0;
//	for(int iscantau=4200;iscantau<=6000;iscantau+=200) // supplemental step 1248
//	for(int iscantau=0;iscantau<=30;iscantau+=4)// original step 2838, 4270, 3541
//	for(int iscantau=30;iscantau<=60;iscantau+=10)// supplemental step 2838, 4270, 3541
	{	
		for(int iscanE=0;iscanE<=0;iscanE++)//No need to change! iscanE<0 means shift the simulated spectrum towards low energy
			//(-8,9) 17 energy points
		{//shift simulated spectrum with different gamma energies 1369 i=34; 2754 i=29; 2870 i=30; 4238 i=27
			//iscanE+1 bin shift +1, energy +2 keV
			double fitrangelow1, fitrangelow2, fitrangehigh1, fitrangehigh2;
			fitrangelow1=minrange; fitrangelow2=startpeak; fitrangehigh1=endpeak; fitrangehigh2=maxrange;


			//readHists();
			Chi2=0;//Chi2 is global varible, use Chi2 to extract the minimized Chi2 obtained in fcn(), and output Chi2 in main(). Must clean Chi2 before each run in main().
			char simurootname[400], outrootname[400], outtxtname[400], outfigname[400], tname[400],hname[400], ename[400], pname[400];

			sprintf(pname,"%d_",(int)((E_CoM[0]+0.000001)*1000));
			for(int ii=1;ii<nBranches;ii++)
			{
				sprintf(ename,"%d_",(int)((E_CoM[ii]+0.000001)*1000));
				strcat(pname,ename);//for the output root file name
				//cout<<pname<<endl;
			}
			//sprintf(tname,"%s%.2f%s%d%s%.1f%s%d%s%s%s%d%s%d%s%d%s%d%s%d%s%d%s%.1f%s","D:/X/out/Thomas/Sim25Si_T",T0_lifetime,"_Step",NumberofSteps,"_Eg",E0_gamma,"_feed",nBranches,"_En",pname,"SP",Flag_SP,"_Det",Flag_Det,"_Ang",Flag_Ang,"_st",Flag_st,"_sig",Flag_sig,"_tau",Flag_tau,"_binwid",0.1,"_P10_Sun_Thomas");
			//		sprintf(tname,"%s%.0f%s%.1f","D:/X/out/DSL/S31_Eg",E0_gamma,"_tau",T0_lifetime+iscantau*Lifetimestep);//read in lots of simulation files for 0.5 fs file
			sprintf(tname,"%s%.0f%s%.0f%s%.2f%s%.1f%s%.2f%s%.1f","D:/X/out/DSL/Eg",E0_gamma,"/S31_Gamma",E0_gamma,"_Eg",flag_Eg,"_Tau",T0_lifetime+iscantau*Lifetimestep,"_SP",flag_SP,"_AC",flag_AC);
			sprintf(simurootname,"%s%s",tname,".root"); //read in lots of simulation root files from Oliver's G4 simulation. Input root file names don't contain flag_bkg, but output fig/dat/root file names contain flag_bkg.
			fin_simu0 = new TFile(simurootname); //Rootfile with Simulation histograms

			sprintf(tname,"%s%.0f%s%.0f%s%.2f%s%.1f%s%.2f%s%.1f","D:/X/out/DSL/Eg",E0_gamma,"/S31_Gamma",E0_gamma,"_Eg",flag_Eg,"_Tau",T0_lifetime+iscantau*Lifetimestep,"_SP",flag_SP,"_AC",flag_AC);
			sprintf(simurootname,"%s%s",tname,"_agate38-40.root"); //read in lots of simulation root files from Oliver's G4 simulation. Input root file names don't contain flag_bkg, but output fig/dat/root file names contain flag_bkg.
			fin_simu1 = new TFile(simurootname); //Rootfile with Simulation histograms

			fin_data = new TFile("D:/X/out/DSL/nice.root"); //Rootfile with Experimental Data Histogram modify
			//fin_data = new TFile("D:/X/out/Si25/100eV-bin-gated_sun_25Si.root"); //Rootfile with Experimental Data Histogram gated

			//Get Histograms from Simulation
			sprintf(hname,"Eg");
			proton_branch[0]=(TH1F*)fin_simu0->Get(hname);
			proton_branch[0]->Rebin(factor_rebin);//Rebin simulated histograms from 1-keV bin to 2-keV bin because data are 2-keV bin
//			sprintf(hname,"Eg");
//			maximum_branch[0]=(TH1F*)fin_simu0->Get(hname);

			sprintf(hname,"Eg");
			proton_branch[1]=(TH1F*)fin_simu1->Get(hname);
			proton_branch[1]->Rebin(factor_rebin);//Rebin simulated histograms from 1-keV bin to 2-keV bin because data are 2-keV bin
//			sprintf(hname,"Eg");
//			maximum_branch[1]=(TH1F*)fin_simu1->Get(hname);

// 			for(int ii=1;ii<nBranches;ii++)
// 			{
// 				maximum_branch[0]->Add(maximum_branch[0],maximum_branch[ii]);
// 			}
//			centroid=maximum_branch[0]->GetMean();//what the original energy for a unbroadened gamma peak would be, mean energy of the simulated gamma peak, useless for DSL
			// 	sprintf(hname,"TotalGamma_AllneutronFeedings");
			// 	totalmax=(TH1F*)fin_simu->Get(hname);
			// 	totalmax->Rebin(10);
			// 	int whichbin=totalmax->GetMaximumBin();
			// 	centroid=totalmax->GetBinCenter(whichbin);
// 			sprintf(hname,"%s%d%s","hGriffinAddback_alpha",Ea,"_t0");//alpha gate
// 			hSega_data= (TH1F*)fin_data->Get(hname); //Get Gamma ray spectrum
			//	hSega_data->Rebin(factor_rebin);//Rebin data histograms already 2-keV bin, no need to rebin
			hSega_data = new TH1F("hSega_data","hSega_data",4000,0,8000);
			for(int ii=0;ii<=12;ii+=2) // ii = 0, 2, 4, 6, 8, 10, 12; Ea = 49, 47, 45, 43, 41, 39, 37
			{
				sprintf(hname,"%s%d%s","hGriffinAddback_alpha",Ea-ii,"_t0");//alpha gate
				hSega_data_one_alpha = (TH1F*)fin_data->Get(hname); //Get Gamma ray spectrum
				hSega_data->Add(hSega_data_one_alpha,1);
			}

			hSega=new TH1F("hSega","hSega",Nbins,minrange,maxrange);
			bkgshort=new TH1F("bkgshort","bkgshort",Nbins,minrange,maxrange);
			background=new TH1F("background","background",Nbinsb,minrangeb,maxrangeb);

			for(int ii=0;ii<nBranches;ii++)
			{
				sprintf(hname,"%s%d","Egcopy",int(ii));
				proton_copy[ii]=new TH1F(hname,hname,Nbins,minrange,maxrange);//for outroot
			}

//*****************************************************************************************************
			//bkg algorithm created by Brent
/*
			TF1 *f1=new TF1("f1","[0]+x*[1]",fitrangelow1,fitrangelow2);//range of left half of hist you want to fit the background to //Brent created, Chris doesn't like 1
			TF1 *f2=new TF1("f2","[0]+x*[1]",fitrangehigh1,fitrangehigh2);//range of right half of hist you want to fit the background to //Brent created, Chris doesn't like 2
			hSega_data->Fit("f1","qR");
			hSega_data->Fit("f2","qR");

			for(int i=1;i<=Nbinsb;i++)
			{
				int bkgpower=1;//change this to whatever you want. 3 for high statistics
				// 		TF1 *f1=new TF1("f1","[0]+x*[1]",minrangeb+50,minrangeb+80);//range of left half of hist you want to fit the background to //useless for now
				// 		TF1 *f2=new TF1("f2","[0]+x*[1]",maxrangeb-80,maxrangeb-50);//range of right half of hist you want to fit the background to //useless for now
				double a = (f1->GetParameter(0)+(minrangeb+(double)i/Binsperkev)*f1->GetParameter(1));//Lower BG //Brent created, Chris doesn't like 3
				double c = (f2->GetParameter(0)+(minrangeb+(double)i/Binsperkev)*f2->GetParameter(1));//Upper BG //Brent created, Chris doesn't like 4
				double f=((double)Nbinsb-(double)i)/(double)Nbinsb;
				double g =(double)i/(double)Nbinsb;
				double h =pow(f,bkgpower)+pow(g,bkgpower);//if bkgpower = 1, h = 1
				double b=0;
				//if(minrangeb+(double)i/(double)Binsperkev<minrangeb+80) b = a; //test for fun
				//else if(minrangeb+(double)i/(double)Binsperkev>maxrangeb-80) b = c; //test for fun
				//else b = (a*pow(f,bkgpower)+c*pow(g,bkgpower))/h; //test for fun
				b = (a*pow(f,bkgpower)+c*pow(g,bkgpower))/h; //Brent created, it's good for constant pol0 f1f2 or linear pol1 f1f2
				//b=(a+c)/2; //lower bkg and higher bkg will be the same constant, not very precise
*/
//*****************************************************************************************************
			// or constant bkg algorithm. Chris prefers this for DSL
/*
			TF1 *f1=new TF1("f1","[0]",fitrangelow1,fitrangelow2);//range of left half of hist you want to fit the background to //Chris likes 1
			TF1 *f2=new TF1("f2","[0]",fitrangehigh1,fitrangehigh2);//range of right half of hist you want to fit the background to //Chris likes 2
			hSega_data->Fit("f1","qR");
			hSega_data->Fit("f2","qR");

			for(int i=1;i<=Nbinsb;i++)
			{
				int bkgpower=1;//change this to whatever you want. 3 for high statistics
				// 		TF1 *f1=new TF1("f1","[0]+x*[1]",minrangeb+50,minrangeb+80);//range of left half of hist you want to fit the background to //useless for now
				// 		TF1 *f2=new TF1("f2","[0]+x*[1]",maxrangeb-80,maxrangeb-50);//range of right half of hist you want to fit the background to //useless for now
 				double a = (f1->GetParameter(0));//Lower BG //Chris likes 3
 				double c = (f2->GetParameter(0));//Upper BG //Chris likes 4
				double f=((double)Nbinsb-(double)i)/(double)Nbinsb;
				double g =(double)i/(double)Nbinsb;
				double h =pow(f,bkgpower)+pow(g,bkgpower);//if bkgpower = 1, h = 1
				double b=0;
				//if(minrangeb+(double)i/(double)Binsperkev<minrangeb+80) b = a; //test for fun
				//else if(minrangeb+(double)i/(double)Binsperkev>maxrangeb-80) b = c; //test for fun
				//else b = (a*pow(f,bkgpower)+c*pow(g,bkgpower))/h; //test for fun
				b = (a*pow(f,bkgpower)+c*pow(g,bkgpower))/h; //Brent created, it's good for either constant pol0 f1f2 or linear pol1 f1f2
				//b=(a+c)/2; //lower bkg and higher bkg will be the same constant, not very precise
*/

//*****************************************************************************************************
			// or manual bkg algorithm. Only for DSL low stat peaks

			double bincounts_left = 0, bincounts_left_ave = 0;
			double bincounts_right = 0, bincounts_right_ave = 0;
			int jc = 0;
			for(int i=1;i<=Nbins;i++)
			{
				hSega->SetBinContent(i,hSega_data->GetBinContent((minrange-0)/binwidth+i));//part of the experimental histogram is stored in hSega//NSCL
				//bkgshort->SetBinContent(i,background->GetBinContent((minrange-minrangeb)/binwidth+i));//part of the background histogram is stored in bkgshort
			}
			//low energy bkg
// 			for(int i=1;i<=(startpeak-minrange)/binwidth;i++)//only for DSL low stat peaks
// 			{
// 				bkgshort->SetBinContent(i,hSega->GetBinContent(i));//Directly adopt the hSega histogram as the low energy bkgshort
// 				//bincounts_left+=hSega->GetBinContent(i);
// 			}
// 			bincounts_left_ave=bincounts_left/((startpeak-minrange)/binwidth-1+1);
// 			bincounts_left=0;
// 
// 			// if (flag_bkg==1)//execute high bkg counts, which is the default average
// 			if (flag_bkg==0)//execute middle bkg counts
// 			{
// 				for(int i=1;i<=(startpeak-minrange)/binwidth;i++)//only for DSL low stat peaks
// 				{
// 					if(hSega->GetBinContent(i)<=bincounts_left_ave+sqrt(bincounts_left_ave))
// 					{
// 						bincounts_left+=hSega->GetBinContent(i); jc++;
// 					}
// 				}
// 				if (jc>0) bincounts_left_ave=bincounts_left/jc;
// 				if (bincounts_left_ave<0.2)	{		bincounts_left_ave=0.2;		}
// 			}
// 
// 			if (flag_bkg==-1)//execute low bkg counts
// 			{
// 				for(int i=1;i<=(startpeak-minrange)/binwidth;i++)//only for DSL low stat peaks
// 				{
// 					if(hSega->GetBinContent(i)<=bincounts_left_ave-sqrt(bincounts_left_ave))
// 					{
// 						bincounts_left+=hSega->GetBinContent(i); jc++;
// 					}
// 				}
// 				if (jc>0) bincounts_left_ave=bincounts_left/jc;
// 				if (bincounts_left_ave<0.2)	{		bincounts_left_ave=0.2;		}
// 			}
//
//
//			//high energy bkg
//			jc=0;
// 			for(int i=(endpeak-minrange)/binwidth+1;i<=Nbins;i++)//only for DSL low stat peaks
// 			{
// 				bkgshort->SetBinContent(i,hSega->GetBinContent(i));//Directly adopt the hSega histogram as the high energy bkgshort
// 				//bincounts_right+=hSega->GetBinContent(i);
// 			}
// 			bincounts_right_ave=bincounts_right/(Nbins-((endpeak-minrange)/binwidth+1)+1);
// 			bincounts_right=0;
// 			
// 			// if (flag_bkg==1)//execute high bkg counts, which is the default average
// 			if (flag_bkg==0)//execute middle bkg counts
// 			{
// 				for(int i=(endpeak-minrange)/binwidth+1;i<=Nbins;i++)//only for DSL low stat peaks
// 				{
// 					if(hSega->GetBinContent(i)<=bincounts_right_ave+sqrt(bincounts_right_ave))
// 					{
// 						bincounts_right+=hSega->GetBinContent(i); jc++;
// 					}
// 				}
// 				if (jc>0) bincounts_right_ave=bincounts_right/jc;
// 				if (bincounts_right_ave<0.2)	{		bincounts_right_ave=0.2;		}
// 			}
// 			
// 			if (flag_bkg==-1)//execute low bkg counts
// 			{
// 				for(int i=(endpeak-minrange)/binwidth+1;i<=Nbins;i++)//only for DSL low stat peaks
// 				{
// 					if(hSega->GetBinContent(i)<=bincounts_right_ave-sqrt(bincounts_right_ave))
// 					{
// 						bincounts_right+=hSega->GetBinContent(i); jc++;
// 					}
// 				}
// 				if (jc>0) bincounts_right_ave=bincounts_right/jc;
// 				if (bincounts_right_ave<0.2)	{		bincounts_right_ave=0.2;		}
// 			}
			


//			cout<<"lalalalala="<<bincounts_left_ave<<"	"<<bincounts_right_ave<<"	"<<jc<<endl;


// 			TF1 *f1=new TF1("f1","[0]+x*[1]",fitrangelow1,fitrangelow2);
// 			TF1 *f2=new TF1("f2","[0]+x*[1]",fitrangehigh1,fitrangehigh2);
			TF1 *fl = new TF1("fl",fline,fitrangelow1,fitrangehigh2,2); //npar=2 is the number of free parameters used by the function pol1
			//TF1 *fl = new TF1("fl",fline,fitrangelow1,fitrangehigh2,5); //npar=5 is the number of free parameters used by the function pol1+gaus
			//we want to fit only the linear background excluding the signal area
			reject = kTRUE;
// 			fl->SetParameter(2, 29);
// 			fl->SetParameter(3, 2308);
// 			fl->SetParameter(4, 2);
			hSega->Fit(fl,"L0"); // "1" will draw a full curve, "0" will not draw the fitted curve in the excluded region
			// option "L" a likelihood fit is used instead of the default chi2 square fit. 
			reject = kFALSE;
			
			//store 3 separate functions for visualization (three bands) can be scaled
			TH1D *hint_tau1 = new TH1D("hint_tau1", "Fitted func with conf.band", (fitrangelow2-fitrangelow1)*10, fitrangelow1,fitrangelow2);
			(TVirtualFitter::GetFitter())->GetConfidenceIntervals(hint_tau1, 0.900);
			hint_tau1->SetStats(kFALSE);
			hint_tau1->SetFillColor(kGreen-7);

			TH1D *hint_tau2 = new TH1D("hint_tau2", "Fitted func with conf.band", (fitrangehigh2-fitrangehigh1)*10, fitrangehigh1,fitrangehigh2);
			(TVirtualFitter::GetFitter())->GetConfidenceIntervals(hint_tau2, 0.900);
			hint_tau2->SetStats(kFALSE);
			hint_tau2->SetFillColor(kGreen-7);

			TH1D *hint_tau3 = new TH1D("hint_tau3", "Fitted func with conf.band", (fitrangehigh2-fitrangelow1)*10, fitrangelow1,fitrangehigh2);
			(TVirtualFitter::GetFitter())->GetConfidenceIntervals(hint_tau3, 0.900);
			hint_tau3->SetStats(kFALSE);
			hint_tau3->SetFillColor(kGreen-7);

			//store 2 separate functions for visualization (two lines) cannot be scaled
// 			TF1 *fleft = new TF1("fleft",fline,fitrangelow1,fitrangelow2,2);
// 			fleft->SetParameters(fl->GetParameters());
// 			hSega->GetListOfFunctions()->Add(fleft);
// 			fleft->SetLineColor(4);
// 			TF1 *fright = new TF1("fright",fline,fitrangehigh1,fitrangehigh2,2);
// 			fright->SetParameters(fl->GetParameters());
// 			hSega->GetListOfFunctions()->Add(fright);
// 			fright->SetLineColor(4);

//			for(int i=(startpeak-minrange)/binwidth+1;i<=(endpeak-minrange)/binwidth;i++)
			for(int i=1;i<=Nbins;i++)
			{
				int bkgpower=1;//change this to whatever you want. 3 for high statistics
				double a = (fl->GetParameter(0)+(minrangeb+(double)i/Binsperkev)*fl->GetParameter(1));//LE BG pol1
//				double a = (fl->GetParameter(0)+(minrangeb+(double)i/Binsperkev)*fl->GetParameter(1))+fl->GetParameter(2)*exp(-((minrangeb+(double)i/Binsperkev)-fl->GetParameter(3))*((minrangeb+(double)i/Binsperkev)-fl->GetParameter(3))/(2*fl->GetParameter(4)*fl->GetParameter(4)));//LE BG pol1+gaus
//				double c = (fl->GetParameter(0)+(minrangeb+(double)i/Binsperkev)*fl->GetParameter(1));//HE BG
//  				double a = (f1->GetParameter(0));//Lower BG //Chris likes 3
//  				double c = (f2->GetParameter(0));//Upper BG //Chris likes 4
// 				double a = bincounts_left_ave;//Lower BG //I liked 3
// 				double c = bincounts_right_ave;//Upper BG //I liked 4
// 				double f=((double)Nbinsb-(double)i)/(double)Nbinsb;
// 				double g =(double)i/(double)Nbinsb;
// 				double h =pow(f,bkgpower)+pow(g,bkgpower);//if bkgpower = 1, h = 1
// 				double b=0;
				//if(minrangeb+(double)i/(double)Binsperkev<minrangeb+80) b = a; //test for fun
				//else if(minrangeb+(double)i/(double)Binsperkev>maxrangeb-80) b = c; //test for fun
				//else b = (a*pow(f,bkgpower)+c*pow(g,bkgpower))/h; //test for fun
//				b = (a*pow(f,bkgpower)+c*pow(g,bkgpower))/h; //Brent created, it's good for either constant pol0 f1f2 or linear pol1 f1f2
				b = a; //DSL final, a or c is the same for entire range, I didn't use Brent weighted average method here
				//b=(a+c)/2; //lower bkg and higher bkg will be the same constant, not very precise
//*****************************************************************************************************
				if (flag_bkg==1.1)		{		b=b*ULbkg;		}
				if (flag_bkg==0.9)		{		b=b*LLbkg;		}
				bkgshort->SetBinContent(i,b);//fit data, get background histogram
			}

			for(int i=1;i<=Nbins;i++)
			{
				for(int ii=0; ii<nBranches; ii++)
				{
					sprintf(hname,"%s%d","Egcopy",int(ii));
					proton_copy[ii]->SetBinContent(i,proton_branch[ii]->GetBinContent(((minrange-0)+Ebin_gamma-E0_gamma)/binwidth-iscanE+i));//part of the simulated proton-γ histogram is stored in proton_copy
				}//proton_copy first bin is minrange, proton_branch first bin is 0
				//Baseline E0_gamma =readin simuroot, iscanE is used to shift the simulated spectrum. How many bins shift towards high energy? iscanE. How many keVs shift towards high energy? iscanE*binwidth.
				//iscanE<0 means shift the simulated spectrum towards low energy
			}
			// 	for(int i=1;i<=iscanE;i++)
			// 	{
			// 		for(int ii=0;ii<nBranches;ii++)
			// 		{
			// 			sprintf(hname,"p%df2",int(ii));
			// 			proton_copy[ii]->SetBinContent(i,0);
			// 		}
			// 	}

			for(int ii=0; ii<nBranches; ii++)
			{
				proton_copy[ii]->SetStats(0);
				proton_copy[ii]->SetLineColor((Color_t)(ii+3));
				if(ii>=2) proton_copy[ii]->SetLineColor((Color_t)(ii+4));
				if(ii>=6) proton_copy[ii]->SetLineColor((Color_t)(56-ii));
				proton_copy[ii]->SetMarkerColor((Color_t)(ii+3));
				//proton_branch[ii]->GetXaxis()->SetRangeUser(minrange,maxrange);
			}
			bkgshort->SetStats(0);
			hSega->SetStats(0);

			bkgshort->SetLineColor(4);
			hSega->SetLineColor(1);

			bkgshort->SetMarkerColor(3);
			hSega->SetMarkerColor(1);
			hSega->SetLineWidth(1);
			hSega->Sumw2(kFALSE);
			hSega->SetBinErrorOption(TH1::kPoisson);//TH1::kNormal or TH1::kPoisson

			//minuit**************************************************************************************************************************//
			const int nParams=2; //number of parameters that are to be found 
			TMinuit *gMin = new TMinuit(nParams);  //initialize TMinuit with a maximum of n parameters
			gMin->SetFCN(fcn);//set the address of the minimization function// fcn里会调用CompareHists子函数
			Double_t arglist[10];
			Int_t ierflg = 0;//command executed normally
			arglist[0] = 1;
			Double_t vstart[2] = {0.04,1.0};//initial guess par_Fit, par_bkg
			Double_t step[2] = {0.0001,0.01};                //fitting step

			gMin->mnparm(0, "a1", vstart[0], step[0], 0.00005,10.0,ierflg);// par for simu
			//parID, parName, initialGuess, step, lowLimit, highLimit, irrelevant //0.00几，因为是从百万高统计模拟的histo scale下来的
			gMin->mnparm(1, "a2", vstart[1], step[1], bkgDown,bkgUp,ierflg);// par for bkg //bkg 基本0.9-1.2之间，不会跟拟合的本底水平差别太大
			//gMin->mnparm(2, "a3", vstart[2], step[2], 0,25,ierflg);  

			arglist[0] = 2000;  //max number of fitting iterations
			arglist[1] = 1.;  //tolerance 1 means one sigma
			gMin->mnexcm("MIGRAD", arglist ,2,ierflg);  //run the minimization using MIGRAD			//minuit.mnexcm(command, argument list, argument number, error flag);			//command is a character array holding the command
			//argument list a double array as the arguments
			//argument number is the number of arguments			//error flag return value (!= 0)			//Precisely find the minimum, but you need to be careful since it assumes the inputs are close to the minimum
			Double_t amin,edm,errdef;
			Int_t nvpar,nparx,icstat;
			gMin->mnstat(amin,edm,errdef,nvpar,nparx,icstat);

			double pars[nParams],errs[nParams];
			for (int i=0;i<nParams;i++)
			{
				gMin->GetParameter(i,pars[i],errs[i]);
				printf("pars=	%f	errs=	%f\n",pars[i],errs[i]);
			}

			//visualization

			hFit->SetLineColor(2);
			hFit->SetMarkerColor(2);
			hFit->SetLineWidth(2);

			for(int ii=0;ii<nBranches;ii++)
			{
				proton_copy[ii]->Scale(pars[0]);//Multiply this histogram by a constant pars[0].
				proton_copy[ii]->Add(bkgshort,pars[1]);
			}
			// 	p1->Scale(pars[0]);
			// 	p2->Scale(2.9/6.73*pars[0]);
			// 	p3->Scale(1.53/6.73*pars[0]);
			//  p5->Scale(0.0392*pars[0]);
			//  p6->Scale(0.0196*pars[0]);
			// 	p1->Add(background,pars[1]);
			// 	p2->Add(background,pars[1]);
			// 	p3->Add(background,pars[1]);
			//  p5->Add(background,pars[1]);
			//  p6->Add(background,pars[1]);

			TCanvas* c1   = new TCanvas("c1","c1",1000,700);
			c1->cd();
			TPad *pad1 = new TPad("pad1", "The pad 70% of the height",0.0,0.3,1.0,1.0);// Double_t xlow, Double_t ylow, Double_t xup, Double_t yup,
			TPad *pad2 = new TPad("pad2", "The pad 30% of the height",0.0,0.0,1.0,0.3);
			pad1->SetTopMargin(0.09);
			pad1->SetRightMargin(0.02);
			pad1->SetLeftMargin(0.08);
			pad1->SetBottomMargin(0.08);
			//pad1->SetBorderMode(0);

			pad2->SetTopMargin(0.03);
			pad2->SetRightMargin(0.02);
			pad2->SetLeftMargin(0.08);
			pad2->SetBottomMargin(0.21);
			//pad2->SetBorderMode(0);

			pad1->Draw();
			pad2->Draw();
			pad1->cd();
			gStyle->SetOptTitle(0);
			
			hSega->Draw("e");
			sprintf(tname,"%s%d%s","Counts per ",binwidth," keV");
			hSega->GetXaxis()->SetTitle("Energy (keV)");
			hSega->GetYaxis()->SetTitle(tname);
			//hSega->GetXaxis()->SetTitle("Energy (keV)");
			//hSega->GetYaxis()->SetTitle("counts/0.1 keV");
			//hSega->GetXaxis()->SetLabelSize(0.08);
			hSega->GetXaxis()->CenterTitle();
			hSega->GetYaxis()->CenterTitle();
			//hSega->GetXaxis()->SetTitleOffset(0.8);
			//hSega->GetYaxis()->SetTitleOffset(0.9);
			//hSega->GetXaxis()->SetTitleSize(0.06);
			//hSega->GetYaxis()->SetTitleSize(0.1);
			//hSega->GetYaxis()->SetLabelSize(0.07);
			hSega->GetXaxis()->SetNdivisions(510);//n = n1 + 100*n2 + 10000*n3
			hSega->GetYaxis()->SetNdivisions(505);//n = n1 + 100*n2 + 10000*n3
			hSega->GetXaxis()->SetRangeUser(minrangezoom,maxrangezoom);
			
			bkgshort->Scale(pars[1]);//bkgshort scaled up/down by up to 10%
// 			hint_tau1->Scale(pars[1]);
//  			hint_tau1->Draw("e3 same");
// 			hint_tau2->Scale(pars[1]);
//  			hint_tau2->Draw("e3 same");
			hint_tau3->Scale(pars[1]);
			hint_tau3->Draw("e3 same");
			bkgshort->SetLineWidth(2);
			bkgshort->Draw("sameh");

			for (int ibin = 1; ibin <= Nbins; ibin++)
			{
				hFit->SetBinError(ibin,errs[0]/pars[0]*hFit->GetBinContent(ibin));
			}

			
			for(int ii=0;ii<nBranches;ii++)
			{
				proton_copy[ii]->SetLineWidth(2);
				proton_copy[ii]->Draw("sameh");//proton_copy[0]=p1_scale+bkgshort_scale
			}
			hFit->Draw("samehist");//hFit=p1_scale+p2_scale+p3_scale+bkgshort_scale

			char paraprint[100];
			TPaveText *textchi = new TPaveText(0.08,0.92,0.98,0.98,"brNDC");//left, down, right, up
			textchi->SetBorderSize(1);
			textchi->SetFillColor(0);
			textchi->SetTextAlign(12);//align = 10*HorizontalAlign + VerticalAlign, 12 means水平左对齐、垂直居中对齐
			textchi->SetTextFont(132);//font = 10 * fontID + precision, 12+2,12 means Symbol; 13+2, 13 means News Times Roman
			//sprintf(paraprint,"#chi^{2} = %.2f / %d = %.5f   E_{#gamma} = %.1f keV   #tau = %.1f fs   N = %.1f +/- %.1f",Chi2,(Nbinsz-2),Chi2/(Nbinsz-2),centroid+iscanE*binwidth,T0_lifetime+iscantau*Lifetimestep,pars[0]*20000,errs[0]*20000);
			sprintf(paraprint,"#chi^{2} = %.2f / %d = %.5f   E_{#gamma} = %.2f keV   #tau = %.1f fs   N = %.1f +/- %.1f", Chi2,(Nbinspeak-2),Chi2/(Nbinspeak-2),flag_Eg+iscanE*binwidth,T0_lifetime+iscantau*Lifetimestep,pars[0]*20000,errs[0]*20000);//only for DSL low stat peaks
			textchi->AddText(paraprint);
			textchi->Draw();

			pad2->cd();
			int rebin=1;//1,2,3,5,6,9,10,15,18,25,30,45,50,75,90,150,225,450// why on earth do you need to rebin here?

			TH1F *hResid=new TH1F("hResid","hResid",Nbins,minrange,maxrange);
			const int Nbinsr=Nbins;//Nbins替换的文本是个变量，还不能用于初始化，need a new const
			double Resids[Nbinsr]={0};
			double ResidsUncert[Nbinsr]={0};
			double Xvals[Nbinsr]={0};
			double XvalsUncert[Nbinsr]={0};
			int count=0;
			for(int i=1;i<=Nbins;i++)
			{
				Xvals[i-1]=minrange+binwidth*i-0.5*binwidth;//X value is at the center of bin
				// cout<<Xvals[i]<<endl;
				ResidsUncert[i-1]=0;
				XvalsUncert[i-1]=0.0;

				Resids[i-1] = hFit->GetBinContent(i)-hSega->GetBinContent(i);//simu-data //residual
				//if(hSega->GetBinContent(i)!=0) Resids[i-1] = (hFit->GetBinContent(i)-hSega->GetBinContent(i))/hFit->GetBinContent(i);//(simu-data)/simu //relative residual

				hResid->AddBinContent(i,Resids[i-1]);//simu-data
				//ResidsUncert[i-1]= pow(hSega->GetBinContent(i),0.5);
				ResidsUncert[i-1]= hSega->GetBinErrorUp(i);
				hResid->SetBinError(i,ResidsUncert[i-1]);//simu-data
				//hResid->Sumw2(kFALSE);
				if(abs(Resids[i-1])<ResidsUncert[i-1]) count++;
				//cout<< abs(Resids[i])<<"         "<<ResidsUncert[i]<<endl;
			}
			//	cout<<endl<<count<<" of "<<Nbins<<" error bars go through the fit."<<endl;//count: the larger the better
			//TGraphErrors *gr;
			//gr=new TGraphErrors(Nbins,Xvals,Resids,XvalsUncert,ResidsUncert); 
			TGraph *gr=new TGraph(Nbinsr,Xvals,Resids); //TGraph *gr1=new TGraph(n,x,y);建立曲线图、散点图//Xvals and Resids are data array, start from [0], you cannot store in Xvals and Resids starting from [1]
			hResid->GetXaxis()->SetTitle("Energy (keV)");

			sprintf(tname,"#splitline{   Residual}{counts/%d keV}",binwidth);//residual
			//sprintf(tname,"#splitline{Relative residual}{counts/%0.1f keV}",binwidth);//relative residual

			hResid->GetYaxis()->SetTitle(tname);
			hResid->GetXaxis()->CenterTitle();
			hResid->GetYaxis()->CenterTitle();
			hResid->GetXaxis()->SetTitleOffset(1.1);
			hResid->GetYaxis()->SetTitleOffset(0.4);
			hResid->GetXaxis()->SetTitleSize(0.08);
			hResid->GetXaxis()->SetLabelOffset(0.015);
			hResid->GetXaxis()->SetLabelSize(0.08);
			hResid->GetYaxis()->SetLabelSize(0.08);
			hResid->GetYaxis()->SetTitleSize(0.08);
			hResid->GetXaxis()->SetNdivisions(510);//n = n1 + 100*n2 + 10000*n3
			// 	gr->GetXaxis()->SetNdivisions(10,10,1);
			hResid->GetYaxis()->SetNdivisions(505);

			hResid->GetXaxis()->SetRangeUser(minrangezoom,maxrangezoom);
			hResid->GetYaxis()->SetRangeUser(-Yaxisrange,Yaxisrange);
			//gr->SetMarkerStyle(3);
			hResid->SetStats(0);

			//hResid->SetBinErrorOption(TH1::kPoisson);//TH1::kNormal or TH1::kPoisson
			hResid->Draw("e");
			TLine *T1=new TLine(minrangezoom,0,maxrangezoom,0);
			cout<<"=======	Eg=	"<<flag_Eg+iscanE*binwidth<<"	tau=	"<<T0_lifetime+iscantau*Lifetimestep<<" fs"<<endl;
			T1->Draw("R");
			//save text
			//fprintf(outfile2, "%d	%e", T0_lifetime+iscantau*Lifetimestep, LBayesian);//LB is global varible, use LB to extract the LB corresponding to the minimum Chi2 iterated by fcn(), and output LB in main().//DSLmodify use this sentence to output 2D LB matrix
			//fprintf(outfile, "%d	%.3f", T0_lifetime+iscantau*Lifetimestep, Chi2);//Chi2 is global varible, use Chi2 to extract the minimized Chi2 obtained in fcn(), and output Chi2 in main().//DSLmodify use this sentence to output 2D chisquare matrix
			//fprintf(outfile, "Ea=	%d	tau=	%.1f	E0=	%.3f	Emax=	%.3f	-2lnL =	%.3f	-2lnL/NDF =	%.5f	%d\n", Ea,T0_lifetime+iscantau*Lifetimestep,centroid,centroid+iscanE*binwidth,Chi2,Chi2/(Nbinsz-2),Nbinsz-2);//Chi2 is global varible, use Chi2 to extract the minimized Chi2 obtained in fcn(), and output Chi2 in main(). Use this sentence to output 1D chisquare.

			//save figure
			sprintf(tname,"%s%.0f%s%.0f%s%.2f%s%.1f%s%.1f%s%.2f%s%.1f%s%.1f","D:/X/out/DSL_Comparison/Eg",E0_gamma,"test/S31_Gamma",E0_gamma,"_Eg",flag_Eg,"_Tau",T0_lifetime+iscantau*Lifetimestep,"_",T0_lifetime+iscantau*Lifetimestep,"_SP",flag_SP,"_AC",flag_AC,"_bkg",flag_bkg);
			sprintf(outfigname,"%s%s",tname,".png");
			c1->SaveAs(outfigname);
			
			//save root file
			sprintf(tname,"%s%.0f%s%.0f%s%.2f%s%.1f%s%.1f%s%.2f%s%.1f%s%.1f","D:/X/out/DSL_Comparison/Eg",E0_gamma,"test/S31_Gamma",E0_gamma,"_Eg",flag_Eg,"_Tau",T0_lifetime+iscantau*Lifetimestep,"_",T0_lifetime+iscantau*Lifetimestep,"_SP",flag_SP,"_AC",flag_AC,"_bkg",flag_bkg);
			sprintf(outrootname,"%s%s",tname,".root");
			TFile *fout = new TFile(outrootname,"RECREATE");

			hSega->Write("hSega");
			hFit->Write("hFit");
			gr->Write("Residual Plot");
			hResid->Write("hResidual");
			for(int ii=0;ii<nBranches;ii++)
			{
				proton_copy[ii]->Write();
			}
			bkgshort->Write("background");


// 			//save histogram txt for MADAI "model results"
// 			sprintf(tname,"%s%.0f%s%.0f%s%.2f%s%.1f%s%.2f%s%.1f%s%.1f","D:/X/out/DSL_Comparison/Eg",E0_gamma,"/S31_Gamma",E0_gamma,"_Eg",flag_Eg,"_Tau",T0_lifetime+iscantau*Lifetimestep,"_SP",flag_SP,"_AC",flag_AC,"_bkg",flag_bkg);
// 			sprintf(outtxtname,"%s%s",tname,"results.dat"); // dat files for backup
// 			FILE *outfilehist = fopen (outtxtname,"w"); //ofstream doesn't work
// 			for (int ibin = (startpeak - minrange) / binwidth + 1; ibin <= (endpeak - minrange) / binwidth; ibin++)
// 			{
// 				//cout << hFit->GetBinCenter(ibin) << "	" << hFit->GetBinContent(ibin) << endl;
// 				if (ibin < (endpeak - minrange) / binwidth) fprintf(outfilehist, "%.0f	%.4f	%.4f\n", hFit->GetBinCenter(ibin), hFit->GetBinContent(ibin), errs[0]/pars[0]*hFit->GetBinContent(ibin));
// 				if (ibin == (endpeak - minrange) / binwidth) fprintf(outfilehist, "%.0f	%.4f	%.4f", hFit->GetBinCenter(ibin), hFit->GetBinContent(ibin), errs[0]/pars[0]*hFit->GetBinContent(ibin));
// 			}
// 			fclose(outfilehist);
// 
// 			sprintf(tname,"%s%.0f%s%04d","D:/X/out/DSL_MADAI/Eg",E0_gamma,"/run",ifolder);
// 			sprintf(outtxtname,"%s%s",tname,"/results.dat"); // dat files to MADAI
// 			FILE *outfilehistMADAI = fopen (outtxtname,"w"); //ofstream doesn't work
// 			for (int ibin = (startpeak - minrange) / binwidth + 1; ibin <= (endpeak - minrange) / binwidth; ibin++)
// 			{
// 				//cout << hFit->GetBinCenter(ibin) << "	" << hFit->GetBinContent(ibin) << endl;
// 				if (ibin < (endpeak - minrange) / binwidth) fprintf(outfilehistMADAI, "%.0f	%.4f\n", hFit->GetBinCenter(ibin), hFit->GetBinContent(ibin));//, errs[0]/pars[0]*hFit->GetBinContent(ibin));
// 				if (ibin == (endpeak - minrange) / binwidth) fprintf(outfilehistMADAI, "%.0f	%.4f", hFit->GetBinCenter(ibin), hFit->GetBinContent(ibin));//, errs[0]/pars[0]*hFit->GetBinContent(ibin));
// 			}
// 			fclose(outfilehistMADAI);
// 
// 			sprintf(tname,"%s%.0f%s%04d","D:/X/out/DSL_MADAI_model_uncertainty/Eg",E0_gamma,"/run",ifolder);
// 			sprintf(outtxtname,"%s%s",tname,"/results.dat"); // dat files to MADAI with model uncertainty
// 			FILE *outfilehistMADAIu = fopen (outtxtname,"w"); //ofstream doesn't work
// 			for (int ibin = (startpeak - minrange) / binwidth + 1; ibin <= (endpeak - minrange) / binwidth; ibin++)
// 			{
// 				//cout << hFit->GetBinCenter(ibin) << "	" << hFit->GetBinContent(ibin) << endl;
// 				if (ibin < (endpeak - minrange) / binwidth) fprintf(outfilehistMADAIu, "%.0f	%.4f	%.4f\n", hFit->GetBinCenter(ibin), hFit->GetBinContent(ibin), errs[0]/pars[0]*hFit->GetBinContent(ibin));
// 				if (ibin == (endpeak - minrange) / binwidth) fprintf(outfilehistMADAIu, "%.0f	%.4f	%.4f", hFit->GetBinCenter(ibin), hFit->GetBinContent(ibin), errs[0]/pars[0]*hFit->GetBinContent(ibin));
// 			}
// 			fclose(outfilehistMADAIu);
// 
// 
// 			//save parameters txt for MADAI "model input parameters"
// 			sprintf(tname,"%s%.0f%s%.0f%s%.2f%s%.1f%s%.2f%s%.1f%s%.1f","D:/X/out/DSL_Comparison/Eg",E0_gamma,"/S31_Gamma",E0_gamma,"_Eg",flag_Eg,"_Tau",T0_lifetime+iscantau*Lifetimestep,"_SP",flag_SP,"_AC",flag_AC,"_bkg",flag_bkg);
// 			sprintf(outtxtname,"%s%s",tname,"parameters.dat"); // backup dat files
// 			FILE *outfilepara = fopen (outtxtname,"w"); //ofstream doesn't work
// 				//cout << hFit->GetBinCenter(ibin) << "	" << hFit->GetBinContent(ibin) << endl;
// 			fprintf(outfilepara, "Tau	%.1f\nEg	%.2f\nbkg	%.1f\nSP	%.2f\nAC	%.1f", T0_lifetime+iscantau*Lifetimestep, flag_Eg, flag_bkg, flag_SP, flag_AC);
// 			fclose(outfilepara);
// 
// 			sprintf(tname,"%s%.0f%s%04d","D:/X/out/DSL_MADAI/Eg",E0_gamma,"/run",ifolder);
// 			sprintf(outtxtname,"%s%s",tname,"/parameters.dat"); // dat files to MADAI
// 			FILE *outfileparaMADAI = fopen (outtxtname,"w"); //ofstream doesn't work
// 			//cout << hFit->GetBinCenter(ibin) << "	" << hFit->GetBinContent(ibin) << endl;
// 			fprintf(outfileparaMADAI, "Tau	%.1f\nEg	%.2f\nbkg	%.1f\nSP	%.2f\nAC	%.1f", T0_lifetime+iscantau*Lifetimestep, flag_Eg, flag_bkg, flag_SP, flag_AC);
// 			fclose(outfileparaMADAI);
// 
// 			sprintf(tname,"%s%.0f%s%04d","D:/X/out/DSL_MADAI_model_uncertainty/Eg",E0_gamma,"/run",ifolder);
// 			sprintf(outtxtname,"%s%s",tname,"/parameters.dat"); // dat files to MADAI with model uncertainty
// 			FILE *outfileparaMADAIu = fopen (outtxtname,"w"); //ofstream doesn't work
// 			//cout << hFit->GetBinCenter(ibin) << "	" << hFit->GetBinContent(ibin) << endl;
// 			fprintf(outfileparaMADAIu, "Tau	%.1f\nEg	%.2f\nbkg	%.1f\nSP	%.2f\nAC	%.1f", T0_lifetime+iscantau*Lifetimestep, flag_Eg, flag_bkg, flag_SP, flag_AC);
// 			fclose(outfileparaMADAIu);

			cout<<"ifolder=	"<<ifolder<<endl;
			ifolder++;



			
			delete fout;
			//fin_data->Close();
			//fin_simu0->Close();
			//fin_simu1->Close();
		}//for(int iscanE=0;iscanE<1;iscanE++)
		//fprintf(outfile, "\n");//DSLmodify use this sentence to output 2D chisquare matrix
		//fprintf(outfile2, "\n");//DSLmodify use this sentence to output 2D chisquare matrix
	}//for(int iscantau=0;iscantau<1;iscantau++)
	//	return 0;
	//fprintf(outfile, "\n");//DSLmodify use this sentence to output 2D chisquare matrix
	//fprintf(outfile2, "\n");//DSLmodify use this sentence to output 2D chisquare matrix
				}
			}
		}
	}// for (flag_Eg)

}//main_readhists; fcn(); comparehists(); main_draw_save;








void fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
	//npar number of free parameters involved in minimization
	//gin partial derivatives (return values) computed gradient values (optional)
	//f the function value itself (return value)
	//par parameter values
	//iflag flag word to switch between several actions of FCN	//Usually only f and par are essential
{
	if (hFit) delete hFit;
	hFit = (TH1F*)proton_copy[0]->Clone("hFit");//Making a copy of histogram p1
	//hFit->Sumw2();//histogram is already filled, the sum of squares of weights is filled with the existing bin contents
	//This function is automatically called when the histogram is created

	//Set your different feeding intensities here
	for(int ii=1;ii<nBranches;ii++)
	{
		hFit->Add(proton_copy[ii],0.075);//p2 p3的相对intensity都用文献值比
	}
	hFit->Scale(par[0]);  //Multiply this histogram by a constant
	hFit->Add(bkgshort,par[1]);//hFit=p1+p2+p3+bkgshort
//	f = CompareHists(hFit,hSega,(minrangezoom-minrange)/binwidth+1,(maxrangezoom-minrange)/binwidth);//first bin -> last bin entire region
	f = CompareHists(hFit,hSega,(startpeak-minrange)/binwidth+1,(endpeak-minrange)/binwidth);//first bin -> last bin peak region
//	f = CompareHists(hFit,hSega,1,Nbinsz);//first bin -> last bin the same as the sentence above
	Chi2=f;//Chi2 is global varible, use Chi2 to extract the minimized Chi2 obtained in fcn(), and output Chi2 in main().
	//cout<<Chi2<<endl;
	//ofstream outfile("D:/X/out/Compareshow.dat",ios::out);
	//FILE *outfile=fopen ("D:/X/out/Compareshow.dat","w");//ofstream报错
	//fprintf(outfile, "E0_gamma=	%.1f	ChiSquare =	%.3f	ChiSquare/NDF =	%.5f\n", Nbins,f,f/(Nbins-2));//keV
	//cout<<"ChiSquare = " <<f<<", ChiSquare/NDF = "<<f/(Nbins-2)<<endl;
}






double CompareHists(TH1F* his1, TH1F* his2, int binMin, int binMax)
	//compare two histograms (data and simulation). returns chi-square
	//Two histograms must be with the same length
{
	double chisquareNeyman=0,chisquarePearson=0,likelihood=0;     //chi square
	long double ydata_factorial=1;
	LBayesian=1;
	int binN = his1->GetNbinsX();                         //get number of bins in his1
	//   return 0;

	//if not identical to his2, exit with error.
	if (binN!=his2->GetNbinsX())
	{
		printf("ERROR <CompareHists>: Different number of bins in histograms!%d\t%d\n",binN,his2->GetNbinsX());
		exit(1);
	}
	double err_ydata=0,residual,yfit,ydata;
	for (int i=binMin;i<=binMax;i++)                         //loop over all bins   
	{
		ydata=his2->GetBinContent(i); //hSega
		err_ydata=his2->GetBinError(i);
		yfit=his1->GetBinContent(i); //hFit
		residual=ydata-yfit; //ML method

		if(ydata!=0) {chisquareNeyman += residual*residual/ydata;}//Neyman Chi2 method
		if(yfit!=0) {chisquarePearson += residual*residual/yfit;}//Pearson Chi2 method
		if(yfit!=0&&ydata!=0) {likelihood+=yfit-ydata+ydata*log(ydata/yfit);} //ML method
		if(yfit!=0&&ydata==0) {likelihood+=yfit-ydata;} //ML method for bins with ydata = 0, the log term is zero.
		//if(y[0]>0)	dy[0] = (his2->GetBinErrorLow(i)+his2->GetBinErrorUp(i))/2;
		//if(y[0]==0)	dy[0] = his2->GetBinErrorUp(i);//empty bin ErrorUp=1.8, ErrorLow=0, Error=0;

		ydata_factorial=1;
		for (long ii=1; ii<=ydata; ii++) ydata_factorial = ydata_factorial*ii; //Bayesian
		//cout<<"yfit= "<<yfit<<" ydata= "<<ydata<<" ydata_factorial= "<<ydata_factorial<<" pow(yfit,ydata)= "<<pow(yfit,ydata)<<" LB= "<<(pow(yfit,ydata)/ydata_factorial)*exp(-yfit)<<endl; //Bayesian
		if(yfit!=0&&ydata_factorial!=0) LBayesian*=(pow(yfit,ydata)/ydata_factorial)*exp(-yfit); //Bayesian
	}
	
	likelihood=likelihood*2; //ML method
	return likelihood; //ML method
}

Double_t fline(Double_t *x, Double_t *par) // pol1
{ // Double_t fcn(Double_t *x, Double_t *params)
	if (reject && x[0] > startpeak && x[0] < endpeak) {
		TF1::RejectPoint();
		return 0;
	}
	return par[0] + par[1]*x[0];
}

// Double_t fline(Double_t *x, Double_t *par) // pol1+gaus
// { // Double_t fcn(Double_t *x, Double_t *params)
// 	if (reject && x[0] > startpeak && x[0] < endpeak) {
// 		TF1::RejectPoint();
// 		return 0;
// 	}
// 	return par[0] + par[1]*x[0]+par[2]*exp(-(x[0]-par[3])*(x[0]-par[3])/(2*par[4]*par[4]));
// }
