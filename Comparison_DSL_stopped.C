#include "TH1F.h"
#include <cmath>
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

TFile *fin_simu,*fin_data;
TH1F *hFit,*hSega,*hSega_data,*bkgshort,*background;//

const double Binsperkev=1; //Number of bins per keV, the only place to change binwidth, all other variables is related to Binsperkev // 10 for NSCL, 2 for RIBLL, useless for DSL
const int binwidth=1; //binwidths in units of keV
const int factor_rebin=1; //simu and data Rebin factor
const double E0_gamma=2814; //DSLmodify //E0+i*binwidth, a baseline E0 determined by the file name of the simuroot.root, 2752 4236 2868 1366, decide which simuroot to use// if you want to use for() loop and E0_gamma is related to iscanE, E0_gamma cannot be global variable out of the main(), as loop has to be in the main().
//E0-Ebin to shift the simulated spectrum.h E0-Ebin>0 is moving simulated spectrum towards higher energy
const double Ebin_gamma=E0_gamma; //modify center for drawing figures: 2754 4238 2870 1369, don't change to other values // if you want to use for() loop and E0_gamma is related to iscanE, E0_gamma cannot be global variable out of the main(), as loop has to be in the main().

//Enable one of the energies and others should be commented out.
// if(E0_gamma==1248)
// 	const int up=170;//Ex=1248.43, Eg=1248.40
// 	const int down=84;
// 	const int Ea=49;
// 	double T0_lifetime = 200;// read in simulation files   // This is Excited State Mean Lifetime (fs), Mean Lifetime τ=T/ln(2); half-life T=τ*ln(2)
// 	int Lifetimestep = 200; // read in simulation files

// if(E0_gamma==2234)
// 	const int up=316;//Ex=2234.06, Eg=2233.97
// 	const int down=134;
// 	const int Ea=47;
// 	double T0_lifetime = 80;
// 	int Lifetimestep = 10;

//if(E0_gamma==3076)//DSLmodify
// 	const int up=328;//Ex=3076.40, Eg=3076.24
// 	const int down=-120;
// 	const int Ea=46;
// 	double T0_lifetime = 0;
// 	int Lifetimestep = 1;

// if(E0_gamma==4971)
// 	const int up=470;//Ex=4970.7, Eg=4970.2
// 	const int down=-235;
// 	const int Ea=42;
// 	double T0_lifetime = 0;
// 	int Lifetimestep = 1;

// if(E0_gamma==5156)
// 	const int up=480;//Ex=5156.1, Eg=5155.7
// 	const int down=-226;
// 	const int Ea=42;
// 	double T0_lifetime = 0;
// 	int Lifetimestep = 1;

//	if(E0_gamma==6129)
// 	const int up=16;
// 	const int down=40;
// 	const int Ea=0;
// 	double T0_lifetime = 10000000;
// 	int Lifetimestep = 1;

//	if(E0_gamma==2814)
	const int up=10;
	const int down=18;
	const int Ea=0;
	double T0_lifetime = 10000000;
	int Lifetimestep = 1;

	const double minrange=Ebin_gamma-down; //always 100-keV window
	const double maxrange=Ebin_gamma+up; //always 100-keV window
	const double minrangeb=Ebin_gamma-down; //longer range to fit the background correctly
	const double maxrangeb=Ebin_gamma+up; //longer range to fit the background correctly
	const double minrangezoom=Ebin_gamma-down; //shorter range to calculate chi2 25 30 30 40
	const double maxrangezoom=Ebin_gamma+up; //shorter range to calculate chi2 25 30 30 40

	const int Nbins=(maxrange-minrange)/binwidth;//make sure this is an integer // there is another global Nbins at the beginning because fcn() also needs an Nbins.
	const int Nbinsb=(maxrangeb-minrangeb)/binwidth;//make sure this is an integer // there is another global Nbins at the beginning, you'd better keep them consistent
	const int Nbinsz=(maxrangezoom-minrangezoom)/binwidth;//make sure this is an integer // there is another global Nbins at the beginning, you'd better keep them consistent


//double E0_gamma=2753.0+double(iscanE)*2.0/10.0;                  // Energy of Gamma-ray (keV) useless
const int nBranches=1;                  //13,19;2,3 Number of proton Branches feeding daughter modify
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

int Flag_SP = 1;                         //stopping power default=1
int Flag_Det = 1;                      //Detector resolution. No f1 sampling when Flag_Det is false, will be very fast. default=1
int Flag_Ang = 1;                       //Tests random cos vs random between -1 and 1 for angular dist (false tests cos and true tests random X(-1,1) default=1
int Flag_Escape = 0;					//If take recoil escaping into consideration default=0
int Flag_sig = 0;					//If take σ error into consideration default=0, -1 or 1, minus err or plus err
int Flag_tau = 0;					//If take τ error into consideration default=0, -1 or 1, minus err or plus err
int Flag_st = 0;					//If take stopping power error into consideration default=0, -1 or 1, minus err or plus err
int NumberofSteps=200;	//Divide the entire lifetime into several steps
double Chi2=0;//Chi2 is global varible, use Chi2 to extract the minimized Chi2 obtained in fcn(), and output Chi2 in main().
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

void Comparison_DSL_stopped()
{
	FILE *outfile=fopen ("D:/X/out/Compare.dat","a");//ofstream报错
	for(int iscantau=0;iscantau<1;iscantau++)//27 for Eg1248, 43 for Eg2234, 25 for Eg3076,4971,5156
	{
		for(int iscanE=-3;iscanE<4;iscanE++)//iscanE<0 means shift the simulated spectrum towards low energy
			//(-8,9) 17 energy points
		{//DSLportal, shift simulated spectrum with different gamma energies 1369 i=34; 2754 i=29; 2870 i=30; 4238 i=27
			//iscanE+1 bin shift +1, energy +2 keV
			double fitrangelow1, fitrangelow2, fitrangehigh1, fitrangehigh2;
			if(Ebin_gamma==1248)	{fitrangelow1=minrangeb; fitrangelow2=minrangeb+60; fitrangehigh1=maxrangeb-60; fitrangehigh2=maxrangeb;}
			if(Ebin_gamma==2234)	{fitrangelow1=minrangeb; fitrangelow2=minrangeb+120; fitrangehigh1=maxrangeb-120; fitrangehigh2=maxrangeb;}
			if(Ebin_gamma==3076)	{fitrangelow1=minrangeb; fitrangelow2=minrangeb+50; fitrangehigh1=maxrangeb-60; fitrangehigh2=maxrangeb;}
			if(Ebin_gamma==4971)	{fitrangelow1=minrangeb; fitrangelow2=minrangeb+50; fitrangehigh1=maxrangeb-60; fitrangehigh2=maxrangeb;}
			if(Ebin_gamma==5156)	{fitrangelow1=minrangeb; fitrangelow2=minrangeb+60; fitrangehigh1=maxrangeb-60; fitrangehigh2=maxrangeb;}
			if(Ebin_gamma==6129)	{fitrangelow1=minrangeb; fitrangelow2=minrangeb+10; fitrangehigh1=maxrangeb-8; fitrangehigh2=maxrangeb;}
			if(Ebin_gamma==2814)	{fitrangelow1=minrangeb; fitrangelow2=minrangeb+6; fitrangehigh1=maxrangeb-5; fitrangehigh2=maxrangeb;}

			//readHists();
			Chi2=0;//Chi2 is global varible, use Chi2 to extract the minimized Chi2 obtained in fcn(), and output Chi2 in main(). Must clean Chi2 before each run in main().
			char simurootname[400],outtxtname[400],tname[400],hname[400], ename[400], pname[400];
			sprintf(pname,"%d_",(int)((E_CoM[0]+0.000001)*1000));
			for(int ii=1;ii<nBranches;ii++)
			{
				sprintf(ename,"%d_",(int)((E_CoM[ii]+0.000001)*1000));
				strcat(pname,ename);//for the output root file name
				//cout<<pname<<endl;
			}
			//Thomas don't change binwidth 0.1, all simulated root files are 0.1-keV bin. The png name is associated with binwidth. modify
			//sprintf(tname,"%s%.2f%s%d%s%.1f%s%d%s%s%s%d%s%d%s%d%s%d%s%d%s%d%s%.1f%s","D:/X/out/Thomas/Sim25Si_T",T0_lifetime,"_Step",NumberofSteps,"_Eg",E0_gamma,"_feed",nBranches,"_En",pname,"SP",Flag_SP,"_Det",Flag_Det,"_Ang",Flag_Ang,"_st",Flag_st,"_sig",Flag_sig,"_tau",Flag_tau,"_binwid",0.1,"_P10_Sun_Thomas");
			//		sprintf(tname,"%s%.0f%s%.1f","D:/X/out/DSL/S31_Eg",E0_gamma,"_tau",T0_lifetime+iscantau*Lifetimestep);//read in lots of simulation files for 0.5 fs file
			sprintf(tname,"%s%.0f%s","D:/X/out/DSL/K39_Eg",E0_gamma,"_tauinfinity");//read in lots of simulation files DSLmodify for others
			sprintf(simurootname,"%s%s",tname,"_EMG.root");
			//Lit don't change binwidth 0.1, all simulated root files are 0.1-keV bin. The png name is associated with binwidth. RIBLL
			//	sprintf(simurootname,"%s%.2f%s%d%s%.1f%s%d%s%s%s%d%s%d%s%d%s%.1f%s","D:/X/out/Lit/Sim25Si_T",T0_lifetime,"_Step",NumberofSteps,"_Eg",E0_gamma,"_feed",nBranches,"_En",pname,"SP",Flag_SP,"_Det",Flag_Det,"_Ang",Flag_Ang,"_binwid",0.1,"_Si_Sun_Lit.root");
			fin_simu = new TFile(simurootname); //Rootfile with Simulation histograms
			//fin_data = new TFile("D:/X/out/DSL/nice.root"); //Rootfile with Experimental Data Histogram modify
			fin_data = new TFile("D:/X/out/DSL/moreAlphaEgates.root"); //Rootfile with Experimental Data Histogram modify
			//fin_data = new TFile("D:/X/out/Si25/100eV-bin-gated_sun_25Si.root"); //Rootfile with Experimental Data Histogram gated
			double centroid=0;

			if(Ebin_gamma==1248)	{centroid=1248;}
			if(Ebin_gamma==2234)	{centroid=2234;}
			if(Ebin_gamma==3076)	{centroid=3076;}
			if(Ebin_gamma==4971)	{centroid=4971;}
			if(Ebin_gamma==5156)	{centroid=5156;}
			if(Ebin_gamma==6129)	{centroid=6128.63;}
			if(Ebin_gamma==2814)	{centroid=2814.06;}

			//Get Histograms from Simulation
			for(int ii=0;ii<nBranches;ii++)
			{
				sprintf(hname,"Eg");
				proton_branch[ii]=(TH1F*)fin_simu->Get(hname);
				proton_branch[ii]->Rebin(factor_rebin);//Rebin simulated histograms from 1-keV bin to 2-keV bin because data are 2-keV bin
				sprintf(hname,"Eg");
				maximum_branch[ii]=(TH1F*)fin_simu->Get(hname);
			}
			for(int ii=1;ii<nBranches;ii++)
			{
				maximum_branch[0]->Add(maximum_branch[0],maximum_branch[ii]);
			}
//			centroid=maximum_branch[0]->GetMean();//what the original energy for a unbroadened gamma peak would be, mean energy of the simulated gamma peak, useless for DSL
			// 	sprintf(hname,"TotalGamma_AllneutronFeedings");
			// 	totalmax=(TH1F*)fin_simu->Get(hname);
			// 	totalmax->Rebin(10);
			// 	int whichbin=totalmax->GetMaximumBin();
			// 	centroid=totalmax->GetBinCenter(whichbin);
			//sprintf(hname,"%s%d%s","hGriffinAddback_alpha",Ea,"_t0");//alpha gate DSLmodify
			sprintf(hname,"%s","hGriffinAddback_0");//
			hSega_data= (TH1F*)fin_data->Get(hname); //Get Gamma ray spectrum
			//	hSega_data->Rebin(factor_rebin);//Rebin data histograms already 2-keV bin, no need to rebin

			hSega=new TH1F("hSega","hSega",Nbins,minrange,maxrange);
			bkgshort=new TH1F("bkgshort","bkgshort",Nbins,minrange,maxrange);
			background=new TH1F("background","background",Nbinsb,minrangeb,maxrangeb);

			for(int ii=0;ii<nBranches;ii++)
			{
				sprintf(hname,"Egcopy");
				proton_copy[ii]=new TH1F(hname,hname,Nbins,minrange,maxrange);//for outroot
			}

//*****************************************************************************************************
			//bkg algorithm created by Brent

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

//*****************************************************************************************************
			// or constant bkg algorithm. Chris prefers this for DSL

// 			TF1 *f1=new TF1("f1","[0]",fitrangelow1,fitrangelow2);//range of left half of hist you want to fit the background to //Chris likes 1
// 			TF1 *f2=new TF1("f2","[0]",fitrangehigh1,fitrangehigh2);//range of right half of hist you want to fit the background to //Chris likes 2
// 			hSega_data->Fit("f1","qR");
// 			hSega_data->Fit("f2","qR");
// 
// 			for(int i=1;i<=Nbinsb;i++)
// 			{
// 				int bkgpower=1;//change this to whatever you want. 3 for high statistics
// 				// 		TF1 *f1=new TF1("f1","[0]+x*[1]",minrangeb+50,minrangeb+80);//range of left half of hist you want to fit the background to //useless for now
// 				// 		TF1 *f2=new TF1("f2","[0]+x*[1]",maxrangeb-80,maxrangeb-50);//range of right half of hist you want to fit the background to //useless for now
// 				double a = (f1->GetParameter(0));//Lower BG //Chris likes 3
// 				double c = (f2->GetParameter(0));//Upper BG //Chris likes 4
// 				double f=((double)Nbinsb-(double)i)/(double)Nbinsb;
// 				double g =(double)i/(double)Nbinsb;
// 				double h =pow(f,bkgpower)+pow(g,bkgpower);//if bkgpower = 1, h = 1
// 				double b=0;
// 				//if(minrangeb+(double)i/(double)Binsperkev<minrangeb+80) b = a; //test for fun
// 				//else if(minrangeb+(double)i/(double)Binsperkev>maxrangeb-80) b = c; //test for fun
// 				//else b = (a*pow(f,bkgpower)+c*pow(g,bkgpower))/h; //test for fun
// 				b = (a*pow(f,bkgpower)+c*pow(g,bkgpower))/h; //Brent created, it's good for either constant pol0 f1f2 or linear pol1 f1f2
// 				//b=(a+c)/2; //lower bkg and higher bkg will be the same constant, not very precise
//*****************************************************************************************************
				background->SetBinContent(i,b);//fit data, get background histogram
			}

			for(int i=1;i<=Nbins;i++)
			{
				//hSega->SetBinContent(i,hSega_data->GetBinContent((minrange-100)*Binsperkev+i));//part of the experimental histogram is stored in hSega//RIBLL
				hSega->SetBinContent(i,hSega_data->GetBinContent((minrange-0)/binwidth+i));//part of the experimental histogram is stored in hSega//NSCL
				bkgshort->SetBinContent(i,background->GetBinContent((minrange-minrangeb)/binwidth+i));//part of the background histogram is stored in bkgshort
			}
			for(int i=1;i<=Nbins;i++)
			{
				for(int ii=0;ii<nBranches;ii++)
				{
					sprintf(hname,"Egcopy");
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

			for(int ii=0;ii<nBranches;ii++)
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

			bkgshort->SetLineColor(7);
			hSega->SetLineColor(1);

			bkgshort->SetMarkerColor(7);
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

			gMin->mnparm(0, "a1", vstart[0], step[0], 0,40,ierflg);// par for simu
			//parID, parName, initialGuess, step, lowLimit, highLimit, irrelevant //0.00几，因为是从百万高统计模拟的histo scale下来的
			gMin->mnparm(1, "a2", vstart[1], step[1], 0.83,1.2,ierflg);// par for bkg //bkg 基本0.9-1.2之间，不会跟拟合的本底水平差别太大
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
			bkgshort->Scale(pars[1]);

			TCanvas* c1   = new TCanvas("c1","c1",1000,700);
			c1->cd();
			TPad *pad1 = new TPad("pad1", "The pad 70% of the height",0.0,0.3,1.0,1.0);// Double_t xlow, Double_t ylow, Double_t xup, Double_t yup,
			TPad *pad2 = new TPad("pad2", "The pad 30% of the height",0.0,0.0,1.0,0.3);
			pad1->SetTopMargin(0.02);
			pad1->SetRightMargin(0.02);
			pad1->SetLeftMargin(0.08);
			pad1->SetBottomMargin(0.08);
			//pad1->SetBorderMode(0);

			pad2->SetTopMargin(0.01);
			pad2->SetRightMargin(0.02);
			pad2->SetLeftMargin(0.08);
			pad2->SetBottomMargin(0.2);
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
			hFit->Draw("samehist");//hFit=p1_scale+p2_scale+p3_scale+bkgshort_scale
			for(int ii=0;ii<nBranches;ii++)
			{
				proton_copy[ii]->SetLineWidth(2);
				//proton_copy[ii]->Draw("sameh");//proton_copy[0]=p1_scale+bkgshort_scale
			}
			// 	p1->Draw("sameh");
			// 	p2->Draw("sameh");
			// 	p3->Draw("sameh");
			//  p5->Draw("sameh");
			//  p6->Draw("sameh");
			bkgshort->SetLineWidth(2);
			bkgshort->Draw("sameh");
			// 	background->SetLineWidth(2);
			// 	background->Draw("sameh");

			char paraprint[100];
			TPaveText *textchi = new TPaveText(0.48,0.88,0.97,0.97,"brNDC");//left, down, right, up
			textchi->SetBorderSize(1);
			textchi->SetFillColor(0);
			textchi->SetTextAlign(12);//align = 10*HorizontalAlign + VerticalAlign, 12 means水平左对齐、垂直居中对齐
			textchi->SetTextFont(132);//font = 10 * fontID + precision, 12+2,12 means Symbol; 13+2, 13 means News Times Roman
			sprintf(paraprint,"-2lnL/ndf = %.2f / %d = %.5f     t = infinity",Chi2,(Nbinsz-2),Chi2/(Nbinsz-2));
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
				ResidsUncert[i-1]= pow(hSega->GetBinContent(i),0.5);
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
			hResid->GetXaxis()->SetTitleOffset(1.0);
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
			//gr->SetMarkerStyle(3);
			hResid->SetStats(0);

			//hResid->SetBinErrorOption(TH1::kPoisson);//TH1::kNormal or TH1::kPoisson
			hResid->Draw("e");
			TLine *T1=new TLine(minrangezoom,0,maxrangezoom,0);
			cout<<"=======	Eg=	"<<E0_gamma+iscanE*binwidth<<"	tau=	"<<T0_lifetime+iscantau*Lifetimestep<<" fs"<<endl;
			T1->Draw("R");
			//save text

			//fprintf(outfile, "%.3f	", Chi2);//Chi2 is global varible, use Chi2 to extract the minimized Chi2 obtained in fcn(), and output Chi2 in main().//DSLmodify use this sentence to output 2D chisquare matrix
			fprintf(outfile, "Ea=	%d	tau=	%.1f	E0=	%.3f	Emax=	%.3f	-2lnL =	%.3f	-2lnL/NDF =	%.5f	%d\n", Ea,T0_lifetime+iscantau*Lifetimestep,centroid,centroid+iscanE*binwidth,Chi2,Chi2/(Nbinsz-2),Nbinsz-2);//Chi2 is global varible, use Chi2 to extract the minimized Chi2 obtained in fcn(), and output Chi2 in main(). Use this sentence to output 1D chisquare.

			//Thomas save figure modify
			//sprintf(tname,"%s%.2f%s%d%s%.1f%s%.1f%s%d%s%s%s%d%s%d%s%d%s%d%s%d%s%d%s%.1f%s","D:/X/out/Figure_ungated/Sim25Si_T",T0_lifetime,"_Step",NumberofSteps,"_Emax",centroid+iscanE*binwidth,"_Emean",E0_gamma+iscanE*binwidth,"_feed",nBranches,"_En",pname,"SP",Flag_SP,"_Det",Flag_Det,"_Ang",Flag_Ang,"_st",Flag_st,"_sig",Flag_sig,"_tau",Flag_tau,"_binwid",binwidth,"_P10_Sun_Thomas");
			//save figure
			sprintf(tname,"%s%d%s%.1f%s%.1f","D:/X/out/Figure_DSL/SimO16_Eg",E0_gamma+iscanE*binwidth,"_tau",T0_lifetime+iscantau*Lifetimestep);
			sprintf(outtxtname,"%s%s",tname,".png");
			c1->SaveAs(outtxtname);
			//save root file
			TFile *fout = new TFile("D:/X/out/Compare.root","RECREATE");

			hSega->Write("hSega");
			hFit->Write("hFit");
			gr->Write("Residual Plot");
			hResid->Write("hResidual");
			for(int ii=0;ii<nBranches;ii++)
			{
				proton_copy[ii]->Write();
			}
			// 	p1->Write("p_1804");
			// 	p2->Write("p_1501");
			// 	p3->Write("p_1039");
			//  p4->Write("p_1056");
			//  p5->Write("p_3194");
			//  p6->Write("p_3714");
			bkgshort->Write("background");
			delete fout;
		}//for(int iscanE=0;iscanE<1;iscanE++)
		//fprintf(outfile, "\n");//DSLmodify use this sentence to output 2D chisquare matrix
	}//for(int iscantau=0;iscantau<1;iscantau++)
	//	return 0;
}//main_readhists; fcn(); comparehists(); main_draw_save;








void fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
	//npar number of free parameters involved in minimization
	//gin partial derivatives (return values) computed gradient values (optional)
	//f the function value itself (return value)
	//par parameter values
	//iflag flag word to switch between several actions of FCN	//Usually only f and par are important
{
	if (hFit) delete hFit;
	hFit = (TH1F*)proton_copy[0]->Clone("hFit");//Making a copy of histogram p1
	//hFit->Sumw2();//histogram is already filled, the sum of squares of weights is filled with the existing bin contents
	//This function is automatically called when the histogram is created

	//Set your different feeding intensities here
	for(int ii=1;ii<nBranches;ii++)
	{
		hFit->Add(proton_copy[ii],1);//p2 p3的相对intensity都用文献值比
	}
	hFit->Scale(par[0]);  //Multiply this histogram by a constant
	hFit->Add(bkgshort,par[1]);//hFit=p1+p2+p3+bkgshort
	f = CompareHists(hFit,hSega,(minrangezoom-minrange)/binwidth+1,(maxrangezoom-minrange)/binwidth);//first bin -> last bin
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
	double chisquareNeyman=0,chisquarePearson=0,likelihood=0;                                   //chi square 
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
		residual=ydata-yfit;

		if(ydata!=0) {chisquareNeyman += residual*residual/ydata;}//Neyman Chi2 method
		if(yfit!=0) {chisquarePearson += residual*residual/yfit;}//Pearson Chi2 method
		if(yfit!=0&&ydata!=0) {likelihood+=yfit-ydata+ydata*log(ydata/yfit);}
		if(yfit!=0&&ydata==0) {likelihood+=yfit-ydata;}
		//if(y[0]>0)	dy[0] = (his2->GetBinErrorLow(i)+his2->GetBinErrorUp(i))/2;
		//if(y[0]==0)	dy[0] = his2->GetBinErrorUp(i);//empty bin ErrorUp=1.8, ErrorLow=0, Error=0;
	}
	likelihood=likelihood*2;
	return likelihood;
}
