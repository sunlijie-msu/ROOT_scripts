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
//main_readhists; fcn(); comparehists(); main_draw_save;//search m o d i f y to change peak, search R I B L L to change data set

TFile *fin_simu,*fin_data;
TH1F *hFit,*hSega,*hSega_data,*bkgshort,*background;//

const int Binsperkev=2; //Number of bins per keV, the only place to change binwidth, all other variables is related to Binsperkev //modify 10 for NSCL, 2 for RIBLL
const double binwidth=1.0/Binsperkev; //binwidths in units of keV
const int factor_rebin=10/Binsperkev; //simu and data Rebin factor
const double E0_gamma=4236; //E0+i*binwidth modify, a baseline E0 determined by the file name of the simuroot.root, 2752 4236 2868 1366, decide which simuroot to use// if you want to use for() loop and E0_gamma is related to iscan, E0_gamma cannot be global variable out of the main(), as loop has to be in the main().
const double Ebin_gamma=4238; //modify center for drawing figures: 2754 4238 2870 1369, don't change to other values // if you want to use for() loop and E0_gamma is related to iscan, E0_gamma cannot be global variable out of the main(), as loop has to be in the main().
const double minrange=Ebin_gamma-50; //always 100-keV window
const double maxrange=Ebin_gamma+50; //always 100-keV window
const double minrangeb=Ebin_gamma-100; //longer range to fit the background correctly
const double maxrangeb=Ebin_gamma+100; //longer range to fit the background correctly
const double minrangezoom=Ebin_gamma-45; //shorter range to calculate chi2 25 30 30 40 modify RIBLL
const double maxrangezoom=Ebin_gamma+35; //shorter range to calculate chi2 25 30 30 40
const int Nbins=(maxrange-minrange)*Binsperkev;//make sure this is an integer // there is another global Nbins at the beginning because fcn() also needs an Nbins.
const int Nbinsb=(maxrangeb-minrangeb)*Binsperkev;//make sure this is an integer // there is another global Nbins at the beginning, you'd better keep them consistent
const int Nbinsz=(maxrangezoom-minrangezoom)*Binsperkev;//make sure this is an integer // there is another global Nbins at the beginning, you'd better keep them consistent


double T_Lifetime = 59.15;//31.74;59.15;1918.78                   // This is Excited State Mean Lifetime (fs), Mean Lifetime τ=T/ln(2); half-life T=τ*ln(2) modify
//double E0_gamma=2753.0+double(iscan)*2.0/10.0;                  // Energy of Gamma-ray (keV) useless
const int nBranches=2;                  //13,19;2,3 Number of proton Branches feeding daughter modify
TH1F *proton_branch[nBranches],*proton_copy[nBranches],*maximum_branch[nBranches],*totalmax;
//Robertson 1369-keV γ 19
//	double E_CoM[nBranches] = {4.556,4.303,4.2583,3.597,3.466,3.342,3.237,3.021,2.486,2.453,2.165,1.272,0.9434,0.55,1.805,1.501,1.04,1.685,1.396};         //Center of Mass energy (MeV) in descending order
// double feed[nBranches]={1.28/216.71,3.32/216.71,100/216.71,10.86/216.71,34.5/216.71,6.57/216.71,4.15/216.71,3.74/216.71,0.96/216.71,0.4/216.71,17.2/216.71,2.26/216.71,17/216.71,2.5/216.71,6.73/216.71,2.9/216.71,1.53/216.71,0.2/216.71,0.61/216.71};		//beta feeding intensity of each proton, if the sum of all the feeding intensities is less than 1, the extra events will assume a En=0.
//Robertson 2754-keV γ 3
// double E_CoM[nBranches] = {1.805,1.501,1.040};         //Center of Mass energy (MeV) in descending order
// double feed[nBranches]={6.73/(6.73+2.90+1.53),2.90/(6.73+2.90+1.53),1.53/(6.73+2.90+1.53)};		//beta feeding intensity of each proton, if the sum of all the feeding intensities is less than 1, the extra events will assume a En=0.
//Robertson 4238-keV γ 2
 	double E_CoM[nBranches] = {1.685,1.396};         //Center of Mass energy (MeV) in descending order
 	double feed[nBranches]={0.93/(2.89+0.93),2.89/(2.89+0.93)};		//beta feeding intensity of each proton, if the sum of all the feeding intensities is less than 1, the extra events will assume a En=0.
//Robertson 2870-keV γ 2
// 	double E_CoM[nBranches] = {1.685,1.396};         //Center of Mass energy (MeV) in descending order
// 	double feed[nBranches]={0.93/(2.89+0.93),2.89/(2.89+0.93)};		//beta feeding intensity of each proton, if the sum of all the feeding intensities is less than 1, the extra events will assume a En=0.
//modify
//Thomas 1369-keV γ 13
// 	double E_CoM[nBranches] = {4.545,4.252,3.610,3.463,3.231,2.98,2.162,1.268,0.943,0.555,1.804,1.489,1.377};         //Center of Mass energy (MeV) in descending order
// 	double feed[nBranches]={6.6/208.21,100/208.21,5.9/208.21,28.1/208.21,5.4/208.21,1.7/208.21,18.1/208.21,6.1/208.21,17.1/208.21,7.2/208.21,6.1/208.21,5.0/208.21,0.91/208.21};		//beta feeding intensity of each proton, if the sum of all the feeding intensities is less than 1, the extra events will assume a En=0.
//Thomas 2754-keV γ 2
//  	double E_CoM[nBranches] = {1.804,1.489};         //Center of Mass energy (MeV) in descending order
//  	double feed[nBranches]={6.1/(5+6.1),5.0/(5+6.1)};		//beta feeding intensity of each proton, if the sum of all the feeding intensities is less than 1, the extra events will assume a En=0.
//Thomas 4238-keV γ 1
//  	double E_CoM[nBranches] = {1.377};         //Center of Mass energy (MeV) in descending order
//  	double feed[nBranches]={1};		//beta feeding intensity of each proton, if the sum of all the feeding intensities is less than 1, the extra events will assume a En=0.
//Thomas 2870-keV γ 1
// 	double E_CoM[nBranches] = {1.377};         //Center of Mass energy (MeV) in descending order
// 	double feed[nBranches]={1};		//beta feeding intensity of each proton, if the sum of all the feeding intensities is less than 1, the extra events will assume a En=0.

int Flag_SP = 1;                         //Stopping power default=1
int Flag_Det = 1;                      //Detector resolution. No f1 sampling when Flag_Det is false, will be very fast. default=1
int Flag_Ang = 1;                       //Tests random cos vs random between -1 and 1 for angular dist (false tests cos and true tests random X(-1,1) default=1
int Flag_Escape = 0;					//If take recoil escaping into consideration default=0
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

void Comparison()
{
	for(int iscan=0;iscan<13;iscan++){//portal, scan many simulation rootfile with different gamma energies 4238 i=35 RIBLL

		double fitrangelow1, fitrangelow2, fitrangehigh1, fitrangehigh2;
		if(Ebin_gamma==2754)	{fitrangelow1=minrangeb+30; fitrangelow2=minrangeb+85; fitrangehigh1=maxrangeb-85; fitrangehigh2=maxrangeb;}
		if(Ebin_gamma==4238)	{fitrangelow1=minrangeb; fitrangelow2=minrangeb+75; fitrangehigh1=maxrangeb-75; fitrangehigh2=maxrangeb;}
		if(Ebin_gamma==2870)	{fitrangelow1=minrangeb+40; fitrangelow2=minrangeb+83; fitrangehigh1=maxrangeb-80; fitrangehigh2=maxrangeb-5;}
		if(Ebin_gamma==1369)	{fitrangelow1=minrangeb+20; fitrangelow2=minrangeb+80; fitrangehigh1=maxrangeb-80; fitrangehigh2=maxrangeb-20;}

		//readHists();
		Chi2=0;//Chi2 is global varible, use Chi2 to extract the minimized Chi2 obtained in fcn(), and output Chi2 in main(). Must clean Chi2 before each run in main().
		char simurootname[300], hname[300], ename[300], pname[300], tname[300];//tname for temporary use
		sprintf(pname,"%d_",(int)((E_CoM[0]+0.000001)*1000));
		for(int ii=1;ii<nBranches;ii++)
		{
			sprintf(ename,"%d_",(int)((E_CoM[ii]+0.000001)*1000));
			strcat(pname,ename);//for the output root file name
			//cout<<pname<<endl;
		}
		//Thomas don't change binwidth 0.1, all simulated root files are 0.1-keV bin. The png name is associated with binwidth. modify
		//	sprintf(simurootname,"%s%.2f%s%d%s%.1f%s%d%s%s%s%d%s%d%s%d%s%.1f%s","D:/X/out/Thomas/Sim25Si_T",T_Lifetime,"_Step",NumberofSteps,"_Eg",E0_gamma,"_feed",nBranches,"_En",pname,"SP",Flag_SP,"_Det",Flag_Det,"_Ang",Flag_Ang,"_binwid",0.1,"_P10_Sun_Thomas.root");
		//Robertson don't change binwidth 0.1, all simulated root files are 0.1-keV bin. The png name is associated with binwidth.
		//sprintf(simurootname,"%s%.2f%s%d%s%.1f%s%d%s%s%s%d%s%d%s%d%s%.1f%s","D:/X/out/Robertson/Sim25Si_T",T_Lifetime,"_Step",NumberofSteps,"_Eg",E0_gamma,"_feed",nBranches,"_En",pname,"SP",Flag_SP,"_Det",Flag_Det,"_Ang",Flag_Ang,"_binwid",0.1,"_P10_Sun_Robertson.root");
		//Robertson don't change binwidth 0.1, all simulated root files are 0.1-keV bin. The png name is associated with binwidth. RIBLL
		sprintf(simurootname,"%s%.2f%s%d%s%.1f%s%d%s%s%s%d%s%d%s%d%s%.1f%s","D:/X/out/Robertson/Sim25Si_T",T_Lifetime,"_Step",NumberofSteps,"_Eg",E0_gamma,"_feed",nBranches,"_En",pname,"SP",Flag_SP,"_Det",Flag_Det,"_Ang",Flag_Ang,"_binwid",0.1,"_Si_Sun_Robertson.root");
		fin_simu = new TFile(simurootname); //Rootfile with Simulation histograms
		//fin_data = new TFile("D:/X/out/Si25/100eV-bin-ungated_sun.root"); //Rootfile with Experimental Data Histogram ungated modify
		//fin_data = new TFile("D:/X/out/Si25/100eV-bin-gated_sun.root"); //Rootfile with Experimental Data Histogram gated
		fin_data = new TFile("X:/T999/Si25_0154_0345_p450_c250_t6T_pg_Doppler.root"); //Rootfile with Experimental Data Histogram gated//RIBLL
		double centroid=0;
		//Get Histograms from Simulation
		for(int ii=0;ii<nBranches;ii++)
		{
			sprintf(hname,"n%df2",int(ii));
			proton_branch[ii]=(TH1F*)fin_simu->Get(hname);
			proton_branch[ii]->Rebin(factor_rebin);//Rebin simulated histograms
			sprintf(hname,"n%df1",int(ii));
			maximum_branch[ii]=(TH1F*)fin_simu->Get(hname);
		}
		for(int ii=1;ii<nBranches;ii++)
		{
			maximum_branch[0]->Add(maximum_branch[0],maximum_branch[ii]);
		}
		centroid=maximum_branch[0]->GetMean();//what the original energy for a unbroadened gamma peak would be
		// 	sprintf(hname,"TotalGamma_AllneutronFeedings");
		// 	totalmax=(TH1F*)fin_simu->Get(hname);
		// 	totalmax->Rebin(10);
		// 	int whichbin=totalmax->GetMaximumBin();
		// 	centroid=totalmax->GetBinCenter(whichbin);
		//hSega_data= (TH1F*)fin_data->Get("SeGA_00"); //Get Gamma ray spectrum ungated modify
		//hSega_data= (TH1F*)fin_data->Get("gbCoincidences"); //Get Gamma ray spectrum gated
		hSega_data= (TH1F*)fin_data->Get("hpgtot"); //Get Gamma ray spectrum gated RIBLL
		//	hSega_data->Rebin(factor_rebin);//Rebin data histograms

		hSega=new TH1F("hSega","hSega",Nbins,minrange,maxrange);
		bkgshort=new TH1F("bkgshort","bkgshort",Nbins,minrange,maxrange);
		background=new TH1F("background","background",Nbinsb,minrangeb,maxrangeb);

		for(int ii=0;ii<nBranches;ii++)
		{
			sprintf(hname,"p%df2",int(ii));
			proton_copy[ii]=new TH1F(hname,hname,Nbins,minrange,maxrange);//for outroot
		}
		TF1 *f1=new TF1("f1","[0]+x*[1]",fitrangelow1,fitrangelow2);//range of left half of hist you want to fit the background to
		TF1 *f2=new TF1("f2","[0]+x*[1]",fitrangehigh1,fitrangehigh2);//range of right half of hist you want to fit the background to
		hSega_data->Fit("f1","qR");
		hSega_data->Fit("f2","qR");

		for(int i=1;i<=Nbinsb;i++)
		{
			int power=1;//change this to whatever you want. 3 for high statistics
			// 		TF1 *f1=new TF1("f1","[0]+x*[1]",minrangeb+50,minrangeb+80);//range of left half of hist you want to fit the background to
			// 		TF1 *f2=new TF1("f2","[0]+x*[1]",maxrangeb-80,maxrangeb-50);//range of right half of hist you want to fit the background to

			double a = (f1->GetParameter(0)+(minrangeb+(double)i/(double)Binsperkev)*f1->GetParameter(1));//Lower BG
			double c = (f2->GetParameter(0)+(minrangeb+(double)i/(double)Binsperkev)*f2->GetParameter(1));//Upper BG
			double f=((double)Nbinsb-(double)i)/(double)Nbinsb;
			double g =(double)i/(double)Nbinsb;
			double h =pow(f,power)+pow(g,power);//if power = 1, h = 1
			double b=0;
			// 		if(minrangeb+(double)i/(double)Binsperkev<minrangeb+80) b = a;
			// 		else if(minrangeb+(double)i/(double)Binsperkev>maxrangeb-80) b = c;
			// 		else b = (a*pow(f,power)+c*pow(g,power))/h;
			b = (a*pow(f,power)+c*pow(g,power))/h;

			//hSega->AddBinContent(i+1,hSega_data->GetBinContent(minrange*Binsperkev+i));//part of the original data are stored in hSega
			background->SetBinContent(i,b);//fit data, get background histogram
		}
		//Baseline E0=readin simuroot, How many bins shift towards high energy? iscan. How many keVs shift towards high energy? iscan*binwidth.
		for(int i=1;i<=Nbins;i++)
		{
			hSega->SetBinContent(i,hSega_data->GetBinContent((minrange-100)*Binsperkev+i));//part of the experimental histogram is stored in hSega//RIBLL
			//hSega->SetBinContent(i,hSega_data->GetBinContent((minrange-0)*Binsperkev+i));//part of the experimental histogram is stored in hSega//NSCL
			bkgshort->SetBinContent(i,background->GetBinContent((minrange-minrangeb)*Binsperkev+i));//part of the background histogram is stored in bkgshort
		}
		for(int i=1;i<=Nbins;i++)
		{
			for(int ii=0;ii<nBranches;ii++)
			{
				sprintf(hname,"p%df2",int(ii));
				proton_copy[ii]->SetBinContent(i+iscan,proton_branch[ii]->GetBinContent((Ebin_gamma-E0_gamma)*Binsperkev+i));//part of the simulated proton-γ histogram is stored in proton_copy
			}
		}
		// 	for(int i=1;i<=iscan;i++)
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
		Double_t vstart[2] = {0.0004,1.0};//initial guess par_Fit, par_bkg
		Double_t step[2] = {0.000001,0.01};                //fitting step

		gMin->mnparm(0, "a1", vstart[0], step[0], 0,10,ierflg);// par for simu
		//parID, parName, initialGuess, step, lowLimit, highLimit, irrelevant //0.00几，因为是从百万高统计模拟的histo scale下来的
		gMin->mnparm(1, "a2", vstart[1], step[1], 0.01,100,ierflg);// par for bkg //bkg 基本0.9-1.2之间，不会跟拟合的本底水平差别太大
		//gMin->mnparm(2, "a3", vstart[2], step[2], 0,25,ierflg);  

		arglist[0] = 2000;  //max number of fitting iterations
		arglist[1] = 1.;  //tolerance 1 means one sigma
		gMin->mnexcm("MIGRAD", arglist ,2,ierflg);  //run the minimization using MIGRAD		//minuit.mnexcm(command, argument list, argument number, error flag);		//command is a character array holding the command
		//argument list a double array as the arguments
		//argument number is the number of arguments		//error flag return value (!= 0)		//Precisely find the minimum, but you need to be careful since it assumes the inputs are close to the minimum
		Double_t amin,edm,errdef;
		Int_t nvpar,nparx,icstat;
		gMin->mnstat(amin,edm,errdef,nvpar,nparx,icstat);

		double pars[nParams],errs[nParams];
		for (int i=0;i<nParams;i++)
		{
			gMin->GetParameter(i,pars[i],errs[i]);
			printf("pars=%f errs=%f\n",pars[i],errs[i]);    
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
		pad1->SetRightMargin(0.01);
		pad1->SetLeftMargin(0.08);
		pad1->SetBottomMargin(0.08);
		//pad1->SetBorderMode(0);

		pad2->SetTopMargin(0.01);
		pad2->SetRightMargin(0.01);
		pad2->SetLeftMargin(0.08);
		pad2->SetBottomMargin(0.2);
		//pad2->SetBorderMode(0);

		pad1->Draw();
		pad2->Draw();
		pad1->cd();
		gStyle->SetOptTitle(0);

		hSega->Draw("e");
		sprintf(tname,"%s%.1f%s","Counts per ",binwidth," keV");
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
			proton_copy[ii]->Draw("sameh");//proton_copy[0]=p1_scale+bkgshort_scale
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
		TPaveText *textchi = new TPaveText(0.10,0.87,0.28,0.97,"brNDC");//left, down, right, up
		textchi->SetBorderSize(1);
		textchi->SetFillColor(0);
		textchi->SetTextAlign(12);//align = 10*HorizontalAlign + VerticalAlign, 12 means水平左对齐、垂直居中对齐
		textchi->SetTextFont(132);//font = 10 * fontID + precision, 12+2,12 means Symbol; 13+2, 13 means News Times Roman
		sprintf(paraprint,"Chi2/ndf = %.5f",Chi2/(Nbinsz-2));
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
			//		if(hSega->GetBinContent(i)!=0) Resids[i-1] = (hFit->GetBinContent(i)-hSega->GetBinContent(i))/hSega->GetBinContent(i);//(simu-data)/data //relative residual
			hResid->AddBinContent(i,Resids[i-1]);//simu-data
			ResidsUncert[i-1]= pow(hSega->GetBinContent(i),0.5);
			if(abs(Resids[i-1])<ResidsUncert[i-1]) count++;
			//cout<< abs(Resids[i])<<"         "<<ResidsUncert[i]<<endl;
		}
		//	cout<<endl<<count<<" of "<<Nbins<<" error bars go through the fit."<<endl;//count: the larger the better
		//TGraphErrors *gr;
		//gr=new TGraphErrors(Nbins,Xvals,Resids,XvalsUncert,ResidsUncert); 
		TGraph *gr=new TGraph(Nbinsr,Xvals,Resids); //TGraph *gr1=new TGraph(n,x,y);建立曲线图、散点图//Xvals and Resids are data array, start from [0], you cannot store in Xvals and Resids starting from [1]
		gr->GetXaxis()->SetTitle("Energy (keV)");
		sprintf(tname,"#splitline{   Residual}{counts/%0.1f keV}",binwidth);//residual
		//	sprintf(tname,"#splitline{Relative residual}{counts/%0.1f keV}",binwidth);//relative residual
		gr->GetYaxis()->SetTitle(tname);
		gr->GetXaxis()->CenterTitle();
		gr->GetYaxis()->CenterTitle();
		gr->GetXaxis()->SetTitleOffset(1.0);
		gr->GetYaxis()->SetTitleOffset(0.4);
		gr->GetXaxis()->SetTitleSize(0.08);
		gr->GetXaxis()->SetLabelOffset(0.015);
		gr->GetXaxis()->SetLabelSize(0.08);
		gr->GetYaxis()->SetLabelSize(0.08);
		gr->GetYaxis()->SetTitleSize(0.08);
		gr->GetXaxis()->SetNdivisions(520);//n = n1 + 100*n2 + 10000*n3
		// 	gr->GetXaxis()->SetNdivisions(10,10,1);
		gr->GetYaxis()->SetNdivisions(505);
		gr->GetYaxis()->SetLabelSize(0.06);
		gr->GetXaxis()->SetLabelSize(0.06);
		gr->GetXaxis()->SetRangeUser(minrangezoom,maxrangezoom);
		gr->SetMarkerStyle(3);
		gr->Draw("AP same");
		TLine *T1=new TLine(minrangezoom,0,maxrangezoom,0);
		cout<<"*********	"<<E0_gamma+iscan*binwidth<<"	"<<centroid+iscan*binwidth<<endl;
		T1->Draw("R");
		//save text
		FILE *outfile=fopen ("D:/X/out/Compare.dat","a");//ofstream报错

		fprintf(outfile, "Emax=	%.2f	Emean=	%.1f	Chi2 =	%.3f	Chi2/NDF =	%.5f	%d\n", centroid+iscan*binwidth,E0_gamma+iscan*binwidth,Chi2,Chi2/(Nbinsz-2),Nbinsz-2);//Chi2 is global varible, use Chi2 to extract the minimized Chi2 obtained in fcn(), and output Chi2 in main().
		//save figure
		//Thomas modify
		//	sprintf(tname,"%s%.2f%s%d%s%.1f%s%.1f%s%d%s%s%s%d%s%d%s%d%s%.1f%s","D:/X/out/Figure_gated/Sim25Si_T",T_Lifetime,"_Step",NumberofSteps,"_Emax",centroid+iscan*binwidth,"_Emean",E0_gamma+iscan*binwidth,"_feed",nBranches,"_En",pname,"SP",Flag_SP,"_Det",Flag_Det,"_Ang",Flag_Ang,"_binwid",binwidth,"_P10_Sun_Thomas.png");
		//Robertson
		//sprintf(tname,"%s%.2f%s%d%s%.1f%s%.1f%s%d%s%s%s%d%s%d%s%d%s%.1f%s","D:/X/out/Figure_gated/Sim25Si_T",T_Lifetime,"_Step",NumberofSteps,"_Emax",centroid+iscan*binwidth,"_Emean",E0_gamma+iscan*binwidth,"_feed",nBranches,"_En",pname,"SP",Flag_SP,"_Det",Flag_Det,"_Ang",Flag_Ang,"_binwid",binwidth,"_P10_Sun_Robertson.png");
		//Robertson RIBLL
		sprintf(tname,"%s%.2f%s%d%s%.1f%s%.1f%s%d%s%s%s%d%s%d%s%d%s%.1f%s","D:/X/out/Figure_gated_RIBLL/Sim25Si_T",T_Lifetime,"_Step",NumberofSteps,"_Emax",centroid+iscan*binwidth,"_Emean",E0_gamma+iscan*binwidth,"_feed",nBranches,"_En",pname,"SP",Flag_SP,"_Det",Flag_Det,"_Ang",Flag_Ang,"_binwid",binwidth,"_Si_Sun_Robertson.png");
		c1->SaveAs(tname);
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
	}//for(int iscan=0;iscan<1;iscan++)
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
	hFit->Sumw2();//histogram is already filled, the sum of squares of weights is filled with the existing bin contents
	//This function is automatically called when the histogram is created

	//Set your different feeding intensities here
	for(int ii=1;ii<nBranches;ii++)
	{
		hFit->Add(proton_copy[ii],1);//p2 p3的相对intensity都用文献值比
	}
	hFit->Scale(par[0]);  //Multiply this histogram by a constant
	hFit->Add(bkgshort,par[1]);//hFit=p1+p2+p3+bkgshort
	f = CompareHists(hFit,hSega,(minrangezoom-minrange)*Binsperkev+1,(maxrangezoom-minrange)*Binsperkev);//first bin -> last bin
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
	double chiSquare = 0;                                   //chi square 
	int binN = his1->GetNbinsX();                         //get number of bins in his1
	//   return 0;

	//if not identical to his2, exit with error.
	if (binN!=his2->GetNbinsX())
	{
		printf("ERROR <CompareHists>: Different number of bins in histograms!%d\t%d\n",binN,his2->GetNbinsX());
		exit(1);
	}
	double y[2],dy[2];                                         //bin information
	double delta,err;                                            //bin difference and relevant uncertainty
	for (int i=binMin;i<=binMax;i++)                         //loop over all bins   
	{
		//if(i<251||i>301){
		//get bin information
		y[0] = his2->GetBinContent(i); //hSega
		//if(y[0]>50) dy[0] = his2->GetBinError(i);
		if(y[0]>0)	dy[0] = (his2->GetBinErrorLow(i)+his2->GetBinErrorUp(i))/2;
		if(y[0]==0)	dy[0] = his2->GetBinErrorUp(i);//empty bin ErrorUp=1.8, ErrorLow=0, Error=0;

		y[1] = his1->GetBinContent(i); //hFit
		//dy[1] = his2->GetBinError(i);//useless
		delta = y[1]-y[0];
		//err = sqrt(dy[0]*dy[0]+dy[1]*dy[1]);   //common uncertainty
		err = dy[0];   //common uncertainty
		if(err!=0)chiSquare += delta*delta/err/err;
		//}      //add to chi square
	} //for 
	return chiSquare;
}
