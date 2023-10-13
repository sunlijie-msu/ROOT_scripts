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
#include "TSpline.h"
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
void Concatenatetxt_DSL()// concatenate two txt files
	// To generate 4156 and 5141 joint results txt files, copy 4156's parameter files for 41565141, this script only concatenates results files.
{
	int Eg,Egbinwidth;
	int Taufirst,Taubinwidth;
	int Nbinsx_Tau, Nbinsy_E;
	double Xlow, Xup, Ylow, Yup;
	double input_value[80][80];
	double Sum_input_X[80]={0};
	double Sum_input_Y[80]={0};
	int i,j,k, ;
	double x, y, ye;
	double hight90, lowt, hight, centralt;
	TCanvas *hcanvas[30];
	TCanvas *gcanvas[30];
	TH1D *h1DX [30];
	TGraph *g1DX[30];
	double scale_Y[6000]={0};
	double scale_X[6000]={0};
	char paraprint[300], h_name[300],g_name[300],b_name[300];
	char results_dat_name1[300], results_dat_name2[300], results_joint_dat_name[300];
	char parameters_dat_name1[300], parameters_dat_name2[300], parameters_joint_dat_name[300];
	int number_run;

	int Eg_flag=41565141;
//	int Eg_flag=21863435; //modify choose one Eg
	int un_flag=1; //modify choose 0 or 1 without/with uncertainty
	
	if (Eg_flag==41565141)	{	number_run = 890; }
	if (Eg_flag==21863435)	{	number_run = 1295; }
	for (int irun=0; irun<=number_run; irun++)
	{
		if (Eg_flag==41565141&&un_flag==0)
		{
			sprintf(results_dat_name1,"%s%04d%s","D:/X/out/DSL_MADAI/Eg4156/run",irun,"/results.dat");
			sprintf(results_dat_name2,"%s%04d%s","D:/X/out/DSL_MADAI/Eg5141/run",irun,"/results.dat");
			sprintf(results_joint_dat_name,"%s%04d%s","D:/X/out/DSL_MADAI/Eg41565141/run",irun,"/results.dat");
		}
		if (Eg_flag==41565141&&un_flag==1)
		{
			sprintf(results_dat_name1,"%s%04d%s","D:/X/out/DSL_MADAI_model_uncertainty/Eg4156/run",irun,"/results.dat");
			sprintf(results_dat_name2,"%s%04d%s","D:/X/out/DSL_MADAI_model_uncertainty/Eg5141/run",irun,"/results.dat");
			sprintf(results_joint_dat_name,"%s%04d%s","D:/X/out/DSL_MADAI_model_uncertainty/Eg41565141/run",irun,"/results.dat");
		}
		if (Eg_flag==21863435&&un_flag==0)
		{
			sprintf(results_dat_name1,"%s%04d%s","D:/X/out/DSL_MADAI/Eg2186/run",irun,"/results.dat");
			sprintf(results_dat_name2,"%s%04d%s","D:/X/out/DSL_MADAI/Eg3435/run",irun,"/results.dat");
			sprintf(results_joint_dat_name,"%s%04d%s","D:/X/out/DSL_MADAI/Eg21863435/run",irun,"/results.dat");
		}
		if (Eg_flag==21863435&&un_flag==1)
		{
			sprintf(results_dat_name1,"%s%04d%s","D:/X/out/DSL_MADAI_model_uncertainty/Eg2186/run",irun,"/results.dat");
			sprintf(results_dat_name2,"%s%04d%s","D:/X/out/DSL_MADAI_model_uncertainty/Eg3435/run",irun,"/results.dat");
			sprintf(results_joint_dat_name,"%s%04d%s","D:/X/out/DSL_MADAI_model_uncertainty/Eg21863435/run",irun,"/results.dat");
		}


		ofstream outfiler(results_joint_dat_name,ios::out);

		int iline=0;
		fstream InputByChar_infiler1;
		InputByChar_infiler1.open(results_dat_name1,ios::in);
		while(!InputByChar_infiler1.eof())
		{
			iline++; //cout<<iline<<endl;
			if (un_flag==0)
			{
				InputByChar_infiler1>>x>>y;  // for no uncertainty
				outfiler<<x<<'	'<<y<<endl;  // for no uncertainty
			}
			if (un_flag==1)
			{
				InputByChar_infiler1>>x>>y>>ye;  // for uncertainty
				outfiler<<x<<'	'<<y<<'	'<<ye<<endl;  // for uncertainty
			}
		}
		fstream InputByChar_infiler2;
		InputByChar_infiler2.open(results_dat_name2,ios::in);
		while(!InputByChar_infiler2.eof())
		{
			iline++; //cout<<iline<<endl;
			if (Eg_flag==41565141&&un_flag==0)
			{
				InputByChar_infiler2>>x>>y;  // for no uncertainty
				if (iline<90) outfiler<<x<<'	'<<y<<endl;  // for41565141 no uncertainty
				if (iline==90) outfiler<<x<<'	'<<y;  // for41565141 no uncertainty
			}
			if (Eg_flag==41565141&&un_flag==1)
			{
				InputByChar_infiler2>>x>>y>>ye;  // for uncertainty
				if (iline<90) outfiler<<x<<'	'<<y<<'	'<<ye<<endl;  // for41565141 uncertainty
				if (iline==90) outfiler<<x<<'	'<<y<<'	'<<ye;  // for41565141 uncertainty
			}
			if (Eg_flag==21863435&&un_flag==0)
			{
				InputByChar_infiler2>>x>>y;  // for no uncertainty
				if (iline<65) outfiler<<x<<'	'<<y<<endl;  // for21863435 no uncertainty
				if (iline==65) outfiler<<x<<'	'<<y;  // for21863435 no uncertainty
			}
			if (Eg_flag==21863435&&un_flag==1)
			{
				InputByChar_infiler2>>x>>y>>ye;  // for uncertainty
				if (iline<65) outfiler<<x<<'	'<<y<<'	'<<ye<<endl;  // for21863435 uncertainty
				if (iline==65) outfiler<<x<<'	'<<y<<'	'<<ye;  // for21863435 uncertainty
			}
		}
		cout<<irun<<endl;
// 		sprintf(parameters_dat_name1,"%s%04d%s","D:/X/out/DSL_MADAI_model_uncertainty/Eg4156/run",irun,"/parameters.dat");
// 		sprintf(parameters_dat_name2,"%s%04d%s","D:/X/out/DSL_MADAI_model_uncertainty/Eg5141/run",irun,"/parameters.dat");
// 		sprintf(parameters_joint_dat_name,"%s%04d%s","D:/X/out/DSL_MADAI_model_uncertainty/Eg41565141/run",irun,"/parameters.dat");
// 		ofstream outfilep(parameters_joint_dat_name,ios::out);
	}
}