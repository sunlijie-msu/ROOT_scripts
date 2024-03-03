#include <iostream>
#include <fstream>
#include <iomanip>
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
void peakfit_gausnerfcpol1_band_pxct() // get histogram and EMG fit one peak
{
	const int ID1 = 0;// i=ID1//which detector
	const int ID2 = 0;// i<=ID2// which detector
	double binwidth;
	double fitrange_min = 0, fitrange_max = 0;
	double histomin = 0, histomax = 0, histoNbins = 0;
	int minbin = 0, maxbin = 0;
	float Eg, Eg2nd, gaplow = 70., gaphigh = 70.;//fitting range随分辨不同调整
	int i, ii, jj, ibin;
	char paraprint[300], histo_name[300], hfit_name[300];
	TH1D* histo[ID2 + 1];//TH1D peak search+gauss fit, create histograms
	int fit_Nbins;
	const int peaknum = 50;//search peak numbers
	TGraph* graph_residual[ID2 + 1];//create graphs
	TF1* fEMG[peaknum];//create function
	TF1* p[peaknum], * g[peaknum], * b[peaknum], * p2[peaknum];
	Double_t peakx[peaknum], peakxerr[peaknum];
	Double_t peaky[peaknum], peakyerr[peaknum];
	Double_t constant[peaknum], constant_err[peaknum];
	Double_t sig[peaknum], sig_err[peaknum];
	Double_t tau[peaknum], tau_err[peaknum];
	Double_t mean[peaknum], mean_err[peaknum];
	Double_t FWHM[peaknum], FWHM_err[peaknum];
	double inflation_factor = 1.0;
	TH1D* h_confidence_interval[ID2 + 1][peaknum];
	TCanvas* canvaspeak[ID2 + 1][peaknum];
	char pathname[300];
	char filename[300];
	sprintf(pathname, "%s", "F:/e21010/pxct/");
	sprintf(filename, "%s%s", pathname, "lower_bounds.dat");
	ifstream infilelowerbounds(filename, ios::in);
	sprintf(filename, "%s%s", pathname, "upper_bounds.dat");
	ifstream infileupperbounds(filename, ios::in);
	double par[peaknum][10], par_err[peaknum][10];
	double parChi[peaknum], parNDF[peaknum], p_value[peaknum];
	double lower_search_bound[ID2 + 1][peaknum], upper_search_bound[ID2 + 1][peaknum];

	for (i = ID1; i <= ID2; i++)//which detector no need to change
	{
		for (ii = 0; ii < peaknum; ii++)//which peak
		{
			lower_search_bound[i][ii] = 0;
			upper_search_bound[i][ii] = 30000;
		}
	}
	lower_search_bound[0][0] = 5;	upper_search_bound[0][0] = 7; //55Fe 5.9 and 6.5 LEGe
	lower_search_bound[0][1] = 6;	upper_search_bound[0][1] = 7; //55Fe 6.5 LEGe
	lower_search_bound[0][2] = 119;	upper_search_bound[0][2] = 124; //152Eu 122 LEGe
	lower_search_bound[0][3] = 44;	upper_search_bound[0][3] = 46; //152Eu 45.3 and 45.4 LEGe
	lower_search_bound[0][4] = 58;	upper_search_bound[0][4] = 61; //241Am 59.5 LEGe
	lower_search_bound[0][5] = 25;	upper_search_bound[0][5] = 27; //241Am 26.3 LEGe
	lower_search_bound[0][6] = 38;	upper_search_bound[0][6] = 43; //152Eu 40.1 and 39.5 LEGe
	lower_search_bound[0][7] = 32;	upper_search_bound[0][7] = 34; //241Am 33.2 LEGe
	lower_search_bound[0][8] = 12;	upper_search_bound[0][8] = 14; //241Am 13.9 and 13.8 LEGe
	lower_search_bound[0][9] = 119;	upper_search_bound[0][9] = 125; //152Eu 122 XtRa
	lower_search_bound[0][10] = 335;	upper_search_bound[0][10] = 350; //152Eu 344 XtRa
	lower_search_bound[0][11] = 1296;	upper_search_bound[0][11] = 1302; //152Eu 1299 XtRa
	lower_search_bound[0][12] = 5440;	upper_search_bound[0][12] = 5500; //241Am 5486 MSD26
	lower_search_bound[0][13] = 11.6;	upper_search_bound[0][13] = 12.2; //241Am 11.89 LEGe
	lower_search_bound[0][14] = 1405;	upper_search_bound[0][14] = 1411; //152Eu 1408 XtRa
	lower_search_bound[0][15] = 3100;	upper_search_bound[0][15] = 3200; //148Gd 3183 MSD26
	lower_search_bound[0][16] = 241;	upper_search_bound[0][16] = 246; //152Eu 245 XtRa
	lower_search_bound[0][17] = 409;	upper_search_bound[0][17] = 415; //152Eu 411 XtRa
	lower_search_bound[0][18] = 441;	upper_search_bound[0][18] = 447; //152Eu 444 XtRa
	lower_search_bound[0][19] = 776;	upper_search_bound[0][19] = 782; //152Eu 779 XtRa
	lower_search_bound[0][20] = 863;	upper_search_bound[0][20] = 870; //152Eu 867 XtRa
	lower_search_bound[0][21] = 960;	upper_search_bound[0][21] = 968; //152Eu 964 XtRa
	lower_search_bound[0][22] = 1109;	upper_search_bound[0][22] = 1115; //152Eu 1112 XtRa
	lower_search_bound[0][23] = -1000;	upper_search_bound[0][23] = 1000; //241Am Timing LEGe-MSD
	lower_search_bound[0][24] = 483;	upper_search_bound[0][24] = 493; //152Eu 489 XtRa
	lower_search_bound[0][25] = 560;	upper_search_bound[0][25] = 568; //152Eu 564 XtRa
	lower_search_bound[0][26] = 585;	upper_search_bound[0][26] = 588; //152Eu 586 XtRa
	lower_search_bound[0][27] = 683;	upper_search_bound[0][27] = 693; //152Eu 689 XtRa
	lower_search_bound[0][28] = 1002;	upper_search_bound[0][28] = 1008; //152Eu 1005 XtRa
	lower_search_bound[0][29] = 1209;	upper_search_bound[0][29] = 1216; //152Eu 1213 XtRa
	lower_search_bound[0][30] = 1524;	upper_search_bound[0][30] = 1532; //152Eu 1528 XtRa
	lower_search_bound[0][31] = 1169;	upper_search_bound[0][31] = 1177; //60Co 1173.228 XtRa
	lower_search_bound[0][32] = 1327;	upper_search_bound[0][32] = 1337; //60Co 1332.501 XtRa
	lower_search_bound[0][33] = 2500;	upper_search_bound[0][33] = 2510; //60Co sum peak XtRa


	for (i = ID1; i <= ID2; i++)//which detector no need to change
	{
		for (ii = 0; ii < peaknum; ii++)//which peak
		{
			cout << lower_search_bound[i][ii] << "	";
			cout << upper_search_bound[i][ii] << "	";
		}
		cout << endl;
	}

	sprintf(filename, "%s%s", pathname, "peakpara.dat");
	ofstream outfile(filename, ios::out);
	sprintf(filename, "%s%s", pathname, "run0216_0217_0218_LEGe_XtRa_MSD26_60Co_I7281_inChamber_vacuum_XtRa_12mm_away_window1us_CFDdelay_adjusted_AnalogGain1.0_cal.root"); // modify
	TFile* fin = new TFile(filename);//after this statement, you can use any ROOT command1 for this rootfile
	cout << filename << endl;

	for (i = ID1; i <= ID2; i++)//which detector no need to change
	{
		sprintf(histo_name, "%s", "hsouth_e"); // modify
		histo[i] = (TH1D*)fin->Get(histo_name); //Get spectrum
		histo[i]->Rebin(1);
		histo[i]->Sumw2(kFALSE);
		histo[i]->SetBinErrorOption(TH1::kPoisson);
		binwidth = histo[i]->GetBinWidth(1);
	}

	for (i = ID1; i <= ID2; i++) // no need to change
	{
		for (ii = 33; ii <= 33; ii++)// modify which peak in one detector =0<=13
		{
			sprintf(hfit_name, "%s%d%s%d", histo_name, i, "_peak", ii);
			canvaspeak[i][ii] = new TCanvas(hfit_name, hfit_name, 1400, 740);//建立画布
			canvaspeak[i][ii]->cd();//进入画布
			// 			canvaspeak[i][ii]->SetTopMargin(0.02);
			// 			canvaspeak[i][ii]->SetRightMargin(0.03);
			// 			canvaspeak[i][ii]->SetLeftMargin(0.09);
			// 			canvaspeak[i][ii]->SetBottomMargin(0.13);
			TPad* pad1 = new TPad("pad1", "The pad 70% of the height", 0.0, 0.3, 1.0, 1.0);
			// xlow, ylow, xup, yup
			TPad* pad2 = new TPad("pad2", "The pad 30% of the height", 0.0, 0.0, 1.0, 0.3);

			// Set margins for pad1
			pad1->SetTopMargin(0.055);  // relative to pad1 height
			pad1->SetBottomMargin(0.13); // no bottom margin
			pad1->SetLeftMargin(0.09);  // relative to pad1 width
			pad1->SetRightMargin(0.025); // relative to pad1 width

			// Set margins for pad2
			pad2->SetTopMargin(0.03);    // no top margin
			pad2->SetBottomMargin(0.35); // relative to pad2 height
			pad2->SetLeftMargin(0.09);    // relative to pad2 width
			pad2->SetRightMargin(0.025);   // relative to pad2 width

			pad1->Draw();
			pad1->cd();

			//gStyle->SetOptTitle(0);

			histo[i]->SetTitle("");//图名
			histo[i]->GetXaxis()->SetTitle("Energy (keV)");//轴名
			histo[i]->GetYaxis()->SetTitle("Counts per 146 eV");// modify
			histo[i]->GetXaxis()->CenterTitle();//居中
			histo[i]->GetYaxis()->CenterTitle();//居中
			histo[i]->GetXaxis()->SetLabelFont(132);//坐标字体
			histo[i]->GetYaxis()->SetLabelFont(132);//坐标字体
			histo[i]->GetXaxis()->SetLabelSize(0.05);
			histo[i]->GetYaxis()->SetLabelSize(0.05);
			histo[i]->GetXaxis()->SetTitleFont(132);//轴名字体
			histo[i]->GetYaxis()->SetTitleFont(132);//轴名字体
			histo[i]->GetXaxis()->SetTitleOffset(1.0);//轴名偏移
			histo[i]->GetYaxis()->SetTitleOffset(0.7);//轴名偏移
			histo[i]->GetXaxis()->SetTitleSize(0.06);
			histo[i]->GetYaxis()->SetTitleSize(0.06);
			histo[i]->GetYaxis()->SetNdivisions(505);
			histo[i]->GetYaxis()->SetTickLength(0.02);
			histo[i]->SetLineWidth(2);
			histo[i]->SetStats(0);
			histo[i]->Draw("e");

			peaky[ii] = 0; peakx[ii] = 0; sig[ii] = 0; tau[ii] = 0; mean[ii] = 0; FWHM[ii] = 0;
			histo[i]->GetXaxis()->SetRangeUser(lower_search_bound[i][ii], upper_search_bound[i][ii]);
			peaky[ii] = histo[i]->GetMaximum();
			peakx[ii] = histo[i]->GetBinCenter(histo[i]->GetMaximumBin());
			cout << "	MaxX = " << peakx[ii] << "	MaxY = " << peaky[ii] << endl;
			if (ii == 0) { gaplow = 0.5; gaphigh = 0.4; Eg = 5895; Eg2nd = 6490; }//55Fe 5.9 and 6.5 LEGe
			if (ii == 2) { gaplow = 1.4; gaphigh = 1.4; Eg = 121.7817; }//152Eu 122 LEGe
			if (ii == 4) { gaplow = 1.0; gaphigh = 1.0; Eg = 59.541; }//241Am	59.5 LEGe
			if (ii == 5) { gaplow = 0.8; gaphigh = 0.8; Eg = 26.345; }//241Am 26.3 LEGe
			if (ii == 7) { gaplow = 0.7; gaphigh = 0.7; Eg = 33.196; }//241Am 33.2 LEGe
			if (ii == 9) { gaplow = 2.1; gaphigh = 2.1; Eg = 121.7817; }//152Eu 122 XtRa
			if (ii == 10) { gaplow = 2.5; gaphigh = 2.5; Eg = 344.2785; }//152Eu 344 XtRa
			if (ii == 11) { gaplow = 3.0; gaphigh = 2.9; Eg = 1299.1420; }//152Eu 1299 XtRa
			if (ii == 12) { gaplow = 120; gaphigh = 90; Eg = 5486; }//241Am 5486 MSD26
			if (ii == 13) { gaplow = 0.4; gaphigh = 0.4; Eg = 11.89; }//241Am 11.89 LEGe
			if (ii == 14) { gaplow = 3.5; gaphigh = 3.5; Eg = 1408.013; }//152Eu 1408 XtRa
			if (ii == 15) { gaplow = 100; gaphigh = 60; Eg = 3182.69; }//148Gd 3183 MSD26
			if (ii == 16) { gaplow = 2.0; gaphigh = 2.1; Eg = 244.6974; }//152Eu 245 XtRa
			if (ii == 17) { gaplow = 2.5; gaphigh = 2.5; Eg = 411.1165; }//152Eu 411 XtRa
			if (ii == 18) { gaplow = 2.5; gaphigh = 2.5; Eg = 443.9606; }//152Eu 444 XtRa
			if (ii == 19) { gaplow = 2.7; gaphigh = 2.7; Eg = 778.9045; }//152Eu 779 XtRa
			if (ii == 20) { gaplow = 2.7; gaphigh = 2.7; Eg = 867.380; }//152Eu 867 XtRa
			if (ii == 21) { gaplow = 2.9; gaphigh = 2.9; Eg = 964.057; }//152Eu 964 XtRa
			if (ii == 22) { gaplow = 3.2; gaphigh = 3.2; Eg = 1112.076; }//152Eu 1112 XtRa
			if (ii == 23) { gaplow = 300; gaphigh = 1000; Eg = 59.541; }//241Am Timing MSD-LEGe
			if (ii == 24) { gaplow = 2.9; gaphigh = 2.9; Eg = 488.6792; }//152Eu 489 XtRa
			if (ii == 25) { gaplow = 4.0; gaphigh = 6.0; Eg = 563.986; }//152Eu 564 XtRa double peak needed
			if (ii == 26) { gaplow = 6.0; gaphigh = 2.5; Eg = 586.2648; }//152Eu 586 XtRa double peak needed
			if (ii == 27) { gaplow = 2.7; gaphigh = 2.7; Eg = 688.670; }//152Eu 689 XtRa
			if (ii == 28) { gaplow = 2.9; gaphigh = 3.0; Eg = 1005.27; }//152Eu 1005 XtRa
			if (ii == 29) { gaplow = 3.2; gaphigh = 3.3; Eg = 1212.948; }//152Eu 1213 XtRa
			if (ii == 30) { gaplow = 4.0; gaphigh = 5.3; Eg = 1528.10; }//152Eu 1528 XtRa
			if (ii == 31) { gaplow = 4.0; gaphigh = 4.0; Eg = 1173.228; } //60Co 1173.228 XtRa
			if (ii == 32) { gaplow = 4.0; gaphigh = 4.0; Eg = 1332.501; } //60Co 1332.501 XtRa
			if (ii == 33) { gaplow = 4.3; gaphigh = 4.3; Eg = 2505.692; } //60Co sum peak XtRa

			histomin = histo[i]->GetXaxis()->GetXmin();
			histomax = histo[i]->GetXaxis()->GetXmax();
			histoNbins = histo[i]->GetNbinsX();
			fitrange_min = peakx[ii] - gaplow;
			fitrange_max = peakx[ii] + gaphigh;
			// cout << histomin << "	" << histomax << "	" << histoNbins << endl;
			histo[i]->GetXaxis()->SetRangeUser(fitrange_min, fitrange_max);//zoom the axis
			//cout<<"************"<<peakx[ii]<<"	"<<peaky[ii]<<endl;

			minbin = histo[i]->FindBin(fitrange_min);
			maxbin = histo[i]->FindBin(fitrange_max);
			fit_Nbins = maxbin - minbin + 1; // Don't forget the +1, or you lose the last bin!
			// cout << "minbin = " << minbin << "	maxbin = " << maxbin << "	fit_Nbins = " << fit_Nbins << " fitrange_min = " << fitrange_min << "	fitrange_max = " << fitrange_max << endl;

			ibin = minbin;
			while (histo[i]->GetBinContent(ibin) < (peaky[ii] / 2)) { ibin++;				if (ibin >= maxbin)break; }
			double sigmaguess = 2 * (peakx[ii] - histo[i]->GetBinCenter(ibin)) / 2.355;
			float highcounts = 0, lowcounts = 0;
			for (jj = 0; jj < 10; jj++)
			{
				highcounts += histo[i]->GetBinContent(maxbin - jj);
				lowcounts += histo[i]->GetBinContent(minbin + jj);
			}
			highcounts = highcounts / 10; lowcounts = lowcounts / 10;
			double aguess = (highcounts - lowcounts) / (fitrange_max - fitrange_min);
			double bguess = lowcounts - fitrange_min * aguess;
			cout << "initial guesses= " << aguess << ",	" << bguess << ",	" << sigmaguess << ",	" << peakx[ii] << ",	" << peaky[ii] << endl;
			fEMG[ii] = new TF1("fEMG", "[0]*x+[1]+[2]/2/[3]*exp(0.5*([4]*[4]/([3]*[3]))+(x-[5])/[3])*ROOT::Math::erfc(1/sqrt(2)*([4]/[3]+(x-[5])/[4]))", histomin, histomax);// Sun PRC2021 low-energy tail
			//fEMG[ii] = new TF1("fEMG", "[0]*x+[1]+[2]/2/[3]*exp(0.5*([4]*[4]/([3]*[3]))-(x-[5])/[3])*ROOT::Math::erfc(1/sqrt(2)*([4]/[3]-(x-[5])/[4]))", histomin, histomax);// Sun PRC2021 high-energy tail
			g[ii] = new TF1("g", "gausn", histomin, histomax);// The [2]-N parameter in total is equivalent to the Constant in gausn
			p[ii] = new TF1("p", "[0]*x+[1]-[0]*x-[1]+[2]/2/[3]*exp(0.5*([4]*[4]/([3]*[3]))+(x-[5])/[3])*ROOT::Math::erfc(1/sqrt(2)*([4]/[3]+(x-[5])/[4]))", histomin, histomax);//pure peak low-energy tail
			//p2[ii] = new TF1("p2", "[0]*exp(-0.5*((x-[1])/[2])^2) / (sqrt(2*3.141592654)*[2])", histomin, histomax);//pure peak2
			b[ii] = new TF1("b", "[0]*x+[1]", histomin, histomax);//pure bkg

			// 			fEMG[ii]=new TF1("total","[0]*x+[1]+sqrt(3.141592654/2)*[2]/[3]*[4]*exp(0.5*([4]*[4]/([3]*[3]))+(x-[5])/[3])*ROOT::Math::erfc(1/sqrt(2)*([4]/[3]+(x-[5])/[4]))",histomin, histomax);// Glassman PRC2019 low-energy tail
			// 			g[ii]=new TF1("g","gausn",histomin, histomax);// The [2]-N parameter in total is equivalent to the Constant in gausn
			// 			p[ii]=new TF1("p","[0]*x+[1]-[0]*x-[1]+sqrt(3.141592654/2)*[2]/[3]*[4]*exp(0.5*([4]*[4]/([3]*[3]))+(x-[5])/[3])*ROOT::Math::erfc(1/sqrt(2)*([4]/[3]+(x-[5])/[4]))",histomin, histomax);//pure peak
			// 			b[ii]=new TF1("b","[0]*x+[1]",histomin, histomax);//pure bkg
			fEMG[ii]->SetNpx(histoNbins * 10);
			g[ii]->SetNpx(histoNbins * 10);
			p[ii]->SetNpx(histoNbins * 10);
			// 			p2[ii]->SetNpx(histoNbins);
			b[ii]->SetNpx(histoNbins * 10);
			//fEMG[ii]->SetParameters(0.1,15,peaky[ii],10,10,peakx[ii]);//initial value [0]-A, [1]-B, [2]-N, [3]-τ, [4]-σ, [5]-μ
			fEMG[ii]->SetParameters(aguess, bguess, peaky[ii], 0.2, sigmaguess, peakx[ii]);//initial value [0]-A, [1]-B, [2]-N, [3]-τ, [4]-σ, [5]-μ
			// 			fEMG[ii]->SetParLimits(0,-500,500);//Bkg A
			// 			fEMG[ii]->SetParLimits(1,-50000,300000);//Bkg B
			fEMG[ii]->SetParLimits(2, 1e3, 2e4);//Constant,min,max
			fEMG[ii]->SetParLimits(3, 0.1, 1.1);//Tau
			fEMG[ii]->SetParLimits(4, 0.1, 1.1);//Sigma
			fEMG[ii]->SetParLimits(5, peakx[ii] - gaplow / 3, peakx[ii] + gaphigh / 2);//Mean
			//fEMG[ii]->SetParLimits(5, 1298.8, 1299.9);//Mean
			fEMG[ii]->SetParNames("BkgA", "BkgB", "Const*bin", "Tau", "Sigma", "Mean");
			histo[i]->Fit("fEMG", "MLE", "", fitrange_min, fitrange_max);
			TFitResultPtr Fit_result_pointer = histo[i]->Fit("fEMG", "MLES", "", fitrange_min, fitrange_max);
			//"S" means the result of the fit is returned in the TFitResultPtr
			//“E” Perform better errors estimation using the Minos technique.
			//“M” Improve fit results, by using the IMPROVE algorithm of TMinuit.
			//If the option "S" is instead used, TFitResultPtr contains the TFitResult and behaves as a smart pointer to it.
			//The fit parameters, error and chi2 (but not covariance matrix) can be retrieved also from the fitted function.

			fEMG[ii]->GetParameters(par[ii]);//二维数组的par[ii]是地址,pointer to the TF1, GetParameters的数组得是double类型Obtaining the value of parameters and saving them to par[]; 
			par_err[ii][0] = fEMG[ii]->GetParError(0);//Obtaining the error of the 1st parameter
			par_err[ii][1] = fEMG[ii]->GetParError(1);//Obtaining the error of the 2nd parameter
			par_err[ii][2] = fEMG[ii]->GetParError(2);//Obtaining the error of the 3rd parameter
			par_err[ii][3] = fEMG[ii]->GetParError(3);//Obtaining the error of the 4th parameter
			par_err[ii][4] = fEMG[ii]->GetParError(4);//Obtaining the error of the 5th parameter
			par_err[ii][5] = fEMG[ii]->GetParError(5);//Obtaining the error of the 6th parameter
			parChi[ii] = fEMG[ii]->GetChisquare();
			parNDF[ii] = fEMG[ii]->GetNDF();
			p_value[ii] = fEMG[ii]->GetProb();//This probability is not the “probability that your fit is good.” If you did many fake experiments (draw many random samples of data points from the assumed distribution (your fit function)), this is the percentage of experiments that would give χ2 values ≥ to the one you got in this experiment.
			g[ii]->SetParameters(par[ii][2], par[ii][5], par[ii][4]);//set parameters for drawing gausn
			g[ii]->SetLineColor(kGreen + 1);
			p[ii]->SetParameters(par[ii][0], par[ii][1], par[ii][2], par[ii][3], par[ii][4], par[ii][5]);//set parameters for drawing peak
			p[ii]->SetLineColor(kGreen + 1);
			// 			p2[ii]->SetParameters(par[ii][6], par[ii][7], par[ii][8]);//set parameters for drawing peak
			// 			p2[ii]->SetLineColor(4);
			b[ii]->SetParameters(par[ii][0], par[ii][1]);//set parameters for drawing bkg
			b[ii]->SetLineColor(8);
			fEMG[ii]->SetLineWidth(2);
			// The confidence band is not always properly displayed.

			// Uncertainty Band
			sprintf(filename, "%s%s%d%s%d", "h_confidence_interval", "_", i, "_", ii);
			h_confidence_interval[i][ii] = (TH1D*)histo[i]->Clone(filename);//Create a histogram to hold the confidence intervals
			TVirtualFitter* fitter = TVirtualFitter::GetFitter();//The method TVirtualFitter::GetFitter())->Get the parameters of your fitting function after having it fitted to an histogram.
			fitter->GetConfidenceIntervals(h_confidence_interval[i][ii], 0.95);//By default the intervals are inflated using the chi2/ndf value of the fit if a chi2 fit is performed
			//confidence interval for the colored band: 1σ confidence interval: P=0.683, 1σ confidence interval: P=0.95, 3σ confidence interval: P=0.997
			//h_confidence_interval will contain the CL result that you can draw on top of your fitted graph.
			//where h_confidence_interval will hold the errors and could superimpose it on the same canvas where you plot central values.
			h_confidence_interval[i][ii]->SetStats(kFALSE);
			h_confidence_interval[i][ii]->SetFillColor(kRed - 10);
			//histo[i]->SetLineColor(kGreen+1);
			//histo[i]->SetMarkerColor(kGreen+1);
			h_confidence_interval[i][ii]->Draw("e3 same"); // plot the uncertainty band
			fEMG[ii]->Draw("same");
			//g[ii]->Draw("same");
			//p[ii]->Draw("same");
			//p2[ii]->Draw("same");
			b[ii]->Draw("same");
			histo[i]->Draw("e same");

			TMatrixD cov = Fit_result_pointer->GetCovarianceMatrix();//error matrix
			TMatrixD cor = Fit_result_pointer->GetCorrelationMatrix();//parameter correlation coefficients
			cov.Print();
			cor.Print();
			double Utau_Utau = fitter->GetCovarianceMatrixElement(3, 3);//(Utau)^2
			double Usigma_Usigma = fitter->GetCovarianceMatrixElement(4, 4);//(Usigma)^2
			double Umean_Umean = fitter->GetCovarianceMatrixElement(5, 5);//(Umean)^2
			double rou_Utau_Umean = fitter->GetCovarianceMatrixElement(3, 5);//(ρ*Utau*Umean), (3,5) (5,3) doesn't matter.
			double rou_Usigma_Umean = fitter->GetCovarianceMatrixElement(4, 5);//
			// cout << rou_Usigma_Umean << endl;

			peaky[ii] = fEMG[ii]->GetMaximum(fitrange_min, fitrange_max);
			peakx[ii] = fEMG[ii]->GetMaximumX(fitrange_min, fitrange_max);
			peakxerr[ii] = sqrt(Utau_Utau + Umean_Umean + 2 * rou_Utau_Umean);

			//FWHM[ii] = par[ii][4] / par[ii][5] * 2.355 * Eg;
			//FWHM_err[ii] = 2.355 * Eg * sqrt(Usigma_Usigma / (par[ii][5] * par[ii][5]) + Umean_Umean * (-par[ii][4] / par[ii][5] / par[ii][5]) * (-par[ii][4] / par[ii][5] / par[ii][5]) + 2 * rou_Usigma_Umean / par[ii][5] * (-par[ii][4] / par[ii][5] / par[ii][5]));

// 			inflation_factor = sqrt(parChi[ii] / parNDF[ii]);
// 			if (inflation_factor < 1)
			inflation_factor = 1;
			constant[ii] = par[ii][2]; constant_err[ii] = par_err[ii][2] * inflation_factor;
			tau[ii] = par[ii][3]; tau_err[ii] = par_err[ii][3] * inflation_factor;
			sig[ii] = par[ii][4]; sig_err[ii] = par_err[ii][4] * inflation_factor;
			mean[ii] = par[ii][5]; mean_err[ii] = par_err[ii][5] * inflation_factor;
			FWHM[ii] = par[ii][4] * 2.355; FWHM_err[ii] = par_err[ii][4] * 2.355 * inflation_factor;

			Double_t topy, topx, lower_half, higher_half;
			topy = p[ii]->GetMaximum(fitrange_min, fitrange_max);
			topx = p[ii]->GetMaximumX(fitrange_min, fitrange_max);
			lower_half = p[ii]->GetX(topy / 2.0, fitrange_min, topx, 1E-12);
			higher_half = p[ii]->GetX(topy / 2.0, topx, fitrange_max, 1E-12);
			FWHM[ii] = higher_half - lower_half;

			outfile << histo_name << i << "	Constant*binsize" << ii << "=	" << constant[ii] << "	+/-	" << constant_err[ii] << "	Mean" << ii << "=	" << mean[ii] << "	+/-	" << mean_err[ii] << "	Maximum" << ii << "=	" << peakx[ii] << "	+/-	" << peakxerr[ii] << "	Sigma" << ii << "=	" << sig[ii] << "	+/-	" << sig_err[ii] << "	Tau" << ii << "=	" << tau[ii] << "	+/-	" << tau_err[ii] << "	A" << ii << "=	" << par[ii][0] << "	+/-	" << par_err[ii][0] << "	B" << ii << "=	" << par[ii][1] << "	+/-	" << par_err[ii][1] << "	Chi2" << ii << "=	" << parChi[ii] << "	NDF" << ii << "=	" << parNDF[ii] << "	Area" << ii << "=	" << par[ii][2] / binwidth << "	FWHM" << ii << "=	" << FWHM[ii] << "	+/-	" << FWHM_err[ii] << endl;

			TPaveText* textgaus = new TPaveText(0.7, 0.4, 0.975, 0.945, "brNDC");//加标注left, down, right, up
			textgaus->SetBorderSize(1);//边框宽度
			textgaus->SetFillColor(0);//填充颜色
			textgaus->SetTextAlign(12);//align = 10*HorizontalAlign + VerticalAlign, 12 means水平左对齐、垂直居中对齐
			textgaus->SetTextFont(132);//font = 10 * fontID + precision, 12+2,12 means Symbol; 13+2, 13 means News Times Roman
			//text->SetTextColor(2);//文本颜色
			sprintf(paraprint, "Constant*binsize%d=%.4f%s%.4f", ii, par[ii][2], "+/-", par_err[ii][2]);//par数组还保持着刚才的参数
			textgaus->AddText(paraprint);
			sprintf(paraprint, "Mean%d=%.6f%s%.6f", ii, par[ii][5], "+/-", par_err[ii][5]);
			textgaus->AddText(paraprint);
			sprintf(paraprint, "Sigma%d=%.6f%s%.6f", ii, par[ii][4], "+/-", par_err[ii][4]);
			textgaus->AddText(paraprint);
			sprintf(paraprint, "FWHM%d=%.6f%s%.6f", ii, FWHM[ii], "+/-", FWHM_err[ii]);
			textgaus->AddText(paraprint);
			sprintf(paraprint, "Res%d=%.4f%%", ii, FWHM[ii] / par[ii][5] * 100);
			textgaus->AddText(paraprint);
			sprintf(paraprint, "Tau%d=%.6f%s%.6f", ii, par[ii][3], "+/-", par_err[ii][3]);
			textgaus->AddText(paraprint);
			sprintf(paraprint, "A%d(%.4f%s%.4f)*x+B%d(%.3f%s%.3f)", ii, par[ii][0], "+/-", par_err[ii][0], ii, par[ii][1], "+/-", par_err[ii][1]);
			textgaus->AddText(paraprint);
			sprintf(paraprint, "Chisquare%d=%.6f", ii, parChi[ii]);
			textgaus->AddText(paraprint);
			sprintf(paraprint, "NDF%d=%.1f", ii, parNDF[ii]);
			textgaus->AddText(paraprint);
			sprintf(paraprint, "p-val=%e", p_value[ii]);
			textgaus->AddText(paraprint);
			sprintf(paraprint, "Maximum%d=%.6f%s%.6f", ii, peakx[ii], "+/-", peakxerr[ii]);
			textgaus->AddText(paraprint);
			textgaus->Draw();

			// for residuals plot
			double x_values[fit_Nbins], y_values[fit_Nbins], y_fit_values[fit_Nbins], y_residuals[fit_Nbins], x_errors[fit_Nbins], y_errors[fit_Nbins];

			for (ibin = 1; ibin <= fit_Nbins; ibin++)
			{ // data array starts from [0], bin starts from [1]
				x_values[ibin - 1] = histo[i]->GetBinCenter(minbin + ibin - 1);
				// cout << x_values[ibin - 1] << " " << minbin+ibin-1 << endl;
				x_errors[ibin - 1] = binwidth / 2.0;
				y_values[ibin - 1] = histo[i]->GetBinContent(minbin + ibin - 1);
				y_errors[ibin - 1] = (histo[i]->GetBinErrorUp(minbin + ibin - 1) + histo[i]->GetBinErrorLow(minbin + ibin - 1)) / 2.0;

				y_fit_values[ibin - 1] = fEMG[ii]->Eval(x_values[ibin - 1]);
				y_residuals[ibin - 1] = y_fit_values[ibin - 1] - y_values[ibin - 1];
			}

			graph_residual[i] = new TGraphErrors(fit_Nbins, x_values, y_residuals, x_errors, y_errors); //TGraph(n,x,y,ex,ey);
			graph_residual[i]->SetTitle("");//图名
			graph_residual[i]->GetXaxis()->SetTitle("Energy (keV)");//轴名
			graph_residual[i]->GetYaxis()->SetTitle("Fit - Data");//轴名
			graph_residual[i]->GetXaxis()->CenterTitle();//居中
			graph_residual[i]->GetYaxis()->CenterTitle();//居中
			graph_residual[i]->GetXaxis()->SetLabelFont(132);//坐标字体
			graph_residual[i]->GetYaxis()->SetLabelFont(132);//坐标字体
			graph_residual[i]->GetXaxis()->SetLabelSize(0.11);
			graph_residual[i]->GetYaxis()->SetLabelSize(0.11);
			graph_residual[i]->GetXaxis()->SetTitleFont(132);//轴名字体
			graph_residual[i]->GetYaxis()->SetTitleFont(132);//轴名字体
			graph_residual[i]->GetXaxis()->SetTitleOffset(1.0);//轴名偏移
			graph_residual[i]->GetYaxis()->SetTitleOffset(0.3);//轴名偏移
			graph_residual[i]->GetXaxis()->SetTitleSize(0.14);
			graph_residual[i]->GetYaxis()->SetTitleSize(0.14);
			graph_residual[i]->GetYaxis()->SetNdivisions(505);
			graph_residual[i]->GetYaxis()->SetTickLength(0.02);
			graph_residual[i]->SetStats(0);
			graph_residual[i]->GetXaxis()->SetRangeUser(fitrange_min, fitrange_max);
			// graph_residual[i]->GetYaxis()->SetRangeUser(-50, 50); 
			graph_residual[i]->SetLineWidth(2);
			graph_residual[i]->SetLineColor(kBlue + 2);
			graph_residual[i]->SetMarkerStyle(6);
			graph_residual[i]->SetMarkerColor(kBlue + 2);
			canvaspeak[i][ii]->cd();//进入画布
			pad2->Draw();
			pad2->cd();
			graph_residual[i]->Draw("APZ");//"A": Axis are drawn around the graph, "P": The current marker is plotted at each point, "Z": Do not draw small horizontal and vertical lines the end of the error bars. Without "Z", the default is to draw these.
			TLine* T1 = new TLine(fitrange_min, 0, fitrange_max, 0);
			T1->Draw("R");//"R" means the line is drawn with the current line attributes

			sprintf(filename, "%s%s%s", pathname, hfit_name, ".png");
			canvaspeak[i][ii]->SaveAs(filename);

		}//for(ii=0;ii<peaknum;ii++)
		outfile << "\n\n" << endl;
	}//for (i=0;i<ID;i++)
}//peakcali main


