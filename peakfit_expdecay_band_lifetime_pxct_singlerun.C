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
void peakfit_expdecay_band_lifetime_pxct_singlerun() // get histogram and exponential fit a time histogram fit Run0079 for figure 237Np half-life PXCT PRC paper
{
	double binwidth = 1;
	double fitrange_min = 0, fitrange_max = 0;
	double histomin = 0, histomax = 0, histoNbins = 0;
	int minbin = 0, maxbin = 0;
	float Eg, Eg2nd, gaplow = 70., gaphigh = 70.;//fitting range随分辨不同调整
	int i, ii, jj, ibin;
	char paraprint[100], histo_name[100], hfit_name[100];
	TH1D* histo[100];//TH1D peak search+gauss fit, create histograms
	int fit_Nbins;
	const int peaknum = 50;//search peak numbers
	TGraph* graph_residual[100];//create graphs
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
	TH1D* h_confidence_interval[200][peaknum];
	TCanvas* canvaspeak[100][peaknum];
	char pathname[150];
	char filename[150];
	sprintf(pathname, "%s", "F:/e21010/pxct/");
	// 	sprintf(filename, "%s%s", pathname, "lower_bounds.dat");
	// 	ifstream infilelowerbounds(filename, ios::in);
	// 	sprintf(filename, "%s%s", pathname, "upper_bounds.dat");
	// 	ifstream infileupperbounds(filename, ios::in);
	double par[peaknum][10], par_err[peaknum][10];
	double parChi[peaknum], parNDF[peaknum], p_value[peaknum];
	TFile* fin;

	sprintf(filename, "%s%s", pathname, "peakpara.dat");
	ofstream outfile(filename, ios::app);
	int Ea_central = 0, msd_e_cut_low = 0, msd_e_cut_high = 0, Ea_gate_start = 0, Ea_gate_end = 0;
	int bin_start_low = 0;
	int bin_start_high = 0;
	int Which_Dataset = 1; //  1 for MSDtotal; 2 for MSD26;
	int Which_MSD;
	if (Which_Dataset == 1)
	{
		Which_MSD = 26; // Modify: 12 for MSD12; 26 for MSD26;
		Ea_central = 5421; // 5421 for MSDtotal, based on LISE++ calculation
		if (Which_MSD == 12)	bin_start_low = 239; // Don't change
		if (Which_MSD == 26)	bin_start_low = 230; // Don't change
		Ea_gate_start = 80; // 3; 3 means +/-3 keV = 6 keV; 20 means +/-20 keV = 40 keV, which is good for MSD26
		Ea_gate_end = 80; // 30; Keep end - start <= 4, due to Windows OS limitation
	}
	if (Which_Dataset == 2)
	{
		Which_MSD = 26; // 26 for MSD26;
		Ea_central = 5479; // 5479 for MSD26; 5421 for MSDtotal, based on LISE++ calculation
		Ea_gate_start = 20; // 3; 3 means +/-3 keV = 6 keV; 20 means +/-20 keV = 40 keV, which is good for MSD26
		Ea_gate_end = 20; // 30; Keep end - start <= 4, due to Windows OS limitation
	}

	for (i = Ea_gate_start; i <= Ea_gate_end; i+=10) // read in many root files with different Ea gates
	{
		//if (i >= 96 && i <= 99) continue;
		msd_e_cut_low = Ea_central - i;
		msd_e_cut_high = Ea_central + i;
		//if (Which_Dataset == 1) sprintf(filename, "%s%s%04d%s%04d%s", pathname, "timing_msdtotal_e_", msd_e_cut_low, "_", msd_e_cut_high, "_msdtotal_tootoot");
		//if (Which_Dataset == 2) sprintf(filename, "%s%s%04d%s%04d%s", pathname, "timing_msd26_e_", msd_e_cut_low, "_", msd_e_cut_high, "_msd26_t.root");
		//if (Which_Dataset == 1) sprintf(filename, "%s%s", pathname, "run0079_LEGe_MSD_241Am_inChamber_window1.5us_CFDdelay_adjusted_cal_5400_5560_forPRC.root");
		if (Which_Dataset == 1) sprintf(filename, "%s%s", pathname, "run0330_LEGe_MSD_241Am_Z7117_ChamberCenter_window1.5us_TrigRise0.064_0.016_0.016us_TrigGap0.952_1.000_1.000us_Th350_2700_1000_CFDDelay0.304us_Scale7_cal.root");
		fin = new TFile(filename);//after this statement, you can use any ROOT command1 for this rootfile
		cout << filename << endl;
		//sprintf(histo_name, "%s%d%s", "htiming_lege_msd", Which_MSD, "_bin1ns");
		sprintf(histo_name, "%s%d", "htiming_lege_msd", Which_MSD);
		//sprintf(histo_name, "%s", "h_decay_exponential_GetRandom"); // Fake test
		histo[i] = (TH1D*)fin->Get(histo_name); //Get spectrum
		histo[i]->Rebin(1);
		histo[i]->Sumw2(kFALSE);
		histo[i]->SetBinErrorOption(TH1::kPoisson);
		for (ii = 0; ii <= 0; ii++) // modify: ii <=39, fit one histogram with many different fit start points = ii * 20 + bin_start_low
		{
			sprintf(hfit_name, "%s%s%d%s%d", histo_name, "_Ea", i, "_Fitrange", ii);
			canvaspeak[i][ii] = new TCanvas(hfit_name, hfit_name, 1300, 700);//建立画布
			canvaspeak[i][ii]->cd();//进入画布
			// 			canvaspeak[i][ii]->SetTopMargin(0.02);
			// 			canvaspeak[i][ii]->SetRightMargin(0.03);
			// 			canvaspeak[i][ii]->SetLeftMargin(0.09);
			// 			canvaspeak[i][ii]->SetBottomMargin(0.13);
			TPad* pad1 = new TPad("pad1", "The pad 70% of the height", 0.0, 0.42, 1.0, 1.0);
			// xlow, ylow, xup, yup
			TPad* pad2 = new TPad("pad2", "The pad 30% of the height", 0.0, 0.0, 1.0, 0.42);

			// Set margins for pad1
			pad1->SetTopMargin(0.04);  // relative to pad1 height
			pad1->SetBottomMargin(0.05); // no bottom margin
			pad1->SetLeftMargin(0.12);  // relative to pad1 width
			pad1->SetRightMargin(0.02); // relative to pad1 width

			// Set margins for pad2
			pad2->SetTopMargin(0.05);    // no top margin
			pad2->SetBottomMargin(0.4); // relative to pad2 height
			pad2->SetLeftMargin(0.12);    // relative to pad2 width
			pad2->SetRightMargin(0.02);   // relative to pad2 width

			pad1->Draw();
			pad1->cd();
			pad1->SetFrameLineWidth(3);
			//gStyle->SetOptTitle(0);

			histo[i]->SetTitle("");//图名
			histo[i]->GetXaxis()->SetTitle("");// Time difference LEGe - MSD26 (ns)
			histo[i]->GetYaxis()->SetTitle("Counts per 1 ns");// modify
			histo[i]->GetXaxis()->CenterTitle();//居中
			histo[i]->GetYaxis()->CenterTitle();//居中
			histo[i]->GetXaxis()->SetLabelFont(132);//坐标字体
			histo[i]->GetYaxis()->SetLabelFont(132);//坐标字体
			histo[i]->GetXaxis()->SetLabelSize(0.135);
			histo[i]->GetYaxis()->SetLabelSize(0.135);
			histo[i]->GetXaxis()->SetLabelOffset(1.0);//轴名偏移
			histo[i]->GetXaxis()->SetTitleFont(132);//轴名字体
			histo[i]->GetYaxis()->SetTitleFont(132);//轴名字体
			histo[i]->GetXaxis()->SetTitleOffset(1.0);//轴名偏移
			histo[i]->GetYaxis()->SetTitleOffset(0.45);//轴名偏移
			histo[i]->GetXaxis()->SetTitleSize(0.135);
			histo[i]->GetYaxis()->SetTitleSize(0.135);
			histo[i]->GetYaxis()->SetNdivisions(505);
			histo[i]->GetYaxis()->SetTickLength(0.010);
			//histo[i]->GetYaxis()->SetRangeUser(0.25, 400);
			histo[i]->SetLineWidth(1);
			histo[i]->SetStats(0);
			histo[i]->Draw("e");

			peaky[ii] = 0; peakx[ii] = 0; sig[ii] = 0; tau[ii] = 0; mean[ii] = 0; FWHM[ii] = 0;
			histo[i]->GetXaxis()->SetRangeUser(0, 1500);
			//peaky[ii] = histo[i]->GetMaximum();
			//peakx[ii] = histo[i]->GetBinCenter(histo[i]->GetMaximumBin());
			//cout << "	MaxX = " << peakx[ii] << "	MaxY = " << peaky[ii] << endl;

			histomin = histo[i]->GetXaxis()->GetXmin();
			histomax = histo[i]->GetXaxis()->GetXmax();
			histoNbins = histo[i]->GetNbinsX();
			fitrange_min = bin_start_low + ii * 20; fitrange_max = bin_start_low+680; // modify // 10*T1/2 = 10*68 = 680
			//fitrange_min = 0 + ii * 20; fitrange_max = 1420; // Fake test

			// cout << histomin << "	" << histomax << "	" << histoNbins << endl;
			histo[i]->GetXaxis()->SetRangeUser(fitrange_min, fitrange_max);//zoom the axis
			//cout<<"************"<<peakx[ii]<<"	"<<peaky[ii]<<endl;

			minbin = histo[i]->FindBin(fitrange_min);
			maxbin = histo[i]->FindBin(fitrange_max);
			fit_Nbins = maxbin - minbin + 1; // Don't forget the +1, or you lose the last bin!
			// cout << "minbin = " << minbin << "	maxbin = " << maxbin << "	fit_Nbins = " << fit_Nbins << " fitrange_min = " << fitrange_min << "	fitrange_max = " << fitrange_max << endl;

			//ibin = minbin;
			//while (histo[i]->GetBinContent(ibin) < (peaky[ii] / 2)) { ibin++;				if (ibin >= maxbin)break; }
			//double sigmaguess = 2 * (peakx[ii] - histo[i]->GetBinCenter(ibin)) / 2.355;
			//float highcounts = 0, lowcounts = 0;
			//for (jj = 0; jj < 10; jj++)
			//{
			//	highcounts += histo[i]->GetBinContent(maxbin - jj);
			//	lowcounts += histo[i]->GetBinContent(minbin + jj);
			//}
			//highcounts = highcounts / 10; lowcounts = lowcounts / 10;
			//double aguess = (highcounts - lowcounts) / (fitrange_max - fitrange_min);
			//double bguess = lowcounts - fitrange_min * aguess;
			//cout << "initial guesses= " << aguess << ",	" << bguess << ",	" << sigmaguess << ",	" << peakx[ii] << ",	" << peaky[ii] << endl;

			fEMG[ii] = new TF1("fEMG", "[0]*0.693147/[1]*exp(x/(-[1]/0.693147))+[2]", histomin, histomax); // exponential decay (N, T, B)
			//fEMG[ii] = new TF1("fEMG", "[0]*exp(x/(-[1]/0.693147))+[2]", histomin, histomax); // exponential decay (A, T, B)
			// 			fEMG[ii] = new TF1("fEMG", "[0]*x+[1]+[2]/2/[3]*exp(0.5*([4]*[4]/([3]*[3]))+(x-[5])/[3])*ROOT::Math::erfc(1/sqrt(2)*([4]/[3]+(x-[5])/[4]))", histomin, histomax);// Sun PRC2021 low-energy tail
			// 			fEMG[ii] = new TF1("fEMG", "[0]*x+[1]+[2]/2/[3]*exp(0.5*([4]*[4]/([3]*[3]))-(x-[5])/[3])*ROOT::Math::erfc(1/sqrt(2)*([4]/[3]-(x-[5])/[4]))", histomin, histomax);// Sun PRC2021 high-energy tail
			//g[ii] = new TF1("g", "gausn", histomin, histomax);// The [2]-N parameter in total is equivalent to the Constant in gausn
			//p[ii] = new TF1("p", "[0]*x+[1]-[0]*x-[1]+[2]/2/[3]*exp(0.5*([4]*[4]/([3]*[3]))-(x-[5])/[3])*ROOT::Math::erfc(1/sqrt(2)*([4]/[3]-(x-[5])/[4]))", histomin, histomax);//pure peak low-energy tail
			//p2[ii] = new TF1("p2", "[0]*exp(-0.5*((x-[1])/[2])^2) / (sqrt(2*3.141592654)*[2])", histomin, histomax);//pure peak2
			//b[ii] = new TF1("b", "[0]*x+[1]", histomin, histomax);//pure bkg

			// 			fEMG[ii]=new TF1("total","[0]*x+[1]+sqrt(3.141592654/2)*[2]/[3]*[4]*exp(0.5*([4]*[4]/([3]*[3]))+(x-[5])/[3])*ROOT::Math::erfc(1/sqrt(2)*([4]/[3]+(x-[5])/[4]))",histomin, histomax);// Glassman PRC2019 low-energy tail
			// 			g[ii]=new TF1("g","gausn",histomin, histomax);// The [2]-N parameter in total is equivalent to the Constant in gausn
			// 			p[ii]=new TF1("p","[0]*x+[1]-[0]*x-[1]+sqrt(3.141592654/2)*[2]/[3]*[4]*exp(0.5*([4]*[4]/([3]*[3]))+(x-[5])/[3])*ROOT::Math::erfc(1/sqrt(2)*([4]/[3]+(x-[5])/[4]))",histomin, histomax);//pure peak
			// 			b[ii]=new TF1("b","[0]*x+[1]",histomin, histomax);//pure bkg
			fEMG[ii]->SetNpx(histoNbins * 10);
			//g[ii]->SetNpx(histoNbins * 10);
			//p[ii]->SetNpx(histoNbins * 10);
			// 			p2[ii]->SetNpx(histoNbins);
			//b[ii]->SetNpx(histoNbins * 10);
			//fEMG[ii]->SetParameters(0.1,15,peaky[ii],10,10,peakx[ii]);//initial value [0]-A, [1]-B, [2]-N, [3]-τ, [4]-σ, [5]-μ
			int Total_decays_guess = 2000 * i;
			fEMG[ii]->SetParameters(6e6, 67.8, 1);//initial value
			//fEMG[ii]->SetParameters(2e6, 68, 2);//initial value Fake test
			fEMG[ii]->SetParLimits(0, 1e4, 1e6);//N
			fEMG[ii]->SetParLimits(1, 40, 140);//T
			fEMG[ii]->SetParLimits(2, 0, 4);//B
			// 			fEMG[ii]->SetParLimits(3, 80, 120);//Tau
			// 			fEMG[ii]->SetParLimits(4, 10, 120);//Sigma
			// 			fEMG[ii]->SetParLimits(5, peakx[ii] - gaplow / 4, peakx[ii] + gaphigh / 2);//Mean
			// 			fEMG[ii]->SetParLimits(5, -100, 200);//Mean
			fEMG[ii]->SetParNames("Total_decays", "Half_life", "Background");
			//fEMG[ii]->SetParNames("BkgA", "BkgB", "Const*bin", "Tau", "Sigma", "Mean");
			histo[i]->Fit("fEMG", "MLE", "", fitrange_min, fitrange_max);
			//histo[i]->Fit("fEMG", "MLEN", "", fitrange_min, fitrange_max);
			TFitResultPtr Fit_result_pointer = histo[i]->Fit("fEMG", "MLES", "0", fitrange_min, fitrange_max);
			//"S" means the result of the fit is returned in the TFitResultPtr
			//“E” Perform better errors estimation using the Minos technique.
			//“M” Improve fit results, by using the IMPROVE algorithm of TMinuit.
			//If the option "S" is instead used, TFitResultPtr contains the TFitResult and behaves as a smart pointer to it.
			//The fit parameters, error and chi2 (but not covariance matrix) can be retrieved also from the fitted function.

			fEMG[ii]->GetParameters(par[ii]);//二维数组的par[ii]是地址,pointer to the TF1, GetParameters的数组得是double类型Obtaining the value of parameters and saving them to par[]; 
			par_err[ii][0] = fEMG[ii]->GetParError(0);//Obtaining the error of the 1st parameter
			par_err[ii][1] = fEMG[ii]->GetParError(1);//Obtaining the error of the 2nd parameter
			par_err[ii][2] = fEMG[ii]->GetParError(2);//Obtaining the error of the 3rd parameter
			parChi[ii] = fEMG[ii]->GetChisquare();
			parNDF[ii] = fEMG[ii]->GetNDF();
			p_value[ii] = fEMG[ii]->GetProb();//This probability is not the “probability that your fit is good.” If you did many fake experiments (draw many random samples of data points from the assumed distribution (your fit function)), this is the percentage of experiments that would give χ2 values ≥ to the one you got in this experiment.
			//g[ii]->SetParameters(par[ii][2], par[ii][5], par[ii][4]);//set parameters for drawing gausn
			//g[ii]->SetLineColor(kViolet+1);
			//p[ii]->SetParameters(par[ii][0], par[ii][1], par[ii][2], par[ii][3], par[ii][4], par[ii][5]);//set parameters for drawing peak
			//p[ii]->SetLineColor(kViolet+1);
			// 			p2[ii]->SetParameters(par[ii][6], par[ii][7], par[ii][8]);//set parameters for drawing peak
			// 			p2[ii]->SetLineColor(4);
			//b[ii]->SetParameters(par[ii][0], par[ii][1]);//set parameters for drawing bkg
			//b[ii]->SetLineColor(8);
			fEMG[ii]->SetLineWidth(3);
			fEMG[ii]->SetLineColor(kRed);
			// The confidence band is not always properly displayed.

			// Uncertainty Band
			sprintf(filename, "%s%s%d%s%d", "h_confidence_interval", "_", i, "_", ii);
			h_confidence_interval[i][ii] = (TH1D*)histo[i]->Clone(filename);//Create a histogram to hold the confidence intervals
			TVirtualFitter* fitter = TVirtualFitter::GetFitter();//The method TVirtualFitter::GetFitter())->Get the parameters of your fitting function after having it fitted to an histogram.
			fitter->GetConfidenceIntervals(h_confidence_interval[i][ii], 0.683);//By default the intervals are inflated using the chi2/ndf value of the fit if a chi2 fit is performed
			//confidence interval for the colored band: 1σ confidence interval: P=0.683, 1σ confidence interval: P=0.95, 3σ confidence interval: P=0.997
			//h_confidence_interval will contain the CL result that you can draw on top of your fitted graph.
			//where h_confidence_interval will hold the errors and could superimpose it on the same canvas where you plot central values.
			h_confidence_interval[i][ii]->SetStats(kFALSE);
			h_confidence_interval[i][ii]->SetFillColor(kRed - 9);
			histo[i]->SetLineColor(kBlack);
			histo[i]->SetMarkerColor(kBlack);
			h_confidence_interval[i][ii]->Draw("e3 same"); // plot the uncertainty band
			fEMG[ii]->Draw("same");
			pad1->RedrawAxis();
			//g[ii]->Draw("same");
			//p[ii]->Draw("same");
			//p2[ii]->Draw("same");
			//b[ii]->Draw("same");
			//histo[i]->Draw("e same");

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

// 			peaky[ii] = fEMG[ii]->GetMaximum(fitrange_min, fitrange_max);
// 			peakx[ii] = fEMG[ii]->GetMaximumX(fitrange_min, fitrange_max);
// 			peakxerr[ii] = sqrt(Utau_Utau + Umean_Umean + 2 * rou_Utau_Umean);

			//FWHM[ii] = par[ii][4] / par[ii][5] * 2.355 * Eg;
			//FWHM_err[ii] = 2.355 * Eg * sqrt(Usigma_Usigma / (par[ii][5] * par[ii][5]) + Umean_Umean * (-par[ii][4] / par[ii][5] / par[ii][5]) * (-par[ii][4] / par[ii][5] / par[ii][5]) + 2 * rou_Usigma_Umean / par[ii][5] * (-par[ii][4] / par[ii][5] / par[ii][5]));

// 			inflation_factor = sqrt(parChi[ii] / parNDF[ii]);
// 			if (inflation_factor < 1)
// 			inflation_factor = 1;
// 			constant[ii] = par[ii][2]; constant_err[ii] = par_err[ii][2] * inflation_factor;
// 			tau[ii] = par[ii][3]; tau_err[ii] = par_err[ii][3] * inflation_factor;
// 			sig[ii] = par[ii][4]; sig_err[ii] = par_err[ii][4] * inflation_factor;
// 			mean[ii] = par[ii][5]; mean_err[ii] = par_err[ii][5] * inflation_factor;
// 			FWHM[ii] = par[ii][4] * 2.355; FWHM_err[ii] = par_err[ii][4] * 2.355 * inflation_factor;
// 
// 			Double_t topy, topx, lower_half, higher_half;
// 			topy = p[ii]->GetMaximum(fitrange_min, fitrange_max);
// 			topx = p[ii]->GetMaximumX(fitrange_min, fitrange_max);
// 			lower_half = p[ii]->GetX(topy / 2.0, fitrange_min, topx, 1E-12);
// 			higher_half = p[ii]->GetX(topy / 2.0, topx, fitrange_max, 1E-12);
// 			FWHM[ii] = higher_half - lower_half;
// 
// 			outfile << histo_name << i << "	Constant*binsize" << ii << "=	" << constant[ii] << "	+/-	" << constant_err[ii] << "	Mean" << ii << "=	" << mean[ii] << "	+/-	" << mean_err[ii] << "	Maximum" << ii << "=	" << peakx[ii] << "	+/-	" << peakxerr[ii] << "	Sigma" << ii << "=	" << sig[ii] << "	+/-	" << sig_err[ii] << "	Tau" << ii << "=	" << tau[ii] << "	+/-	" << tau_err[ii] << "	A" << ii << "=	" << par[ii][0] << "	+/-	" << par_err[ii][0] << "	B" << ii << "=	" << par[ii][1] << "	+/-	" << par_err[ii][1] << "	Chi2" << ii << "=	" << parChi[ii] << "	NDF" << ii << "=	" << parNDF[ii] << "	Area" << ii << "=	" << par[ii][2] / binwidth << "	FWHM" << ii << "=	" << FWHM[ii] << "	+/-	" << FWHM_err[ii] << endl;

			outfile << std::scientific << std::setprecision(10) << "msd_e_cut_low	" << msd_e_cut_low << "	msd_e_cut_high	" << msd_e_cut_high << "	fitrange_min	" << fitrange_min << "	fitrange_max	" << fitrange_max << "	" << hfit_name << "	Total_decays_" << ii << "=	" << par[ii][0] << "	+/-	" << par_err[ii][0] << "	Half-life_" << ii << "=	" << par[ii][1] << "	+/-	" << par_err[ii][1] << "	Background_" << ii << "=	" << par[ii][2] << "	+/-	" << par_err[ii][2] << "	Chi2_" << ii << "=	" << parChi[ii] << "	NDF_" << ii << "=	" << parNDF[ii] << "	p-val_" << ii << "=	" << p_value[ii] << endl;

			TPaveText* textgaus = new TPaveText(0.60, 0.30, 0.975, 0.95, "brNDC");//加标注left, down, right, up
			textgaus->SetBorderSize(1);//边框宽度
			textgaus->SetFillColor(0);//填充颜色
			textgaus->SetTextAlign(12);//align = 10*HorizontalAlign + VerticalAlign, 12 means水平左对齐、垂直居中对齐
			textgaus->SetTextFont(132);//font = 10 * fontID + precision, 12+2,12 means Symbol; 13+2, 13 means News Times Roman
			//text->SetTextColor(2);//文本颜色
			sprintf(paraprint, "Total decays%d=%.1f%s%.1f", ii, par[ii][0], "+/-", par_err[ii][0]);//par数组还保持着刚才的参数
			textgaus->AddText(paraprint);
			sprintf(paraprint, "Half-life%d=%.2f%s%.2f", ii, par[ii][1], "+/-", par_err[ii][1]);
			textgaus->AddText(paraprint);
			sprintf(paraprint, "Background%d=%.3f%s%.3f", ii, par[ii][2], "+/-", par_err[ii][2]);
			textgaus->AddText(paraprint);
			sprintf(paraprint, "Chisquare%d=%.1f", ii, parChi[ii]);
			textgaus->AddText(paraprint);
			sprintf(paraprint, "NDF%d=%.1f", ii, parNDF[ii]);
			textgaus->AddText(paraprint);
			sprintf(paraprint, "p-val=%.2f", p_value[ii]);
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
				y_residuals[ibin - 1] = y_values[ibin - 1] - y_fit_values[ibin - 1];
			}

			graph_residual[i] = new TGraphErrors(fit_Nbins, x_values, y_residuals, x_errors, y_errors); //TGraph(n,x,y,ex,ey);
			graph_residual[i]->SetTitle("");//图名
			sprintf(paraprint, "Time difference LEGe #minus MSD%d (ns)", Which_MSD);
			graph_residual[i]->GetXaxis()->SetTitle(paraprint);
			graph_residual[i]->GetYaxis()->SetTitle("Data #minus Fit");
			graph_residual[i]->GetXaxis()->CenterTitle();//居中
			//graph_residual[i]->GetYaxis()->CenterTitle();//居中
			graph_residual[i]->GetXaxis()->SetLabelFont(132);//坐标字体
			graph_residual[i]->GetYaxis()->SetLabelFont(132);//坐标字体
			graph_residual[i]->GetXaxis()->SetLabelSize(0.18);
			graph_residual[i]->GetYaxis()->SetLabelSize(0.18);
			graph_residual[i]->GetXaxis()->SetTitleFont(132);//轴名字体
			graph_residual[i]->GetYaxis()->SetTitleFont(132);//轴名字体
			graph_residual[i]->GetXaxis()->SetTitleOffset(1.11);//轴名偏移
			graph_residual[i]->GetYaxis()->SetTitleOffset(0.35);//轴名偏移
			graph_residual[i]->GetXaxis()->SetTitleSize(0.18);
			graph_residual[i]->GetYaxis()->SetTitleSize(0.18);
			graph_residual[i]->GetYaxis()->SetNdivisions(105);
			graph_residual[i]->GetYaxis()->SetTickLength(0.010);
			graph_residual[i]->SetStats(0);
			graph_residual[i]->GetXaxis()->SetRangeUser(fitrange_min, fitrange_max);
			graph_residual[i]->GetYaxis()->SetRangeUser(-53, 53); 
			graph_residual[i]->SetLineWidth(1);
			graph_residual[i]->SetLineColor(kBlack);
			graph_residual[i]->SetMarkerStyle(6);
			graph_residual[i]->SetMarkerColor(kBlack);
			canvaspeak[i][ii]->cd();//进入画布
			pad2->SetFrameLineWidth(3);
			pad2->Draw();
			pad2->RedrawAxis();
			pad2->cd();
			graph_residual[i]->Draw("APZ");//"A": Axis are drawn around the graph, "P": The current marker is plotted at each point, "Z": Do not draw small horizontal and vertical lines the end of the error bars. Without "Z", the default is to draw these.
			TLine* T1 = new TLine(fitrange_min, 0, fitrange_max, 0);
			T1->SetLineColor(kRed);
			T1->SetLineWidth(3);
			T1->Draw("R");//"R" means the line is drawn with the current line attributes

			pad1->SetLogy(0); // residuals are wrong if logy is turn on earlier
			sprintf(filename, "%s%s%s%s", pathname, "png/", hfit_name, ".png");
			canvaspeak[i][ii]->SaveAs(filename);

		}//for(ii=0;ii<peaknum;ii++)
		//outfile << "\n\n" << endl;
	}//for (i=0;i<ID;i++)
}//peakcali main


