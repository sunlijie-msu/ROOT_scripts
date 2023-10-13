#include <iostream>
#include <fstream>
//#include <iomanip.h>
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
void peakfit_gausnerfcpol1_band_doublepeak_dsl2() // get histogram and EMG fit two peaks
{
	const int ID1 = 0;// i=ID1//which detector
	const int ID2 = 0;// i<=ID2// modify which detector
	double binwidth = 2;
	double fitrange_min = 0, fitrange_max = 0;
	double histomin = 0, histomax = 0, histoNbins = 0;
	int minbin = 0, maxbin = 0;
	float Eg, Eg2nd, gaplow = 70., gaphigh = 70.;//fitting range随分辨不同调整
	int i, ii, jj, ibin;
	char paraprint[100], histo_name[200], hfit_name[200];
	TH1D* histo[ID2 + 1];//TH1D peak search+gauss fit, create histograms
	int fit_Nbins;
	const int peaknum = 1;//search peak numbers, modify
	TGraph* graph_residual[ID2 + 1];//create graphs
	TF1* fEMG[peaknum];//create function
	TF1* p[peaknum], * p2[peaknum], * g[peaknum], * b[peaknum], * g2[peaknum];
	Double_t peakx[peaknum], peakxerr[peaknum];
	Double_t peakx2nd[peaknum], peakx2nderr[peaknum];
	Double_t peaky[peaknum], peakyerr[peaknum];
	Double_t constant[peaknum], constant_err[peaknum];
	Double_t sig[peaknum], sig_err[peaknum];
	Double_t tau[peaknum], tau_err[peaknum];
	Double_t mean[peaknum], mean_err[peaknum];
	Double_t FWHM[peaknum], FWHM_err[peaknum];
	Double_t constant2nd[peaknum], constant2nd_err[peaknum];
	Double_t sig2nd[peaknum], sig2nd_err[peaknum];
	Double_t tau2nd[peaknum], tau2nd_err[peaknum];
	Double_t mean2nd[peaknum], mean2nd_err[peaknum];
	Double_t FWHM2nd[peaknum], FWHM2nd_err[peaknum];
	double inflation_factor = 1.0;
	TH1D* h_confidence_interval[ID2 + 1][peaknum];
	TCanvas* canvaspeak[ID2 + 1][peaknum];
	char pathname[150];
	char filename[150];
	sprintf(pathname, "%s", "F:/mid_files/");
	sprintf(filename, "%s%s", pathname, "lower_bounds.dat");
	ifstream infilelowerbounds(filename, ios::in);
	sprintf(filename, "%s%s", pathname, "upper_bounds.dat");
	ifstream infileupperbounds(filename, ios::in);
	double par[peaknum][10], par_err[peaknum][10];
	double parChi[peaknum], parNDF[peaknum], p_value[peaknum];
	int lower_search_bound[ID2 + 1][peaknum], upper_search_bound[ID2 + 1][peaknum];
	for (i = ID1; i <= ID2; i++)//which detector no need to modify
	{
		for (ii = 0; ii < peaknum; ii++)//which peak
		{
			infilelowerbounds >> lower_search_bound[i][ii];
			infileupperbounds >> upper_search_bound[i][ii];
		}
	}
	for (i = ID1; i <= ID2; i++)//which detector no need to modify
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
	sprintf(filename, "%s%s", pathname, "Sum-tree_00370.root");
	TFile* fin = new TFile(filename);//after this statement, you can use any ROOT command1 for this rootfile
	cout << filename << endl;

	for (i = ID1; i <= ID2; i++)//which detector no need to modify
	{
		sprintf(histo_name, "%s", "h9");
		histo[i] = (TH1D*)fin->Get(histo_name); //Get spectrum
		histo[i]->Rebin(1);
		histo[i]->Sumw2(kFALSE);
		histo[i]->SetBinErrorOption(TH1::kPoisson);
	}

	for (i = ID1; i <= ID2; i++) // no need to modify
	{
		for (ii = 0; ii <= 0; ii++)// modify which peak in one detector =0<=13
		{
			sprintf(hfit_name, "%s%d%s%d", "LEGe", i, "_peak", ii);
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
			pad1->SetTopMargin(0.02);  // relative to pad1 height
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
			histo[i]->GetXaxis()->SetTitle("ADC channel");//轴名
			histo[i]->GetYaxis()->SetTitle("Counts per 2 channels");//轴名
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
			histo[i]->SetLineWidth(2);
			histo[i]->SetStats(0);
			histo[i]->Draw("e");

			peaky[ii] = 0; peakx[ii] = 0; sig[ii] = 0; tau[ii] = 0; mean[ii] = 0; FWHM[ii] = 0; FWHM2nd[ii] = 0;
			histo[i]->GetXaxis()->SetRangeUser(lower_search_bound[i][ii], upper_search_bound[i][ii]);
			peaky[ii] = histo[i]->GetMaximum();
			peakx[ii] = histo[i]->GetBinCenter(histo[i]->GetMaximumBin());
			cout << peakx[ii] << "	" << peaky[ii] << endl;
			if (ii == 0) { gaplow = 70; gaphigh = 30; Eg = 5763; Eg2nd = 5806;}//55Fe 5.9 and 6.5
			if (ii == 3) {
				gaplow = 1.05; gaphigh = 0.75; Eg = 45294; Eg2nd = 45413; }//152Eu 45.3 and 45.4
			if (ii == 6) { gaplow = 1.5; gaphigh = 0.9; Eg = 40118; Eg2nd = 39522; }//152Eu 40.1 and 39.5
			if (ii == 8) { gaplow = 0.6; gaphigh = 0.7; Eg = 13946; Eg2nd = 13761; }//241Am 13.9 and 13.8


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
			while (histo[i]->GetBinContent(ibin) < (peaky[ii] / 2))
			{
				ibin++;
				if (ibin >= maxbin)break;
			}
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
			cout << aguess << "	" << bguess << "	" << sigmaguess << "	" << peakx[ii] << "	" << peaky[ii] << endl;
			fEMG[ii] = new TF1("fEMG", "[0]*x+[1]+[2]/2/[3]*exp(0.5*([4]*[4]/([3]*[3]))+(x-[5])/[3])*ROOT::Math::erfc(1/sqrt(2)*([4]/[3]+(x-[5])/[4]))+[6]/2/[7]*exp(0.5*([8]*[8]/([7]*[7]))+(x-[9])/[7])*ROOT::Math::erfc(1/sqrt(2)*([8]/[7]+(x-[9])/[8]))", histomin, histomax);// Sun PRC2021 low-energy tail

			g[ii] = new TF1("g", "gausn", histomin, histomax);// The [2]-N parameter in total is equivalent to the Constant in gausn
			g2[ii] = new TF1("g2", "gausn", histomin, histomax);// The [2]-N parameter in total is equivalent to the Constant in gausn
			p[ii]=new TF1("p","[0]*x+[1]+[2]/2/[3]*exp(0.5*([4]*[4]/([3]*[3]))+(x-[5])/[3])*ROOT::Math::erfc(1/sqrt(2)*([4]/[3]+(x-[5])/[4]))", histomin, histomax);//pure peak
			p2[ii] = new TF1("p2", "[0]*x+[1]+[2]/2/[3]*exp(0.5*([4]*[4]/([3]*[3]))+(x-[5])/[3])*ROOT::Math::erfc(1/sqrt(2)*([4]/[3]+(x-[5])/[4]))", histomin, histomax);//pure peak
			b[ii] = new TF1("b", "[0]*x+[1]", histomin, histomax);//pure bkg

			// 			fEMG[ii]=new TF1("total","[0]*x+[1]+sqrt(3.141592654/2)*[2]/[3]*[4]*exp(0.5*([4]*[4]/([3]*[3]))+(x-[5])/[3])*ROOT::Math::erfc(1/sqrt(2)*([4]/[3]+(x-[5])/[4]))",histomin, histomax);// Glassman PRC2019 low-energy tail
			// 			g[ii]=new TF1("g","gausn",histomin, histomax);// The [2]-N parameter in total is equivalent to the Constant in gausn
			// 			p[ii]=new TF1("p","[0]*x+[1]-[0]*x-[1]+sqrt(3.141592654/2)*[2]/[3]*[4]*exp(0.5*([4]*[4]/([3]*[3]))+(x-[5])/[3])*ROOT::Math::erfc(1/sqrt(2)*([4]/[3]+(x-[5])/[4]))",histomin, histomax);//pure peak
			// 			b[ii]=new TF1("b","[0]*x+[1]",histomin, histomax);//pure bkg
			fEMG[ii]->SetNpx(histoNbins);
			g[ii]->SetNpx(histoNbins);
			p[ii]->SetNpx(histoNbins);
			// 			p2[ii]->SetNpx(histoNbins);
			b[ii]->SetNpx(histoNbins);
			//fEMG[ii]->SetParameters(0.1,15,peaky[ii],10,10,peakx[ii]);//initial value [0]-A, [1]-B, [2]-N, [3]-τ, [4]-σ, [5]-μ
			fEMG[ii]->SetParameters(aguess, bguess, peaky[ii], 1.0, sigmaguess, peakx[ii], peaky[ii], 1.0, sigmaguess, peakx[ii] + gaplow/2);//initial value [0]-A, [1]-B, [2]-N, [3]-τ, [4]-σ, [5]-μ, [6]-N2, [7]-τ2, [8]-σ2, [9]-μ2
			//fEMG[ii]->SetParLimits(0,0,700);//Bkg A
			// 	fEMG[ii]->SetParLimits(1,-50000,300000);//Bkg B
			fEMG[ii]->SetParLimits(2, 42000, 51000);//Constant,min,max
			fEMG[ii]->SetParLimits(3, 0.02, 18);//Tau
			fEMG[ii]->SetParLimits(4, 0.01, 18);//Sigma
			//fEMG[ii]->SetParLimits(5, peakx[ii] - gaplow, peakx[ii] + gaphigh);//Mean
			fEMG[ii]->SetParLimits(5, 2060, 2080);//Mean
			fEMG[ii]->SetParLimits(6, 12500, 15300);//Constant_2nd,min,max
			fEMG[ii]->SetParLimits(7, 0.02, 18);//Tau_2nd
			fEMG[ii]->SetParLimits(8, 0.1, 18);//Sigma_2nd
			//fEMG[ii]->SetParLimits(9, peakx[ii] + gaplow, peakx[ii] + gaplow);//Mean_2nd
			fEMG[ii]->SetParLimits(9, 2040, 2060);//Mean_2nd

			fEMG[ii]->SetParNames("BkgA", "BkgB", "Const*bin", "Tau", "Sigma", "Mean", "Const_2nd*bin", "Tau_2nd", "Sigma_2nd", "Mean_2nd");
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
			par_err[ii][6] = fEMG[ii]->GetParError(6);//Obtaining the error of the 7th parameter
			par_err[ii][7] = fEMG[ii]->GetParError(7);//Obtaining the error of the 8th parameter
			par_err[ii][8] = fEMG[ii]->GetParError(8);//Obtaining the error of the 9th parameter
			par_err[ii][9] = fEMG[ii]->GetParError(9);//Obtaining the error of the 10th parameter
			parChi[ii] = fEMG[ii]->GetChisquare();
			parNDF[ii] = fEMG[ii]->GetNDF();
			p_value[ii] = fEMG[ii]->GetProb();//This probability is not the “probability that your fit is good.” If you did many fake experiments (draw many random samples of data points from the assumed distribution (your fit function)), this is the percentage of experiments that would give χ2 values ≥ to the one you got in this experiment.
			g[ii]->SetParameters(par[ii][2], par[ii][5], par[ii][4]);//set parameters for drawing gausn
			g[ii]->SetLineColor(7); 
			g2[ii]->SetParameters(par[ii][6], par[ii][9], par[ii][8]);//set parameters for drawing gausn
			g2[ii]->SetLineColor(7);
			p[ii]->SetParameters(par[ii][0],par[ii][1],par[ii][2],par[ii][3],par[ii][4],par[ii][5]);//set parameters for drawing peak
			p[ii]->SetLineColor(4);
			p2[ii]->SetParameters(par[ii][0], par[ii][1], par[ii][6], par[ii][7], par[ii][8], par[ii][9]);//set parameters for drawing peak
			p2[ii]->SetLineColor(4);
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
			h_confidence_interval[i][ii]->Draw("e3 same"); // plot the uncertainty band
			fEMG[ii]->Draw("same");
// 			g[ii]->Draw("same"); // plot the gaussian
// 			g2[ii]->Draw("same"); // plot the gaussian_2nd
			p[ii]->Draw("same");
 			p2[ii]->Draw("same");
			b[ii]->Draw("same");
			histo[i]->Draw("e same");

			TMatrixD cov = Fit_result_pointer->GetCovarianceMatrix();//error matrix
			TMatrixD cor = Fit_result_pointer->GetCorrelationMatrix();//parameter correlation coefficients
			cov.Print();
			cor.Print();
			double Utau_Utau = fitter->GetCovarianceMatrixElement(3, 3);//(Utau)^2
			double Usigma_Usigma = fitter->GetCovarianceMatrixElement(4, 4);//(Usigma)^2
			double Umean_Umean = fitter->GetCovarianceMatrixElement(5, 5);//(Umean)^2
			double Usigma2nd_Usigma2nd = fitter->GetCovarianceMatrixElement(8, 8);//(Usigma2nd)^2
			double Umean2nd_Umean2nd = fitter->GetCovarianceMatrixElement(9, 9);//(Umean2nd)^2
			double rou_Utau_Umean = fitter->GetCovarianceMatrixElement(3, 5);//(ρ*Utau*Umean), (3,5) (5,3) doesn't matter.
			double rou_Usigma_Umean = fitter->GetCovarianceMatrixElement(4, 5);//
			cout << rou_Usigma_Umean << endl; 
			double rou_Usigma2nd_Umean2nd = fitter->GetCovarianceMatrixElement(8, 9);//
			cout << rou_Usigma2nd_Umean2nd << endl;

			peaky[ii] = fEMG[ii]->GetMaximum();
			peakx[ii] = fEMG[ii]->GetMaximumX();
			peakxerr[ii] = sqrt(Utau_Utau + Umean_Umean + 2 * rou_Utau_Umean);

// 			inflation_factor = sqrt(parChi[ii] / parNDF[ii]);
// 			if (inflation_factor < 1)
				inflation_factor = 1;
			constant[ii] = par[ii][2]; constant_err[ii] = par_err[ii][2] * inflation_factor;
			tau[ii] = par[ii][3]; tau_err[ii] = par_err[ii][3] * inflation_factor;
			sig[ii] = par[ii][4]; sig_err[ii] = par_err[ii][4] * inflation_factor;
			mean[ii] = par[ii][5]; mean_err[ii] = par_err[ii][5] * inflation_factor;
			FWHM[ii] = par[ii][4] * 2.355; FWHM_err[ii] = par_err[ii][4] * 2.355 * inflation_factor;
			constant2nd[ii] = par[ii][6]; constant2nd_err[ii] = par_err[ii][6] * inflation_factor;
			tau2nd[ii] = par[ii][7]; tau2nd_err[ii] = par_err[ii][7] * inflation_factor;
			sig2nd[ii] = par[ii][8]; sig2nd_err[ii] = par_err[ii][8] * inflation_factor;
			mean2nd[ii] = par[ii][9]; mean2nd_err[ii] = par_err[ii][9] * inflation_factor;
			FWHM2nd[ii] = par[ii][8] * 2.355; FWHM2nd_err[ii] = par_err[ii][8] * 2.355 * inflation_factor;

			outfile << "LEGe" << i << "	Constant*binsize" << ii << "=	" << constant[ii] << "	+/-	" << constant_err[ii] << "	Mean" << ii << "=	" << mean[ii] << "	+/-	" << mean_err[ii] << "	Maximum" << ii << "=	" << peakx[ii] << "	+/-	" << peakxerr[ii] << "	Sigma" << ii << "=	" << sig[ii] << "	+/-	" << sig_err[ii] << "	Tau" << ii << "=	" << tau[ii] << "	+/-	" << tau_err[ii] << "	A" << ii << "=	" << par[ii][0] << "	+/-	" << par_err[ii][0] << "	B" << ii << "=	" << par[ii][1] << "	+/-	" << par_err[ii][1] << "	Chi2" << ii << "=	" << parChi[ii] << "	NDF" << ii << "=	" << parNDF[ii] << "	Area" << ii << "=	" << par[ii][2] / binwidth << "	FWHM" << ii << "=	" << FWHM[ii] << "	+/-	" << FWHM_err[ii] << "	Constant2*binsize" << ii << "=	" << constant2nd[ii] << "	+/-	" << constant2nd_err[ii] << "	Mean2" << ii << "=	" << mean2nd[ii] << "	+/-	" << mean2nd_err[ii] << "	Sigma2" << ii << "=	" << sig2nd[ii] << "	+/-	" << sig2nd_err[ii] << "	Tau2" << ii << "=	" << tau2nd[ii] << "	+/-	" << tau2nd_err[ii] << "	Area2" << ii << "=	" << par[ii][6] / binwidth << "	FWHM2" << ii << "=	" << FWHM2nd[ii] << "	+/-	" << FWHM2nd_err[ii] << endl;


			TPaveText* textgaus = new TPaveText(0.7, 0.4, 0.975, 0.98, "brNDC");//加标注left, down, right, up
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
			sprintf(paraprint, "Res%d=%.4f%%", ii, par[ii][4] / par[ii][5] * 2.355 * 100);
			textgaus->AddText(paraprint);
			sprintf(paraprint, "Tau%d=%.6f%s%.6f", ii, par[ii][3], "+/-", par_err[ii][3]);
			textgaus->AddText(paraprint);
			sprintf(paraprint, "Constant2nd*binsize%d=%.4f%s%.4f", ii, par[ii][6], "+/-", par_err[ii][6]);//par数组还保持着刚才的参数
			textgaus->AddText(paraprint);
			sprintf(paraprint, "Mean2nd%d=%.6f%s%.6f", ii, par[ii][9], "+/-", par_err[ii][9]);
			textgaus->AddText(paraprint);
			sprintf(paraprint, "Sigma2nd%d=%.6f%s%.6f", ii, par[ii][8], "+/-", par_err[ii][8]);
			textgaus->AddText(paraprint);
			sprintf(paraprint, "FWHM2nd%d=%.6f%s%.6f", ii, FWHM2nd[ii], "+/-", FWHM2nd_err[ii]);
			textgaus->AddText(paraprint);
			sprintf(paraprint, "Res2nd%d=%.4f%%", ii, par[ii][8] / par[ii][9] * 2.355 * 100);
			textgaus->AddText(paraprint);
			sprintf(paraprint, "Tau2nd%d=%.6f%s%.6f", ii, par[ii][7], "+/-", par_err[ii][7]);
			textgaus->AddText(paraprint);
			sprintf(paraprint, "A%d(%.6f%s%.6f)*x+B%d(%.5f%s%.5f)", ii, par[ii][0], "+/-", par_err[ii][0], ii, par[ii][1], "+/-", par_err[ii][1]);
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
			graph_residual[i]->GetXaxis()->SetTitle("ADC channel");//轴名
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
