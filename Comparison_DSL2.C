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
//Run with ROOT5. Somehow ROOT6 fits one spectrum and stops.
// ROOTSYS put in C:\root5 or C:\root6. simple.
//main_readhists; fcn(); comparehists(); main_draw_save;//search m o d i f y to change something
// one peak + linear fline two-part bkg. Recommended for final results!
// could subtract 2234's low-energy Compton bkg. // It could output some final figures.
TFile* fin_simu, * fin_data;
TH1F* h_simulated_spec, * h_measured_spec, * h_fit_background;
TH1F* h_simulated_spec_scaled_plus_fit_background_scaled;

const double Binsperkev = 1; //Number of bins per keV, the only place to change binwidth, all other variables is related to Binsperkev // 10 for NSCL, 2 for RIBLL, useless for DSL
const int binwidth = 1; //binwidths in units of keV
const int factor_rebin = 1; //simu and data Rebin factor
float bkgDown = 0.50;
float bkgUp = 0.51; // if you set bkgDown and bkgUp the same value, the bkg will be set to zero in minimization, which is wrong

const double E0_gamma = 7333; //Ex=7784.7
const int peakrange_min = 7740; // bin = 4411; bin center = 4410.5
const int peakrange_max = 7870; // bin = 4520; bin center = 4523.5
const double fitrange_min = 7700; // bin = 4281; bin center = 4280.5
const double fitrange_max = 7910; // bin = 4660; bin center = 4659.5
const int num_bins_peak = peakrange_max - peakrange_min;
double Low_bkg = 0.92;
double High_bkg = 1.08;

// const double E0_gamma = 4156; //Ex=6390
// const int peakrange_min = 4410; // bin = 4411; bin center = 4410.5
// const int peakrange_max = 4524; // bin = 4520; bin center = 4523.5
// const int Ea = 39;
// double T0_lifetime = 0;
// int Lifetimestep = 1;
// double centroid = 4155.84;
// double Eg_uncertainty = 0.31;
// const double fitrange_min = 4280; // bin = 4281; bin center = 4280.5
// const double fitrange_max = 4660; // bin = 4660; bin center = 4659.5
// const int num_bins_fullrange = 380;
// const int num_bins_peak = 114;
// double Low_bkg = 0.97;
// double High_bkg = 1.03;

// double Tau_values[] = { 0.0, 5.0, 10.0, 15.0, 20.0 };
// double Eg_values[] = { 7331.20, 7333.20, 7335.20 };
// double Bkg_values[] = { 0.90, 1.00, 1.10 };
// double SP_values[] = { 0.90, 1.00, 1.10 };
// double AC_values[] = { 0.0 };

double Tau_values[] = { 0.0, 5.0, 10.0, 15.0, 20.0 };
double Eg_values[] = { 7333.20 };
double Bkg_values[] = { 1.00 };
double SP_values[] = { 1.00 };
double AC_values[] = { 0.0 };

Bool_t reject;

double Chi2 = 0;//Chi2 is global varible, use Chi2 to extract the minimized Chi2 obtained in fcn(), and output Chi2 in main().
long double LBayesian = 1;

double CompareHists(TH1F* his_fit, TH1F* his_data, int first_bin, int last_bin)
//compare two histograms (data and simulation). returns chi-square
//Two histograms must be with the same length
{
	double chisquareNeyman = 0, chisquarePearson = 0, likelihood = 0;     //chi square
	long double ydata_factorial = 1;
	LBayesian = 1;

	double err_ydata = 0, residual, yfit, ydata;
	for (int i = first_bin; i <= last_bin; i++)                         //loop over all bins   
	{
		ydata = his_data->GetBinContent(i); //h_measured_spec
		err_ydata = his_data->GetBinError(i);
		yfit = his_fit->GetBinContent(i); //h_simulated_gammaspec_and_fit_background_combined
		residual = ydata - yfit; //ML method

		if (ydata != 0) { chisquareNeyman += residual * residual / ydata; }//Neyman Chi2 method
		if (yfit != 0) { chisquarePearson += residual * residual / yfit; }//Pearson Chi2 method
		if (yfit != 0 && ydata != 0) { likelihood += yfit - ydata + ydata * log(ydata / yfit); } //ML method, highly recommended!
		if (yfit != 0 && ydata == 0) { likelihood += yfit - ydata; } //ML method for bins with ydata = 0, the log term is -infinity
		//if(y[0]>0)	dy[0] = (his_data->GetBinErrorLow(i)+his_data->GetBinErrorUp(i))/2;
		//if(y[0]==0)	dy[0] = his_data->GetBinErrorUp(i);//empty bin ErrorUp=1.8, ErrorLow=0, Error=0;

		ydata_factorial = 1;
		for (long ii = 1; ii <= ydata; ii++) ydata_factorial = ydata_factorial * ii; //Bayesian
		//cout<<"yfit= "<<yfit<<" ydata= "<<ydata<<" ydata_factorial= "<<ydata_factorial<<" pow(yfit,ydata)= "<<pow(yfit,ydata)<<" LB= "<<(pow(yfit,ydata)/ydata_factorial)*exp(-yfit)<<endl;
		if (yfit != 0 && ydata_factorial != 0) LBayesian *= (pow(yfit, ydata) / ydata_factorial) * exp(-yfit); //Bayesian Naomi Galinski, equivalent to the ML method above.
	}

	likelihood = likelihood * 2; //ML method
	return likelihood; //ML method
}

Double_t fline(Double_t* x, Double_t* par) // pol1
{ // Double_t fcn(Double_t *x, Double_t *params)
	if (reject && x[0] > peakrange_min && x[0] < peakrange_max)
	{
		TF1::RejectPoint();
		return 0;
	}
	return par[0] + par[1] * x[0];
}


void fcn(Int_t& npar, Double_t* gin, Double_t& f, Double_t* par, Int_t iflag)
//npar number of free parameters involved in minimization
//gin partial derivatives (return values) computed gradient values (optional)
//f the function value itself (return value)
//par parameter values
//iflag flag word to switch between several actions of FCN//Usually only f and par are essential
{
	if (h_simulated_spec_scaled_plus_fit_background_scaled) delete h_simulated_spec_scaled_plus_fit_background_scaled;
	h_simulated_spec_scaled_plus_fit_background_scaled = new TH1F("h_simulated_spec_scaled_plus_fit_background_scaled", "h_simulated_spec_scaled_plus_fit_background_scaled", 10000, 0, 10000);
	//h_simulated_spec_scaled_plus_fit_background_scaled->Scale(par[0]);  //Scale this histogram by a constant (relative intensity), determined by the minuit above
	//h_simulated_gammaspec_and_fit_background_combined->Sumw2();//histogram is already filled, the sum of squares of weights is filled with the existing bin contents
	//This function is automatically called when the histogram is created
	h_simulated_spec_scaled_plus_fit_background_scaled->Add(h_simulated_spec, par[0]);
	h_simulated_spec_scaled_plus_fit_background_scaled->Add(h_fit_background, par[1]);//h_simulated_gammaspec_and_fit_background_combined=p1+p2+p3+h_fit_background
	//	f = CompareHists
	f = CompareHists(h_simulated_spec_scaled_plus_fit_background_scaled, h_measured_spec, peakrange_min + 1, peakrange_max);//first bin -> last bin peak region
	//	f = CompareHists(h_simulated_gammaspec_and_fit_background_combined,h_measured_spec,1,num_bins_fullrange);//first bin -> last bin the same as the sentence above
	Chi2 = f;//Chi2 is global varible, use Chi2 to extract the minimized Chi2 obtained in fcn(), and output Chi2 in main().
}








void Comparison_DSL2()
{
	char simurootname[400], outrootname[400], outtxtname[400], outfigname[400], tname[400], hname[400];
	sprintf(outtxtname, "%s%.0f%s%.0f%s", "D:/X/out/DSL2_Comparison/Eg", E0_gamma, "/DSL_23Mg", E0_gamma, "_model_y_values.dat");
	FILE* outfilehist1 = fopen(outtxtname, "w"); // "w" stands for write, and it will create a new file if it doesn't exist or clear the file to zero length if it already exists.

	sprintf(outtxtname, "%s%.0f%s%.0f%s", "D:/X/out/DSL2_Comparison/Eg", E0_gamma, "/DSL_23Mg", E0_gamma, "_model_y_values_var.dat");
	FILE* outfilehist2 = fopen(outtxtname, "w"); // "w" stands for write, and it will create a new file if it doesn't exist or clear the file to zero length if it already exists.

	sprintf(outtxtname, "%s%.0f%s%.0f%s", "D:/X/out/DSL2_Comparison/Eg", E0_gamma, "/DSL_23Mg", E0_gamma, "_model_parameter_values.dat");
	FILE* outfilepara = fopen(outtxtname, "w"); // "w" stands for write, and it will create a new file if it doesn't exist or clear the file to zero length if it already exists.

	int i_model_run = 0;
	for (int iEg = 0; iEg < sizeof(Eg_values) / sizeof(Eg_values[0]); iEg++)
	{
		for (int iBkg = 0; iBkg < sizeof(Bkg_values) / sizeof(Bkg_values[0]); iBkg++)
		{
			for (int iSP = 0; iSP < sizeof(SP_values) / sizeof(SP_values[0]); iSP++)
			{
				for (int iAC = 0; iAC < sizeof(AC_values) / sizeof(AC_values[0]); iAC++)
				{
// 					sprintf(b_name, "%s%.0f%s", "D:/X/out/DSL/Chi2/Chi2_Gamma", E0_gamma, "slice_vTau_G4.dat");
// 					FILE* outChi2file = fopen(b_name, "a"); //ofstream doesn't work
					for (int iTau = 0; iTau < sizeof(Tau_values) / sizeof(Tau_values[0]); iTau++)
					{
						//readHists();
						Chi2 = 0;//Chi2 is global varible, use Chi2 to extract the minimized Chi2 obtained in fcn(), and output Chi2 in main(). Must clean Chi2 before each run in main().

						sprintf(simurootname, "%s%.0f%s%.0f%s%.2f%s%.1f%s%.2f%s%.1f%s", "F:/out/G4_rootfiles_with_tree_Eg", E0_gamma, "/Mg23_Gamma", E0_gamma, "_Eg", Eg_values[iEg], "_Tau", Tau_values[iTau], "_SP", SP_values[iSP], "_AC", AC_values[iAC], ".root");
						//read in lots of simulation root files. Input root file names don't contain Bkg_values[iBkg], but output fig/dat/root file names contain Bkg_values[iBkg].
						fin_simu = new TFile(simurootname); //Rootfile with Simulation histograms
						cout << simurootname << endl;
						fin_data = new TFile("F:/out/testadd.root"); //Rootfile with Experimental Data Histogram modify

						//Get Histograms from Simulation
						sprintf(hname, "Eg");
						h_simulated_spec = (TH1F*)fin_simu->Get(hname);
						h_simulated_spec->Rebin(factor_rebin); //Rebin the simulation histogram

						int G4bkg = 1; // 1 means keep G4 bkg in the spectrum, 0 means subtract G4 bkg from the spectrum (preferable).
						if (E0_gamma == 7333 && G4bkg == 0) // get rid of the low-energy Geant bkg so that pure peak component can be used to do data-to-simu comparison
						{
							TF1* G4pol1 = new TF1("G4pol1", "[0]*x+[1]", 0, 10000);//G4 bkg
							h_simulated_spec->Fit("G4pol1", "MLE0", "", fitrange_min, peakrange_min - 10);
							h_simulated_spec->Fit("G4pol1", "MLE0", "", fitrange_min, peakrange_min - 10);
							double aG4 = G4pol1->GetParameter(0);
							double bG4 = G4pol1->GetParameter(1);
							for (int ibin = fitrange_min+1; ibin <= fitrange_max; ibin++)
							{
								double ycount = h_simulated_spec->GetBinContent(ibin);
								double xvalue = h_simulated_spec->GetBinCenter(ibin);
								double ybkg = bG4 + aG4 * xvalue;
								if (ybkg < 0) break;
								h_simulated_spec->SetBinContent(ibin, ycount - ybkg); //subtract fitted Compton bkg from G4 simulated spectrum
								//cout << ibin << "	" << ycount << "	" << ybkg << endl;
							}
						}

						sprintf(hname, "%s", "centersum2");//data histogram name
						h_measured_spec = (TH1F*)fin_data->Get(hname); //Get Gamma ray spectrum from data
						h_measured_spec->Rebin(factor_rebin); //Rebin the data histogram
						h_fit_background = new TH1F("h_fit_background", "h_fit_background", 10000, 0, 10000);

						// 			TF1 *f1=new TF1("f1","[0]+x*[1]",fitrange_min,fitrangelow2);
						// 			TF1 *f2=new TF1("f2","[0]+x*[1]",fitrangehigh1,fitrangehigh2);
						TF1* fl = new TF1("fl", fline, fitrange_min, fitrange_max, 2); //npar=2 is the number of free parameters used by the function pol1
						//we want to fit only the linear background excluding the peak part
						reject = kTRUE;
						//fl->SetParLimits(0, 0, 100); //set limits for intercept
						//fl->SetParLimits(1, -0.02, 0.0001); //set limits for slope
						h_measured_spec->Fit(fl, "M0"); // "1" will draw a full curve, "0" will not draw the fitted curve in the excluded region
						// option "L" a likelihood fit is used instead of the default chi2 square fit. Sometimes L doesn't work.
						reject = kFALSE;

						//store 3 separate functions for visualization (three bands) can be scaled
						TH1D* h_confidence_interval1 = new TH1D("h_confidence_interval1", "Fitted func with conf.band", (peakrange_min - fitrange_min) * 10, fitrange_min, peakrange_min);
						(TVirtualFitter::GetFitter())->GetConfidenceIntervals(h_confidence_interval1, 0.683);
						h_confidence_interval1->SetStats(kFALSE);
						h_confidence_interval1->SetFillColor(kGreen - 7);

						TH1D* h_confidence_interval2 = new TH1D("h_confidence_interval2", "Fitted func with conf.band", (fitrange_max - peakrange_max) * 10, peakrange_max, fitrange_max);
						(TVirtualFitter::GetFitter())->GetConfidenceIntervals(h_confidence_interval2, 0.683);
						h_confidence_interval2->SetStats(kFALSE);
						h_confidence_interval2->SetFillColor(kGreen - 7);

						TH1D* h_confidence_interval3 = new TH1D("h_confidence_interval3", "Fitted func with conf.band", (fitrange_max - fitrange_min) * 10, fitrange_min, fitrange_max);
						(TVirtualFitter::GetFitter())->GetConfidenceIntervals(h_confidence_interval3, 0.683);
						h_confidence_interval3->SetStats(kFALSE);
						h_confidence_interval3->SetFillColor(kGreen - 7);

						for (int ibin = fitrange_min+1; ibin <= fitrange_max; ibin++)
						{
							double b = fl->GetParameter(0) + h_fit_background->GetBinCenter(ibin) * fl->GetParameter(1);
							if (Bkg_values[iBkg] == 1.1) { b = b * High_bkg; }
							if (Bkg_values[iBkg] == 0.9) { b = b * Low_bkg; }
							h_fit_background->SetBinContent(ibin, b);//fit data, get background histogram
						}

						//minuit**************************************************************************************************************************//
			//			const int nParams=4; //number of parameters that are to be found
						const int nParams = 2; //number of parameters that are to be found
						TMinuit* gMin = new TMinuit(nParams);  //initialize TMinuit with a maximum of n parameters
						gMin->SetFCN(fcn);//set the address of the minimization function// fcn里会调用CompareHists子函数
						Double_t arglist[10];
						Int_t ierflg = 0;//command executed normally
						arglist[0] = 1;
						// 			Double_t vstart[4] = {0.04,1.0,0.004,0.004};//initial guess par_Fit, par_bkg
						// 			Double_t step[4] = {0.0001,0.01,0.0001,0.0001};                //fitting step
						Double_t vstart[2] = { 0.001,1.0 };//initial guess par_Fit, par_bkg
						Double_t step[2] = { 0.000001,0.0001 };                //fitting step

						gMin->mnparm(0, "a1", vstart[0], step[0], 0.0001, 1.0, ierflg);// par for simu
						//parID, parName, initialGuess, step, lowLimit, highLimit, irrelevant //0.00几，因为是从百万高统计模拟的histo scale下来的
						gMin->mnparm(1, "a2", vstart[1], step[1], bkgDown, bkgUp, ierflg);// par for bkg //bkg 基本0.9-1.2之间，不会跟拟合的本底水平差别太大

						arglist[0] = 2000;  //max number of fitting iterations
						arglist[1] = 1.;  //tolerance 1 means one sigma
						gMin->mnexcm("MIGRAD", arglist, 2, ierflg);  //run the minimization using MIGRAD // no need to change the number "2".						//minuit.mnexcm(command, argument list, argument number, error flag);						//command is a character array holding the command
						//argument list a double array as the arguments
						//argument number is the number of arguments						//error flag return value (!= 0)						//Precisely 						//  the minimum, but you need to be careful since it assumes the inputs are close to the minimum
						Double_t amin, edm, errdef;
						Int_t nvpar, nparx, icstat;
						gMin->mnstat(amin, edm, errdef, nvpar, nparx, icstat);

						double pars[nParams], errs[nParams];
						for (int i = 0; i < nParams; i++)
						{
							gMin->GetParameter(i, pars[i], errs[i]);
							printf("pars=	%f	errs=	%f\n", pars[i], errs[i]);
						}

						//visualization

						TCanvas* c1 = new TCanvas("c1", "c1", 1100, 700);
						c1->cd();
						TPad* pad1 = new TPad("pad1", "The pad 70% of the height", 0.0, 0.3, 1.0, 1.0);// Double_t xlow, Double_t ylow, Double_t xup, Double_t yup,
						TPad* pad2 = new TPad("pad2", "The pad 30% of the height", 0.0, 0.0, 1.0, 0.3);
						pad1->SetTopMargin(0.08);
						pad1->SetRightMargin(0.035);
						pad1->SetLeftMargin(0.11);
						pad1->SetBottomMargin(0.06);
						pad1->SetFrameLineWidth(2);
						//pad1->SetBorderMode(0);

						pad2->SetTopMargin(0.03);
						pad2->SetRightMargin(0.035);
						pad2->SetLeftMargin(0.11);
						pad2->SetBottomMargin(0.35);
						pad2->SetFrameLineWidth(2);
						//pad2->SetBorderMode(0);

						pad1->Draw();
						pad2->Draw();
						pad1->cd();
						gStyle->SetOptTitle(0);

						h_measured_spec->SetStats(0);
						h_measured_spec->SetLineColor(1);
						h_measured_spec->SetMarkerColor(1);
						h_measured_spec->SetLineWidth(1);
						h_measured_spec->Sumw2(kFALSE);
						h_measured_spec->SetBinErrorOption(TH1::kPoisson);//TH1::kNormal or TH1::kPoisson
						h_measured_spec->Draw("e");
						sprintf(tname, "%s%d%s", "Counts per ", binwidth, " keV");
						//h_measured_spec->GetXaxis()->SetTitle("Energy (keV)");
						h_measured_spec->GetYaxis()->SetTitle(tname);
						h_measured_spec->GetXaxis()->SetLabelFont(132);
						h_measured_spec->GetYaxis()->SetLabelFont(132);
						h_measured_spec->GetXaxis()->SetTitleFont(132);
						h_measured_spec->GetYaxis()->SetTitleFont(132);
						h_measured_spec->GetXaxis()->SetLabelSize(0.06);
						h_measured_spec->GetYaxis()->SetLabelSize(0.06);
						h_measured_spec->GetXaxis()->CenterTitle();
						h_measured_spec->GetYaxis()->CenterTitle();
						//h_measured_spec->GetXaxis()->SetTitleOffset(0.8);
						h_measured_spec->GetYaxis()->SetTitleOffset(0.84);
						//h_measured_spec->GetXaxis()->SetTitleSize(0.05);
						h_measured_spec->GetYaxis()->SetTitleSize(0.06);
						h_measured_spec->GetXaxis()->SetNdivisions(505);//n = n1 + 100*n2 + 10000*n3
						h_measured_spec->GetYaxis()->SetNdivisions(505);//n = n1 + 100*n2 + 10000*n3
						h_measured_spec->GetYaxis()->SetTickLength(0.015);
						h_measured_spec->GetXaxis()->SetRangeUser(fitrange_min, fitrange_max);

						h_fit_background->Scale(pars[1]);//h_fit_background scaled up/down by up to 10%
						// 			h_confidence_interval1->Scale(pars[1]);
						//  			h_confidence_interval1->Draw("e3 same");
						// 			h_confidence_interval2->Scale(pars[1]);
						//  			h_confidence_interval2->Draw("e3 same");
						h_confidence_interval3->Scale(pars[1]);
						h_confidence_interval3->Draw("e3 same"); // background uncertainty band
						h_fit_background->SetLineWidth(2);
						h_simulated_spec->SetStats(0);
						h_simulated_spec->SetLineColor(3);
						h_simulated_spec->SetMarkerColor(3);
						h_fit_background->SetStats(0);
						h_fit_background->SetLineColor(4);
						h_fit_background->SetMarkerColor(4);
						h_simulated_spec_scaled_plus_fit_background_scaled->SetLineColor(kRed);
						h_simulated_spec_scaled_plus_fit_background_scaled->SetMarkerColor(kRed);
						h_simulated_spec_scaled_plus_fit_background_scaled->SetLineWidth(2);
						h_simulated_spec_scaled_plus_fit_background_scaled->Draw("samec"); //“][”: Draw histogram without the vertical lines for the first and the last bin. Use it when superposing many histograms on the same picture.
						pad1->RedrawAxis();

						char paraprint[100];
						TPaveText* textchi = new TPaveText(0.11, 0.93, 0.965, 0.99, "brNDC");//left, down, right, up
						textchi->SetBorderSize(1);
						textchi->SetFillColor(0);
						textchi->SetTextAlign(12);//align = 10*HorizontalAlign + VerticalAlign, 12 means水平左对齐、垂直居中对齐
						textchi->SetTextFont(132);//font = 10 * fontID + precision, 12+2,12 means Symbol; 13+2, 13 means News Times Roman
						sprintf(paraprint, "#chi^{2} = %.2f / %d = %.5f   E_{#gamma} = %.2f keV   #tau = %.1f fs   N = %.1f +/- %.1f", Chi2, (num_bins_peak - 2), Chi2 / (num_bins_peak - 2), Eg_values[iEg], Tau_values[iTau], pars[0] * 14000, errs[0] * 14000);
						textchi->AddText(paraprint);
						textchi->Draw();

						// for residuals plot

						pad2->cd();
						double x_values[10000] = { 0 }, y_residuals[10000] = { 0 }, x_errors[10000] = { 0 }, y_errors[10000] = { 0 };

						for (int ibin = fitrange_min + 1; ibin <= fitrange_max; ibin++)
						{ // data array starts from [0], bin starts from [1]
							x_values[ibin - 1] = h_measured_spec->GetBinCenter(ibin);
							x_errors[ibin - 1] = binwidth / 2.0;
							y_residuals[ibin - 1] = h_measured_spec->GetBinContent(ibin) - h_simulated_spec_scaled_plus_fit_background_scaled->GetBinContent(ibin);
							y_errors[ibin - 1] = (h_measured_spec->GetBinErrorUp(ibin) + h_measured_spec->GetBinErrorLow(ibin)) / 2.0;
						}

						TGraph* graph_residual = new TGraphErrors(10000, x_values, y_residuals, x_errors, y_errors); //TGraph(n,x,y,ex,ey);
						graph_residual->SetTitle("");//图名
						graph_residual->GetXaxis()->SetTitle("Energy (keV)");
						graph_residual->GetYaxis()->SetTitle("Data #minus Fit");
						graph_residual->GetXaxis()->CenterTitle();
						graph_residual->GetYaxis()->CenterTitle();
						graph_residual->GetXaxis()->SetLabelFont(132);
						graph_residual->GetYaxis()->SetLabelFont(132);
						graph_residual->GetXaxis()->SetLabelSize(0.15);
						graph_residual->GetYaxis()->SetLabelSize(0.15);
						graph_residual->GetXaxis()->SetTitleFont(132);
						graph_residual->GetYaxis()->SetTitleFont(132);
						graph_residual->GetXaxis()->SetTitleOffset(1.1);
						graph_residual->GetYaxis()->SetTitleOffset(0.35);
						graph_residual->GetXaxis()->SetTitleSize(0.15);
						graph_residual->GetYaxis()->SetTitleSize(0.15);
						graph_residual->GetXaxis()->SetNdivisions(505);
						graph_residual->GetYaxis()->SetNdivisions(505);
						graph_residual->GetYaxis()->SetTickLength(0.015);
						//graph_residual->SetStats(0);
						graph_residual->GetXaxis()->SetRangeUser(fitrange_min, fitrange_max);
						// graph_residual->GetYaxis()->SetRangeUser(-50, 50); 
						graph_residual->SetLineWidth(1);
						graph_residual->SetLineColor(1);
						graph_residual->SetMarkerStyle(6);
						graph_residual->SetMarkerColor(1);

						graph_residual->Draw("APZ");//"A": Axis are drawn around the graph, "P": The current marker is plotted at each point, "Z": Do not draw small horizontal and vertical lines the end of the error bars. Without "Z", the default is to draw these.
						TLine* T1 = new TLine(fitrange_min, 0, fitrange_max, 0);
						T1->Draw("R");//"R" means the line is drawn with the current line attributes

						cout << "=======	Eg=	" << Eg_values[iEg] << "	tau=	" << Tau_values[iTau] << " fs" << endl;

						//save text
						//fprintf(outfile2, "%d	%e", Tau_values[iTau], LBayesian);//LB is global varible, use LB to extract the LB corresponding to the minimum Chi2 iterated by fcn(), and output LB in main().// use this sentence to output 2D LB matrix
						//fprintf(outChi2file, "%.0f	%.3f", Tau_values[iTau], Chi2);//Chi2 is global varible, use Chi2 to extract the minimized Chi2 obtained in fcn(), and output Chi2 in main(). // use this sentence to output 1D2D chisquare matrix

						sprintf(tname, "%s%.0f%s%.0f%s%.2f%s%.1f%s%.2f%s%.1f%s%.2f", "D:/X/out/DSL2_Comparison/Eg", E0_gamma, "/Mg23_Gamma", E0_gamma, "_Eg", Eg_values[iEg], "_Tau", Tau_values[iTau], "_SP", SP_values[iSP], "_AC", AC_values[iAC], "_Bkg", Bkg_values[iBkg]);
						//save figure
						sprintf(outfigname, "%s%s", tname, ".png");
						c1->SaveAs(outfigname);

						//save root file
						sprintf(outrootname, "%s%s", tname, ".root");
						TFile* outfileroot = new TFile(outrootname, "RECREATE");

						h_measured_spec->Write("h_measured_spec");
						h_simulated_spec_scaled_plus_fit_background_scaled->Write("h_simulated_gammaspec_and_fit_background_combined");
						graph_residual->Write("graph_residual");
						h_simulated_spec->Write("h_simulated_spec");
						h_fit_background->Write("h_fit_background");


						//save histogram bin counts for Surmise "model y values full range"
			 			sprintf(outtxtname,"%s%.0f%s%.0f%s","D:/X/out/DSL2_Comparison/Eg",E0_gamma,"/DSL_23Mg",E0_gamma,"_model_y_values.dat");
			 			outfilehist1 = fopen (outtxtname,"a"); //ofstream doesn't work
						fprintf(outfilehist1, "%d", i_model_run);
						for (int ibin = fitrange_min + 1; ibin <= fitrange_max; ibin++)
			 			{
			 				fprintf(outfilehist1, "	%.4f", h_simulated_spec_scaled_plus_fit_background_scaled->GetBinContent(ibin));
			 			}
						fprintf(outfilehist1, "\n");
						fclose(outfilehist1);

						//save histogram bin counts for Surmise "model y values errors full range" (this may be less important)
						sprintf(outtxtname, "%s%.0f%s%.0f%s", "D:/X/out/DSL2_Comparison/Eg", E0_gamma, "/DSL_23Mg", E0_gamma, "_model_y_values_var.dat");
						outfilehist2 = fopen(outtxtname, "a"); //ofstream doesn't work
						fprintf(outfilehist2, "%d", i_model_run);
						for (int ibin = fitrange_min + 1; ibin <= fitrange_max; ibin++)
						{
							fprintf(outfilehist2, "	%.6f", errs[0] / pars[0] * h_simulated_spec_scaled_plus_fit_background_scaled->GetBinContent(ibin));
						}
						fprintf(outfilehist2, "\n");
						fclose(outfilehist2);

						// save parameters txt for Surmise "model parameters"
						sprintf(outtxtname, "%s%.0f%s%.0f%s", "D:/X/out/DSL2_Comparison/Eg", E0_gamma, "/DSL_23Mg", E0_gamma, "_model_parameter_values.dat");
			 			outfilepara = fopen (outtxtname,"a"); //ofstream doesn't work
						fprintf(outfilepara, "%d	%.1f	%.2f	%.2f	%.2f\n", i_model_run, Tau_values[iTau], Eg_values[iEg], Bkg_values[iBkg], SP_values[iSP]);
						fclose(outfilepara);

						cout << "i_model_run=	" << i_model_run << endl;
						i_model_run++;

						delete outfileroot;
					} // for (iTau)
				} // for (iAC)
			} // for (iSP)
		} // for (iBkg)
	} // for(iEg)
		
}//main_readhists; fcn(); comparehists(); main_draw_save;

