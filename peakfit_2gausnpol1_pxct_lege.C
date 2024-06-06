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
void peakfit_2gausnpol1_pxct_lege() // get histogram and Gausn fit some peaks
{
	TFile* _file0 = TFile::Open("F:/e21010/pxct/run0228_0229_0230_LEGe_XtRa_MSD26_152Eu_Z2707_inChamber_vacuum_XtRa_12mm_away_window1us_XtRaCFDdelay_0.2us_for_efficiency_cal.root");
	//TFile* _file0 = TFile::Open("F:/e21010/pxct/run0334_LEGe_241Am_Z7117_ChamberCenter_window1.5us_TrigRise0.064us_TrigGap0.952us_Th350_CFDDelay0.304us_Scale7_for_X_efficiency_cal.root");
	//double fitrange_min = 44.6, fitrange_max = 46.05;
	double fitrange_min = 46.1, fitrange_max = 47.3;
	TH1D* histo = (TH1D*)_file0->Get("hlege_e");
	TCanvas* canvaspeak = new TCanvas("LEGe", "LEGe", 1400, 553);//
	canvaspeak->cd();//
	canvaspeak->SetTopMargin(0.065);
	canvaspeak->SetRightMargin(0.035);
	canvaspeak->SetLeftMargin(0.10);
	canvaspeak->SetBottomMargin(0.16);
	canvaspeak->SetFrameLineWidth(2);

	histo->SetLineWidth(3);
	histo->Rebin(1);
	histo->SetLineColor(1);
	//histo->SetStats(0);
	histo->SetTitle("");
	histo->GetXaxis()->SetTitle("Energy (keV)");
	histo->GetYaxis()->SetTitle("Counts per 7 eV");
	histo->GetXaxis()->CenterTitle();
	histo->GetYaxis()->CenterTitle();
	histo->GetXaxis()->SetLabelFont(132);
	histo->GetYaxis()->SetLabelFont(132);
	histo->GetXaxis()->SetLabelSize(0.06);
	histo->GetYaxis()->SetLabelSize(0.06);
	histo->GetXaxis()->SetTitleFont(132);
	histo->GetYaxis()->SetTitleFont(132);
	histo->GetXaxis()->SetTitleOffset(1.0);
	histo->GetYaxis()->SetTitleOffset(0.72);
	histo->GetXaxis()->SetTitleSize(0.07);
	histo->GetYaxis()->SetTitleSize(0.07);
	histo->GetYaxis()->SetTickLength(0.02);
	histo->GetXaxis()->SetRangeUser(fitrange_min, fitrange_max);
	histo->GetYaxis()->SetNdivisions(505);
	//histo->GetYaxis()->SetRangeUser(1,20000);
	histo->Draw("hist");

	// Define gtotal with the combined functions
	TF1* gtotal = new TF1("gtotal", "pol1(0) + gausn(2) + gausn(5)", 0, 100);
	gtotal->SetNpx(50000);
	// Create an array to hold the parameters
	double params[] = { 100, 0.1, 20000, 46.5, 0.12, 20000, 46.7, 0.12 };

	// Set the parameters for gtotal using the array
	gtotal->SetParameters(params);

	// Set the parameter names individually
	gtotal->SetParName(0, "Bkg_Intercept");
	gtotal->SetParName(1, "Bkg_Slope");
	gtotal->SetParName(2, "Const1*bin");
	gtotal->SetParName(3, "Mean1");
	gtotal->SetParName(4, "Sigma1");
	gtotal->SetParName(5, "Const2*bin");
	gtotal->SetParName(6, "Mean2");
	gtotal->SetParName(7, "Sigma2");
	gtotal->SetParLimits(2, 10, 100000);//Const1*bin
	gtotal->SetParLimits(3, 46.3, 46.6);//Mean_1
 	gtotal->SetParLimits(4, 0.1, 0.20);//Sigma_1
	gtotal->SetParLimits(7, 0.1, 0.20);//Sigma_2
	gtotal->SetParLimits(5, 1, 100000);//Const2*bin
	gtotal->SetParLimits(6, 46.6, 46.8);//Mean_2
	// Perform the fit
	histo->Fit("gtotal", "MLE", "", fitrange_min, fitrange_max);
	histo->Fit("gtotal", "MLE", "", fitrange_min, fitrange_max);
	histo->Fit("gtotal", "MLE", "", fitrange_min, fitrange_max);
	gtotal->Draw("same");
	// print all parameters names values with errors
	for (int i = 0; i < gtotal->GetNpar(); i++)
	{
		cout << "Par" << i << " " << gtotal->GetParName(i) << "	" << gtotal->GetParameter(i) << "	" << gtotal->GetParError(i) << endl;
	}
	cout << "Chi2/NDF	" << gtotal->GetChisquare() / gtotal->GetNDF() << endl;
	cout << "Prob	" << gtotal->GetProb() << endl;
	double binwidth = histo->GetBinWidth(1);
	cout << "binwidth	" << binwidth << endl;
}
	