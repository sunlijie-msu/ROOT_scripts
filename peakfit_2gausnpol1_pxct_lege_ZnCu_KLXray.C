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
void peakfit_2gausnpol1_pxct_lege_ZnCu_KLXray() // get histogram and Gausn fit some peaks
{
	TFile* _file0 = TFile::Open("F:/e21010/pxct/Adams_5days_9000pps.root");
	TTree* tree = (TTree*)_file0->Get("tree");
	//TFile* _file0 = TFile::Open("F:/e21010/pxct/run0334_LEGe_241Am_Z7117_ChamberCenter_window1.5us_TrigRise0.064us_TrigGap0.952us_Th350_CFDDelay0.304us_Scale7_for_X_efficiency_cal.root");
	//double fitrange_min = 44.6, fitrange_max = 46.05;
	double fitrange_min = 7.2, fitrange_max = 9.2;
	TH1D* histo = new TH1D("histo", "histo", 1000, 0, 50);

	TCanvas* canvaspeak = new TCanvas("LEGe", "LEGe", 1400, 500);
	canvaspeak->cd();
	canvaspeak->SetTopMargin(0.035);
	canvaspeak->SetRightMargin(0.03);
	canvaspeak->SetLeftMargin(0.09);
	canvaspeak->SetBottomMargin(0.19);
	canvaspeak->SetFrameLineWidth(3);
	gStyle->SetFrameLineWidth(3);

	tree->Draw("LEGe_e>>histo", "LEGe_e>0.11&&LEGe_e<50&&MSD12_e>121&&MSD12_e<853&&(MSD12_e+MSD26_e)>1500&&(MSD12_e+MSD26_e)<2500", "");

	histo->SetLineWidth(2);
	histo->Rebin(1);
	histo->SetLineColor(1);
	histo->SetStats(0);
	histo->SetTitle("");
	histo->GetXaxis()->SetTitle("Energy (keV)");
	histo->GetYaxis()->SetTitle("Counts per 0.05 keV");
	histo->GetXaxis()->CenterTitle();
	histo->GetYaxis()->CenterTitle();
	histo->GetXaxis()->SetLabelFont(132);
	histo->GetYaxis()->SetLabelFont(132);
	histo->GetXaxis()->SetLabelSize(0.07);
	histo->GetYaxis()->SetLabelSize(0.07);
	histo->GetXaxis()->SetTitleFont(132);
	histo->GetYaxis()->SetTitleFont(132);
	histo->GetXaxis()->SetTitleOffset(1.1);
	histo->GetYaxis()->SetTitleOffset(0.55);
	histo->GetXaxis()->SetTitleSize(0.08);
	histo->GetYaxis()->SetTitleSize(0.08);
	histo->GetYaxis()->SetTickLength(0.015);
	histo->GetXaxis()->SetRangeUser(6, 10);
	//histo->GetYaxis()->SetRangeUser(0,4);
	histo->GetYaxis()->SetNdivisions(505);
	histo->SetBinErrorOption(TH1::kPoisson);
	histo->Draw("e");

	// Define gtotal with the combined functions
	TF1* gtotal = new TF1("gtotal", "pol1(0) + gausn(2) + gausn(5)", fitrange_min, fitrange_max);
	gtotal->SetNpx(50000);
	// Create an array to hold the parameters
	double params[] = { 0.3, 0.3, 10, 8.05, 0.102, 10, 8.64, 0.102 };

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
	gtotal->SetParLimits(2, 0.5, 50);//Const1*bin
	gtotal->SetParLimits(3, 7.9, 8.1);//Mean_1
 	gtotal->SetParLimits(4, 0.1, 0.11);//Sigma_1
	gtotal->SetParLimits(5, 0.5, 50);//Const2*bin
	gtotal->SetParLimits(6, 8.5, 8.7);//Mean_2
	gtotal->SetParLimits(7, 0.1, 0.11);//Sigma_2
	// Perform the fit
	histo->Fit("gtotal", "MLE", "", fitrange_min, fitrange_max);
	histo->Fit("gtotal", "MLE", "", fitrange_min, fitrange_max);
	histo->Fit("gtotal", "MLE", "", fitrange_min, fitrange_max);
	gtotal->SetLineWidth(3);
	gtotal->Draw("same");
	canvaspeak->SaveAs("F:/e21010/pxct/60Ga_Geant4_px_spec.png");

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
	