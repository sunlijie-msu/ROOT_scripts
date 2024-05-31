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
void peakfit_6gausnpol1() // get histogram and Gausn fit one peak
{
	//TFile* _file0 = TFile::Open("F:/e21010/pxct/run0010_16_LEGe_152Eu_inChamber_vacuum_window0.5us_CFDdelay_adjusted_for_LEGe_efficiency_1440min_cal.root");
	TFile* _file0 = TFile::Open("F:/e21010/pxct/run0228_0229_0230_LEGe_XtRa_MSD26_152Eu_Z2707_inChamber_vacuum_XtRa_12mm_away_window1us_XtRaCFDdelay_0.2us_for_efficiency_cal.root");
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
	histo->GetXaxis()->SetRangeUser(4, 8.3);
	histo->GetYaxis()->SetNdivisions(505);
	//histo->GetYaxis()->SetRangeUser(1,20000);
	histo->Draw("hist");

	// Define gtotal with the combined functions
	TF1* gtotal = new TF1("gtotal", "pol1(0) + gausn(2) + gausn(5) + gausn(8) + gausn(11) + gausn(14) + gausn(17)", 0, 10);
	gtotal->SetNpx(5000);
	// Create an array to hold the parameters
	double params[] = { 2000, -1, 900, 4.93, 0.12, 15000, 5.58, 0.12, 14000, 6.17, 0.12, 4000, 6.58, 0.12, 3000, 7.10, 0.12, 200, 7.5, 0.12 };

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
	gtotal->SetParName(8, "Const3*bin");
	gtotal->SetParName(9, "Mean3");
	gtotal->SetParName(10, "Sigma3");
	gtotal->SetParName(11, "Const4*bin");
	gtotal->SetParName(12, "Mean4");
	gtotal->SetParName(13, "Sigma4");
	gtotal->SetParName(14, "Const5*bin");
	gtotal->SetParName(15, "Mean5");
	gtotal->SetParName(16, "Sigma5");
	gtotal->SetParName(17, "Const6*bin");
	gtotal->SetParName(18, "Mean6");
	gtotal->SetParName(19, "Sigma6");
	gtotal->SetParLimits(4, 0.1, 0.15);//Sigma_1
	gtotal->SetParLimits(15, 7.10, 7.15);//Mean_5
	gtotal->SetParLimits(16, 0.1, 0.2);//Sigma_5
	// Perform the fit
	histo->Fit("gtotal", "MLE", "", 4, 8.3);
	histo->Fit("gtotal", "MLE", "", 4, 8.3);
	histo->Fit("gtotal", "MLE", "", 4, 8.3);
	gtotal->Draw("same");

	double binwidth = histo->GetBinWidth(1);
	cout << "binwidth = " << binwidth << endl;
}
	