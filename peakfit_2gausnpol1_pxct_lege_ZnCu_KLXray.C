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
	double fitrange_min = 7.4, fitrange_max = 9.1;
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
	double params[] = { 0.3, 0.3, 1, 8.05, 0.102, 1, 8.64, 0.102 };

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
	gtotal->SetParLimits(2, 0.1, 20);//Const1*bin
	gtotal->SetParLimits(3, 7.9, 8.15);//Mean_1
 	gtotal->SetParLimits(4, 0.100, 0.106);//Sigma_1
	gtotal->SetParLimits(5, 0.1, 20);//Const2*bin
	gtotal->SetParLimits(6, 8.5, 8.7);//Mean_2
	gtotal->SetParLimits(7, 0.100, 0.106);//Sigma_2
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
	cout << "Cu counts	" << gtotal->GetParameter(2) / binwidth << "	" << gtotal->GetParError(2) / binwidth * sqrt(gtotal->GetChisquare() / gtotal->GetNDF()) << endl;
	cout << "Zn counts	" << gtotal->GetParameter(5) / binwidth << "	" << gtotal->GetParError(5) / binwidth * sqrt(gtotal->GetChisquare() / gtotal->GetNDF()) << endl;

	// Define the constants and parameters
	double Ratio_Ka_Emission_Zn_Cu = 1.0707;
	double Ratio_Detection_Efficiency_Zn_Cu = 1.0541;
	double lifetime_Zn_K_vacancy = 0.422; // in fs

	// Retrieve the peak counts and uncertainties from the fit
	double Cu_Ka_counts = gtotal->GetParameter(2) / binwidth;
	double Zn_Ka_counts = gtotal->GetParameter(5) / binwidth;
	double uncertainty_Cu_Ka_counts = (gtotal->GetParError(2) / binwidth) * sqrt(gtotal->GetChisquare() / gtotal->GetNDF());
	double uncertainty_Zn_Ka_counts = (gtotal->GetParError(5) / binwidth) * sqrt(gtotal->GetChisquare() / gtotal->GetNDF());

	// Calculate the real peak counts
	double real_Cu_Ka_counts = Cu_Ka_counts * Ratio_Ka_Emission_Zn_Cu * Ratio_Detection_Efficiency_Zn_Cu;
	double real_Zn_Ka_counts = Zn_Ka_counts;

	// Propagate the uncertainties for the real peak counts
	double uncertainty_real_Cu_Ka_counts = uncertainty_Cu_Ka_counts * Ratio_Ka_Emission_Zn_Cu * Ratio_Detection_Efficiency_Zn_Cu;

	double uncertainty_real_Zn_Ka_counts = uncertainty_Zn_Ka_counts;

	// Calculate the ratio and its uncertainty
	double Ratio_Zn_Cu = real_Zn_Ka_counts / real_Cu_Ka_counts;
	double uncertainty_Ratio_Zn_Cu = Ratio_Zn_Cu * sqrt(
		pow(uncertainty_real_Zn_Ka_counts / real_Zn_Ka_counts, 2) +
		pow(uncertainty_real_Cu_Ka_counts / real_Cu_Ka_counts, 2));

	// Calculate the lifetime of the proton emitting state and its uncertainty
	double Lifetime_proton_emitting_state = Ratio_Zn_Cu * lifetime_Zn_K_vacancy;
	double Uncertainty_lifetime_proton_emitting_state = Lifetime_proton_emitting_state * (uncertainty_Ratio_Zn_Cu / Ratio_Zn_Cu);
	
	// Output the results
	cout << "Real Cu Ka peak counts: " << real_Cu_Ka_counts << "	" << uncertainty_real_Cu_Ka_counts << endl;
	cout << "Real Zn Ka peak counts: " << real_Zn_Ka_counts << "	" << uncertainty_real_Zn_Ka_counts << endl;
	cout << "Ratio: " << Ratio_Zn_Cu << "	" << uncertainty_Ratio_Zn_Cu << endl;
	cout << "Lifetime of proton emitting state: " << Lifetime_proton_emitting_state << "	" << Uncertainty_lifetime_proton_emitting_state << " fs" << endl;

}
	