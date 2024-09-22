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
void peakfit_2gausnpol1_pxct_lege_ZnCu_KXray() // get histogram and Gausn fit some peaks
{
	TFile* _file0 = TFile::Open("F:/e21010/pxct/ExG4_60Ga_all.root");
	TTree* tree = (TTree*)_file0->Get("tree");
	//TFile* _file0 = TFile::Open("F:/e21010/pxct/run0334_LEGe_241Am_Z7117_ChamberCenter_window1.5us_TrigRise0.064us_TrigGap0.952us_Th350_CFDDelay0.304us_Scale7_for_X_efficiency_cal.root");
	//double fitrange_min = 44.6, fitrange_max = 46.05;
	double fitrange_min = 7.3, fitrange_max = 9.3;
	TH1D* histo = new TH1D("histo", "histo", 1000, 0, 50); // bin width = 50/1000 = 0.05 keV
	double binwidth = histo->GetBinWidth(1);
	TCanvas* canvaspeak = new TCanvas("LEGe", "LEGe", 1300, 700);
	canvaspeak->cd();
	canvaspeak->SetTopMargin(0.035);
	canvaspeak->SetRightMargin(0.04);
	canvaspeak->SetLeftMargin(0.14);
	canvaspeak->SetBottomMargin(0.20);
	canvaspeak->SetFrameLineWidth(3);
	gStyle->SetFrameLineWidth(3);

	tree->Draw("LEGe_e>>histo", "LEGe_e>0.11&&LEGe_e<50&&MSD12_e>30&&MSD12_e<1000&&(MSD12_e+MSD26_e)>30&&(MSD12_e+MSD26_e)<9000", ""); // All proton

	//tree->Draw("LEGe_e>>histo", "LEGe_e>0.11&&LEGe_e<50&&MSD12_e>1200&&MSD12_e<4000&&(MSD12_e+MSD26_e)>1200&&(MSD12_e+MSD26_e)<8000", ""); // All alpha

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
	histo->GetXaxis()->SetTitleOffset(1.2);
	histo->GetYaxis()->SetTitleOffset(0.9);
	histo->GetXaxis()->SetTitleSize(0.08);
	histo->GetYaxis()->SetTitleSize(0.08);
	histo->GetYaxis()->SetTickLength(0.015);
	histo->GetXaxis()->SetRangeUser(6, 10);
	int ymax = histo->GetMaximum();
	//histo->GetYaxis()->SetRangeUser(0, ymax * 1.24);
	histo->GetYaxis()->SetNdivisions(505);
	histo->SetBinErrorOption(TH1::kPoisson);
	histo->Draw("e");

	TPaveText* textgaus = new TPaveText(0.52, 0.88, 0.58, 0.95, "brNDC");//left, down, right, up
	textgaus->SetBorderSize(0);
	textgaus->SetFillColor(0);
	textgaus->SetTextAlign(12);
	textgaus->SetTextFont(132);
	textgaus->SetTextSize(0.08);
	//text->SetTextColor(2);
	textgaus->AddText("Cu");
	textgaus->Draw();
	textgaus = new TPaveText(0.65, 0.52, 0.71, 0.59, "brNDC");//left, down, right, up
	textgaus->SetBorderSize(0);
	textgaus->SetFillColor(0);
	textgaus->SetTextAlign(12);
	textgaus->SetTextFont(132);
	textgaus->SetTextSize(0.08);
	textgaus->AddText("Zn");
	textgaus->Draw();

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
	gtotal->SetParLimits(2, 0.1, 90);//Const1*bin
	gtotal->SetParLimits(3, 7.9, 8.1);//Mean_1
 	gtotal->SetParLimits(4, 0.103, 0.112);//Sigma_1
	gtotal->SetParLimits(5, 0.1, 90);//Const2*bin
	gtotal->SetParLimits(6, 8.5, 8.7);//Mean_2
	gtotal->SetParLimits(7, 0.103, 0.119);//Sigma_2
	// Perform the fit
	histo->Fit("gtotal", "MLE", "", fitrange_min, fitrange_max);
	histo->Fit("gtotal", "MLE", "", fitrange_min, fitrange_max);
	histo->Fit("gtotal", "MLE", "", fitrange_min, fitrange_max);
	gtotal->SetLineWidth(3);
	gtotal->Draw("same");
	//canvaspeak->SaveAs("F:/e21010/pxct/Fig_PXCT_60Ga_Geant4_Total_PX_Fit.eps");
	//canvaspeak->SaveAs("F:/e21010/pxct/Fig_PXCT_60Ga_Geant4_Total_PX_Fit.png");

	// print all parameters names values with errors
	for (int i = 0; i < gtotal->GetNpar(); i++)
	{
		cout << "Par" << i << " " << gtotal->GetParName(i) << "	" << gtotal->GetParameter(i) << "	" << gtotal->GetParError(i) << endl;
	}
	double Chi2_Inflated = sqrt(gtotal->GetChisquare() / gtotal->GetNDF());
	cout << "Chi2/NDF	" << pow(Chi2_Inflated,2) << endl;
	cout << "Prob	" << gtotal->GetProb() << endl;
	cout << "binwidth	" << binwidth << endl;


	if (Chi2_Inflated < 1.0)
	{
		Chi2_Inflated = 1.0;
	}
	// Extract the peak counts and uncertainties from the fit
	double Cu_Ka_counts = gtotal->GetParameter(2) / binwidth;
	double Zn_Ka_counts = gtotal->GetParameter(5) / binwidth;
	double Cu_Ka_counts_Uncertainty = (gtotal->GetParError(2) / binwidth) * Chi2_Inflated;
	double Zn_Ka_counts_Uncertainty = (gtotal->GetParError(5) / binwidth) * Chi2_Inflated;

	cout << "Cu_Ka_counts	" << Cu_Ka_counts << "	" << Cu_Ka_counts_Uncertainty << endl;
	cout << "Zn_Ka_counts	" << Zn_Ka_counts << "	" << Zn_Ka_counts_Uncertainty << endl;

	// Define the atomic parameters
	double Ratio_Ka_Emission_ZnCu = 1.0715;
	double Ratio_Ka_Emission_ZnCu_Uncertainty = Ratio_Ka_Emission_ZnCu * 0.0000010;
	double Ratio_Detection_Efficiency_ZnCu = 1.0458;
	double Ratio_Detection_Efficiency_ZnCu_Uncertainty = Ratio_Detection_Efficiency_ZnCu * 0.00000010;
	double Lifetime_Zn_K_shell_vacancy = 0.406; // in fs
	double Lifetime_Zn_K_shell_vacancy_Uncertainty = Lifetime_Zn_K_shell_vacancy * 0.000000025; // in fs

	// Calculate the real peak counts
	double Real_Cu_Ka_counts = Cu_Ka_counts * Ratio_Ka_Emission_ZnCu * Ratio_Detection_Efficiency_ZnCu;
	double Real_Zn_Ka_counts = Zn_Ka_counts;

	// Propagate the uncertainties for the real peak counts
	double Real_Cu_Ka_counts_Uncertainty = sqrt(
		pow(Cu_Ka_counts_Uncertainty / Cu_Ka_counts, 2) +
		pow(Ratio_Ka_Emission_ZnCu_Uncertainty / Ratio_Ka_Emission_ZnCu, 2) +
		pow(Ratio_Detection_Efficiency_ZnCu_Uncertainty / Ratio_Detection_Efficiency_ZnCu, 2)) * Real_Cu_Ka_counts;
	double Real_Zn_Ka_counts_Uncertainty = Zn_Ka_counts_Uncertainty;

	// Calculate the ratio and its uncertainty
	double Ratio_Cu_Zn = Real_Cu_Ka_counts / Real_Zn_Ka_counts;
	double Ratio_Cu_Zn_Uncertainty = sqrt(
		pow(Real_Cu_Ka_counts_Uncertainty / Real_Cu_Ka_counts, 2) +
		pow(Real_Zn_Ka_counts_Uncertainty / Real_Zn_Ka_counts, 2)) * Ratio_Cu_Zn;

	// Calculate the lifetime of the proton emitting state and its uncertainty
	double Lifetime_proton_emitting_state = Lifetime_Zn_K_shell_vacancy / Ratio_Cu_Zn;
	double Lifetime_proton_emitting_state_Uncertainty = sqrt(
		pow(Ratio_Cu_Zn_Uncertainty / Ratio_Cu_Zn, 2) +
		pow(Lifetime_Zn_K_shell_vacancy_Uncertainty / Lifetime_Zn_K_shell_vacancy, 2)) * Lifetime_proton_emitting_state;
	
	// Output the results
	cout << "Real_Cu_Ka_counts	" << Real_Cu_Ka_counts << "	" << Real_Cu_Ka_counts_Uncertainty << endl;
	cout << "Real_Zn_Ka_counts	" << Real_Zn_Ka_counts << "	" << Real_Zn_Ka_counts_Uncertainty << endl;
	cout << "Ratio_Zn_Cu	" << Ratio_Cu_Zn << "	" << Ratio_Cu_Zn_Uncertainty << endl;
	cout << "Lifetime_proton_emitting_state	" << Lifetime_proton_emitting_state << "	" << Lifetime_proton_emitting_state_Uncertainty << endl;

}
	