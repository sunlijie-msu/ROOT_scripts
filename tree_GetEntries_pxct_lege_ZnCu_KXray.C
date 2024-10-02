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
void tree_GetEntries_pxct_lege_ZnCu_KXray() // tree GetEntries for PXCT LEGe Zn Cu K X-ray Ratio Plot
{
	TFile* _file0 = TFile::Open("F:/e21010/pxct/ExG4_60Ga_20days_9000pps.root");
	TTree* tree = (TTree*)_file0->Get("tree");

	// Define atomic parameters
	double Ratio_Ka_Emission_CuZn = 1.0715;
	double Ratio_Ka_Emission_CuZn_Uncertainty = Ratio_Ka_Emission_CuZn * 1e-10;
	double Ratio_Detection_Efficiency_CuZn = 1.0458;
	double Ratio_Detection_Efficiency_CuZn_Uncertainty = Ratio_Detection_Efficiency_CuZn * 1e-10;
	double Lifetime_Zn_K_shell_vacancy = 0.406; // in fs
	double Lifetime_Zn_K_shell_vacancy_Uncertainty = Lifetime_Zn_K_shell_vacancy * 2.5e-10; // in fs

	// Define Ep distribution flexibly
	vector<double> Ep;
	vector<double> DEp;
	double Ep_start = 6600;
	double Ep_end = 8000;
	double val = Ep_start;

	while (val <= Ep_end)
	{
		double Ep_step;

		// Adjust Ep_step based on current Ep value
		if (val >= 700 && val < 950)
		{
			Ep_step = 600; // Larger step when Ep is small
		}
		else if (val >= 950 && val <= 6000)
		{
			Ep_step = 400; // Smaller step
		}
		else if (val > 6000 && val <= 9000)
		{
			Ep_step = 600;
		}

		Ep.push_back(val);
		DEp.push_back(Ep_step / 2.0);

		val += Ep_step;
	}

	// Prepare vectors for results
	size_t N = Ep.size();
	vector<double> Cu_counts(N), Cu_counts_Uncertainty(N);
	vector<double> Zn_counts(N), Zn_counts_Uncertainty(N);
	vector<double> Real_Cu_counts(N), Real_Cu_counts_Uncertainty(N);
	vector<double> Real_Zn_counts(N), Real_Zn_counts_Uncertainty(N);
	vector<double> Ratio_Cu_Zn(N), Ratio_Cu_Zn_Uncertainty(N);
	vector<double> Lifetime_proton_emitting_state(N), Lifetime_proton_emitting_state_Uncertainty(N);

	cout << "Ep_min\tEp_max\t"
		<< "Cu_counts_integral\tZn_counts_integral\t"
		<< "bg_left_counts\tbg_right_counts\t"
		<< "background_counts\tbackground_uncertainty\t"
		<< "Cu_counts\tCu_counts_Uncertainty\t"
		<< "Zn_counts\tZn_counts_Uncertainty\t"
		<< "Real_Cu_counts\tReal_Cu_counts_Uncertainty\t"
		<< "Real_Zn_counts\tReal_Zn_counts_Uncertainty\t"
		<< "Ratio_Cu_Zn\tRatio_Cu_Zn_Uncertainty\t"
		<< "Lifetime_proton_emitting_state\tLifetime_proton_emitting_state_Uncertainty"
		<< "\n";

	for (size_t i = 0; i < N; ++i)
	{
		double Ep_min = Ep[i] - DEp[i];
		double Ep_max = Ep[i] + DEp[i];

		// Conditions
		TString condition_Cu = Form("LEGe_e > 7.6 && LEGe_e < 8.3 && MSD12_e > 10 && MSD12_e < 1000 && "
			"(MSD12_e + MSD26_e) > %f && (MSD12_e + MSD26_e) < %f", Ep_min, Ep_max);
		TString condition_Zn = Form("LEGe_e > 8.3 && LEGe_e < 9.0 && MSD12_e > 10 && MSD12_e < 1000 && "
			"(MSD12_e + MSD26_e) > %f && (MSD12_e + MSD26_e) < %f", Ep_min, Ep_max);
		TString condition_bg_left = Form("LEGe_e > 4.0 && LEGe_e < 7.0 && MSD12_e > 10 && MSD12_e < 1000 && "
			"(MSD12_e + MSD26_e) > %f && (MSD12_e + MSD26_e) < %f", Ep_min, Ep_max);
		TString condition_bg_right = Form("LEGe_e > 9.5 && LEGe_e < 12.5 && MSD12_e > 10 && MSD12_e < 1000 && "
			"(MSD12_e + MSD26_e) > %f && (MSD12_e + MSD26_e) < %f", Ep_min, Ep_max);

		// Counts
		double Cu_counts_integral = tree->GetEntries(condition_Cu);
		double Zn_counts_integral = tree->GetEntries(condition_Zn);
		double bg_left_counts = tree->GetEntries(condition_bg_left);
		double bg_right_counts = tree->GetEntries(condition_bg_right);
		double background_counts = (bg_left_counts + bg_right_counts) / 6.0 * 0.7;
		if (background_counts <= 0) background_counts = 1e-10;
		double background_uncertainty = sqrt(background_counts);

		cout << Ep_min << "\t" << Ep_max << "\t"
			<< Cu_counts_integral << "\t" << Zn_counts_integral << "\t"
			<< bg_left_counts << "\t" << bg_right_counts << "\t"
			<< background_counts << "\t" << background_uncertainty << "\t";

		// Subtract background
		Cu_counts[i] = Cu_counts_integral - background_counts;
		if (Cu_counts[i] <= 0) Cu_counts[i] = 1e-10;
		Cu_counts_Uncertainty[i] = sqrt(pow(sqrt(Cu_counts_integral), 2) + pow(background_uncertainty, 2));

		Zn_counts[i] = Zn_counts_integral - background_counts;
		if (Zn_counts[i] <= 0) Zn_counts[i] = 1e-10;
		Zn_counts_Uncertainty[i] = sqrt(pow(sqrt(Zn_counts_integral), 2) + pow(background_uncertainty, 2));

		// Real counts
		Real_Cu_counts[i] = Cu_counts[i] * Ratio_Ka_Emission_CuZn * Ratio_Detection_Efficiency_CuZn;
		Real_Cu_counts_Uncertainty[i] = Real_Cu_counts[i] * sqrt(
			pow(Cu_counts_Uncertainty[i] / Cu_counts[i], 2) +
			pow(Ratio_Ka_Emission_CuZn_Uncertainty / Ratio_Ka_Emission_CuZn, 2) +
			pow(Ratio_Detection_Efficiency_CuZn_Uncertainty / Ratio_Detection_Efficiency_CuZn, 2)
		);
		Real_Zn_counts[i] = Zn_counts[i];
		Real_Zn_counts_Uncertainty[i] = Zn_counts_Uncertainty[i];

		// Ratio and lifetime
		Ratio_Cu_Zn[i] = Real_Cu_counts[i] / Real_Zn_counts[i];
		Ratio_Cu_Zn_Uncertainty[i] = Ratio_Cu_Zn[i] * sqrt(	pow(Real_Cu_counts_Uncertainty[i] / Real_Cu_counts[i], 2) + pow(Real_Zn_counts_Uncertainty[i] / Real_Zn_counts[i], 2));
		Lifetime_proton_emitting_state[i] = Lifetime_Zn_K_shell_vacancy / Ratio_Cu_Zn[i];
		Lifetime_proton_emitting_state_Uncertainty[i] = Lifetime_proton_emitting_state[i] * sqrt(pow(Ratio_Cu_Zn_Uncertainty[i] / Ratio_Cu_Zn[i], 2) + pow(Lifetime_Zn_K_shell_vacancy_Uncertainty / Lifetime_Zn_K_shell_vacancy, 2));


		// Output
		cout << Cu_counts[i] << "\t" << Cu_counts_Uncertainty[i] << "\t"
			   << Zn_counts[i] << "\t" << Zn_counts_Uncertainty[i] << "\t"
			   << Real_Cu_counts[i] << "\t" << Real_Cu_counts_Uncertainty[i] << "\t"
			   << Real_Zn_counts[i] << "\t" << Real_Zn_counts_Uncertainty[i] << "\t"
			   << Ratio_Cu_Zn[i] << "\t" << Ratio_Cu_Zn_Uncertainty[i] << "\t"
			   << Lifetime_proton_emitting_state[i] << "\t" << Lifetime_proton_emitting_state_Uncertainty[i]
			   << "\n";
	}

	TCanvas* canvaspeak = new TCanvas("RX", "RX", 1300, 700);
	canvaspeak->cd();
	canvaspeak->SetTopMargin(0.035);
	canvaspeak->SetRightMargin(0.04);
	canvaspeak->SetLeftMargin(0.14);
	canvaspeak->SetBottomMargin(0.20);
	canvaspeak->SetFrameLineWidth(3);
	gStyle->SetFrameLineWidth(3);

	TGraphErrors* graph = new TGraphErrors(N, &Ep[0], &Ratio_Cu_Zn[0], 0, &Ratio_Cu_Zn_Uncertainty[0]);

	graph->SetTitle("");
	graph->SetStats(0);
	graph->GetXaxis()->SetTitle("Proton Energy (keV)");
	graph->GetYaxis()->SetTitle("#it{R}_{Cu/Zn}");
	graph->GetXaxis()->CenterTitle();
	graph->GetYaxis()->CenterTitle();
	graph->GetXaxis()->SetLabelFont(132);
	graph->GetYaxis()->SetLabelFont(132);
	graph->GetXaxis()->SetLabelSize(0.07);
	graph->GetYaxis()->SetLabelSize(0.07);
	graph->GetXaxis()->SetTitleFont(132);
	graph->GetYaxis()->SetTitleFont(132);
	graph->GetXaxis()->SetTitleOffset(1.2);
	graph->GetYaxis()->SetTitleOffset(0.9);
	graph->GetXaxis()->SetTitleSize(0.08);
	graph->GetYaxis()->SetTitleSize(0.08);
	//graph->GetXaxis()->SetTickLength(0.015);
	graph->GetYaxis()->SetTickLength(0.015);
	graph->GetXaxis()->SetNdivisions(505);
	graph->GetYaxis()->SetNdivisions(505);
	//graph->GetXaxis()->SetLimits(0, 7000);
	graph->GetYaxis()->SetRangeUser(0, 50);

	// Draw graph
	graph->SetMarkerStyle(21);
	graph->SetMarkerSize(1);
	graph->SetMarkerColor(kBlack);
	graph->SetLineColor(kBlack);
	graph->SetLineWidth(3);
	
	graph->Draw("AP"); // "A" Axis are drawn around the graph, "P" to draw only markers, "L" to draw only lines, "PL" to draw both

	canvaspeak->SaveAs("F:/e21010/pxct/Fig_PXCT_60Ga_Xratio.png");
	canvaspeak->SaveAs("F:/e21010/pxct/Fig_PXCT_60Ga_Xratio.eps");
}

