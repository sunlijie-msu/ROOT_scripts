#include <iostream>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <time.h>
#include <map>
#include <TROOT.h>
#include <TSystem.h>
#include <TStyle.h>
#include <TCutG.h>
#include "TChain.h"
#include "TStyle.h"
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
#include "TLatex.h"
using namespace std;
void average_half_life_pxct_237Np()// graph-pol0 fit several points to get the weighted average
{
	TCanvas* canvas;
	TGraph* graph, * graph_uncertainty_low, * graph_uncertainty_high;//TGraph
	double x_value[100], x_error[100];
	double y_value[100], y_error[100];
	int num_values;
	double x_point[1000];// set location of point for uncertainty line drawing, can be denser than the x_value
	for (int ii = 0; ii < 1000; ii++)
	{
		x_point[ii] = ii / 10.0;
	}
	double y_point_error[1000];  // error on the function at point x0 for uncertainty line drawing
	double y_point_lower_bound[1000], y_point_higher_bound[1000];
	double p0, p0err, parChi2, parNDF, p_value;

	// Alternative (disabled): read the same values from a text file instead of the hard-coded arrays.
	// string line;
	// char h_name[100], infile_name[100], paraprint[100];
	// sprintf(infile_name, "%s", "D:/X/out/weightedaverage.dat");
	// cout << infile_name << endl;
	// ifstream infile(infile_name, ios::in);//The data that need to be fitted
	// int number = 0;
	// while (getline(infile, line))
	// {
	//     stringstream(line) >> y_value[number] >> y_error[number];
	//     x_value[number] = number + 1;
	//     x_error[number] = 0;
	//     cout << x_value[number] << '\t' << y_value[number] << '\t' << y_error[number] << endl;
	//     number++;
	// }

	char h_name[100], paraprint[100];
	
//63 ns {I5} (Beling_PR1952), \\
//62 ns {I3} (Vartapetian_AP1958), \\
//63 ns {I5} (Unik_Thesis1960), \\
//64.2 ns {I20} (Oberhofer_Thesis1968), \\
//66.7 ns {I7} (Garg_ZP1971), \\
//66.9 ns {I10} (McBeth_NIM1972), \\
//67.8 ns {I10} (McBeth_NIM1972), \\
//68.3 ns {I2} (Miller_NIM1972), \\
//68.5 ns {I4} (Bishop_NIM1974), \\
//70 ns {I4} (Finck_NIMA1985), \\
//67.7 ns {I1} (Vretenar_AJP2019), \\
//67.60 ns {I25} (Dutsov_ARI2021), \\
//67.86 ns {I9} (Takacs_ARI2021), \\
//67.60 ns {I20} (Santos_NPA2023), \\
//68.03 ns {I7} (Present_work)

	const int number = 15; // Half-life (ns); uncertainty (ns) = ENSDF last-digit uncertainty {In}
	double y_data[number] = { // half-life (ns)
		63.0,  //  1  Beling_PR1952
		62.0,  //  2  Vartapetian_AP1958
		63.0,  //  3  Unik_Thesis1960
		64.2,  //  4  Oberhofer_Thesis1968
		66.7,  //  5  Garg_ZP1971
		66.9,  //  6  McBeth_NIM1972 (66.9)
		67.8,  //  7  McBeth_NIM1972 (67.8)
		68.3,  //  8  Miller_NIM1972
		68.5,  //  9  Bishop_NIM1974
		70.0,  // 10  Finck_NIMA1985
		67.7,  // 11  Vretenar_AJP2019
		67.60, // 12  Dutsov_ARI2021
		67.86, // 13  Takacs_ARI2021
		67.60, // 14  Santos_NPA2023
		68.03  // 15  Present_work
	};
	double y_data_err[number] = { // uncertainty (ns)
		5.0,   //  1  Beling_PR1952
		3.0,   //  2  Vartapetian_AP1958
		5.0,   //  3  Unik_Thesis1960
		2.0,   //  4  Oberhofer_Thesis1968
		0.7,   //  5  Garg_ZP1971
		1.0,   //  6  McBeth_NIM1972 (66.9)
		1.0,   //  7  McBeth_NIM1972 (67.8)
		0.2,   //  8  Miller_NIM1972
		0.4,   //  9  Bishop_NIM1974
		4.0,   // 10  Finck_NIMA1985
		0.1,   // 11  Vretenar_AJP2019
		0.25,  // 12  Dutsov_ARI2021
		0.09,  // 13  Takacs_ARI2021
		0.20,  // 14  Santos_NPA2023
		0.07   // 15  Present_work
	};
	for (int ii = 0; ii < number; ii++)
	{
		x_value[ii] = ii + 1;
		x_error[ii] = 0;
		y_value[ii] = y_data[ii];
		y_error[ii] = y_data_err[ii];
	}

	sprintf(h_name, "%s", "weighted_average");
	canvas = new TCanvas(h_name, h_name, 1300, 800);//create canvas
	canvas->cd();//enter canvas
	canvas->SetTopMargin(0.03);
	canvas->SetRightMargin(0.04);
	canvas->SetLeftMargin(0.11);
	canvas->SetBottomMargin(0.15);
	canvas->SetFrameLineWidth(3);
	//canvascali[i]->SetGrid();//show grid
	//graph[i]=new TGraph(peaknum,Number[i],input_value[i]);//TGraph *gr1=new TGraph(n,x,y);
	graph = new TGraphErrors(number, x_value, y_value, x_error, y_error);//error bars TGraph(n,x,y,ex,ey);
	graph->SetTitle(h_name);
	graph->GetXaxis()->SetTitle("Measurements");//axis title
	graph->GetYaxis()->SetTitle("Half-life (ns)");//axis title
	graph->GetXaxis()->CenterTitle();//center title
	graph->GetYaxis()->CenterTitle();//center title
	graph->SetTitle("");
	graph->SetStats(0);
	graph->GetXaxis()->SetLabelFont(132);//axis label font
	graph->GetYaxis()->SetLabelFont(132);//axis label font
	graph->GetXaxis()->SetLabelSize(0.06);
	graph->GetYaxis()->SetLabelSize(0.06);
	graph->GetYaxis()->SetTickLength(0.015);
	graph->GetXaxis()->SetTitleFont(132);//axis title font
	graph->GetYaxis()->SetTitleFont(132);//axis title font
	graph->GetXaxis()->SetTitleOffset(1.0);//axis title offset
	graph->GetYaxis()->SetTitleOffset(0.72);//axis title offset
	graph->GetXaxis()->SetTitleSize(0.07);
	graph->GetYaxis()->SetTitleSize(0.07);
	//		graph[i]->GetYaxis()->SetRangeUser(1,8);
	graph->GetYaxis()->SetRangeUser(54, 78);
	graph->GetYaxis()->SetNdivisions(505);
	graph->GetXaxis()->SetNdivisions(510);//n = n1 + 100*n2 + 10000*n3
	graph->SetMarkerStyle(6);
	graph->SetMarkerColor(kBlack);
	graph->SetLineWidth(2);
	TF1* pol0 = new TF1("pol0", "pol0", 0, 200);//pol0 constant fit
	pol0->SetNpx(2000);
	pol0->SetParNames("p0");
	pol0->SetLineColor(kRed);
	pol0->SetLineWidth(1);

	//graph->Fit("pol1");//pol1 can be used directly without TF1 constructor in CINT
	//graph->Fit("pol0","ME");
	TFitResultPtr Fit_result_pointer = graph->Fit("pol0", "MES");//"S" means the result of the fit is returned in the TFitResultPtr
	//"E" Perform better errors estimation using the Minos technique.
	//"M" Improve fit results, by using the IMPROVE algorithm of TMinuit.
	//If the option "S" is instead used, TFitResultPtr contains the TFitResult and behaves as a smart pointer to it.
	//The fit parameters, error and chi2 (but not covariance matrix) can be retrieved also from the fitted function.

	// Uncertainty Band
	TH1D* h_confidence_interval = new TH1D("h_confidence_interval", "h_confidence_interval", 10000, 0, 100);//Create a histogram to hold the confidence intervals
	// GetConfidenceIntervals fills only the bin ERRORS, so the bin CONTENTS must hold the fit curve;
	// otherwise the band is painted around y=0 and falls outside the plotted range.
	for (int ib = 1; ib <= h_confidence_interval->GetNbinsX(); ib++)
		h_confidence_interval->SetBinContent(ib, pol0->Eval(h_confidence_interval->GetBinCenter(ib)));
	TVirtualFitter* fitter = TVirtualFitter::GetFitter();//The method TVirtualFitter::GetFitter())->Get the parameters of your fitting function after having it fitted to an histogram.
	fitter->GetConfidenceIntervals(h_confidence_interval, 0.683);//By default the intervals are inflated using the chi2/ndf value of the fit if a chi2 fit is performed
	//confidence interval for the colored band: 1Sigma confidence interval: P=0.683, 1Sigma confidence interval: P=0.95, 3Sigma confidence interval: P=0.997
	//h_confidence_interval will contain the CL result that you can draw on top of your fitted graph.
	//where h_confidence_interval will hold the errors and could superimpose it on the same canvas where you plot central values.
	h_confidence_interval->SetStats(kFALSE);
	h_confidence_interval->SetFillColor(kAzure - 9);
	graph->Draw("AP");//"A": Axis are drawn around the graph, "P": The current marker is plotted at each point
	h_confidence_interval->Draw("e3 same"); // plot the uncertainty band

	// Uncertainty Lines
	Fit_result_pointer->GetConfidenceIntervals(1000, 1, 1, x_point, y_point_error, 0.683, false);//get the error of y at x points of interest, this command is independent of the colored band
	//(Number of x points, 1, 1, x, err, confidence level, false); norm is a flag to control if the intervals need to be inflated by the chi2/ndf value. true is inflated = external uncertainty, false is not inflated = internal uncertainty.
	for (int ii = 0; ii < 1000; ii++)
	{
		y_point_lower_bound[ii] = pol0->Eval(x_point[ii]) - y_point_error[ii];
		y_point_higher_bound[ii] = pol0->Eval(x_point[ii]) + y_point_error[ii];
	}
	graph_uncertainty_low = new TGraph(1000, x_point, y_point_lower_bound);//error bars TGraph(n,x,y,ex,ey);
	graph_uncertainty_low->SetLineColor(4);
	graph_uncertainty_low->SetLineWidth(1);
	//graph_uncertainty_low->Draw("C");//"C" A smooth line is drawn
	graph_uncertainty_high = new TGraph(1000, x_point, y_point_higher_bound);//error bars TGraph(n,x,y,ex,ey);
	graph_uncertainty_high->SetLineColor(4);
	graph_uncertainty_high->SetLineWidth(1);
	//graph_uncertainty_high->Draw("C");//"C" A smooth line is drawn
	h_confidence_interval->Draw("e3 same");//uncertainty band from the fit confidence intervals
	pol0->Draw("C same");//fit line
	graph->Draw("PE same");//data points and error bars above the band

	// Write Fit Results on Canvas
	p0 = pol0->GetParameter(0);
	p0err = pol0->GetParError(0);
	parChi2 = pol0->GetChisquare();
	parNDF = pol0->GetNDF();
	p_value = pol0->GetProb();//This probability is not the "probability that your fit is good." If you did many fake experiments (draw many random samples of data points from the assumed distribution (your fit function)), this is the percentage of experiments that would give chi2 values >= to the one you got in this experiment.

	// Fit results: plain text, no box, matching peakfit style
	double para_textsize = 0.05; // font size, fraction of pad height
	double para_dy = 1.2 * para_textsize; // line spacing, fraction of pad height

	// Block 1 (p0, averages, inflation factor): left-aligned. fit_para_x = left edge of text, fit_para_y = vertical center of FIRST line
	double fit_para_x = 0.15, fit_para_y = 0.92;
	TLatex* textgaus = new TLatex();
	textgaus->SetNDC();
	textgaus->SetTextFont(132);
	textgaus->SetTextSize(para_textsize);
	textgaus->SetTextAlign(12); // left, vertically centered
	sprintf(paraprint, "Weighted average = %.2f #pm %.2f", p0, p0err);
	textgaus->DrawLatex(fit_para_x, fit_para_y, paraprint);
	sprintf(paraprint, "Internal = %.2f #pm %.2f", p0, p0err);
	textgaus->DrawLatex(fit_para_x, fit_para_y - para_dy, paraprint);
	sprintf(paraprint, "sqrt(#chi^{2}/#it{NDF}) = %.2f", sqrt(parChi2 / parNDF));
	textgaus->DrawLatex(fit_para_x, fit_para_y - 2 * para_dy, paraprint);
	if (sqrt(parChi2 / parNDF) > 1)
	{
		sprintf(paraprint, "External = %.2f #pm %.2f", p0, p0err * sqrt(parChi2 / parNDF));
		textgaus->DrawLatex(fit_para_x, fit_para_y - 3 * para_dy, paraprint);
	}

	// Block 2 (chi2/NDF, p-value): right-aligned. stat_x = right edge of text, stat_y = vertical center of FIRST line
	double stat_x = 0.95, stat_y = 0.31;
	TLatex* textstat = new TLatex();
	textstat->SetNDC();
	textstat->SetTextFont(132);
	textstat->SetTextSize(para_textsize);
	textstat->SetTextAlign(32); // right, vertically centered
	sprintf(paraprint, "#chi^{2}/#it{NDF} = %.1f / %.0f = %.2f", parChi2, parNDF, parChi2 / parNDF);
	textstat->DrawLatex(stat_x, stat_y, paraprint);
	sprintf(paraprint, "#it{p}-value = %.3f", p_value);
	textstat->DrawLatex(stat_x, stat_y - para_dy, paraprint);
	TString outFile = h_name;
	if (!outFile.Contains(":/") && !outFile.BeginsWith("/"))
	{
		outFile = TString("D:/X/out/") + outFile;
	}
	outFile += ".png";
	gPad->RedrawAxis();//repaint frame and axes above the confidence band; also required so the pad is painted before saving
	canvas->SaveAs(outFile);//save figure
}// main