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
using namespace std;
void average()// graph-pol0 fit several points to get the weighted average
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

	string line;
	char h_name[100], infile_name[100], paraprint[100];
	sprintf(infile_name, "%s", "D:/X/out/weightedaverage.dat");
	cout << infile_name << endl;
	ifstream infile(infile_name, ios::in);//The data that need to be fitted
	int number = 0;
	while (getline(infile, line))
	{
		stringstream(line) >> y_value[number] >> y_error[number];
		x_value[number] = number + 1;
		x_error[number] = 0;
		cout << x_value[number] << '	' << y_value[number] << '	' << y_error[number] << endl;
		number++;
	}

	sprintf(h_name, "%s", "weighted_average");
	canvas = new TCanvas(h_name, h_name, 1300, 800);//建立画布
	canvas->cd();//进入画布
	canvas->SetTopMargin(0.03);
	canvas->SetRightMargin(0.02);
	canvas->SetLeftMargin(0.13);
	canvas->SetBottomMargin(0.13);
	canvas->SetFrameLineWidth(2);
	//canvascali[i]->SetGrid();//显示网格
	//graph[i]=new TGraph(peaknum,Number[i],input_value[i]);//TGraph *gr1=new TGraph(n,x,y);
	graph = new TGraphErrors(number, x_value, y_value, x_error, y_error);//画error bars TGraph(n,x,y,ex,ey);
	graph->SetTitle(h_name);
	graph->GetXaxis()->SetTitle("Measurements");//轴名
	graph->GetYaxis()->SetTitle("Half-life (ns)");//轴名
	graph->GetXaxis()->CenterTitle();//居中
	graph->GetYaxis()->CenterTitle();//居中
	graph->SetTitle("");
	graph->GetXaxis()->SetLabelFont(132);//坐标字体
	graph->GetYaxis()->SetLabelFont(132);//坐标字体
	graph->GetXaxis()->SetLabelSize(0.06);
	graph->GetYaxis()->SetLabelSize(0.06);
	graph->GetXaxis()->SetTitleFont(132);//轴名字体
	graph->GetYaxis()->SetTitleFont(132);//轴名字体
	graph->GetXaxis()->SetTitleOffset(0.8);//轴名偏移
	graph->GetYaxis()->SetTitleOffset(0.9);//轴名偏移
	graph->GetXaxis()->SetTitleSize(0.07);
	graph->GetYaxis()->SetTitleSize(0.07);
	//		graph[i]->GetYaxis()->SetRangeUser(1,8);
	graph->GetYaxis()->SetNdivisions(505);
	graph->GetXaxis()->SetNdivisions(100);//n = n1 + 100*n2 + 10000*n3
	graph->SetMarkerStyle(21);
	graph->SetMarkerColor(1);
	graph->SetLineWidth(2);
	TF1* pol0 = new TF1("pol0", "pol0", 0, 200);//多项式拟合
	pol0->SetNpx(2000);
	pol0->SetParNames("p0");
	pol0->SetLineWidth(3);

	//graph->Fit("pol1");//pol1 can be used directly without TF1 constructor in CINT
	//graph->Fit("pol0","ME");
	TFitResultPtr Fit_result_pointer = graph->Fit("pol0", "MES");//"S" means the result of the fit is returned in the TFitResultPtr
	//“E” Perform better errors estimation using the Minos technique.
	//“M” Improve fit results, by using the IMPROVE algorithm of TMinuit.
	//If the option "S" is instead used, TFitResultPtr contains the TFitResult and behaves as a smart pointer to it.
	//The fit parameters, error and chi2 (but not covariance matrix) can be retrieved also from the fitted function.

	// Uncertainty Band
	TH1D* h_confidence_interval = new TH1D("h_confidence_interval", "h_confidence_interval", 10000, 0, 100);//Create a histogram to hold the confidence intervals
	TVirtualFitter* fitter = TVirtualFitter::GetFitter();//The method TVirtualFitter::GetFitter())->Get the parameters of your fitting function after having it fitted to an histogram.
	fitter->GetConfidenceIntervals(h_confidence_interval, 0.683);//By default the intervals are inflated using the chi2/ndf value of the fit if a chi2 fit is performed
	//confidence interval for the colored band: 1σ confidence interval: P=0.683, 1σ confidence interval: P=0.95, 3σ confidence interval: P=0.997
	//h_confidence_interval will contain the CL result that you can draw on top of your fitted graph.
	//where h_confidence_interval will hold the errors and could superimpose it on the same canvas where you plot central values.
	h_confidence_interval->SetStats(kFALSE);
	h_confidence_interval->SetFillColor(kRed - 10);
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

	graph->Draw("P same");//draw the points again above error band for better visibility

	// Write Fit Results on Canvas
	p0 = pol0->GetParameter(0);
	p0err = pol0->GetParError(0);
	parChi2 = pol0->GetChisquare();
	parNDF = pol0->GetNDF();
	p_value = pol0->GetProb();//This probability is not the “probability that your fit is good.” If you did many fake experiments (draw many random samples of data points from the assumed distribution (your fit function)), this is the percentage of experiments that would give χ2 values ≥ to the one you got in this experiment.

	TPaveText* textpol2 = new TPaveText(0.32, 0.10, 0.90, 0.35, "brNDC");//left, down, right, up
	textpol2->SetBorderSize(1);
	textpol2->SetFillColor(0);
	textpol2->SetTextAlign(12);//align = 10*HorizontalAlign + VerticalAlign, 12 means水平左对齐、垂直居中对齐
	textpol2->SetTextFont(132);//font = 10 * fontID + precision, 12+2,12 means Symbol; 13+2, 13 means News Times Roman
	sprintf(paraprint, "Chi2 = %.2f", parChi2);
	textpol2->AddText(paraprint);
	sprintf(paraprint, "NDF = %.0f", parNDF);
	textpol2->AddText(paraprint);
	sprintf(paraprint, "p-val = %.2f", p_value);
	textpol2->AddText(paraprint);
	sprintf(paraprint, "Average with internal uncertainty = %.4f+/-%.4f", p0, p0err);
	textpol2->AddText(paraprint);
	sprintf(paraprint, "Inflation factor = sqrt(%.2f/%.0f)=%.2f", parChi2, parNDF, sqrt(parChi2 / parNDF));
	if (sqrt(parChi2 / parNDF) > 1)
	{
		textpol2->AddText(paraprint);
		sprintf(paraprint, "Average with external uncertainty = %.4f+/-%.4f", p0, p0err * sqrt(parChi2 / parNDF));
	}
	textpol2->AddText(paraprint);
	textpol2->Draw();
	sprintf(h_name, "%s%s%s", "D:/X/out/Chi2/", h_name, ".png");
	canvas->SaveAs(h_name);//存图
}// main