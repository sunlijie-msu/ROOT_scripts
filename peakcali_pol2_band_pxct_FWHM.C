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
void peakcali_pol2_band_pxct_FWHM()//fit 12 X-ray/Gamma ray peaks' FWHM values measured by LEGe detector with pol2 with a credible interval band
{
	const int ID1 = 0;// i=ID1//which detector
	const int ID2 = 0;// i<=ID2// modify which detector
	int i, ii;
	char paraprint[100], graph_name[200], hcali_name[200];
	TCanvas* canvascali[ID2 + 1];
	TGraph* graph[ID2 + 1];
	TGraph* graph_residual[ID2 + 1];//creat graphs
	TH1D* h_confidence_interval[ID2 + 1];
	const int num_datapoints = 11;//search peak numbers, modify
	Double_t x_value[ID2 + 1][num_datapoints], x_error[ID2 + 1][num_datapoints];
	Double_t y_value[ID2 + 1][num_datapoints], y_error[ID2 + 1][num_datapoints];
	Double_t residual[ID2 + 1][num_datapoints], residual_err[ID2 + 1][num_datapoints];
	double par[ID2 + 1][3], par_err[ID2 + 1][3];
	double parChi[ID2 + 1], parNDF[ID2 + 1], p_value[ID2 + 1];

	const int Npoint = 2;
	double energy_point[Npoint] = { 8.05, 8.64};// set location of point for single value. modify
	double err_point[Npoint];  // error on the function at point x0 for single value

	char pathname[150], filename[150];

	sprintf(pathname, "%s", "F:/e21010/pxct/");
	for (i = ID1; i <= ID2; i++)
	{
		sprintf(filename, "%s%s%d%s", pathname, "y_values_LEGe", i, ".dat");
		ifstream infile_y_values(filename, ios::in);
		sprintf(filename, "%s%s%d%s", pathname, "x_values_LEGe", i, ".dat");
		ifstream infile_x_values(filename, ios::in);
		for (ii = 0; ii < num_datapoints; ii++)
		{
			infile_x_values >> x_value[i][ii] >> x_error[i][ii];
			infile_y_values >> y_value[i][ii] >> y_error[i][ii];
			x_error[i][ii] = x_error[i][ii]*1; y_error[i][ii] = y_error[i][ii]*1;
			cout << "LEGe" << i << '	' << ii << '	' << x_value[i][ii] << '	' << x_error[i][ii] << '	' << y_value[i][ii] << '	' << y_error[i][ii] << endl;
		}
	}

	sprintf(filename, "%s%s", pathname, "calipara.dat");
	ofstream outfile(filename, ios::out);

	for (i = ID1; i <= ID2; i++) // no need to modify
	{
		sprintf(hcali_name, "%s%d", "FWHM_calibration_LEGe", i);
		canvascali[i] = new TCanvas(hcali_name, hcali_name, 1000, 700);
		canvascali[i]->cd();//进入画布
		TPad* pad1 = new TPad("pad1", "The pad 70% of the height", 0.0, 0.3, 1.0, 1.0);// Double_t xlow, Double_t ylow, Double_t xup, Double_t yup,
		TPad* pad2 = new TPad("pad2", "The pad 30% of the height", 0.0, 0.0, 1.0, 0.3);
		pad1->SetTopMargin(0.04);
		pad1->SetRightMargin(0.03);
		pad1->SetLeftMargin(0.08);
		pad1->SetBottomMargin(0.10);
		//pad1->SetBorderMode(0);
		pad2->SetTopMargin(0.01);
		pad2->SetRightMargin(0.03);
		pad2->SetLeftMargin(0.08);
		pad2->SetBottomMargin(0.23);
		//pad2->SetBorderMode(0);
		pad1->Draw();
		pad2->Draw();

		pad1->cd();
		gStyle->SetOptTitle(0);

		graph[i] = new TGraphErrors(num_datapoints, x_value[i], y_value[i], x_error[i], y_error[i]);//error bars TGraph(n,x,y,ex,ey);
		graph[i]->SetTitle(hcali_name);
		graph[i]->GetXaxis()->SetTitle("Energy (keV)");
		graph[i]->GetYaxis()->SetTitle("FWHM (keV)");
		graph[i]->GetXaxis()->CenterTitle();
		graph[i]->GetYaxis()->CenterTitle();
		graph[i]->GetXaxis()->SetLabelFont(132);
		graph[i]->GetYaxis()->SetLabelFont(132);
		graph[i]->GetXaxis()->SetTitleFont(132);
		graph[i]->GetYaxis()->SetTitleFont(132);
		//graph[i]->GetYaxis()->SetLabelSize(0.05);
		//graph[i]->GetYaxis()->SetTitleSize(0.05);
		graph[i]->GetXaxis()->SetTitleOffset(1.2);
		graph[i]->GetYaxis()->SetTitleOffset(1.0);
		graph[i]->GetXaxis()->SetRangeUser(0, 200);
		graph[i]->SetMarkerStyle(21);
		graph[i]->SetMarkerSize(0.5);
		graph[i]->SetMarkerColor(1);

		TF1* pol2 = new TF1("pol2", "pol2", 0, 200);
		pol2->SetNpx(40000);
		pol2->SetParNames("p0", "p1", "p2");//y=p0+p1*x+p2*x^2
		//graph[i]->Fit("pol1");//pol1 can be used directly without TF1 constructor in CINT
		TFitResultPtr Fit_result_pointer = graph[i]->Fit("pol2", "MS");
		//"S" means the result of the fit is returned in the TFitResultPtr
		//“E” Perform better errors estimation using the Minos technique.
		//“M” Improve fit results, by using the IMPROVE algorithm of TMinuit.
		//If the option "S" is instead used, TFitResultPtr contains the TFitResult and behaves as a smart pointer to it.
		//The fit parameters, error and chi2 (but not covariance matrix) can be retrieved also from the fitted function.

		// Uncertainty Band
		sprintf(filename, "%s%s%d", "h_confidence_interval", "_", i);
		h_confidence_interval[i] = new TH1D(filename, filename, 20000, 0, 200);//Create a histogram to hold the confidence intervals
		TVirtualFitter* fitter = TVirtualFitter::GetFitter();//The method TVirtualFitter::GetFitter())->Get the parameters of your fitting function after having it fitted to an histogram.
		fitter->GetConfidenceIntervals(h_confidence_interval[i], 0.95);//By default the intervals are inflated using the chi2/ndf value of the fit if a chi2 fit is performed
		//confidence interval for the colored band: 1σ confidence interval: P=0.683, 1σ confidence interval: P=0.95, 3σ confidence interval: P=0.997
		//h_confidence_interval[i] will contain the CL result that you can draw on top of your fitted graph.
		//where h_confidence_interval[i] will hold the errors and could superimpose it on the same canvas where you plot central values.
		h_confidence_interval[i]->SetStats(kFALSE);
		h_confidence_interval[i]->SetFillColor(kRed - 10);
		//Fit_result_pointer->GetConfidenceIntervals(Npoint, 1, 1, energy_point, err_point, 0.683, true);//get the error of y at x of interest, this command is independent of the colored band
		//(Number of x points, 1, 1, x, err, confidence level, false); norm is a flag to control if the intervals need to be inflated by the chi2/ndf value. true is inflated, false is not inflated.

		TMatrixD cov = Fit_result_pointer->GetCovarianceMatrix();//error matrix
		TMatrixD cor = Fit_result_pointer->GetCorrelationMatrix();//parameter correlation coefficients
		cov.Print();
		cor.Print();
		//Fit_result_pointer->Print("V");//print main fit information

		graph[i]->Draw("AP");//"A": Axis are drawn around the graph, "P": The current marker is plotted at each point
		h_confidence_interval[i]->Draw("e3 same");//draw the confidence band [colored] showing 1 or 2 or 3 standard deviation uncertainty
		graph[i]->Draw("P same");//draw the points again above error band

		pol2->GetParameters(par[i]);
		par_err[i][0] = pol2->GetParError(0);
		par_err[i][1] = pol2->GetParError(1);
		par_err[i][2] = pol2->GetParError(2);
		parChi[i] = pol2->GetChisquare();
		parNDF[i] = pol2->GetNDF();
		p_value[i] = pol2->GetProb();//This probability is not the “probability that your fit is good.” If you did many fake experiments (draw many random samples of data points from the assumed distribution (your fit function)), this is the percentage of experiments that would give χ2 values ≥ to the one you got in this experiment.

		// To extract FWHM values at energies of interest
		Fit_result_pointer->GetConfidenceIntervals(Npoint, 1, 1, energy_point, err_point, 0.683, true);//get the error of FWHM at energies of interest, this command is independent of the colored band
		//(Number of x points, 1, 1, x, err, confidence level, false); norm is a flag to control if the intervals need to be inflated by the chi2/ndf value. true is inflated, false is not inflated.
		//err_point contains one side of the error bar, so the full error bar length is 2*err_point
		for (int ipeak = 0; ipeak < Npoint; ipeak++)//get sigma/tau at energies of interest for output file
		{
			cout << "Eg=	" << energy_point[ipeak] << "	FWHM=	" << pol2->Eval(energy_point[ipeak]) << "	err_FWHM=	" << err_point[ipeak] << endl;
		}

		TPaveText* textpol1 = new TPaveText(0.69, 0.14, 0.96, 0.45, "brNDC");//left, down, right, up
		textpol1->SetBorderSize(1);
		textpol1->SetFillColor(0);
		textpol1->SetTextAlign(12);//align = 10*HorizontalAlign + VerticalAlign, 12 means水平左对齐、垂直居中对齐
		textpol1->SetTextFont(132);//font = 10 * fontID + precision, 12+2,12 means Symbol; 13+2, 13 means News Times Roman
		sprintf(paraprint, "y=p2*x^2 + p1*x + p0");
		textpol1->AddText(paraprint);
		sprintf(paraprint, "p2=%.9f+/-%.9f", par[i][2], par_err[i][2]);
		textpol1->AddText(paraprint);
		sprintf(paraprint, "p1=%.9f+/-%.9f", par[i][1], par_err[i][1]);
		textpol1->AddText(paraprint);
		sprintf(paraprint, "p0=%.5f+/-%.5f", par[i][0], par_err[i][0]);
		textpol1->AddText(paraprint);
		sprintf(paraprint, "Chi2=%.2f", parChi[i]);
		textpol1->AddText(paraprint);
		sprintf(paraprint, "NDF=%.0f", parNDF[i]);
		textpol1->AddText(paraprint);
		sprintf(paraprint, "p-val=%e", p_value[i]);
		textpol1->AddText(paraprint);
		textpol1->Draw();

		for (ii = 0; ii < num_datapoints; ii++)// for plotting the residuals
		{
			residual[i][ii] = (pol2->Eval(x_value[i][ii]) - y_value[i][ii]);
			residual_err[i][ii] = y_error[i][ii];
			outfile << "Ch=	" << x_value[i][ii] << "	residual=	" << residual[i][ii] << endl;
		}

		pad2->cd();
		graph_residual[i] = new TGraphErrors(num_datapoints, x_value[i], residual[i], x_error[i], residual_err[i]);//画error bars TGraph(n,x,y,ex,ey);
		graph_residual[i]->GetXaxis()->SetTitle("Energy (keV)");
		graph_residual[i]->GetYaxis()->SetTitle("Fit - Data (keV)");
		graph_residual[i]->GetXaxis()->CenterTitle();//居中
		graph_residual[i]->GetYaxis()->CenterTitle();//居中
		graph_residual[i]->GetXaxis()->SetLabelFont(132);//坐标字体
		graph_residual[i]->GetYaxis()->SetLabelFont(132);//坐标字体
		graph_residual[i]->GetXaxis()->SetTitleFont(132);//轴名字体
		graph_residual[i]->GetYaxis()->SetTitleFont(132);//轴名字体
		//graph_residual[i]->GetYaxis()->SetLabelSize(0.05);//坐标字号
		//graph_residual[i]->GetYaxis()->SetTitleSize(0.05);//轴名字号
		graph_residual[i]->GetXaxis()->SetRangeUser(0, 3700);
		graph_residual[i]->GetXaxis()->SetTitleOffset(1.3);
		graph_residual[i]->GetYaxis()->SetTitleOffset(0.4);
		graph_residual[i]->GetXaxis()->SetTitleSize(0.08);
		graph_residual[i]->GetXaxis()->SetLabelOffset(0.015);
		graph_residual[i]->GetXaxis()->SetLabelSize(0.08);
		graph_residual[i]->GetYaxis()->SetLabelSize(0.08);
		graph_residual[i]->GetYaxis()->SetTitleSize(0.08);
		graph_residual[i]->GetXaxis()->SetNdivisions(520);//n = n1 + 100*n2 + 10000*n3
		graph_residual[i]->GetXaxis()->SetNdivisions(10, 10, 1);
		graph_residual[i]->GetYaxis()->SetNdivisions(505);
		graph_residual[i]->SetLineWidth(2);
		graph_residual[i]->SetLineColor(kBlue + 2);
		graph_residual[i]->SetMarkerStyle(21);
		graph_residual[i]->SetMarkerSize(0.5);
		graph_residual[i]->SetMarkerColor(1);
		graph_residual[i]->Draw("APZ");//"A": Axis are drawn around the graph, "P": The current marker is plotted at each point, "Z": Do not draw small horizontal and vertical lines the end of the error bars. Without "Z", the default is to draw these.
		TLine* T1 = new TLine(0, 0, 65.5, 0); // TLine(x1, y1, x2, y2)
		T1->Draw("R");

		sprintf(filename, "%s%s%s", pathname, hcali_name, ".png");
		canvascali[i]->SaveAs(filename);
		outfile << "G" << i << "	y=p2*x^2+p1*x+p0	p2=	" << setprecision(8) << par[i][2] << "	+/-	" << par_err[i][2] << "	p1 = " << setprecision(8) << par[i][1] << " + / -" << par_err[i][1] <<"	p0 = " << par[i][0] << " + / -" << par_err[i][0] << "	Chi2 = " << parChi[i] << "	ndf = " << parNDF[i] << "	p - value = " << p_value[i] << endl;

	}//for (i=0;i<ID;i++)
}//peakcali main
