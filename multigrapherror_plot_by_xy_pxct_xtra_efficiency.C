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
#include "TMultiGraph.h"
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
#include <sstream>
#include <vector>
#include <string>
#include "stdlib.h"
#include "stdio.h"
#include "TLegend.h"
using namespace std;

void multigrapherror_plot_by_xy_pxct_xtra_efficiency() // reading x and y values from a file where the first column is x values and the subsequent columns are y values
// https://root.cern.ch/doc/master/classTMultiGraph.html
{
	char pathname[300], filename[300], detectorname[300], figurename[300];
	vector<string> detector_name = { "XtRa1", "XtRa2" };
	vector<string> header; // To store header names for the legend
	vector<double> x_values, x_errors, y_values, y_errors;
	/*
	The outer vector holds all the rows.Each element of this vector corresponds to one row of your data file.
		Each inner vector holds the y values for a particular x value.The size of each inner vector corresponds to the number of y values that are associated with an x value.
		Let's go through a simplified example to demonstrate how this works:

		Suppose your data file looks like this:
		1.0 2.1 3.2
		4.0 5.1 6.2 7.3
		8.0 9.1
		Here's how you might read this data into a vector<double> for x values and a vector<vector<double>> for y values:

		xValues vector : Holds the x values(1.0, 4.0, 8.0).
		yValues vector : Holds vectors of y values.For the first row, the inner vector holds(2.1, 3.2).For the second row, it holds(5.1, 6.2, 7.3), and for the third row, it holds just(9.1).
	*/
	// Create a TGraphErrors object and set its properties
	TCanvas* canvaspeak = new TCanvas("canvaspeak", "canvaspeak", 1200, 800);
	canvaspeak->cd();
	canvaspeak->SetTopMargin(0.029);
	canvaspeak->SetRightMargin(0.04);
	canvaspeak->SetLeftMargin(0.13);
	canvaspeak->SetBottomMargin(0.18);
	canvaspeak->SetFrameLineWidth(3);
	gStyle->SetFrameLineWidth(3);
	gStyle->SetEndErrorSize(5);
	
	TPad* pad1 = new TPad("pad1", "The pad 70% of the height", 0.0, 0.39, 1.0, 1.0);// Double_t xlow, Double_t ylow, Double_t xup, Double_t yup,
	TPad* pad2 = new TPad("pad2", "The pad 30% of the height", 0.0, 0.259, 1.0, 0.39);
	TPad* pad3 = new TPad("pad3", "The pad 30% of the height", 0.0, 0.0, 1.0, 0.259);
	pad1->SetTopMargin(0.04);
	pad1->SetRightMargin(0.03);
	pad1->SetLeftMargin(0.12);
	pad1->SetBottomMargin(0.03);
	pad1->SetFillColor(0);
	pad1->SetFillStyle(0);
	//pad1->SetBorderMode(0);
	pad2->SetTopMargin(0.01);
	pad2->SetRightMargin(0.03);
	pad2->SetLeftMargin(0.12);
	pad2->SetBottomMargin(0.00);
	pad2->SetFillColor(0);
	pad2->SetFillStyle(0);
	//pad2->SetBorderMode(0);	
	pad3->SetTopMargin(0.00);
	pad3->SetRightMargin(0.03);
	pad3->SetLeftMargin(0.12);
	pad3->SetBottomMargin(0.50);
	pad3->SetFillColor(0);
	pad3->SetFillStyle(0);
	//pad2->SetBorderMode(0);	
	pad1->Draw();
	pad2->Draw();
	pad3->Draw();

	pad1->cd();

	TMultiGraph* multigraph = new TMultiGraph();
	TGraph* graph[2];
	TGraph* graph_residual[2];
	TH1D* h_confidence_interval[2];
	int colors[] = { kRed + 2, kBlue + 2, kRed - 9, kAzure - 9, kBlack, kGreen + 1, kViolet + 1, kOrange + 7, kMagenta, kCyan + 1, kYellow + 1, kSpring + 7 };
	int styles[] = { 24, 20, 24, 25, 22, 23, 26, 27, 28, 29 };
	TF1* logpol6 = new TF1("logpol6", "exp([0]+[1]*log(x)+[2]*pow(log(x),2)+[3]*pow(log(x),3)+[4]*pow(log(x),4)+[5]*pow(log(x),5)+[6]*pow(log(x),6))", 0, 5000);
	logpol6->SetNpx(50000);
	logpol6->SetParNames("p0", "p1", "p2", "p3", "p4", "p5", "p6");//y=exp(p0+p1*lnx+p2*(lnx)^2+...
	logpol6->SetLineWidth(2);

	// Loop over all detectors
	for (size_t i_detector = 0; i_detector < detector_name.size(); i_detector++)
	{
		// Reset your vectors for each detector
		x_values.clear();
		x_errors.clear();
		y_values.clear();
		y_errors.clear();
		sprintf(pathname, "%s", "F:/e21010/pxct/");
		sprintf(filename, "%s%s%s%s", pathname, "x_y_values_for_pxct_", detector_name[i_detector].c_str(), "_efficiency.txt");
		ifstream infile_xy_values(filename, ios::in); // open the txt file for reading
		string line;

		// Read the header line
		if (getline(infile_xy_values, line))
		{
			stringstream ss(line);
			string name;
			while (ss >> name)
			{ // Read each name in the header
				header.push_back(name);
			}
		}
		while (getline(infile_xy_values, line))
		{
			stringstream ss(line);
			double x, x_err, y, y_err;
			// Read an x value and x error from the line
			if (ss >> x >> x_err >> y >> y_err)
			{
				x_values.push_back(x);
				x_errors.push_back(x_err);
				y_values.push_back(y);
				y_errors.push_back(y_err);
			}
		}

		// Print the data for verification
		cout << "Loaded data from " << filename << endl;
		cout << x_values.size() << " data points" << endl;

		for (size_t i_row = 0; i_row < x_values.size(); i_row++)
		{
			cout << "X: " << x_values[i_row] << " +/- " << x_errors[i_row];
			cout << ", Y: " << y_values[i_row] << " +/- " << y_errors[i_row] << endl;
		}

		graph[i_detector] = new TGraphErrors(x_values.size(), &x_values[0], &y_values[0], &x_errors[0], &y_errors[0]);
		graph[i_detector]->SetMarkerColor(colors[i_detector % (sizeof(colors) / sizeof(colors[0]))]);
		graph[i_detector]->SetMarkerStyle(styles[i_detector % (sizeof(styles) / sizeof(styles[0]))]);
		graph[i_detector]->SetMarkerSize(1.2);
		graph[i_detector]->SetLineColor(colors[i_detector % (sizeof(colors) / sizeof(colors[0]))]);
		graph[i_detector]->SetLineWidth(2);
		graph[i_detector]->SetTitle("");
		graph[i_detector]->SetName(detector_name[i_detector].c_str());
		//graph[i_detector]->GetXaxis()->SetLimits(0, 1536);
		multigraph->Add(graph[i_detector], "P"); // Use "P" to draw only markers, "L" to draw only lines, "PL" to draw both

		logpol6->SetLineColor(colors[(i_detector+2) % (sizeof(colors) / sizeof(colors[0]))]);
		logpol6->SetLineWidth(0); // Set the line width to 0 to hide the fit line because it is always above data points
		logpol6->SetParameter(0, -45);
		logpol6->SetParameter(1, 14);
		logpol6->SetParameter(2, -1.8);
		graph[i_detector]->Fit("logpol6", "MQ");//logpol6 can be used directly without TF1 constructor in CINT
		graph[i_detector]->Fit("logpol6", "MQ");//logpol6 can be used directly without TF1 constructor in CINT
		graph[i_detector]->Fit("logpol6", "MQ");//logpol6 can be used directly without TF1 constructor in CINT
		TFitResultPtr Fit_result_pointer = graph[i_detector]->Fit("logpol6", "MS");
		//"S" means the result of the fit is returned in the TFitResultPtr
		//“E” Perform better errors estimation using the Minos technique.
		//“M” Improve fit results, by using the IMPROVE algorithm of TMinuit.
		//If the option "S" is instead used, TFitResultPtr contains the TFitResult and behaves as a smart pointer to it.
		//The fit parameters, error and chi2 (but not covariance matrix) can be retrieved also from the fitted function.

		// Uncertainty Band
		sprintf(detectorname, "%s%s", "h_confidence_interval_", detector_name[i_detector].c_str());
		h_confidence_interval[i_detector] = new TH1D(detectorname, detectorname, 50000, 0, 5000);//Create a histogram to hold the confidence intervals
		TVirtualFitter* fitter = TVirtualFitter::GetFitter();//The method TVirtualFitter::GetFitter())->Get the parameters of your fitting function after having it fitted to an histogram.
		fitter->GetConfidenceIntervals(h_confidence_interval[i_detector], 0.683);//By default the intervals are inflated using the chi2/ndf value of the fit if a chi2 fit is performed
		//confidence interval for the colored band: 1σ confidence interval: P=0.683, 1σ confidence interval: P=0.95, 3σ confidence interval: P=0.997
		//h_confidence_interval[i] will contain the CL result that you can draw on top of your fitted graph.
		//where h_confidence_interval[i] will hold the errors and could superimpose it on the same canvas where you plot central values.
		h_confidence_interval[i_detector]->SetStats(kFALSE);
		h_confidence_interval[i_detector]->SetFillColor(colors[(i_detector+2) % (sizeof(colors) / sizeof(colors[0]))]);
		//Fit_result_pointer->GetConfidenceIntervals(Npoint, 1, 1, energy_point, err_point, 0.683, true);//get the error of y at x of interest, this command is independent of the colored band
		//(Number of x points, 1, 1, x, err, confidence level, false); norm is a flag to control if the intervals need to be inflated by the chi2/ndf value. true is inflated, false is not inflated.

		TMatrixD cov = Fit_result_pointer->GetCovarianceMatrix();//error matrix
		TMatrixD cor = Fit_result_pointer->GetCorrelationMatrix();//parameter correlation coefficients
		cov.Print();
		cor.Print();
		//Fit_result_pointer->Print("V");//print main fit information

		// Extract the fit parameters
		double p7[2], p6[2], p5[2], p4[2], p3[2], p2[2], p1[2], p0[2];
		double p7err[2], p6err[2], p5err[2], p4err[2], p3err[2], p2err[2], p1err[2], p0err[2];
		double parChi[2], parNDF[2], p_value[2];
		p0[i_detector] = logpol6->GetParameter(0);
		p1[i_detector] = logpol6->GetParameter(1);
		p2[i_detector] = logpol6->GetParameter(2);
		p3[i_detector] = logpol6->GetParameter(3);
		p4[i_detector] = logpol6->GetParameter(4);
		p5[i_detector] = logpol6->GetParameter(5);
		p6[i_detector] = logpol6->GetParameter(6);

		p0err[i_detector] = logpol6->GetParError(0);
		p1err[i_detector] = logpol6->GetParError(1);
		p2err[i_detector] = logpol6->GetParError(2);
		p3err[i_detector] = logpol6->GetParError(3);
		p4err[i_detector] = logpol6->GetParError(4);
		p5err[i_detector] = logpol6->GetParError(5);
		p6err[i_detector] = logpol6->GetParError(6);

		parChi[i_detector] = logpol6->GetChisquare();
		parNDF[i_detector] = logpol6->GetNDF();
		p_value[i_detector] = logpol6->GetProb();//This probability is not the “probability that your fit is good.” If you did many fake experiments (draw many random samples of data points from the assumed distribution (your fit function)), this is the percentage of experiments that would give χ2 values ≥ to the one you got in this experiment.

		// To extract FWHM values at energies of interest
		const int Npoint = 4;
		double energy_point[Npoint] = { 100, 500, 1000, 3000 };// set location of point for single value. modify
		double err_point[Npoint];  // error on the function at point x0 for single value
		Fit_result_pointer->GetConfidenceIntervals(Npoint, 1, 1, energy_point, err_point, 0.683, false);//get the error of FWHM at energies of interest, this command is independent of the colored band
		//(Number of x points, 1, 1, x, err, confidence level, false); norm is a flag to control if the intervals need to be inflated by the chi2/ndf value. true is inflated, false is not inflated.
		//err_point contains one side of the error bar, so the full error bar length is 2*err_point
		for (int ipeak = 0; ipeak < Npoint; ipeak++)//get sigma/tau at energies of interest for output file
		{
			cout << "X_value=	" << energy_point[ipeak] << "	Y_value=	" << logpol6->Eval(energy_point[ipeak]) << "	err_Y_value=	" << err_point[ipeak] << endl;
		}

		TPaveText* textpol6 = new TPaveText(0.66, 0.46, 0.969, 0.96, "brNDC");//left, down, right, up
		textpol6->SetBorderSize(1);
		textpol6->SetFillColor(0);
		textpol6->SetTextAlign(12);//align = 10*HorizontalAlign + VerticalAlign, 12 means水平左对齐、垂直居中对齐
		textpol6->SetTextFont(132);//font = 10 * fontID + precision, 12+2,12 means Symbol; 13+2, 13 means News Times Roman
		char paraprint[100];
		sprintf(paraprint, "y=exp[sum(p_i ln(E)^i)]");
		textpol6->AddText(paraprint);
		sprintf(paraprint, "p6=%e+/-%e", p6[i_detector], p6err[i_detector]);
		textpol6->AddText(paraprint);
		sprintf(paraprint, "p5=%e+/-%e", p5[i_detector], p5err[i_detector]);
		textpol6->AddText(paraprint);
		sprintf(paraprint, "p4=%e+/-%e", p4[i_detector], p4err[i_detector]);
		textpol6->AddText(paraprint);
		sprintf(paraprint, "p3=%e+/-%e", p3[i_detector], p3err[i_detector]);
		textpol6->AddText(paraprint);
		sprintf(paraprint, "p2=%e+/-%e", p2[i_detector], p2err[i_detector]);
		textpol6->AddText(paraprint);
		sprintf(paraprint, "p1=%e+/-%e", p1[i_detector], p1err[i_detector]);
		textpol6->AddText(paraprint);
		sprintf(paraprint, "p0=%e+/-%e", p0[i_detector], p0err[i_detector]);
		textpol6->AddText(paraprint);
		sprintf(paraprint, "Chi2=%.2f", parChi[i_detector]);
		textpol6->AddText(paraprint);
		sprintf(paraprint, "NDF=%.0f", parNDF[i_detector]);
		textpol6->AddText(paraprint);
		sprintf(paraprint, "p-val=%e", p_value[i_detector]);
		textpol6->AddText(paraprint);
		//textpol6->Draw();
		// print all the fit parameters
		cout << "p6=" << p6[i_detector] << " +/- " << p6err[i_detector] << endl;
		cout << "p5=" << p5[i_detector] << " +/- " << p5err[i_detector] << endl;
		cout << "p4=" << p4[i_detector] << " +/- " << p4err[i_detector] << endl;
		cout << "p3=" << p3[i_detector] << " +/- " << p3err[i_detector] << endl;
		cout << "p2=" << p2[i_detector] << " +/- " << p2err[i_detector] << endl;
		cout << "p1=" << p1[i_detector] << " +/- " << p1err[i_detector] << endl;
		cout << "p0=" << p0[i_detector] << " +/- " << p0err[i_detector] << endl;
		cout << "Chi2=" << parChi[i_detector] << endl;
		cout << "NDF=" << parNDF[i_detector] << endl;
		cout << "p-val=" << p_value[i_detector] << endl;


		// Residuals
		std::vector<double> residual(x_values.size());
		std::vector<double> residual_err(x_values.size());
		for (int i = 0; i < x_values.size(); i++)
		{
			residual[i] = y_values[i] - logpol6->Eval(x_values[i]);
			residual_err[i] = y_errors[i];
		}
		graph_residual[i_detector] = new TGraphErrors(x_values.size(), &x_values[0], &residual[0], &x_errors[0], &residual_err[0]);
		graph_residual[i_detector]->SetMarkerColor(colors[i_detector % (sizeof(colors) / sizeof(colors[0]))]);
		graph_residual[i_detector]->SetMarkerStyle(styles[i_detector % (sizeof(styles) / sizeof(styles[0]))]);
		graph_residual[i_detector]->SetMarkerSize(1.2);
		graph_residual[i_detector]->SetLineColor(colors[i_detector % (sizeof(colors) / sizeof(colors[0]))]);
		graph_residual[i_detector]->SetLineWidth(2);
		graph_residual[i_detector]->SetTitle("");
		graph_residual[i_detector]->SetName(detector_name[i_detector].c_str());
	}

	multigraph->Draw("A"); // "A" Axis are drawn around the graph
	logpol6->Draw("C same");
	// Overlay the confidence intervals
	for (size_t i = 0; i < detector_name.size(); i++)
	{
		h_confidence_interval[i]->Draw("E3 same");
	}
	multigraph->Draw("P"); // Draw the markers again to make them visible
	// Setting up axis titles
	multigraph->GetXaxis()->SetTitle("Energy (keV)");
	multigraph->GetYaxis()->SetTitle("Detection Efficiency (%)");
	multigraph->GetXaxis()->CenterTitle();
	multigraph->GetYaxis()->CenterTitle();
	multigraph->GetXaxis()->SetLabelFont(132);
	multigraph->GetYaxis()->SetLabelFont(132);
	multigraph->GetXaxis()->SetTitleFont(132);
	multigraph->GetYaxis()->SetTitleFont(132);
	multigraph->GetYaxis()->SetLabelSize(0.10);
	multigraph->GetYaxis()->SetTitleSize(0.10);
	multigraph->GetXaxis()->SetLabelSize(0.10);
	multigraph->GetXaxis()->SetTitleSize(0.10);
	multigraph->GetYaxis()->SetTickLength(0.02);
	multigraph->GetYaxis()->SetNdivisions(505);
	multigraph->GetXaxis()->SetTitleOffset(1.05);
	multigraph->GetYaxis()->SetTitleOffset(0.613);
	multigraph->GetXaxis()->SetLabelOffset(0.09);
	multigraph->GetXaxis()->SetLimits(0, 1536);
	multigraph->GetXaxis()->SetRangeUser(0, 1536);
	multigraph->GetYaxis()->SetRangeUser(0, 0.85);

	// Adding legend
	TLegend* legend1 = new TLegend(0.79, 0.82, 0.96, 0.95);//left, down, right, up
	legend1->SetBorderSize(0);
	legend1->SetTextFont(132);
	legend1->SetTextSize(0.10);
	legend1->SetTextAlign(12);
	legend1->SetFillStyle(0);
	legend1->AddEntry(multigraph->GetListOfGraphs()->At(0), detector_name[0].c_str(), "EP"); // "EP" draws the marker and the error bar
	legend1->SetTextColor(colors[0]);
	legend1->Draw();
	TLegend* legend2 = new TLegend(0.79, 0.69, 0.96, 0.82);//left, down, right, up
	legend2->SetBorderSize(0);
	legend2->SetTextFont(132);
	legend2->SetTextSize(0.10);
	legend2->SetTextAlign(12);
	legend2->SetFillStyle(0);
	legend2->AddEntry(multigraph->GetListOfGraphs()->At(1), detector_name[1].c_str(), "EP");
	legend2->SetTextColor(colors[1]);
	legend2->Draw();
	pad1->RedrawAxis();

	pad2->cd();
	graph_residual[0]->Draw("AP");
	graph_residual[0]->GetYaxis()->SetTitle("Fit (%)");
	graph_residual[0]->GetYaxis()->CenterTitle();
	graph_residual[0]->GetXaxis()->SetLabelFont(132);
	graph_residual[0]->GetYaxis()->SetLabelFont(132);
	graph_residual[0]->GetXaxis()->SetTitleFont(132);
	graph_residual[0]->GetYaxis()->SetTitleFont(132);
	graph_residual[0]->GetYaxis()->SetLabelSize(0.41);
	graph_residual[0]->GetYaxis()->SetTitleSize(0.385);
	graph_residual[0]->GetXaxis()->SetLabelSize(0.07);
	graph_residual[0]->GetXaxis()->SetTitleSize(0.08);
	graph_residual[0]->GetXaxis()->SetTickLength(0.14);
	graph_residual[0]->GetYaxis()->SetTickLength(0.010);
	graph_residual[0]->GetXaxis()->SetTitleOffset(1.05);
	graph_residual[0]->GetYaxis()->SetTitleOffset(0.162);
	graph_residual[0]->GetYaxis()->SetNdivisions(105);
	graph_residual[0]->GetXaxis()->SetLimits(0, 1536);
	graph_residual[0]->GetXaxis()->SetRangeUser(0, 1536);
	graph_residual[0]->GetYaxis()->SetRangeUser(-0.035, 0.035);
	TLine* T1 = new TLine(0, 0, 1536, 0); // TLine(x1, y1, x2, y2)
	T1->Draw("R");

	pad3->cd();
	graph_residual[1]->Draw("AP");
	graph_residual[1]->GetXaxis()->SetTitle("Energy (keV)");
	graph_residual[1]->GetYaxis()->SetTitle("Data #minus ");
	graph_residual[1]->GetXaxis()->CenterTitle();
	graph_residual[1]->GetXaxis()->SetLabelFont(132);
	graph_residual[1]->GetYaxis()->SetLabelFont(132);
	graph_residual[1]->GetXaxis()->SetTitleFont(132);
	graph_residual[1]->GetYaxis()->SetTitleFont(132);
	graph_residual[1]->GetYaxis()->SetLabelSize(0.21);
	graph_residual[1]->GetYaxis()->SetTitleSize(0.22);
	graph_residual[1]->GetXaxis()->SetLabelSize(0.24);
	graph_residual[1]->GetXaxis()->SetTitleSize(0.24);
	graph_residual[1]->GetYaxis()->SetTickLength(0.02);
	graph_residual[1]->GetXaxis()->SetTickLength(0.07);
	graph_residual[1]->GetXaxis()->SetTitleOffset(1.00);
	graph_residual[1]->GetYaxis()->SetTitleOffset(0.29);
	graph_residual[1]->GetYaxis()->SetNdivisions(105);
	graph_residual[1]->GetXaxis()->SetLimits(0, 1536);
	graph_residual[1]->GetXaxis()->SetRangeUser(0, 1536);
	graph_residual[1]->GetYaxis()->SetRangeUser(-0.035, 0.035);
	T1->Draw("R");


	// Adding text with the same colors as the curves
// 	TPaveText* colorText = new TPaveText(0.23, 0.21, 0.52, 0.57, "NDC");//left, down, right, up
// 	colorText->SetBorderSize(0);
// 	colorText->SetFillStyle(0); // Transparent
// 	colorText->SetTextFont(132);
// 	colorText->SetTextSize(0.07);
// 	colorText->SetTextAlign(12); //align = 10*HorizontalAlign + VerticalAlign, 12 means aligns the text to the left (horizontal alignment) and centered vertically relative to the specified coordinates
// 	// Add text entries with corresponding colors
// 	const char* seriesNames[] = { "XtRa1", "XtRa2" }; // to overlay a text box on top of or beside your legend to display the colors corresponding to each series
// 	for (size_t i = 0; i < sizeof(seriesNames) / sizeof(seriesNames[0]); i++)
// 	{
// 		auto textEntry = colorText->AddText(seriesNames[i]);
// 		textEntry->SetTextColor(colors[i]);
// 	}
// 	colorText->Draw();

	canvaspeak->cd();
	canvaspeak->Update();
	sprintf(filename, "%s%s", pathname, "Fig_PXCT_Gamma_Efficiency_XtRa");
	sprintf(figurename, "%s%s", filename, ".png");
	canvaspeak->SaveAs(figurename);
	sprintf(figurename, "%s%s", filename, ".eps");
	canvaspeak->SaveAs(figurename);

}
