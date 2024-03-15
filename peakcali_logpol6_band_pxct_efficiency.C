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
#include "TGraphPainter.h"
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
void peakcali_logpol6_band_pxct_efficiency()//for efficiency curve fit with a credible interval band
{
	const int ID1 = 0;// i=ID1//which detector
	const int ID2 = 1;// i<=ID2// modify which detector
	int i, ii;
	char paraprint[100], graph_name[200], hcali_name[200];
	TCanvas* canvascali[ID2 + 1];
	TGraph* graph[ID2 + 1];
	TGraph* graph_residual[ID2 + 1];//creat graphs
	TGraph* graph_residual_stat[ID2 + 1];//creat graphs
	TH1D* h_confidence_interval[ID2 + 1];
	TPad* pad1;// Double_t xlow, Double_t ylow, Double_t xup, Double_t yup,
	TPad* pad2;
	int num_datapoints = 22;
	Double_t x_value[ID2 + 1][num_datapoints], x_error[ID2 + 1][num_datapoints];
	Double_t y_value[ID2 + 1][num_datapoints], y_error[ID2 + 1][num_datapoints], y_stat_err[ID2 + 1][num_datapoints];
	Double_t residual[ID2 + 1][num_datapoints], residual_err[ID2 + 1][num_datapoints], residual_stat_err[ID2 + 1][num_datapoints];
	double p7[300], p6[300], p5[300], p4[300], p3[300], p2[300], p1[300], p0[300];
	double p7err[300], p6err[300], p5err[300], p4err[300], p3err[300], p2err[300], p1err[300], p0err[300];
	double parChi[ID2 + 1], parNDF[ID2 + 1], p_value[ID2 + 1];

	const int Npoint = 3;
	double energy_point[Npoint] = { 100, 500, 1000 };// set location of point for single value. modify
	double err_point[Npoint];  // error on the function at point x0 for single value

	char pathname[150], filename[150];

	sprintf(pathname, "%s", "F:/e21010/pxct/");
	for (i = ID1; i <= ID2; i++)
	{
		if (i == 0) num_datapoints = 21;
		if (i == 1) num_datapoints = 22;
		sprintf(filename, "%s%s%d%s", pathname, "y_values_XtRa_", i+1, ".dat");
		ifstream infile_y_values(filename, ios::in);
		sprintf(filename, "%s%s%d%s", pathname, "x_values_XtRa_", i+1, ".dat");
		ifstream infile_x_values(filename, ios::in);
		for (ii = 0; ii < num_datapoints; ii++)
		{
			infile_x_values >> x_value[i][ii] >> x_error[i][ii];
			infile_y_values >> y_value[i][ii] >> y_error[i][ii] >> y_stat_err[i][ii];
			y_value[i][ii] = y_value[i][ii] * 100;//convert to percentage
			y_error[i][ii] = y_error[i][ii] * 100;//convert to percentage
			//x_error[i][ii] = x_error[i][ii]*1; y_error[i][ii] = y_error[i][ii]*1;
			cout << "XtRa_" << i+1 << '	' << ii << '	' << x_value[i][ii] << '	' << x_error[i][ii] << '	' << y_value[i][ii] << '	' << y_error[i][ii] << '	' << y_stat_err[i][ii] << endl;
		}
	}
	string line;
	stringstream ss;
	int irow = 0, icolumn = 0;

// 	while (getline(infile_x_values, line)) // works perfectly for txt files tab delimited. Read in row-by-row
// 	{
// 		//if (irow++ == 0) continue; // skip the row headings
// 		ss.clear(); //clear(): Used to clear the stream
// 		ss.str(line); //str(): To get and set the string object whose content is present in stream
// 		icolumn = 0; // column number 0-4 are parameters; 5-n are parameters/bincounts
// 		while (!ss.fail()) // not the end of a line
// 		{
// 			ss >> x_value[i][ii] >> x_error[i][ii]; //operator >> : This is used to read from stringstream object.
// 		}

	sprintf(filename, "%s%s", pathname, "calipara.dat");
	ofstream outfile(filename, ios::out);

	for (i = ID1; i <= ID2; i++) // no need to modify
	{
		sprintf(hcali_name, "%s%d", "Efficiency_XtRa", i + 1);
		canvascali[i] = new TCanvas(hcali_name, hcali_name, 1000, 700);
		canvascali[i]->cd();//进入画布
		gStyle->SetOptTitle(0);
		gStyle->SetFrameLineWidth(2);
		pad1 = new TPad("pad1", "The pad 70% of the height", 0.0, 0.35, 1.0, 1.0);// Double_t xlow, Double_t ylow, Double_t xup, Double_t yup,
		pad2 = new TPad("pad2", "The pad 30% of the height", 0.0, 0.0, 1.0, 0.35);
		pad1->SetTopMargin(0.04);
		pad1->SetRightMargin(0.03);
		pad1->SetLeftMargin(0.12);
		pad1->SetBottomMargin(0.07);
		//pad1->SetBorderMode(0);
		pad2->SetTopMargin(0.08);
		pad2->SetRightMargin(0.03);
		pad2->SetLeftMargin(0.12);
		pad2->SetBottomMargin(0.41);
		//pad2->SetBorderMode(0);	
		pad1->Draw();
		gStyle->SetFrameLineWidth(2);
		pad2->Draw();
		gStyle->SetFrameLineWidth(2);
		pad1->cd();

		if (i == 0) num_datapoints = 21;
		if (i == 1) num_datapoints = 22;
		//graph[i] = new TGraphErrors(num_datapoints, x_value[i], y_value[i], x_error[i], y_error[i]);//error bars TGraph(n,x,y,ex,ey);
		graph[i] = new TGraphErrors(num_datapoints, x_value[i], y_value[i], x_error[i], y_stat_err[i]);//error bars TGraph(n,x,y,ex,ey);
		sprintf(graph_name, "%s%d", "gEfficiency_XtRa", i + 1);
		graph[i]->SetTitle(graph_name);
		graph[i]->SetName(graph_name);
		graph[i]->GetXaxis()->SetTitle("Energy (keV)");
		graph[i]->GetYaxis()->SetTitle("Detection Efficiency (%)");
		graph[i]->GetXaxis()->CenterTitle();
		graph[i]->GetYaxis()->CenterTitle();
		graph[i]->GetXaxis()->SetLabelFont(132);
		graph[i]->GetYaxis()->SetLabelFont(132);
		graph[i]->GetXaxis()->SetTitleFont(132);
		graph[i]->GetYaxis()->SetTitleFont(132);
		graph[i]->GetYaxis()->SetLabelSize(0.07);
		graph[i]->GetYaxis()->SetTitleSize(0.07);
		graph[i]->GetXaxis()->SetLabelSize(0.07);
		graph[i]->GetXaxis()->SetTitleSize(0.07);
		graph[i]->GetYaxis()->SetTickLength(0.02);
		graph[i]->GetYaxis()->SetNdivisions(505);
		graph[i]->GetXaxis()->SetTitleOffset(1.2);
		graph[i]->GetYaxis()->SetTitleOffset(0.8);
		graph[i]->GetXaxis()->SetRangeUser(0, 1536);
		graph[i]->GetYaxis()->SetRangeUser(0, 0.85);
		graph[i]->SetMarkerStyle(21);
		graph[i]->SetMarkerColor(1);
		graph[i]->SetLineWidth(2);
		gStyle->SetEndErrorSize(5);
		TF1* logpol6 = new TF1("logpol6", "exp([0]+[1]*log(x)+[2]*pow(log(x),2)+[3]*pow(log(x),3)+[4]*pow(log(x),4)+[5]*pow(log(x),5)+[6]*pow(log(x),6))", 0, 10000);
		logpol6->SetNpx(50000);
		logpol6->SetParNames("p0", "p1", "p2", "p3", "p4", "p5", "p6");//y=exp(p0+p1*lnx+p2*(lnx)^2+...
		logpol6->SetParameter(0, -45);
		//logpol6->SetParError(0,0.005598814);
		//logpol6->SetParLimits(0,1.777853,1.777853);
		logpol6->SetParameter(1, 14);
		//logpol6->SetParError(1,0.0009436133);
		//logpol6->SetParLimits(1,-1.175923,-1.175923);
		logpol6->SetParameter(2, -1.8);
		//logpol6->SetParError(2,0.0001292948);
		//logpol6->SetParLimits(2,1.106058,1.106058);
		//logpol6->SetParameter(3, -0.60187579);
		//logpol6->SetParError(3,1.688012e-05);
		//logpol6->SetParLimits(3,-0.6004948,-0.6004948);
		//logpol6->SetParameter(4, 0.12472777);
		//logpol6->SetParError(4,2.128892e-06);
		//logpol6->SetParLimits(4,0.1250527,0.1250527);
		//logpol6->SetParameter(5, -0.011253297);
		//logpol6->SetParError(5,2.523441e-07);
		//logpol6->SetParLimits(5,-0.0112302,-0.0112302);
		//logpol6->SetParameter(6, 0.00037136058);
		//logpol6->SetParError(6,2.741992e-08);
		//logpol6->SetParLimits(6,0.0003661849,0.0003661849);
		logpol6->SetLineColor((i+1)*2);
		logpol6->SetLineWidth(2);
		graph[i]->Fit("logpol6","MQ");//logpol6 can be used directly without TF1 constructor in CINT
		graph[i]->Fit("logpol6","Q");//logpol6 can be used directly without TF1 constructor in CINT
		graph[i]->Fit("logpol6","Q");//logpol6 can be used directly without TF1 constructor in CINT
		
		TFitResultPtr Fit_result_pointer = graph[i]->Fit("logpol6", "MS");
		//"S" means the result of the fit is returned in the TFitResultPtr
		//“E” Perform better errors estimation using the Minos technique.
		//“M” Improve fit results, by using the IMPROVE algorithm of TMinuit.
		//If the option "S" is instead used, TFitResultPtr contains the TFitResult and behaves as a smart pointer to it.
		//The fit parameters, error and chi2 (but not covariance matrix) can be retrieved also from the fitted function.

		// Uncertainty Band
		sprintf(filename, "%s%s%d", "h_confidence_interval", "_", i);
		h_confidence_interval[i] = new TH1D(filename, filename, 16000, 0, 1600);//Create a histogram to hold the confidence intervals
		TVirtualFitter* fitter = TVirtualFitter::GetFitter();//The method TVirtualFitter::GetFitter())->Get the parameters of your fitting function after having it fitted to an histogram.
		fitter->GetConfidenceIntervals(h_confidence_interval[i], 0.683);//By default the intervals are inflated using the chi2/ndf value of the fit if a chi2 fit is performed
		//confidence interval for the colored band: 1σ confidence interval: P=0.683, 1σ confidence interval: P=0.95, 3σ confidence interval: P=0.997
		//h_confidence_interval[i] will contain the CL result that you can draw on top of your fitted graph.
		//where h_confidence_interval[i] will hold the errors and could superimpose it on the same canvas where you plot central values.
		h_confidence_interval[i]->SetStats(kFALSE);
		if (i == 0) h_confidence_interval[i]->SetFillColor(kRed - 9);
		if (i == 1) h_confidence_interval[i]->SetFillColor(kAzure - 9);
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
		pad1->RedrawAxis();

		p0[i] = logpol6->GetParameter(0);
		p1[i] = logpol6->GetParameter(1);
		p2[i] = logpol6->GetParameter(2);
		p3[i] = logpol6->GetParameter(3);
		p4[i] = logpol6->GetParameter(4);
		p5[i] = logpol6->GetParameter(5);
		p6[i] = logpol6->GetParameter(6);

		p0err[i] = logpol6->GetParError(0);
		p1err[i] = logpol6->GetParError(1);
		p2err[i] = logpol6->GetParError(2);
		p3err[i] = logpol6->GetParError(3);
		p4err[i] = logpol6->GetParError(4);
		p5err[i] = logpol6->GetParError(5);
		p6err[i] = logpol6->GetParError(6);

		parChi[i] = logpol6->GetChisquare();
		parNDF[i] = logpol6->GetNDF();
		p_value[i] = logpol6->GetProb();//This probability is not the “probability that your fit is good.” If you did many fake experiments (draw many random samples of data points from the assumed distribution (your fit function)), this is the percentage of experiments that would give χ2 values ≥ to the one you got in this experiment.

		// To extract FWHM values at energies of interest
		Fit_result_pointer->GetConfidenceIntervals(Npoint, 1, 1, energy_point, err_point, 0.683, true);//get the error of FWHM at energies of interest, this command is independent of the colored band
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
		sprintf(paraprint, "y=exp[sum(p_i ln(E)^i)]");
		textpol6->AddText(paraprint);
		sprintf(paraprint, "p6=%e+/-%e", p6[i], p6err[i]);
		textpol6->AddText(paraprint);
		sprintf(paraprint, "p5=%e+/-%e", p5[i], p5err[i]);
		textpol6->AddText(paraprint);
		sprintf(paraprint, "p4=%e+/-%e", p4[i], p4err[i]);
		textpol6->AddText(paraprint);
		sprintf(paraprint, "p3=%e+/-%e", p3[i], p3err[i]);
		textpol6->AddText(paraprint);
		sprintf(paraprint, "p2=%e+/-%e", p2[i], p2err[i]);
		textpol6->AddText(paraprint);
		sprintf(paraprint, "p1=%e+/-%e", p1[i], p1err[i]);
		textpol6->AddText(paraprint);
		sprintf(paraprint, "p0=%e+/-%e", p0[i], p0err[i]);
		textpol6->AddText(paraprint);
		sprintf(paraprint, "Chi2=%.2f", parChi[i]);
		textpol6->AddText(paraprint);
		sprintf(paraprint, "NDF=%.0f", parNDF[i]);
		textpol6->AddText(paraprint);
		sprintf(paraprint, "p-val=%e", p_value[i]);
		textpol6->AddText(paraprint);
		//textpol6->Draw();
		
		for (ii = 0; ii < num_datapoints; ii++)// for plotting the residuals
		{
			residual[i][ii] = (logpol6->Eval(x_value[i][ii]) - y_value[i][ii]);
			residual_err[i][ii] = y_error[i][ii];
			residual_stat_err[i][ii] = y_stat_err[i][ii];
			outfile << "Energy=	" << x_value[i][ii] << "	residual=	" << residual[i][ii] << endl;
		}

		pad2->cd();

		//auto graph_residual = new TGraphMultiErrors(num_datapoints, x_value[i], residual[i], x_error[i], x_error[i], residual_err[i], residual_err[i]);//画error bars TGraph(n,x,y,exl,exh,eyl,eyh);
		//graph_residual->AddYError(num_datapoints, residual_stat_err[i], residual_stat_err[i]);
		//graph_residual[i] = new TGraphErrors(num_datapoints, x_value[i], residual[i], x_error[i], residual_err[i]);//画error bars TGraph(n,x,y,ex,ey);
		graph_residual[i] = new TGraphErrors(num_datapoints, x_value[i], residual[i], x_error[i], residual_stat_err[i]);//画error bars TGraph(n,x,y,ex,ey);
		graph_residual[i]->GetXaxis()->SetTitle("Energy (keV)");
		graph_residual[i]->GetYaxis()->SetTitle("Fit - Data (%)");
		graph_residual[i]->GetXaxis()->CenterTitle();//居中
		graph_residual[i]->GetYaxis()->CenterTitle();//居中
		graph_residual[i]->GetXaxis()->SetLabelFont(132);//坐标字体
		graph_residual[i]->GetYaxis()->SetLabelFont(132);//坐标字体
		graph_residual[i]->GetXaxis()->SetTitleFont(132);//轴名字体
		graph_residual[i]->GetYaxis()->SetTitleFont(132);//轴名字体
		graph_residual[i]->GetXaxis()->SetRangeUser(0, 1536);
		graph_residual[i]->GetYaxis()->SetRangeUser(-0.04, 0.04);
		graph_residual[i]->GetXaxis()->SetTitleSize(0.15);
		graph_residual[i]->GetYaxis()->SetTitleSize(0.12); 
		graph_residual[i]->GetXaxis()->SetTitleOffset(1.2);
		graph_residual[i]->GetYaxis()->SetTitleOffset(0.48); 
		graph_residual[i]->GetXaxis()->SetLabelSize(0.135);
		graph_residual[i]->GetYaxis()->SetLabelSize(0.13); 
		graph_residual[i]->GetXaxis()->SetLabelOffset(0.015);
// 		graph_residual[i]->GetXaxis()->SetNdivisions(520);//n = n1 + 100*n2 + 10000*n3
// 		graph_residual[i]->GetXaxis()->SetNdivisions(10, 10, 1);
		graph_residual[i]->GetYaxis()->SetNdivisions(105);
		graph_residual[i]->GetYaxis()->SetTickLength(0.02);
		graph_residual[i]->SetMarkerStyle(21);
		graph_residual[i]->SetMarkerColor(1);
		graph_residual[i]->SetLineWidth(2);
		graph_residual[i]->Draw("AP");//"A": Axis are drawn around the graph, "P": The current marker is plotted at each point
		TLine* T1 = new TLine(0, 0, 1536, 0); // TLine(x1, y1, x2, y2)
		T1->Draw("R");

		graph_residual_stat[i] = new TGraphErrors(num_datapoints, x_value[i], residual[i], x_error[i], residual_stat_err[i]);//画error bars TGraph(n,x,y,ex,ey);
		graph_residual_stat[i]->SetMarkerStyle(1);
		graph_residual_stat[i]->SetMarkerColor(3);
		graph_residual_stat[i]->SetLineWidth(2);
		graph_residual_stat[i]->SetLineColor(3);
		//graph_residual_stat[i]->Draw("P same");//"P": The current marker is plotted at each point

		sprintf(filename, "%s%s%s", pathname, hcali_name, ".png");
		canvascali[i]->SaveAs(filename);
		//outfile << "XtRa_" << i+1 << "	y=p1*x+p0	p1=	" << setprecision(8) << par[i][1] << "	+/-	" << par_err[i][1] << "	p0=	" << par[i][0] << "	+/-	" << par_err[i][0] << "	Chi2=	" << parChi[i] << "	ndf=	" << parNDF[i] << "	p-value=	" << p_value[i] << endl;

	}//for (i=0;i<ID;i++)
}//peakcali main
