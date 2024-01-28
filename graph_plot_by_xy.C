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
void graph_plot_by_xy() // read in x y values from a txt file and plot a graph
{
	const int row_start = 0;
	const int row_end = 2320;
	const int rows_one_group = 40;
	int i, ii, igroup;
	int Ea_gate_all[row_end] = { 0 }, Ea_gate[rows_one_group] = { 0 };
	double fitrange_min_all[row_end] = { 0 }, Half_life_all[row_end] = { 0 }, Half_life_err_all[row_end] = { 0 };
	double fitrange_min[rows_one_group] = { 0 }, Half_life[rows_one_group] = { 0 }, Half_life_err[rows_one_group] = { 0 }, fitrange_min_err[rows_one_group] = { 0 };
	char graph_name[200], pathname[150], filename[150];

	sprintf(pathname, "%s", "F:/e21010/pxct/");
	sprintf(filename, "%s%s", pathname, "x_T_y_fit_start_point_x_y_values.txt");
	ifstream infile_xy_values(filename, ios::in); // open the txt file for reading
	for (i = row_start; i < row_end; i++)
	{
		infile_xy_values >> Ea_gate_all[i] >> fitrange_min_all[i] >> Half_life_all[i] >> Half_life_err_all[i];  // read in x y values by row
		//cout << "fitrange_min_all[" << i << "]=" << fitrange_min_all[i] << "	Half_life_all[" << i << "]=" << Half_life_all[i] << "	Half_life_err_all[" << i << "]=" << Half_life_err_all[i] << endl;
	}
	igroup = 57 * rows_one_group; // modify this number to plot different groups of data
	for (ii = 0; ii < rows_one_group; ii++)
	{
		Ea_gate[ii] = Ea_gate_all[ii + igroup];
		fitrange_min[ii] = fitrange_min_all[ii + igroup];
		Half_life[ii] = Half_life_all[ii + igroup];
		Half_life_err[ii] = Half_life_err_all[ii + igroup];
		cout << "Ea_gate[" << ii << "]=" << Ea_gate[ii] << "	fitrange_min[" << ii << "]=" << fitrange_min[ii] << "	Half_life[" << ii << "]=" << Half_life[ii] << "	Half_life_err[" << ii << "]=" << Half_life_err[ii] << endl;
	}
	
	TCanvas* canvaspeak = new TCanvas("canvaspeak", "canvaspeak", 1300, 505);
	canvaspeak->cd();
	canvaspeak->SetTopMargin(0.07);
	canvaspeak->SetRightMargin(0.03);
	canvaspeak->SetLeftMargin(0.08);
	canvaspeak->SetBottomMargin(0.17);
	canvaspeak->SetFrameLineWidth(2);

	TGraph* graph = new TGraphErrors(rows_one_group, fitrange_min, Half_life, fitrange_min_err, Half_life_err);//error bars TGraph(n,x,y,ex,ey);
	sprintf(graph_name, "%s%d%s", "Ea_gate_", Ea_gate[0], "_keV");
	graph->SetTitle(graph_name);
	graph->SetName(graph_name);
	//graph->GetHistogram()->SetTitleFont(132, "t"); // "t" is for title
	graph->GetXaxis()->SetTitle("Fit start time (ns)");
	graph->GetYaxis()->SetTitle("Half-life (ns)");
	graph->GetXaxis()->CenterTitle();
	graph->GetYaxis()->CenterTitle();
	graph->GetXaxis()->SetLabelFont(132);
	graph->GetYaxis()->SetLabelFont(132);
	graph->GetXaxis()->SetTitleFont(132);
	graph->GetYaxis()->SetTitleFont(132);
	graph->GetYaxis()->SetLabelSize(0.06);
	graph->GetYaxis()->SetTitleSize(0.07);
	graph->GetXaxis()->SetLabelSize(0.06);
	graph->GetXaxis()->SetTitleSize(0.07);
	graph->GetYaxis()->SetTickLength(0.015);
	graph->GetYaxis()->SetNdivisions(505);
	graph->GetXaxis()->SetTitleOffset(1.2);
	graph->GetYaxis()->SetTitleOffset(0.5);
	graph->GetXaxis()->SetRangeUser(140, 550);
	//graph->GetYaxis()->SetRangeUser(0, 0.85);
	graph->SetMarkerStyle(21);
	graph->SetMarkerColor(1);
	graph->SetLineWidth(2);
	gStyle->SetEndErrorSize(5);

	graph->Draw("AP");//"A": Axis are drawn around the graph, "P": The current marker is plotted at each point

	sprintf(filename, "%s%s%s", pathname, graph_name, ".png");
	canvaspeak->SaveAs(filename);
		
}//peakcali main
