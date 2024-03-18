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
void multigraph_plot_by_xy_ddas_dead_time() // reading x and y values from a file where the first column is x values and the subsequent columns are y values
// https://root.cern.ch/doc/master/classTMultiGraph.html
{
	char pathname[300], filename[300], figurename[300];
	vector<string> header; // To store header names for the legend
	vector<double> x_values; 
	vector<vector<double>> y_values; // A vector of vectors. Each inner vector can represent a series of y values corresponding to a specific x value. 
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

	sprintf(pathname, "%s", "F:/e21010/pxct/");
	sprintf(filename, "%s%s", pathname, "x_y_values_for_ddas_deadtime.txt");
	ifstream infile_xy_values(filename, ios::in); // open the txt file for reading
	string line;

	// Reading data from file
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
		double x, y;
		vector<double> y_a_line;
		ss >> x; // Read an x value from ss
		x_values.push_back(x); //add a new element to the end of a vector,

		while (ss >> y) { // Read a row of y values
			y_a_line.push_back(y);
		}
		y_values.push_back(y_a_line);
	}

	// Print the data for verification
	cout << "Loaded data from " << filename << endl;
	cout << x_values.size() << " x values and each x has " << y_values[0].size() << " y values" << endl;
	// y_values[0].size() means the number of y values in the first row.

	for (size_t i_row = 0; i_row < x_values.size(); i_row++)
	{
		cout << "X: " << x_values[i_row] << " Y:";
		for (size_t i_column = 0; i_column < y_values[i_row].size(); i_column++)
		{
			cout << " " << y_values[i_row][i_column];
		}
		cout << endl;
	}
	
	TCanvas* canvaspeak = new TCanvas("canvaspeak", "canvaspeak", 1300, 800);
	canvaspeak->cd();
	canvaspeak->SetTopMargin(0.029);
	canvaspeak->SetRightMargin(0.04);
	canvaspeak->SetLeftMargin(0.13);
	canvaspeak->SetBottomMargin(0.18);
	canvaspeak->SetFrameLineWidth(3);

	TMultiGraph* multigraph = new TMultiGraph();
	vector<TGraph*> graphs;
	int colors[] = { kBlack, kRed - 4, kAzure - 3, kGreen + 1, kViolet + 1, kOrange + 7, kMagenta, kCyan + 1, kYellow + 1, kSpring + 7 };
	int styles[] = { 20, 21, 24, 25, 22, 23, 26, 27, 28, 29 };

	// Assuming all rows have the same number of y values
	// Number of series is determined by the number of y-values in the first row
	// Creating a graph for each series of y-values
	for (size_t i_column = 0; i_column < y_values[0].size(); i_column++)
	{
		graphs.push_back(new TGraph(x_values.size())); // Initialize each graph with the number of x-values

		// Setting points for this column (series)
		for (size_t i_row = 0; i_row < x_values.size(); i_row++)
		{
			graphs[i_column]->SetPoint(i_row, x_values[i_row], y_values[i_row][i_column]);
		}
		graphs[i_column]->SetMarkerColor(colors[i_column % (sizeof(colors) / sizeof(colors[0]))]);
		graphs[i_column]->SetMarkerStyle(styles[i_column % (sizeof(styles) / sizeof(styles[0]))]);
		graphs[i_column]->SetMarkerSize(1.5);
		graphs[i_column]->SetLineColor(colors[i_column % (sizeof(colors) / sizeof(colors[0]))]);
		graphs[i_column]->SetLineWidth(3);
		//graphs[i_column]->SetTitle(header[i_column+1].c_str());
		graphs[i_column]->SetTitle("");
		graphs[i_column]->SetName("");
		multigraph->Add(graphs[i_column], "P"); // Use "P" to draw only markers, "L" to draw only lines, "PL" to draw both
	}

	// Drawing the multi-graph
	multigraph->Draw("A"); // "A" Axis are drawn around the graph
	for (auto& graph : graphs) // Draw a curve on top of the markers if needed
	{
		graph->Draw("C same");
	}

	// Setting up axis titles
	multigraph->GetXaxis()->SetTitle("Trigger Rate / s");
	multigraph->GetYaxis()->SetTitle("Acceptance (%)");
	multigraph->GetXaxis()->CenterTitle();
	multigraph->GetYaxis()->CenterTitle();
	multigraph->GetXaxis()->SetLabelFont(132);
	multigraph->GetYaxis()->SetLabelFont(132);
	multigraph->GetXaxis()->SetTitleFont(132);
	multigraph->GetYaxis()->SetTitleFont(132);
	multigraph->GetYaxis()->SetLabelSize(0.07);
	multigraph->GetYaxis()->SetTitleSize(0.08);
	multigraph->GetXaxis()->SetLabelSize(0.07);
	multigraph->GetXaxis()->SetTitleSize(0.08);
	multigraph->GetYaxis()->SetTickLength(0.02);
	multigraph->GetYaxis()->SetNdivisions(505);
	multigraph->GetXaxis()->SetTitleOffset(1.05);
	multigraph->GetYaxis()->SetTitleOffset(0.8);
	multigraph->GetXaxis()->SetRangeUser(0, 5200);
	multigraph->GetYaxis()->SetRangeUser(85, 101);

	// Adding legend
	TLegend* legend = new TLegend(0.15, 0.21, 0.52, 0.57);//left, down, right, up
	legend->SetBorderSize(0);
	legend->SetTextFont(132);
	legend->SetTextSize(0.07);
	legend->SetTextAlign(12);
	legend->SetFillStyle(0);

	// Skip the first header entry because it's the label for the x-axis
	for (size_t i = 1; i < header.size(); i++)
	{
		//legend->AddEntry(graphs[i - 1], header[i].c_str(), "lp");
		legend->AddEntry(graphs[i-1], "", "LP");
		//"L" stands for line. This option indicates that a line sample should be drawn in the legend to represent the dataset.
		//"P" stands for marker(point).It means that a marker sample will be drawn in the legend.
		//So, when you use "LP" together, it tells ROOT to draw both a line and a marker in the legend entry for that dataset.
	}
	legend->Draw();

	// Adding text with the same colors as the curves
	TPaveText* colorText = new TPaveText(0.23, 0.21, 0.52, 0.57, "NDC");//left, down, right, up
	colorText->SetBorderSize(0);
	colorText->SetFillStyle(0); // Transparent
	colorText->SetTextFont(132);
	colorText->SetTextSize(0.07);
	colorText->SetTextAlign(12); //align = 10*HorizontalAlign + VerticalAlign, 12 means aligns the text to the left (horizontal alignment) and centered vertically relative to the specified coordinates
	// Add text entries with corresponding colors
	const char* seriesNames[] = { "LEGe", "XtRa1", "XtRa2", "MSD12", "MSD26" }; // to overlay a text box on top of or beside your legend to display the colors corresponding to each series
	for (size_t i = 0; i < sizeof(seriesNames) / sizeof(seriesNames[0]); i++)
	{
		auto textEntry = colorText->AddText(seriesNames[i]);
		textEntry->SetTextColor(colors[i]);
	}
	colorText->Draw();

	sprintf(filename, "%s%s", pathname, "Fig_PXCT_DDAS_Dead_Time");
	sprintf(figurename, "%s%s", filename, ".png");
	canvaspeak->SaveAs(figurename);
	sprintf(figurename, "%s%s", filename, ".eps");
	canvaspeak->SaveAs(figurename);

}
