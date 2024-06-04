#include "TH1F.h"
//#include <cmath> //can't use pow() with this header
#include <stdlib.h>
#include "TMinuit.h"
#include "TFumili.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TStyle.h"
#include "TLine.h"
#include <iostream>
#include <fstream>
#include <math.h>
#include <time.h>
#include <map>
#include <TROOT.h>
#include <TSystem.h>
#include <TStyle.h>
#include <TCutG.h>
#include <TChain.h>
#include <TMinuit.h>
#include "TChain.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TEllipse.h"
#include "TString.h"
#include "TGraph.h"
#include "TSpectrum.h"
#include "TVirtualFitter.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "TRandom.h"
#include <TRandom3.h>
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TF1.h"
#include "TMath.h"
#include "TPaveText.h"
#include "TPad.h"
#include "string.h"
#include "stdlib.h"
#include "stdio.h"
#include "TLegend.h"
#include "TPaletteAxis.h"
using namespace std;
void tree_GetEntries_from_G4_rootfiles_PXCT()//get BinContent and BinError from histograms for Surmise data.csv
{
	//vector<double> Eg_values = { 1, 2, 3, 4, 5, 5.6, 6, 6.2, 7, 7.2, 8, 8.6, 9, 10, 11, 11.1, 11.2, 11.5, 12, 13, 14, 15, 16, 17, 18, 19, 20, 25, 30, 35, 40, 45, 50 }; // 152Eu simulation points
	vector<double> Eg_values = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 30, 35, 40, 45, 50 }; // 241Am simulation points
	vector<double> Eg_thresholds(Eg_values.size());

	for (size_t i = 0; i < Eg_values.size(); ++i)
	{
		Eg_thresholds[i] = Eg_values[i] - 0.6;
	}
	char rootname[300];
	char txtname[300];
	char pathname[300];
	sprintf(pathname, "%s", "F:/out/G4_pxct_root_files/");
	sprintf(txtname, "%sExG4_LowEnergy_all.txt", pathname);
	TFile* fin_data;
	ofstream outfile(txtname, ios::out);

	for (size_t i = 0; i < Eg_values.size(); i++)
	{
		sprintf(rootname, "%sExG4_241Am_LowEnergy_%gkeV_all.root", pathname, Eg_values[i]);
		fin_data = new TFile(rootname);//after this statement, you can use any ROOT command for this rootfile
		cout << rootname << endl;
		TTree* tree = (TTree*)fin_data->Get("tree");
		double threshold = Eg_thresholds[i];
		Long64_t nEntries = tree->GetEntries(Form("LEGe_e>%f", threshold));
		outfile << Eg_values[i] << "\t" << nEntries << endl;
	}
	outfile.close();
}