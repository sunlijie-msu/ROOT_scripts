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
void scan_tree() // output tree scan to a txt file
{
	char pathname[300];
	char filename[300];
	sprintf(pathname, "%s", "F:/e21010/pxct/");
	sprintf(filename, "%s%s", pathname, "run0211_LEGe_152Eu_Z2707_inChamber_vacuum_window10us_CFDdelay_adjusted_cal.root"); // modify
	TFile* fin = new TFile(filename);//after this statement, you can use any ROOT command1 for this rootfile
	cout << filename << endl;
	TTree* tree = (TTree*)fin->Get("tree");
	unsigned long nentries = tree->GetEntries();
	cout << "  Entries=" << nentries << endl;

	sprintf(filename, "%s%s", pathname, "11.txt"); // modify
	ofstream outfile(filename, ios::out);

	double lege_t, lege_e;
	tree->SetBranchAddress("lege_t", &lege_t);
	tree->SetBranchAddress("lege_e", &lege_e);
	//tree->Scan("lege_t:lege_e", "lege_e>100", "precision=10");

	for (unsigned long i = 0; i < 1000000; i++)
	{
		tree->GetEntry(i);
		// Write the values to the file, adjust format as needed
		outfile << std::fixed << std::setprecision(10) << lege_t << "\t" << std::setprecision(1) << lege_e << std::endl;
	}

	fin->Close();
	outfile.close(); // Close the output file stream

}