#include "TH1F.h"
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
//#include <math.h>
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
#include <TFile.h>
#include <TH1F.h>
#include <TTree.h>
#include <TFile.h>
#include <TH1F.h>
#include <TTree.h>

void extract_histograms_from_G4_rootfiles_DSL2()
{
// 	double Tau_values[] = { 0.0, 1.0, 2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 18.0, 30.0 };
// 	double SP_values[] = { 0.60, 1.00, 1.20 };
// 	double AC_values = 0.0;
// 	double Eg_values[] = { 4155.53, 4155.84, 4156.15 };

	double Tau_values[] = { 5.0 };
	double Eg_values[] = { 4155.84 };
	double Bkg_values[] = { 1.00 };
	double SP_values[] = { 1.00 };
	double AC_values = 0.0;

	const char* baseInputFileName = "F:/out/G4_rootfiles_with_tree_Eg4156/S31_Gamma4156_Eg%.2f_Tau%.1f_SP%.2f_AC%.1f_all.root";
	const char* baseOutputFileName = "F:/out/G4_rootfiles_with_tree_Eg4156/S31_Gamma4156_Eg%.2f_Tau%.1f_SP%.2f_AC%.1f.root";

	for (int iEg = 0; iEg < sizeof(Eg_values) / sizeof(Eg_values[0]); ++iEg)
	{
		for (int iSP = 0; iSP < sizeof(SP_values) / sizeof(SP_values[0]); ++iSP)
		{
			for (int iTau = 0; iTau < sizeof(Tau_values) / sizeof(Tau_values[0]); ++iTau)
			{
				char inputFileName[300];
				char outputFileName[300];
				sprintf(inputFileName, baseInputFileName, Eg_values[iEg], Tau_values[iTau], SP_values[iSP], AC_values);
				sprintf(outputFileName, baseOutputFileName, Eg_values[iEg], Tau_values[iTau], SP_values[iSP], AC_values);

				TFile* _file0 = TFile::Open(inputFileName);
				TTree* tree = (TTree*)_file0->Get("tree"); // the tree name is "tree"

				TH1F* Eg = new TH1F("Eg", "Eg", 8000, 0, 8000);
				tree->Draw("Clovere>>Eg", Form("Clovere>0&&Clovere<8000&&DSSD1e+DSSD2e>10000&&DSSD1e>100&&DSSD2e>100"));

				TFile* fout = new TFile(outputFileName, "RECREATE");
				Eg->Write();
				fout->Close();
			}
		}
	}
}