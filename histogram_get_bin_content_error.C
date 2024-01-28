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
void histogram_get_bin_content_error()//get BinContent and BinError from histograms for Surmise data.csv
{
	char rootname[300];
	char histoname[300];
	char txtname[300];
	sprintf(rootname, "%s", "F:/e21010/pxct/run0092_93_94_95_LEGe_MSD26_241Am_inChamber_2mmCollimator_window1.5us_CFDdelay_adjusted_cal.root");
	//sprintf(rootname, "%s", "F:/out/G4_rootfiles_with_tree_Eg4156/Fakedata_S31_Gamma4156_Eg4155.84_Tau0.0_SP1.00_AC0.0_manipulated.root");
	TFile* fin_data = new TFile(rootname);//after this statement, you can use any ROOT command for this rootfile
	TH1F* h1;

	sprintf(txtname, "%s", "F:/e21010/pxct/htiming_lege_msd26.dat");
	//sprintf(txtname, "%s", "F:/out/G4_rootfiles_with_tree_Eg4156/Fakedata_S31_Gamma4156_Eg4155.84_Tau0.0_SP1.00_AC0.0_manipulated.dat");
	ofstream outfile(txtname, ios::out);
	sprintf(histoname, "%s", "htiming_lege_msd26");
	//sprintf(histoname, "%s", "hfakedata");
	h1 = (TH1F*)fin_data->Get(histoname); //Get spectrum
	h1->SetBinErrorOption(TH1::kPoisson);
	for (int i = 1; i <= h1->GetNbinsX(); i++)
	{
		outfile << h1->GetBinCenter(i) << "	";
		outfile << h1->GetBinContent(i) << "	";
		outfile << h1->GetBinErrorLow(i) << "	";
		outfile << h1->GetBinErrorUp(i) << endl;
	}
	outfile.close();

}