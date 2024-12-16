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
void histogram_get_bin_content_error()//get BinContent and BinError from DSL data histograms for Surmise data.csv
{
	char rootname[300];
	char histoname[300];
	char txtname[300];
	char pathname[300];
	char filename[300];
	sprintf(pathname, "%s", "F:/out/");
	sprintf(filename, "%s%s", pathname, "testadd");
	sprintf(rootname, "%s%s", filename, ".root");
	sprintf(txtname, "%s%s", filename, ".csv");
	cout << rootname << endl;
	TFile* fin_data = new TFile(rootname);//after this statement, you can use any ROOT command for this rootfile
	TH1F* h1;
	ofstream outfile(txtname, ios::out);
	sprintf(histoname, "%s", "centersum2");
	h1 = (TH1F*)fin_data->Get(histoname); //Get spectrum
	h1->Rebin(5);
	h1->Sumw2(kFALSE);
	h1->SetBinErrorOption(TH1::kPoisson);
	for (int i = 1; i <= h1->GetNbinsX(); i++)
	{
		outfile << h1->GetBinCenter(i) << ",";
		outfile << h1->GetBinContent(i) << ",";
		outfile << h1->GetBinErrorLow(i) << ",";
		outfile << h1->GetBinErrorUp(i) << endl;
	}
	outfile.close();
}