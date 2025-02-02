#include <iostream>
#include <fstream>
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
void Fill_histogram_random_decay_pxct_237Np()// random sample from a exp decay function and fill a histogram; Test Chris's ideas.
{
	// If gRandom is not set, initialize it as a TRandom3 with a time-dependent seed.
	if (!gRandom) {
		gRandom = new TRandom3(0);
	}
	double true_time, measured_time1, measured_time2, measured_time3;

	TH1D* h_decay_exponential_GetRandom = new TH1D("h_decay_exponential_GetRandom", "h_decay_exponential_GetRandom", 4000, -2000, 2000);

	TH1D* h_decay_exponential_FillRandom = new TH1D("h_decay_exponential_FillRandom", "h_decay_exponential_FillRandom", 4000, -2000, 2000);

	TH1D* h_decay_exponential_TRandom3 = new TH1D("h_decay_exponential_TRandom3", "h_decay_exponential_TRandom3", 4000, -2000, 2000);

	TH1D* h_decay_exponential_Gaus1_TRandom3 = new TH1D("h_decay_exponential_Gaus1_TRandom3", "h_decay_exponential_Gaus1_TRandom3", 4000, -2000, 2000);
	TH1D* h_decay_exponential_Gaus2_TRandom3 = new TH1D("h_decay_exponential_Gaus2_TRandom3", "h_decay_exponential_Gaus2_TRandom3", 4000, -2000, 2000);
	TH1D* h_decay_exponential_Gaus3_TRandom3 = new TH1D("h_decay_exponential_Gaus3_TRandom3", "h_decay_exponential_Gaus3_TRandom3", 4000, -2000, 2000);

	TF1* SiDEC = new TF1("SiDEC", "[0]*0.693147/[1]*exp(x/(-[1]/0.693147))+[2]", -2000, 2000); // exponential decay (N, T, B)
	SiDEC->SetParameters(1e7, 68, 0);

	for (int i = 0; i < 1e7; i++)
	{
		true_time = SiDEC->GetRandom();
		h_decay_exponential_GetRandom->Fill(true_time);//Add one value at a time
	}

	h_decay_exponential_FillRandom->FillRandom("SiDEC", 1e7);

	for (int i = 0; i < 1e7; i++)
	{
		true_time = gRandom->Exp(68 / log(2.0)); // 68 / log(2.0) = 98.103263
		h_decay_exponential_TRandom3->Fill(true_time);//Add one value at a time

		measured_time1 = true_time + gRandom->Gaus(0, 5);
		h_decay_exponential_Gaus1_TRandom3->Fill(measured_time1);

		measured_time2 = true_time + gRandom->Gaus(0, 30);
		h_decay_exponential_Gaus2_TRandom3->Fill(measured_time2);

		measured_time3 = true_time + gRandom->Gaus(0, 200);
		h_decay_exponential_Gaus3_TRandom3->Fill(measured_time3);
	}

	TFile* fout = new TFile("F:/e21010/pxct/Fake_decay_1e7.root", "RECREATE");
	h_decay_exponential_GetRandom->Write();
	h_decay_exponential_FillRandom->Write();
	h_decay_exponential_TRandom3->Write();
	h_decay_exponential_Gaus1_TRandom3->Write();
	h_decay_exponential_Gaus2_TRandom3->Write();
	h_decay_exponential_Gaus3_TRandom3->Write();

	fout->Close();
}