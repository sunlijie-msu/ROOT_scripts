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
void Fill_histogram_random_decay_pxct_237Np()// read data from a random function and fill a TH1F
{
	double y;
	TH1D* h_decay_exponential_GetRandom = new TH1D("h_decay_exponential_GetRandom", "h_decay_exponential_GetRandom", 2000, 0, 2000);
	TH1D* h_decay_exponential_FillRandom = new TH1D("h_decay_exponential_FillRandom", "h_decay_exponential_FillRandom", 2000, 0, 2000);
	TF1* SiDEC = new TF1("SiDEC", "[0]*0.693147/[1]*exp(x/(-[1]/0.693147))+[2]", 0, 2000); // exponential decay (N, T, B)
	SiDEC->SetParameters(2e5, 68, 2);
	for (int i = 0; i < 2e5; i++)
	{
		y = SiDEC->GetRandom();
		h_decay_exponential_GetRandom->Fill(y);//Add one y-value at a time
	}
	h_decay_exponential_FillRandom->FillRandom("SiDEC", 2e5);
	TFile* fout = new TFile("F:/e21010/pxct/Fake_decay_2e5.root", "RECREATE");
	h_decay_exponential_GetRandom->Write();
	h_decay_exponential_FillRandom->Write();
	fout->Close();
}