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
void Fill_histogram_random()// read data from a random function and fill a TH1F
{
	double y;
	TH1F* hbackground_flat = new TH1F("hbackground_flat", "hbackground_flat", 8000, 0, 8000);
	TH1F* hbackground_linear = new TH1F("hbackground_linear", "hbackground_linear", 8000, 0, 8000);
	TF1* pol1 = new TF1("pol1", "[0]*x+[1]", 0, 8000);
	pol1->SetParameters(-0.09, 1000);
	for (int i = 0; i < 3200000; i++)
	{
		y = pol1->GetRandom();
		hbackground_linear->Fill(y);//Add one y-value at a time
		y = gRandom->Uniform(0, 1) * 8000;
		hbackground_flat->Fill(y);//Add one y-value at a time
	}
	TFile* fout = new TFile("F:/out/G4_rootfiles_with_tree_Eg4156/Fakebackground_high.root", "RECREATE");
	hbackground_flat->Write();
	hbackground_linear->Write();
	fout->Close();
}