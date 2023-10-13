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
void remove_Compton_bkg()
{
	TFile *fin = new TFile("D:/X/out/DSL/S31_Gamma1248_Eg1248.40_Tau600.0_SP1.00_AC0.0_all1.root");
	TH1F *hGriffinAddback_alpha49_t0=new TH1F("hGriffinAddback_alpha49_t0","hGriffinAddback_alpha49_t0",8000,0.,8000);
	tree->Draw("Clovere>>hGriffinAddback_alpha49_t0","Clovere>0&&DSSD1e+DSSD2e>47000","");
	
	TF1 *G4pol1=new TF1("G4pol1","[0]*x+[1]",1154,1230);//G4 bkg
	hGriffinAddback_alpha49_t0->Fit("G4pol1","MLE","",1154,1230);
	hGriffinAddback_alpha49_t0->Fit("G4pol1","MLE","",1154,1230);
	hGriffinAddback_alpha49_t0->Fit("G4pol1","MLE","",1154,1230);
	double aG4 = G4pol1->GetParameter(0);
	double bG4 = G4pol1->GetParameter(1);
	for (int x=1154;x<1410;x++)
	{
		int xbin=hGriffinAddback_alpha49_t0->FindBin(x);
		int count = hGriffinAddback_alpha49_t0->GetBinContent(xbin);
		if (bG4+aG4*x<0) break;
		hGriffinAddback_alpha49_t0->SetBinContent(xbin,count-(bG4+aG4*x));
		//cout<<x<<""<<count<<""<<(bG4+aG4*x)<<endl;
	}

	for (int x=1154;x<1410;x++)
	{
		int xbin=hGriffinAddback_alpha49_t0->FindBin(x);
		int count = hGriffinAddback_alpha49_t0->GetBinContent(xbin);
		hGriffinAddback_alpha49_t0->SetBinContent(xbin,count+2); // +2 to mimic the bkg level in data
		//cout<<x<<""<<count<<""<<(bG4+aG4*x)<<endl;
	}
	hGriffinAddback_alpha49_t0->Rebin(2);
	TFile *fout = new TFile("D:/X/out/DSL/nice_test1248_systematic_0.4k.root","RECREATE");
	hGriffinAddback_alpha49_t0->Write();
	fout->Close();
}
