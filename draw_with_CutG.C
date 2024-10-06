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
void draw_with_CutG()
{
	// **1. Open the Cut File and Retrieve ProtonCut**
	TFile* cutFile = TFile::Open("F:/e21010/pxct/ProtonCut.root");
	if (!cutFile || cutFile->IsZombie())
	{
		std::cerr << "Error: Cannot open ProtonCut.root" << std::endl;
		return;
	}

	TCutG* ProtonCut = (TCutG*)cutFile->Get("ProtonCut");
	if (!ProtonCut)
	{
		std::cerr << "Error: Cannot retrieve ProtonCut from ProtonCut.root" << std::endl;
		cutFile->Close();
		return;
	}

	// **2. Open the Data File and Access the Tree**
	TFile* dataFile = TFile::Open("F:/e21010/pxct/ExG4_60Ga_20days_9000pps_added_beta_bkg.root");
	if (!dataFile || dataFile->IsZombie())
	{
		std::cerr << "Error: Cannot open data file" << std::endl;
		cutFile->Close();
		return;
	}

	TTree* tree = (TTree*)dataFile->Get("tree"); // Using the actual tree name "tree"
	if (!tree)
	{
		std::cerr << "Error: Cannot find TTree named 'tree' in data file" << std::endl;
		dataFile->Close();
		cutFile->Close();
		return;
	}

	// **3. Set Branch Addresses**
	double MSD12_e, MSD26_e;
	tree->SetBranchAddress("MSD12_e", &MSD12_e);
	tree->SetBranchAddress("MSD26_e", &MSD26_e);

	// **4. Create and Configure the Histogram**
	TH1D* histo = new TH1D("histo", "MSD12_e + MSD26_e", 9000, 0, 9000);

	// **5. Loop Over Tree Entries and Apply Cuts**
	Long64_t nentries = tree->GetEntries();
	for (Long64_t i = 0; i < nentries; i++)
	{
		tree->GetEntry(i);

		// Apply numerical cuts
		if (MSD12_e > 30 && MSD12_e < 1000 && (MSD12_e + MSD26_e) > 30 && (MSD12_e + MSD26_e) < 9000)
		{
			// Apply graphical cut using TCutG. MSD26_e is x and MSD12_e is y!
			if (ProtonCut->IsInside(MSD26_e, MSD12_e))
			{
				histo->Fill(MSD12_e + MSD26_e);
			}
		}
	}

	// **6. Set Up the Canvas and Draw the Histogram**
	TCanvas* canvaspeak = new TCanvas("canvaspeak", "Proton Cut Histogram", 1300, 700);
	canvaspeak->SetTopMargin(0.035);
	canvaspeak->SetRightMargin(0.04);
	canvaspeak->SetLeftMargin(0.14);
	canvaspeak->SetBottomMargin(0.20);
	canvaspeak->SetFrameLineWidth(3);
	gStyle->SetFrameLineWidth(3);


	histo->SetLineWidth(2);
	histo->Rebin(1);
	histo->SetLineColor(2);
	histo->SetStats(0);
	histo->SetTitle("");
	histo->GetXaxis()->SetTitle("Energy (keV)");
	histo->GetYaxis()->SetTitle("Counts per 1 keV");
	histo->GetXaxis()->CenterTitle();
	histo->GetYaxis()->CenterTitle();
	histo->GetXaxis()->SetLabelFont(132);
	histo->GetYaxis()->SetLabelFont(132);
	histo->GetXaxis()->SetLabelSize(0.07);
	histo->GetYaxis()->SetLabelSize(0.07);
	histo->GetXaxis()->SetTitleFont(132);
	histo->GetYaxis()->SetTitleFont(132);
	histo->GetXaxis()->SetTitleOffset(1.2);
	histo->GetYaxis()->SetTitleOffset(0.9);
	histo->GetXaxis()->SetTitleSize(0.08);
	histo->GetYaxis()->SetTitleSize(0.08);
	histo->GetYaxis()->SetTickLength(0.015);
	histo->GetXaxis()->SetRangeUser(0, 9000);
	int ymax = histo->GetMaximum();
	histo->GetYaxis()->SetRangeUser(1, ymax * 1.4);
	histo->GetYaxis()->SetNdivisions(505);
	gPad->SetLogy();
	histo->SetBinErrorOption(TH1::kPoisson);
	histo->Draw("hist");

	//tree->SetLineWidth(2);
	//tree->SetLineColor(1);
	//tree->Draw("MSD12_e", "MSD12_e>30&&MSD12_e<1000&&MSD26_e<100", "same"); // All low-energy proton

	tree->SetLineWidth(2);
	tree->SetLineColor(kGreen + 2);
	tree->Draw("MSD12_e+MSD26_e", "MSD12_e>1200&&MSD12_e<4000&&(MSD12_e+MSD26_e)>1200&&(MSD12_e+MSD26_e)<9000", "same"); // All alpha 98422


	// **7. Optional: Save the Canvas**
	canvaspeak->SaveAs("F:/e21010/pxct/60Ga_histo_with_proton_cut.png");

	// **8. Clean Up (Optional)**
	// Close the files if you're done with them
	// dataFile->Close();
	// cutFile->Close();
}
