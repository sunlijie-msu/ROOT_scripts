#include <iostream>
#include <fstream>
//#include <iomanip.h>
#include <math.h>
#include <time.h>
#include <map>
#include <TROOT.h>
#include <TSystem.h>
#include <TStyle.h>
#include <TCutG.h>
#include<TChain.h>
#include "TColor.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TEllipse.h"
#include "TString.h"
#include "TGraph.h"
#include "TSpectrum.h"
#include "TVirtualFitter.h"
#include "TRandom.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TF2.h"
#include "TExec.h"
#include "TMath.h"
#include "TPaveText.h"
#include "TPad.h"
#include "string.h"
#include "stdlib.h"
#include "stdio.h"
#include "TLegend.h"
using namespace std;
void Fill_2D_from_txt_Tau_0fs_OmegaGamma()//read bincontent from a txt file and create a histo and save it as a new root file.
{// input bin by bin, not event by event
	char outrootname[200], infilename[200], outfigname[200], tname[200];
	int Entries=300;//modify
	sprintf(tname, "%s", "D:/X/out/Bayesian_VS/Resonance_Strength/Tau_0fs_OmegaGamma_Joint");
	sprintf(infilename, "%s%s", tname, ".dat");
	sprintf(outrootname, "%s%s", tname, ".root");
	sprintf(outfigname, "%s%s", tname, ".png");

	TFile* fout = new TFile(outrootname, "RECREATE");
	TCanvas* c1 = new TCanvas("c1", "c1", 900, 900);
	TH1F* hTau = new TH1F("hTau", "hTau", 90, 0, 6);
	TH1F* hOmegaGamma = new TH1F("hOmegaGamma", "hOmegaGamma", 1000, 0, 1000);
	TH2F* hTauOmegaGamma = new TH2F("hTauOmegaGamma", "hTauOmegaGamma", 180, 0, 6, 1000, 0, 1000); // TH2(name, title, nbinsx, xlow, xup, nbinsy, ylow, yup);

	ifstream infile(infilename, ios::in);//The data that need to be fitted
	string line;
	stringstream ss;
	int irow = 0;
	int icolumn = 0;
	double parameters[2] = { 0 };

	while (getline(infile, line)) // works perfectly for txt files tab delimited. Read in row-by-row
	{
		// if (irow++ == 0) continue; // skip the row headings
		ss.clear(); //clear(): Used to clear the stream
		ss.str(line); //str(): To get and set the string object whose content is present in stream
		icolumn = 0; // column number 0-4 are parameters; 5-n are parameters/bincounts
		while (!ss.fail()) // not the end of a line
		{
			ss >> parameters[icolumn++]; //operator >> : This is used to read from stringstream object.
		}
		hTau->Fill(parameters[0]);
		hOmegaGamma->Fill(parameters[1]);
		hTauOmegaGamma->Fill(parameters[0], parameters[1]);
	}
	TPad* pad1 = new TPad("pad1", "pad1", 0.0, 0.8, 0.8, 1.0); // xlow,  ylow,  xup,  yup,
	TPad* pad2 = new TPad("pad2", "pad2", 0.0, 0.0, 0.8, 0.8);
	TPad* pad3 = new TPad("pad3", "pad3", 0.8, 0.0, 1.0, 0.8);
	pad1->SetTopMargin(0.02);
	pad1->SetRightMargin(0.02);
	pad1->SetLeftMargin(0.18);
	pad1->SetBottomMargin(0.01);
	pad1->SetFrameLineWidth(0);
	
	pad2->SetTopMargin(0.02);
	pad2->SetRightMargin(0.02);
	pad2->SetLeftMargin(0.18);
	pad2->SetBottomMargin(0.17);
	pad2->SetFrameLineWidth(0);

	pad3->SetTopMargin(0.02);
	pad3->SetRightMargin(0.02);
	pad3->SetLeftMargin(0.02);
	pad3->SetBottomMargin(0.17);
 	pad3->SetFrameLineWidth(0);

	pad1->Draw();
	pad2->Draw();
	pad3->Draw();

	pad1->cd();
	hTau->SetFillColor(kCyan - 3);
	hTau->Draw("bar0");
	hTau->GetXaxis()->SetNdivisions(0);
	hTau->GetYaxis()->SetNdivisions(0);

	pad3->cd();
	hOmegaGamma->SetFillColor(kCyan - 3);
	hOmegaGamma->Draw("hbar0");
	//hOmegaGamma->SetBarWidth(0.3);
	hOmegaGamma->GetXaxis()->SetNdivisions(0);
	hOmegaGamma->GetYaxis()->SetNdivisions(0);

	pad2->cd();
	//pad2->SetLogz();
	//hTauOmegaGamma->SetDrawOption("col");
	hTauOmegaGamma->Draw("col");
	hTauOmegaGamma->GetXaxis()->SetTitle("#it{#tau} (fs)");
	hTauOmegaGamma->GetYaxis()->SetTitle("#it{#omega#gamma} (#mueV)");
	hTauOmegaGamma->GetXaxis()->SetLabelFont(132);
	hTauOmegaGamma->GetYaxis()->SetLabelFont(132);
	hTauOmegaGamma->GetXaxis()->SetTitleFont(132);
	hTauOmegaGamma->GetYaxis()->SetTitleFont(132);
	hTauOmegaGamma->GetXaxis()->SetLabelSize(0.05);
	hTauOmegaGamma->GetYaxis()->SetLabelSize(0.05);
	hTauOmegaGamma->GetXaxis()->SetTitleSize(0.05);
	hTauOmegaGamma->GetYaxis()->SetTitleSize(0.05);
	hTauOmegaGamma->GetXaxis()->CenterTitle();
	hTauOmegaGamma->GetYaxis()->CenterTitle();
	hTauOmegaGamma->GetXaxis()->SetTitleOffset(1.4);
	hTauOmegaGamma->GetYaxis()->SetTitleOffset(1.6);
	hTauOmegaGamma->GetXaxis()->SetNdivisions(510);//n = n1 + 100*n2 + 10000*n3
	hTauOmegaGamma->GetYaxis()->SetNdivisions(505);//n = n1 + 100*n2 + 10000*n3
	hTauOmegaGamma->GetXaxis()->SetTickLength(0.01);
	hTauOmegaGamma->GetYaxis()->SetTickLength(0.01);
	hTauOmegaGamma->SetContour(99);

	TLine* LineOmegaGamma = new TLine(0, 43.024477061004276, 6, 43.024477061004276); // TLine(x1,y1,x2,y2)
	LineOmegaGamma->SetLineColor(kRed);
	LineOmegaGamma->SetLineWidth(2);
	LineOmegaGamma->SetLineStyle(7);
	LineOmegaGamma->Draw("");//"R" means the line is drawn with the current line attributes
	TLine* LineTau = new TLine(2.457093962700432, 0, 2.457093962700432, 1000); // TLine(x1,y1,x2,y2)
	LineTau->SetLineColor(kRed);
	LineTau->SetLineWidth(2);
	LineTau->SetLineStyle(7);
	LineTau->Draw("");//"R" means the line is drawn with the current line attributes

	gStyle->SetOptStat(0);
	gStyle->SetOptTitle(0);
	gStyle->SetPalette(kDeepSea);
	gStyle->SetLineWidth(2); // set the width of the axis lines, not frame lines
	TColor::InvertPalette();
	//pad2->RedrawAxis();
	
	c1->Update();
	c1->Draw();
	c1->SaveAs(outfigname);
	fout->Write();
	// fout->Close();
}