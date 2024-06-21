#include <iostream>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <time.h>
#include <map>
#include <TROOT.h>
#include <TSystem.h>
#include <TStyle.h>
#include <TCutG.h>
#include<TChain.h>
#include "TChain.h"
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
#include "TMath.h"
#include "TPaveText.h"
#include "TPad.h"
#include "string.h"
#include "stdlib.h"
#include "stdio.h"
#include "TLegend.h"
using namespace std;
void Fill_histogram_from_txt_byFill_p_width()//read bincontent from a txt file and create a histo and save it as a new root file.
{// input bin by bin, not event by event
	int ii,jj;
	float Channel[500];
	float x,y,xe,ye;
	long Count[500];
	char outrootname[200], datname[200], tname[200];
	int Entries=300;//modify
	sprintf(tname, "%s", "D:/X/out/Alex_Brown_Proton_width");
	sprintf(datname, "%s%s", tname, ".dat");
	sprintf(outrootname, "%s%s", tname, ".root");
	ifstream infile(datname,ios::in);
	TFile *fout = new TFile(outrootname,"RECREATE");//It's better to define histograms and then define fout, in case of draw bugs.
	TCanvas* canvaspeak = new TCanvas("canvaspeak", "canvaspeak", 1300, 800);//
	canvaspeak->cd();//
	canvaspeak->SetTopMargin(0.029);
	canvaspeak->SetRightMargin(0.04);
	canvaspeak->SetLeftMargin(0.11);
	canvaspeak->SetBottomMargin(0.15);
	canvaspeak->SetFrameLineWidth(3);
	TH1F *h_p_width_ratio = new TH1F("h_p_width_ratio","h_p_width_ratio",400,-4,4);//create a histogram
	string line;
	while ( getline(infile, line) )
	{
		stringstream(line) >> x;
		h_p_width_ratio->Fill(log10(x));
		//cout <<x<<endl;
	}
	h_p_width_ratio->Draw("e");
	h_p_width_ratio->GetXaxis()->SetTitle("Log #Gammap Ratio of Brown/Rauscher");
	h_p_width_ratio->GetYaxis()->SetTitle("Counts");
	h_p_width_ratio->SetLineWidth(3);
	h_p_width_ratio->SetLineColor(1);
	//h_p_width_ratio->SetStats(0);
	h_p_width_ratio->SetBinErrorOption(TH1::kPoisson);
	h_p_width_ratio->SetTitle("");
	h_p_width_ratio->GetXaxis()->CenterTitle();
	h_p_width_ratio->GetYaxis()->CenterTitle();
	h_p_width_ratio->GetXaxis()->SetLabelFont(132);
	h_p_width_ratio->GetYaxis()->SetLabelFont(132);
	h_p_width_ratio->GetXaxis()->SetLabelSize(0.06);
	h_p_width_ratio->GetYaxis()->SetLabelSize(0.06);
	h_p_width_ratio->GetXaxis()->SetTitleFont(132);
	h_p_width_ratio->GetYaxis()->SetTitleFont(132);
	h_p_width_ratio->GetXaxis()->SetTitleOffset(1.0);
	h_p_width_ratio->GetYaxis()->SetTitleOffset(0.8);
	h_p_width_ratio->GetXaxis()->SetTitleSize(0.07);
	h_p_width_ratio->GetYaxis()->SetTitleSize(0.07);
	//h_p_width_ratio->GetXaxis()->SetRangeUser(100, 6000);
	//h_p_width_ratio->GetYaxis()->SetRangeUser(0, 30000);
	h_p_width_ratio->GetYaxis()->SetNdivisions(505);
	h_p_width_ratio->GetYaxis()->SetTickLength(0.015);
	gPad->RedrawAxis();
	cout<<"RMS= "<<setiosflags(ios::fixed)<<setprecision(9)<<h_p_width_ratio->GetRMS();//you can also find the RMS in the rootfile
	fout->Write();
	fout->Close();
}