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
void extract_time_interval_between_neighboring_entries() // pxct daq dead time test
{
	char inpathname[300], outpathname[300];
	char calrootname[400];
	char anarootname[400];
	char filename[400];
	int runnumber;
	cout << "input run number: ";
	cin >> runnumber;
	sprintf(inpathname, "%s", "F:/e21010/pxct/");
	sprintf(outpathname, "%s", "F:/e21010/pxct/");
	if (runnumber == 304)	sprintf(filename, "%s", "run0304_Pulser_Ch0_100kHz_10nsWindow");
	if (runnumber == 305)	sprintf(filename, "%s", "run0305_Pulser_Ch0_100kHz_1000nsWindow");
	if (runnumber == 306)	sprintf(filename, "%s", "run0306_Pulser_Ch0_100kHz_10000nsWindow");
	if (runnumber == 307)	sprintf(filename, "%s", "run0307_Pulser_Ch0_100kHz_20000nsWindow");
	if (runnumber == 308)	sprintf(filename, "%s", "run0308_Pulser_Ch0_100kHz_30000nsWindow");
	if (runnumber == 309)	sprintf(filename, "%s", "run0309_Pulser_Ch0_10kHz_10000nsWindow");
	if (runnumber == 310)	sprintf(filename, "%s", "run0310_Pulser_Ch0_99kHz_10nsWindow");
	if (runnumber == 311)	sprintf(filename, "%s", "run0311_Pulser_Ch0_99kHz_Random_10nsWindow");
	if (runnumber == 312)	sprintf(filename, "%s", "run0312_Pulser_Ch0_99kHz_Random_100nsWindow");
	if (runnumber == 313)	sprintf(filename, "%s", "run0313_Pulser_Ch0_99kHz_100nsWindow");
	if (runnumber == 314)	sprintf(filename, "%s", "run0314_Pulser_Ch0_99kHz_Random_1000nsWindow");
	if (runnumber == 315)	sprintf(filename, "%s", "run0315_Pulser_Ch0_99kHz_Random_10000nsWindow");
	if (runnumber == 316)	sprintf(filename, "%s", "run0316_Pulser_Ch0_99kHz_Random_20000nsWindow");
	if (runnumber == 317)	sprintf(filename, "%s", "run0317_Pulser_Ch0_10kHz_1000nsWindow");
	if (runnumber == 318)	sprintf(filename, "%s", "run0318_Pulser_Ch0_10kHz_Random_1000nsWindow");
	if (runnumber == 319)	sprintf(filename, "%s", "run0319_Pulser_Ch0_2kHz_1000nsWindow");
	if (runnumber == 320)	sprintf(filename, "%s", "run0320_Pulser_Ch0_2kHz_Random_1000nsWindow");
	if (runnumber == 321)	sprintf(filename, "%s", "run0321_Pulser_Ch0249_3kHz_Random_100nsWindow");
	if (runnumber == 322)	sprintf(filename, "%s", "run0322_Pulser_Ch0249_3kHz_Random_100nsWindow_pileuprejection");
	if (runnumber == 323)	sprintf(filename, "%s", "run0323_Pulser_Ch0249_1kHz_Random_100nsWindow");

	sprintf(calrootname, "%s%s%s", inpathname, filename, "_cal.root");
	sprintf(anarootname, "%s%s%s", outpathname, filename, "_time_interval.root");
	
	TFile* fin = new TFile(calrootname);//after this statement, you can use any ROOT command1 for this rootfile
	cout << calrootname << endl;
	TTree* tree = (TTree*)fin->Get("tree");
	unsigned long totalentries = tree->GetEntries();
	cout << "  Entries=" << totalentries << endl;

	Double_t lege_e_low;
	Double_t lege_t_low;
	Double_t north_e_low;
	Double_t north_t_low;
	Double_t south_e_low;
	Double_t south_t_low;
	Double_t pulser_e;
	Double_t pulser_t;
	TBranch* b_lege_e_low;
	TBranch* b_lege_t_low;
	TBranch* b_north_e_low;
	TBranch* b_north_t_low;
	TBranch* b_south_e_low;
	TBranch* b_south_t_low;
	TBranch* b_pulser_e;
	TBranch* b_pulser_t;
	tree->SetBranchAddress("lege_e_low", &lege_e_low, &b_lege_e_low);
	tree->SetBranchAddress("lege_t_low", &lege_t_low, &b_lege_t_low);
	tree->SetBranchAddress("north_e_low", &north_e_low, &b_north_e_low);
	tree->SetBranchAddress("north_t_low", &north_t_low, &b_north_t_low);
	tree->SetBranchAddress("south_e_low", &south_e_low, &b_south_e_low);
	tree->SetBranchAddress("south_t_low", &south_t_low, &b_south_t_low);
	tree->SetBranchAddress("pulser_e", &pulser_e, &b_pulser_e);
	tree->SetBranchAddress("pulser_t", &pulser_t, &b_pulser_t);
	
	TFile* fout = new TFile(anarootname, "RECREATE");
	cout << "output file: " << anarootname << endl;
	//TTree* tree2 = new TTree("tree2", "tree2");
	TH1D* hinterval_lege = new TH1D("hinterval_lege", "hinterval_lege", 1000000, 0, 100000);
	TH1D* hinterval_north = new TH1D("hinterval_north", "hinterval_north", 1000000, 0, 100000);
	TH1D* hinterval_south = new TH1D("hinterval_south", "hinterval_south", 1000000, 0, 100000);
	TH1D* hinterval_pulser = new TH1D("hinterval_pulser", "hinterval_pulser", 1000000, 0, 100000);

	//sprintf(filename, "%s%s", pathname, "11.txt"); // modify
	//ofstream outfile(filename, ios::out);
	double t1_lege = 0, t2_lege = 0, t1_north = 0, t2_north = 0, t1_south = 0, t2_south = 0, t1_pulser = 0, t2_pulser = 0;

	//tree->Scan("lege_t:lege_e", "lege_e>100", "precision=10");
	tree->GetEntry(0);
	t1_lege = lege_t_low / 1000;
	t1_north = north_t_low / 1000;
	t1_south = south_t_low / 1000;
	t1_pulser = pulser_t / 1000;

	for (unsigned long i = 1; i < totalentries; i++)
	{
		tree->GetEntry(i);

		t2_lege = lege_t_low / 1000;
		hinterval_lege->Fill(t2_lege - t1_lege);
		t1_lege = t2_lege;

		t2_north = north_t_low / 1000;
		//if (north_e_low>11780&&north_e_low<12800)
		hinterval_north->Fill(t2_north - t1_north);
		t1_north = t2_north;

		t2_south = south_t_low / 1000;
		hinterval_south->Fill(t2_south - t1_south);
		t1_south = t2_south;

		t2_pulser = pulser_t / 1000;
		hinterval_pulser->Fill(t2_pulser - t1_pulser);
		t1_pulser = t2_pulser;
		
		// Write the values to the file, adjust format as needed
		//outfile << std::fixed << std::setprecision(10) << lege_t << "\t" << std::setprecision(1) << lege_e << std::endl;
	}
	hinterval_lege->SetTitle(filename);
	hinterval_lege->GetXaxis()->SetTitle("Time Difference between entries (us)");
	hinterval_lege->GetYaxis()->SetTitle("Counts per 100 ns");
	hinterval_lege->GetXaxis()->CenterTitle();
	hinterval_lege->GetYaxis()->CenterTitle();
	hinterval_north->SetTitle(filename);
	hinterval_north->GetXaxis()->SetTitle("Time Difference between entries (us)");
	hinterval_north->GetYaxis()->SetTitle("Counts per 100 ns");
	hinterval_north->GetXaxis()->CenterTitle();
	hinterval_north->GetYaxis()->CenterTitle();
	hinterval_south->SetTitle(filename);
	hinterval_south->GetXaxis()->SetTitle("Time Difference between entries (us)");
	hinterval_south->GetYaxis()->SetTitle("Counts per 100 ns");
	hinterval_south->GetXaxis()->CenterTitle();
	hinterval_south->GetYaxis()->CenterTitle();
	hinterval_pulser->SetTitle(filename);
	hinterval_pulser->GetXaxis()->SetTitle("Time Difference between entries (us)");
	hinterval_pulser->GetYaxis()->SetTitle("Counts per 100 ns");
	hinterval_pulser->GetXaxis()->CenterTitle();
	hinterval_pulser->GetYaxis()->CenterTitle();
	hinterval_lege->Write();
	hinterval_north->Write();
	hinterval_south->Write();
	hinterval_pulser->Write();
	fin->Close();
	//outfile.close(); // Close the output file stream

}