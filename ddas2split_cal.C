/*
 * FOR PXCT. No Calibration! Used to investigate the gain drift of the detectors and gain matching, and calibration. use this after executing ddas2split.C
 *
 * Last update: Feb 5, 2024. Author: Lijie Sun. Contact: sunli@frib.msu.edu
 *
 * This program will convert DDAS ROOT files from the original form to a more simple form independent of the DDAS software package.
 * Currently, only energies for (up to) 16 channels are saved in coincidence, and histograms are produced.
 * This code will probably not run "as is" under different computer either than spdaq01
 *
 * Input: ROOT file from ddasdumper, standard mode. Not legacy mode!
 * Output: ROOT file with events split_tree, 16 histograms for 16 channels and 1 histogram of the total spectrum.
 *
 * Usage:
 * .L /usr/opt/ddas/6.1-001/lib/libddaschannel.so     //load ddas library
	.x ddas2split.C
 *
 */

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include <DDASEvent.h>
#include "TF1.h"
#include <fstream>
#include <iostream>
#include <vector>
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

void ddas2split_cal()
{
	//step 0: cd /user/pxct/readout/rootfile/ or use godata
	//step 1: . /usr/opt/root/root-6.24.06/bin/thisroot.sh or use ./r.sh
	//step 2: Ctrl+R then type .L /usr/opt/ddas/6.1-001/lib/libddaschannel.so     //load ddas library
	//step 3: .x ddas2split_cal.C
	// https://docs.nscl.msu.edu/daq/newsite/ddas-1.1/The.html
	// https://docs.nscl.msu.edu/daq/newsite/ddas-1.1/ddaschannel_8h_source.html

	long totalentries, i = 0;
	int Nchannels;
	char rawrootname[300];
	char calrootname[300];
	char txtname[300];
	char inpathname[300], outpathname[300];
	char filename[300];
	int runnumber = 200; // modify
	double timestampshift = 0;
	string line;

	//sprintf(inpathname, "%s", "/mnt/daqtesting/pxct/stagearea/");
	sprintf(outpathname, "%s", "/user/pxct/readout/rootfile/");
	sprintf(inpathname, "%s", "/user/pxct/readout/rootfile/");

	if (runnumber == 200)
	{
		sprintf(filename, "%s", "run0092_0095_0100_0105_LEGe_MSD26_241Am_inChamber_2mmCollimator_window1.5us_CFDdelay_adjusted");
	}

	// name both input and output root files accordingly
	sprintf(rawrootname, "%s%s%s", inpathname, filename, ".root");
	sprintf(txtname, "%s%s%s", inpathname, filename, ".dat");
	ofstream outfile(txtname, ios::out);
	sprintf(txtname, "%s%s%s", inpathname, filename, "_cal.dat");
	ifstream infile(txtname, ios::in);
	int iline = 0;
	double slope[400];
	while (getline(infile, line))
	{
		stringstream(line) >> iline >> slope[iline];
		cout << slope[iline] << endl;
		iline++;
	}
	TFile* pFile = new TFile(rawrootname);
	cout << "  input file: " << rawrootname << endl;
	TTree* pTree;
	pFile->GetObject("dchan", pTree);
	totalentries = pTree->GetEntries();
	cout << "  Entries=" << totalentries << endl;

	DDASEvent* pEvent = new DDASEvent;
	std::vector <ddaschannel*> pChan;
	pTree->SetBranchAddress("ddasevent", &pEvent);

	Double_t lege_e;
	Double_t lege_t;
	Double_t north_e;
	Double_t north_t;
	Double_t south_e;
	Double_t south_t;
	Double_t msd12_e;
	Double_t msd12_t;
	Double_t msd26_e;
	Double_t msd26_t;
	Double_t msdtotal_e;

	Double_t lege_e_low;
	Double_t lege_t_low;
	Double_t north_e_low;
	Double_t north_t_low;
	Double_t south_e_low;
	Double_t south_t_low;

	for (int iFile = 0; iFile < 192; iFile++)
	{
		sprintf(calrootname, "%s%s%s%04d%s", outpathname, filename, "_", iFile, "_cal.root");
		TFile* fout = new TFile(calrootname, "RECREATE");
		cout << "  output file: " << calrootname << endl;
		TTree* split_tree = new TTree("split_tree", "split_tree");

		// If you have closed the channels in CSRA, please comment out the corresponding branch lines. No other changes are needed.
		split_tree->Branch("lege_e", &lege_e, "lege_e/D");
		split_tree->Branch("lege_t", &lege_t, "lege_t/D");
		split_tree->Branch("msd26_e", &msd26_e, "msd26_e/D");
		split_tree->Branch("msd26_t", &msd26_t, "msd26_t/D");

		TH1D* hlege_e;
		TH1D* hmsd26_e;
		TH1D* htiming_lege_msd26;

		double binmax = slope[iFile] * 60000;
		hlege_e = new TH1D("hlege_e", "hlege_e", 60000, -0.0540888491631, 435.8137733893930); // 60000 channels, 7.3 eV per channel run0060
		hmsd26_e = new TH1D("hmsd26_e", "hmsd26_e", 60000, 0, binmax); // 60000 channels, 115 eV per channel
		htiming_lege_msd26 = new TH1D("htiming_lege_msd26", "htiming_lege_msd26", 3000, -1500, 1500); // 3000 channels, 1 ns per channel

		for (i = iFile * 1e6; i < (iFile+1) * 1e6; i++)
		{
			lege_e = 0; lege_t = 0; north_e = 0; north_t = 0; south_e = 0; south_t = 0; msd12_e = 0; msd12_t = 0; 	msd26_e = 0; msd26_t = 0;
			lege_e_low = 0; lege_t_low = 0; north_e_low = 0; north_t_low = 0; south_e_low = 0; south_t_low = 0;
			pTree->GetEntry(i);
			pChan = pEvent->GetData();
			Nchannels = (int)pChan.size();  //get active number of channels
			// cout<<"  number of active channels = "<<Nchannels<<endl;
			for (int j = 0; j < Nchannels; j++)
			{
				if (pChan[j]->GetChannelID() == 1)
				{
					//lege_e = pChan[j]->GetEnergy();
					//lege_e = pChan[j]->GetEnergy() * 0.0072712838511 - 0.0001420032512;
					lege_e = pChan[j]->GetEnergy() * 0.0072644643706 - 0.054088849163; // based on run0060
					hlege_e->Fill(lege_e);
					lege_t = pChan[j]->GetTime() + timestampshift;
				}

				if (pChan[j]->GetChannelID() == 8)
				{
					//msd26_e = pChan[j]->GetEnergy();
					//msd26_e = pChan[j]->GetEnergy() * 0.1152143446474 + 22.352486402900; // based on 148Gd and 241Am two peaks
					msd26_e = pChan[j]->GetEnergy() * slope[iFile]; // based on 241Am two peaks Run 100
					hmsd26_e->Fill(msd26_e);
					//hmsd26_e1keVbin->Fill(msd26_e + gRandom->Uniform(-0.5, 0.5));
					msd26_t = pChan[j]->GetTime() + timestampshift;
				}

				if (runnumber >= 92 && lege_e > 59 && lege_e < 60.1 && msd26_e > 5470 && msd26_e < 5510)
				{
					htiming_lege_msd26->Fill(lege_t - msd26_t);
				}
			} // (int j = 0; j < Nchannels; j++)
			split_tree->Fill();
			if ( i>= totalentries) break;
		} // for (i = iFile * 1e6; i < (iFile+1) * 1e6; i++)
		cout << " File: " << iFile << ";  Last entry: " << i << endl;
		outfile << "File:	" << iFile << "	Last entry:	" << i << "	Max:	" << hmsd26_e->GetBinCenter(hmsd26_e->GetMaximumBin()) << endl;
		fout->Write();
		fout->Close();
	} // for (long iFile = 0; iFile < 100; iFile++)
	pFile->Close();
}
