/*
 * For PXCT run0079-0091. No gain drift. No need to gain match.
   This code does calibration. Output a huge cal.root file.
 *
 * Last update: Nov 7, 2024. Author: Lijie Sun. Contact: sunli@frib.msu.edu
 *
 * This program will convert DDAS ROOT files from the original form to a more simple form independent of the DDAS software package.
 * Currently, only energies for (up to) 16 channels are saved in coincidence, and histograms are produced.
 * This code will probably not run "as is" under different computer either than spdaq01
 *
 * Input: ROOT file from ddasdumper, standard mode. Not legacy mode!
 * Output: ROOT file with events tree, 16 histograms for 16 channels and 1 histogram of the total spectrum.
 *
 * Usage:
 * .L /usr/opt/ddas/6.1-001/lib/libddaschannel.so     //load ddas library
	.x ddas2root_run0079_0091.C
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

void run0079_0091_ddas2root()
{
	// For Run 1 241Am, 237Np lifetime root to cal
	// For Run 3 241Am, 237Np lifetime root to cal
	//step 0: cd /user/pxct/readout/rootfile/ or use godata
	//step 1: . /usr/opt/root/root-6.24.06/bin/thisroot.sh or use ./r.sh
	//step 2: Ctrl+R then type .L /usr/opt/ddas/6.1-001/lib/libddaschannel.so     //load ddas library
	//step 3: .x run0079_0091_ddas2root.C
	// https://docs.nscl.msu.edu/daq/newsite/ddas-1.1/The.html
	// https://docs.nscl.msu.edu/daq/newsite/ddas-1.1/ddaschannel_8h_source.html

	long totalentries = 0, i = 0;
	int Nchannels;
	char rawrootname[300];
	char calrootname[300];
	char txtname[300];
	char inpathname[300], outpathname[300];
	char filename[300];
	int Which_data = 1; // modify this number to select different data sets; 1 means MSD12+MSD26; 2 means MSD26 with a 2mm collimator
	double timestampshift = 0;

	sprintf(inpathname, "%s", "/mnt/daqtesting/pxct/stagearea/");
	//sprintf(inpathname, "%s", "/user/pxct/readout/rootfile/");
	
	//sprintf(outpathname, "%s", "/mnt/daqtesting/pxct/stagearea/");
	sprintf(outpathname, "%s", "/user/pxct/readout/rootfile/");

	if (Which_data == 1)
	{
		sprintf(filename, "%s", "run0079_0091_LEGe_MSD_241Am_inChamber_window1.5us_CFDdelay_adjusted");
	}
	if (Which_data == 2)
	{
		sprintf(filename, "%s", "run0092_0095_0100_0109_LEGe_MSD26_241Am_window1.5us");
	}
	if (Which_data == 3)
	{
		sprintf(filename, "%s", "run0330_0332_LEGe_MSD_241Am_Z7117_ChamberCenter_window1.5us_TrigRise0.064_0.016_0.016us_TrigGap0.952_1.000_1.000us_Th350_2700_1000_CFDDelay0.304us_Scale7");
	}
	// name both input and output root files accordingly
	sprintf(rawrootname, "%s%s%s", inpathname, filename, ".root");
	sprintf(txtname, "%s%s%s", inpathname, filename, ".dat");
	ofstream outfile(txtname, ios::out);

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

	for (int iFile = 0; iFile < 1; iFile++) // don't change
	{
		sprintf(calrootname, "%s%s%s%04d%s", outpathname, filename, "_", iFile, "_cal.root");
		TFile* fout = new TFile(calrootname, "RECREATE");
		cout << "  output file: " << calrootname << endl;
		TTree* tree = new TTree("tree", "tree");

		// If you have closed the channels in CSRA, please comment out the corresponding branch lines. No other changes are needed.
		tree->Branch("lege_e", &lege_e, "lege_e/D");
		tree->Branch("lege_t", &lege_t, "lege_t/D");
		tree->Branch("msd12_e", &msd12_e, "msd12_e/D");
		tree->Branch("msd12_t", &msd12_t, "msd12_t/D");
		tree->Branch("msd26_e", &msd26_e, "msd26_e/D");
		tree->Branch("msd26_t", &msd26_t, "msd26_t/D");
		tree->Branch("msdtotal_e", &msdtotal_e, "msdtotal_e/D");

		TH1D* hlege_e;
		TH1D* hmsd12_e;
		TH1D* hmsd26_e;
		TH1D* hmsd12_e_05keVbin;
		TH1D* hmsd26_e_05keVbin;
		TH1D* hmsdtotal_e;
		TH1D* htiming_lege_msd12;
		TH1D* htiming_lege_msd26;
		TH1D* htiming_msd12_msd26;

		double slope_msd12_e, slope_msd26_e;
		if (Which_data == 1)
		{
			slope_msd12_e = 0.1260823653643; // based on run0079-0091 241Am Ea_loss=1791 keV
			slope_msd26_e = 0.1172665825814; // based on run0079-0091 241Am Ea_loss=3626 keV
		}
		if (Which_data == 3)
		{
			slope_msd12_e = 0.1248083623693; // based on run0330-0332 241Am Ea_loss=1791 keV
			slope_msd26_e = 0.1180530685333; // based on run0330-0332 241Am Ea_loss=3626 keV
		}

		double binmax_msd12 = slope_msd12_e * 60000; // based on run0330-0332 241Am Ea_loss=1791 keV
		double binmax_msd26 = slope_msd26_e * 60000; // based on run0330-0332 241Am Ea_loss=3626 keV
		hlege_e = new TH1D("hlege_e", "hlege_e", 60000, -0.0540888491631, 435.8137733893930); // 60000 channels, 7.3 eV per channel run0060
		hmsd12_e = new TH1D("hmsd12_e", "hmsd12_e", 60000, 0, binmax_msd12); // 60000 channels, 126 eV per channel
		hmsd26_e = new TH1D("hmsd26_e", "hmsd26_e", 60000, 0, binmax_msd26); // 60000 channels, 117 eV per channel
		hmsd12_e_05keVbin = new TH1D("hmsd12_e_05keVbin", "hmsd12_e_05keVbin", 14000, 0, 7000); // 7000 channels, 0.5 keV per channel
		hmsd26_e_05keVbin = new TH1D("hmsd26_e_05keVbin", "hmsd26_e_05keVbin", 14000, 0, 7000); // 7000 channels, 0.5 keV per channel
		hmsdtotal_e = new TH1D("hmsdtotal_e", "hmsdtotal_e", 14000, 0, 7000); // 7000 channels, 0.5 keV per channel
		htiming_lege_msd12 = new TH1D("htiming_lege_msd12", "htiming_lege_msd12", 3000, -1500, 1500); // 3000 channels, 1 ns per channel
		htiming_lege_msd26 = new TH1D("htiming_lege_msd26", "htiming_lege_msd26", 3000, -1500, 1500); // 3000 channels, 1 ns per channel
		htiming_msd12_msd26 = new TH1D("htiming_msd12_msd26", "htiming_msd12_msd26", 3000, -1500, 1500); // 3000 channels, 1 ns per channel

		//for (i = iFile * 1e6; i < (iFile+1) * 1e6; i++)
		for (i = 0; i < totalentries; i++)
		{
			lege_e = 0; lege_t = 0; north_e = 0; north_t = 0; south_e = 0; south_t = 0; msd12_e = 0; msd12_t = 0; 	msd26_e = 0; msd26_t = 0;
			//lege_e_low = 0; lege_t_low = 0; north_e_low = 0; north_t_low = 0; south_e_low = 0; south_t_low = 0;
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
				if (pChan[j]->GetChannelID() == 6)
				{
					//msd12_e = pChan[j]->GetEnergy();
					msd12_e = pChan[j]->GetEnergy() * slope_msd12_e + 0; // based on 241Am Ea_loss=1791 keV
					hmsd12_e->Fill(msd12_e + gRandom->Uniform(-0.25, 0.25));
					msd12_t = pChan[j]->GetTime() + timestampshift;
				}
				if (pChan[j]->GetChannelID() == 8)
				{
					//msd26_e = pChan[j]->GetEnergy();
					//msd26_e = pChan[j]->GetEnergy() * 0.1152143446474 + 22.352486402900; // based on 148Gd and 241Am two peaks
					//msd26_e = pChan[j]->GetEnergy() * 0.1155273461893 - 0.001613374918; // based on 241Am two peaks Run 100
					msd26_e = pChan[j]->GetEnergy() * slope_msd26_e + 0; // based on 241Am Ea_loss=3626 keV
					hmsd26_e->Fill(msd26_e + gRandom->Uniform(-0.25, 0.25));
					msd26_t = pChan[j]->GetTime() + timestampshift;
				}

				if (msd12_e > 300 && msd26_e > 300)
				{
					hmsd12_e_05keVbin->Fill(msd12_e + gRandom->Uniform(-0.6, 0.6));
					hmsd26_e_05keVbin->Fill(msd26_e + gRandom->Uniform(-0.6, 0.6));
					msdtotal_e = msd12_e + msd26_e;
					hmsdtotal_e->Fill(msdtotal_e + gRandom->Uniform(-0.4, 0.4));
				}

				//if lege_e > 59.0 && lege_e < 60.1 && msd26_e > 5470 && msd26_e < 5510)
				//{
				//	htiming_lege_msd26->Fill(lege_t - msd26_t);
				//}

				if (lege_e > 59.0 && lege_e < 60.1 && msd12_e > 1700 && msd12_e < 2200 && msd26_e > 3300 && msd26_e < 3800 && msd12_e + msd26_e > 5421 - 60 && msd12_e + msd26_e < 5421 + 60)
				{
					htiming_lege_msd12->Fill(lege_t - msd12_t);
					htiming_lege_msd26->Fill(lege_t - msd26_t);
				}

				if (msd12_e > 1700 && msd12_e < 2200 && msd26_e > 3300 && msd26_e < 3800 && msd12_e + msd26_e > 5421 - 60 && msd12_e + msd26_e < 5421 + 60)
				{
					htiming_msd12_msd26->Fill(msd12_t - msd26_t);
				}

			} // (int j = 0; j < Nchannels; j++)
			tree->Fill();
			//if ( i>= totalentries) break;
			if (i % 10000000 == 0 && i != 0) cout << "entry: " << i << endl;
		} // for (i = 0; i < totalentries; i++)
		//outfile << "File:	" << iFile << "	Last entry:	" << i << "	Max:	" << hmsd26_e->GetBinCenter(hmsd26_e->GetMaximumBin()) << endl;
		fout->Write();
		fout->Close();
	} // for (long iFile = 0; iFile < 100; iFile++)
	pFile->Close();
}
