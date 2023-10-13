/*
 * DDAStoROOR.cxx
 * 
 * CALIBRATION FOR PXCT
 *
 * Created : May , 2017.   Author: Moshe Friedman.  Contact: Friedmam@nscl.msu.edu
 * Update: May 5 , 2018.   Author: David Perez Loureiro.  Contact: dperezlo@utk.edu
 * Last update: Sep 12, 2023. Author: Lijie Sun. Contact: sunli@frib.msu.edu
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
	.x DDAStoROOT.C
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
 //#include <iomanip.h>
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

void ddas2root()
{
	//step 0: cd /user/pxct/readout/rootfile/
	//step 1: . /usr/opt/root/root-6.24.06/bin/thisroot.sh
	//step 2: .L /usr/opt/ddas/6.1-001/lib/libddaschannel.so     //load ddas library
	//step 3: .x ddas2root.C
	// https://docs.nscl.msu.edu/daq/newsite/ddas-1.1/The.html
	// https://docs.nscl.msu.edu/daq/newsite/ddas-1.1/ddaschannel_8h_source.html

	long totalentries;
	int Nchannels;
	char rawrootname[250];
	char calrootname[250];
	char pathname[250];
	char filename[250];
	int runnumber = 40; // modify

	sprintf(pathname, "%s", "/user/pxct/readout/rootfile/");

	if (runnumber == 11)	sprintf(rawrootname, "%s%s", pathname, "run0011-00_LEGe_room_background_24h_sources_in_rm1035.root");
	if (runnumber == 13)	sprintf(rawrootname, "%s%s", pathname, "run0013-00_LEGe_room_background_24h_no_sources_in_rm1035.root");
	if(runnumber == 14)	sprintf(rawrootname, "%s%s", pathname, "run0014-00_LEGe_55Fe_62min.root");
	if (runnumber == 6)	sprintf(rawrootname, "%s%s", pathname, "run0006-00_LEGe_152Eu_88min.root");
	if (runnumber == 36)	sprintf(rawrootname, "%s%s", pathname, "run0036-00_XtRa_152Eu_14cmAir_120min_window1us_baseline500us.root");
	if (runnumber == 26)	sprintf(rawrootname, "%s%s", pathname, "run0026-00_LEGe_MSD26_241Am_nocollimator_90min_window1us.root");
	if (runnumber == 37)	sprintf(rawrootname, "%s%s", pathname, "run0037-00_LEGe_MSD26_241Am_2mmcollimator_180min_window1us.root");
	if (runnumber == 38)	sprintf(rawrootname, "%s%s", pathname, "run0038-00_LEGe_MSD26_241Am_2mmcollimator_40min_window1us_CFDdelay0.5us.root");

	if (runnumber == 11)	sprintf(calrootname, "%s%s", pathname, "run0011-00_LEGe_room_background_24h_sources_in_rm1035_cal.root");
	if (runnumber == 13)	sprintf(calrootname, "%s%s", pathname, "run0013-00_LEGe_room_background_24h_no_sources_in_rm1035_cal.root");
	if (runnumber == 14)	sprintf(calrootname, "%s%s", pathname, "run0014-00_LEGe_55Fe_62min_cal.root");
	if (runnumber == 6)	sprintf(calrootname, "%s%s", pathname, "run0006-00_LEGe_152Eu_88min_cal.root");
	if (runnumber == 36)	sprintf(calrootname, "%s%s", pathname, "run0036-00_XtRa_152Eu_14cmAir_120min_window1us_baseline500us_cal.root");
	if (runnumber == 26)	sprintf(calrootname, "%s%s", pathname, "run0026-00_LEGe_MSD26_241Am_nocollimator_90min_window1us_cal.root");
	if (runnumber == 37)	sprintf(calrootname, "%s%s", pathname, "run0037-00_LEGe_MSD26_241Am_2mmcollimator_180min_window1us_cal.root");
	if (runnumber == 38)	sprintf(calrootname, "%s%s", pathname, "run0038-00_LEGe_MSD26_241Am_2mmcollimator_40min_window1us_CFDdelay0.5us_cal.root");
	if (runnumber == 39)
	{
		sprintf(filename, "%s", "run0039-00_LEGe_MSD26_241Am_2mmcollimator_60min_window1us_CFDdelay0.3us");
		sprintf(rawrootname, "%s%s%s", pathname, filename, ".root");
		sprintf(calrootname, "%s%s%s", pathname, filename, "_cal.root");
	}
	if (runnumber == 40)
	{
		sprintf(filename, "%s", "run0040-00_LEGe_MSD26_241Am_2mmcollimator_20min_window1us_CFDdelay0.05us");
		sprintf(rawrootname, "%s%s%s", pathname, filename, ".root");
		sprintf(calrootname, "%s%s%s", pathname, filename, "_cal.root");
	}

	TFile* pFile = new TFile(rawrootname);
	cout<<"  input file: "<<rawrootname<<endl;
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
	Double_t msd26_e;
	Double_t msd26_t;

	Double_t lege_e_low;
	Double_t lege_t_low;
	Double_t north_e_low;
	Double_t north_t_low;
	Double_t south_e_low;
	Double_t south_t_low;

	TFile* fout = new TFile(calrootname, "RECREATE");
	cout << "  output file: " << calrootname << endl;
	TTree* tree = new TTree("tree", "tree");

	tree->Branch("lege_e", &lege_e, "lege_e/D");
	tree->Branch("lege_t", &lege_t, "lege_t/D");
	tree->Branch("north_e", &north_e, "north_e/D");
	tree->Branch("north_t", &north_t, "north_t/D");
	tree->Branch("south_e", &south_e, "south_e/D");
	tree->Branch("south_t", &south_t, "south_t/D");
	tree->Branch("msd26_e", &msd26_e, "msd26_e/D");
	tree->Branch("msd26_t", &msd26_t, "msd26_t/D");

	tree->Branch("lege_e_low", &lege_e_low, "lege_e_low/D");
	tree->Branch("lege_t_low", &lege_t_low, "lege_t_low/D");
	tree->Branch("north_e_low", &north_e_low, "north_e_low/D");
	tree->Branch("north_t_low", &north_t_low, "north_t_low/D");
	tree->Branch("south_e_low", &south_e_low, "south_e_low/D");
	tree->Branch("south_t_low", &south_t_low, "south_t_low/D");

	TH1D* hlege_e;
	TH1D* hnorth_e;
	TH1D* hsouth_e;
	TH1D* hmsd26_e;
	TH1D* hlege_e_low;
	TH1D* hnorth_e_low;
	TH1D* hsouth_e_low;

// 	if (runnumber == 0)	hlege_e = new TH1D("hlege_e", "hlege_e", 60000, -0.0270924805469, 436.3668916340480); // 60000 channels, 7.3 eV per channel
// 	if (runnumber == 1)	hlege_e = new TH1D("hlege_e", "hlege_e", 60000, -0.0584419300571, 435.7108166612020); // 60000 channels, 7.3 eV per channel
// 	if (runnumber == 55)	hlege_e = new TH1D("hlege_e", "hlege_e", 60000, 0.0080746242334, 437.1574593691120); // 60000 channels, 7.3 eV per channel
// 	if (runnumber == 152)	hlege_e = new TH1D("hlege_e", "hlege_e", 60000, -0.0001420032512, 436.2768890605690); // 60000 channels, 7.3 eV per channel
//	if (runnumber == 241)	hlege_e = new TH1D("hlege_e", "hlege_e", 60000, 0.0541091878755, 435.4667414524260); // 60000 channels, 7.3 eV per channel
	hlege_e = new TH1D("hlege_e", "hlege_e", 60000, 0.0541091878755, 435.4667414524260); // 60000 channels, 7.3 eV per channel
	hnorth_e = new TH1D("hnorth_e", "hnorth_e", 60000, -0.1563415655380, 2286.23648597448); // 60000 channels, 38 eV per channel from Excel calibration
	hsouth_e = new TH1D("hsouth_e", "hsouth_e", 60000, -0.1317718501648, 2277.83634038208); // 60000 channels, 38 eV per channel from Excel calibration
	hmsd26_e = new TH1D("hmsd26_e", "hmsd26_e", 60000, -115.801142747884, 6983.41926747833); // 60000 channels, 119 eV per channel

	hlege_e_low = new TH1D("hlege_e_low", "hlege_e_low", 60000, 0, 60000); // 60000 channels
	hnorth_e_low = new TH1D("hnorth_e_low", "hnorth_e_low", 60000, 0, 60000); // 60000 channels
	hsouth_e_low = new TH1D("hsouth_e_low", "hsouth_e_low", 60000, 0, 60000); // 60000 channels

	for (int i = 0; i < totalentries; i++)
	{
		lege_e = 0; lege_t = 0; north_e = 0; north_t = 0; south_e = 0; south_t = 0; msd26_e = 0; msd26_t = 0;
		lege_e_low = 0; lege_t_low = 0; north_e_low = 0; north_t_low = 0; south_e_low = 0; south_t_low = 0;
		pTree->GetEntry(i);
		pChan = pEvent->GetData();
		Nchannels = (int)pChan.size();  //get active number of channels
		// cout<<"  number of active channels = "<<Nchannels<<endl;
		for (int j = 0; j < Nchannels; j++)
		{
			if (pChan[j]->GetChannelID() == 0)
			{
				lege_e_low = pChan[j]->GetEnergy();
				hlege_e_low->Fill(lege_e_low);
				lege_t_low = pChan[j]->GetTime();
			}
			if (pChan[j]->GetChannelID() == 1)
			{
				//lege_e = pChan[j]->GetEnergy();
				lege_e = pChan[j]->GetEnergy() * 0.0072712838511 - 0.0001420032512;
				hlege_e->Fill(lege_e);
				lege_t = pChan[j]->GetTime();
			}
			if (pChan[j]->GetChannelID() == 2)
			{
				north_e_low = pChan[j]->GetEnergy();
				hnorth_e_low->Fill(north_e_low);
				north_t_low = pChan[j]->GetTime();
			}
			if (pChan[j]->GetChannelID() == 3)
			{
				//north_e = pChan[j]->GetEnergy();
				north_e = pChan[j]->GetEnergy() * 0.0381065471257 - 0.1563415655380;
				hnorth_e->Fill(north_e);
				north_t = pChan[j]->GetTime();
			}
			if (pChan[j]->GetChannelID() == 4)
			{
				south_e_low = pChan[j]->GetEnergy();
				hsouth_e_low->Fill(south_e_low);
				south_t_low = pChan[j]->GetTime();
			}
			if (pChan[j]->GetChannelID() == 5)
			{
				//south_e = pChan[j]->GetEnergy();
				south_e = pChan[j]->GetEnergy() * 0.0379661352039 - 0.1317718501648;
				hsouth_e->Fill(south_e);
				south_t = pChan[j]->GetTime();
			}
			if (pChan[j]->GetChannelID() == 8)
			{
				//msd26_e = pChan[j]->GetEnergy();
				msd26_e = pChan[j]->GetEnergy() * 0.1183203401704 - 115.801142747884;
				hmsd26_e->Fill(msd26_e);
				msd26_t = pChan[j]->GetTime();
			}
		} // (int j = 0; j < Nchannels; j++)
		tree->Fill();
	} // for (int i = 0; i < totalentries; i++)
	fout->Write();
	fout->Close();
	pFile->Close();
}
