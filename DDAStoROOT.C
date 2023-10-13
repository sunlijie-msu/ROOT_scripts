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

void DDAStoROOT()
{
	//root
	//.L /usr/opt/ddas/6.1-001/lib/libddaschannel.so     //load ddas library
	//.x DDAStoROOT.C
	// https://docs.nscl.msu.edu/daq/newsite/ddas-1.1/The.html
	// https://docs.nscl.msu.edu/daq/newsite/ddas-1.1/ddaschannel_8h_source.html

	long totalentries;
	int Nchannels;
	char rawrootname[150];
	char calrootname[150];
	char pathname[150];
	char filename[150];
	int source = 152; // modify

	sprintf(pathname, "%s", "/user/pxct/readout/crate_1/");

	if (source == 0)	sprintf(rawrootname, "%s%s", pathname, "run0011-00_LEGe_room_background_24h_sources_in_rm1035.root");
	if (source == 1)	sprintf(rawrootname, "%s%s", pathname, "run0013-00_LEGe_room_background_24h_no_sources_in_rm1035.root");
	if(source == 55)	sprintf(rawrootname, "%s%s", pathname, "run0014-00_LEGe_55Fe_62min.root");
	if (source == 152)	sprintf(rawrootname, "%s%s", pathname, "run0006-00_LEGe_152Eu_88min.root");
	if (source == 241)	sprintf(rawrootname, "%s%s", pathname, "run0008-00_LEGe_241Am_66min.root");

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

	if (source == 0)	sprintf(calrootname, "%s%s", pathname, "run0011-00_LEGe_room_background_24h_sources_in_rm1035_cal.root");
	if (source == 1)	sprintf(calrootname, "%s%s", pathname, "run0013-00_LEGe_room_background_24h_no_sources_in_rm1035_cal.root");
	if (source == 55)	sprintf(calrootname, "%s%s", pathname, "run0014-00_LEGe_55Fe_62min_cal.root");
	if (source == 152)	sprintf(calrootname, "%s%s", pathname, "run0006-00_LEGe_152Eu_88min_cal.root");
	if (source == 241)	sprintf(calrootname, "%s%s", pathname, "run0008-00_LEGe_241Am_66min_cal.root");

	TFile* fout = new TFile(calrootname, "RECREATE");
	TTree* tree = new TTree("tree", "tree");

	tree->Branch("lege_e", &lege_e, "lege_e/D");
	tree->Branch("lege_t", &lege_t, "lege_t/D");
	TH1D* hlege_e;
	if (source == 0)	hlege_e = new TH1D("hlege_e", "hlege_e", 60000, -0.0270924805469, 436.3668916340480); // 60000 channels, 7.3 eV per channel
	if (source == 1)	hlege_e = new TH1D("hlege_e", "hlege_e", 60000, -0.0584419300571, 435.7108166612020); // 60000 channels, 7.3 eV per channel
	if (source == 55)	hlege_e = new TH1D("hlege_e", "hlege_e", 60000, 0.0080746242334, 437.1574593691120); // 60000 channels, 7.3 eV per channel
	if (source == 152)	hlege_e = new TH1D("hlege_e", "hlege_e", 60000, -0.0001420032512, 436.2768890605690); // 60000 channels, 7.3 eV per channel
	if (source == 241)	hlege_e = new TH1D("hlege_e", "hlege_e", 60000, 0.0541091878755, 435.4667414524260); // 10000 channels, 7.3 eV per channel

	for (int i = 0; i < totalentries; i++)
	{
		lege_e = 0; lege_t = 0;
		pTree->GetEntry(i);
		pChan = pEvent->GetData();
		Nchannels = (int)pChan.size();  //get active number of channels
		// cout<<"  number of active channels = "<<Nchannels<<endl;
		for (int j = 0; j < Nchannels; j++)
		{
			if (pChan[j]->GetChannelID() == 13)
			{
				//lege_e = pChan[j]->GetEnergy();
				if (source == 0)      lege_e = pChan[j]->GetEnergy() * 0.0072732330686 - 0.0270924805469;
				if (source == 1)      lege_e = pChan[j]->GetEnergy() * 0.0072628209765 - 0.0584419300571;
				if (source == 55)	    lege_e = pChan[j]->GetEnergy() * 0.0072858230791 + 0.0080746242334;
				if (source == 152)  lege_e = pChan[j]->GetEnergy() * 0.0072712838511 - 0.0001420032512;
				if (source == 241)  lege_e = pChan[j]->GetEnergy() * 0.0072568772044 + 0.0541091878755;
				hlege_e->Fill(lege_e);
				lege_t = pChan[j]->GetTime();
			}
		}
		tree->Fill();
	}
	fout->Write();
	fout->Close();
	pFile->Close();
}
