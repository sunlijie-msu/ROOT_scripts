/*
 * DDAStoROOR.cxx
 *
 * CALIBRATION FOR PXCT
 *
 * Created : May , 2017.   Author: Moshe Friedman.  Contact: Friedmam@nscl.msu.edu
 * Update: May 5 , 2018.   Author: David Perez Loureiro.  Contact: dperezlo@utk.edu
 * Last update: Jan 15, 2024. Author: Lijie Sun. Contact: sunli@frib.msu.edu
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

void ddas2root()
{
	//step 0: cd /user/pxct/readout/rootfile/ or use godata
	//step 1: . /usr/opt/root/root-6.24.06/bin/thisroot.sh or use ./r.sh
	//step 2: Ctrl+R then type .L /usr/opt/ddas/6.1-001/lib/libddaschannel.so     //load ddas library
	//step 3: .x ddas2root.C
	// https://docs.nscl.msu.edu/daq/newsite/ddas-1.1/The.html
	// https://docs.nscl.msu.edu/daq/newsite/ddas-1.1/ddaschannel_8h_source.html

	long totalentries, i;
	int Nchannels;
	char rawrootname[300];
	char calrootname[300];
	char inpathname[300], outpathname[300];
	char filename[300];
	int runnumber = 202; // modify
	double timestampshift = 0;

	//sprintf(inpathname, "%s", "/mnt/daqtesting/pxct/stagearea/");
	sprintf(inpathname, "%s", "/user/pxct/readout/rootfile/");

	//sprintf(outpathname, "%s", "/mnt/daqtesting/pxct/stagearea/");
	sprintf(outpathname, "%s", "/user/pxct/readout/rootfile/");

	if (runnumber == 7)
	{
		sprintf(filename, "%s", "run0007_LEGe_152Eu_inChamber_vacuum_window1.5us_CFDdelay_adjusted_for_LEGe_efficiency_345min");
	}
	if (runnumber == 9)
	{
		sprintf(filename, "%s", "run0009_LEGe_152Eu_inChamber_vacuum_window1.5us_CFDdelay_adjusted_for_LEGe_efficiency_720min");
	}
	if (runnumber == 10)
	{
		sprintf(filename, "%s", "run0010_LEGe_152Eu_inChamber_vacuum_window0.5us_CFDdelay_adjusted_for_LEGe_efficiency_440min");
	}
	if (runnumber == 16)
	{
		sprintf(filename, "%s", "run0016_LEGe_152Eu_inChamber_vacuum_window0.5us_CFDdelay_adjusted_for_LEGe_efficiency_1000min");
	}
	if (runnumber == 11)	sprintf(rawrootname, "%s%s", inpathname, "run0011-00_LEGe_room_background_24h_sources_in_rm1035.root");
	if (runnumber == 13)	sprintf(rawrootname, "%s%s", inpathname, "run0013-00_LEGe_room_background_24h_no_sources_in_rm1035.root");
	if (runnumber == 14)	sprintf(rawrootname, "%s%s", inpathname, "run0014-00_LEGe_55Fe_62min.root");
	if (runnumber == 6)	sprintf(rawrootname, "%s%s", inpathname, "run0006-00_LEGe_152Eu_88min.root");
	if (runnumber == 36)	sprintf(rawrootname, "%s%s", inpathname, "run0036-00_XtRa_152Eu_14cmAir_120min_window1us_baseline500us.root");
	if (runnumber == 26)	sprintf(rawrootname, "%s%s", inpathname, "run0026-00_LEGe_MSD26_241Am_nocollimator_90min_window1us.root");
	if (runnumber == 37)	sprintf(rawrootname, "%s%s", inpathname, "run0037-00_LEGe_MSD26_241Am_2mmcollimator_180min_window1us.root");
	if (runnumber == 38)	sprintf(rawrootname, "%s%s", inpathname, "run0038-00_LEGe_MSD26_241Am_2mmcollimator_40min_window1us_CFDdelay0.5us.root");

	if (runnumber == 11)	sprintf(calrootname, "%s%s", inpathname, "run0011-00_LEGe_room_background_24h_sources_in_rm1035_cal.root");
	if (runnumber == 13)	sprintf(calrootname, "%s%s", inpathname, "run0013-00_LEGe_room_background_24h_no_sources_in_rm1035_cal.root");
	if (runnumber == 14)	sprintf(calrootname, "%s%s", inpathname, "run0014-00_LEGe_55Fe_62min_cal.root");
	if (runnumber == 6)	sprintf(calrootname, "%s%s", inpathname, "run0006-00_LEGe_152Eu_88min_cal.root");
	if (runnumber == 36)	sprintf(calrootname, "%s%s", inpathname, "run0036-00_XtRa_152Eu_14cmAir_120min_window1us_baseline500us_cal.root");
	if (runnumber == 26)	sprintf(calrootname, "%s%s", inpathname, "run0026-00_LEGe_MSD26_241Am_nocollimator_90min_window1us_cal.root");
	if (runnumber == 37)	sprintf(calrootname, "%s%s", inpathname, "run0037-00_LEGe_MSD26_241Am_2mmcollimator_180min_window1us_cal.root");
	if (runnumber == 38)	sprintf(calrootname, "%s%s", inpathname, "run0038-00_LEGe_MSD26_241Am_2mmcollimator_40min_window1us_CFDdelay0.5us_cal.root");
	if (runnumber == 39)
	{
		sprintf(filename, "%s", "run0039-00_LEGe_MSD26_241Am_2mmcollimator_60min_window1us_CFDdelay0.3us");
		sprintf(rawrootname, "%s%s%s", inpathname, filename, ".root");
		sprintf(calrootname, "%s%s%s", inpathname, filename, "_cal.root");
	}
	if (runnumber == 40)
	{
		sprintf(filename, "%s", "run0040-00_LEGe_MSD26_241Am_2mmcollimator_20min_window1us_CFDdelay0.05us");
		sprintf(rawrootname, "%s%s%s", inpathname, filename, ".root");
		sprintf(calrootname, "%s%s%s", inpathname, filename, "_cal.root");
	}
	if (runnumber == 41)
	{
		sprintf(filename, "%s", "run0041-00_LEGe_MSD26_241Am_2mmcollimator_480min_window1us_CFDdelay0.16us");
		sprintf(rawrootname, "%s%s%s", inpathname, filename, ".root");
		sprintf(calrootname, "%s%s%s", inpathname, filename, "_cal.root");
	}
	if (runnumber == 42)
	{
		sprintf(filename, "%s", "run0042-00_LEGe_MSD26_XtRa_241Am_2mmcollimator_152Eu_11min_window1.5us_CFDdelay0.16us");
		sprintf(rawrootname, "%s%s%s", inpathname, filename, ".root");
		sprintf(calrootname, "%s%s%s", inpathname, filename, "_cal.root");
	}
	if (runnumber == 43)
	{
		sprintf(filename, "%s", "run0043-00_LEGe_MSD26_XtRa_241Am_2mmcollimator_152Eu_160min_window1.5us_CFDdelay0.16us");
		sprintf(rawrootname, "%s%s%s", inpathname, filename, ".root");
		sprintf(calrootname, "%s%s%s", inpathname, filename, "_cal.root");
	}
	if (runnumber == 45)
	{
		sprintf(filename, "%s", "run0045-00_LEGe_MSD_XtRa_241Am_3mmcollimator_152Eu_120min_window1.5us");
		sprintf(rawrootname, "%s%s%s", inpathname, filename, ".root");
		sprintf(calrootname, "%s%s%s", inpathname, filename, "_cal.root");
	}
	if (runnumber == 46)
	{
		sprintf(filename, "%s", "run0046-00_LEGe_MSD_XtRa_241Am_3mmcollimator_152Eu_121min_window1.5us");
		sprintf(rawrootname, "%s%s%s", inpathname, filename, ".root");
		sprintf(calrootname, "%s%s%s", inpathname, filename, "_cal.root");
	}
	if (runnumber == 47)
	{
		sprintf(filename, "%s", "run0047-00_LEGe_MSD_XtRa_Pulser_20min_window1.5us");
		sprintf(rawrootname, "%s%s%s", inpathname, filename, ".root");
		sprintf(calrootname, "%s%s%s", inpathname, filename, "_cal.root");
	}
	if (runnumber == 54)
	{
		sprintf(filename, "%s", "run0054-00_MSD26_241Am_2mmcollimator_60min_window1.5us_positivebias");
		sprintf(rawrootname, "%s%s%s", inpathname, filename, ".root");
		sprintf(calrootname, "%s%s%s", inpathname, filename, "_cal.root");
	}
	if (runnumber == 55)
	{
		sprintf(filename, "%s", "run0055-00_MSD26_148Gd_2mmcollimator_80min_window1.5us");
		sprintf(rawrootname, "%s%s%s", inpathname, filename, ".root");
		sprintf(calrootname, "%s%s%s", inpathname, filename, "_cal.root");
	}
	if (runnumber == 56)
	{
		sprintf(filename, "%s", "run0056-00_MSD26_241Am_2mmcollimator_60min_window1.5us_negativebias");
		sprintf(rawrootname, "%s%s%s", inpathname, filename, ".root");
		sprintf(calrootname, "%s%s%s", inpathname, filename, "_cal.root");
	}
	if (runnumber == 57)
	{
		sprintf(filename, "%s", "run0057-00_MSD_148Gd_2mmcollimator_120min_window1.5us");
		sprintf(rawrootname, "%s%s%s", inpathname, filename, ".root");
		sprintf(calrootname, "%s%s%s", inpathname, filename, "_cal.root");
	}
	if (runnumber == 60)
	{
		sprintf(filename, "%s", "run0060_LEGe_XtRa_152Eu_inChamber_400min_efficiency");
		sprintf(rawrootname, "%s%s%s", inpathname, filename, ".root");
		sprintf(calrootname, "%s%s%s", inpathname, filename, "_cal.root");
	}
	if (runnumber == 63)
	{
		sprintf(filename, "%s", "run0063_MSD_LEGe_XtRa_241Am_inChamber_2mmcollimator_600min");
		sprintf(rawrootname, "%s%s%s", inpathname, filename, ".root");
		sprintf(calrootname, "%s%s%s", inpathname, filename, "_cal.root");
	}
	if (runnumber == 64)
	{
		sprintf(filename, "%s", "run0064_LEGe_MSD_XtRa_241Am_inChamber_2mmcollimator_152Eu_outChamber_357min");
		sprintf(rawrootname, "%s%s%s", inpathname, filename, ".root");
		sprintf(calrootname, "%s%s%s", inpathname, filename, "_cal.root");
	}
	if (runnumber == 65)
	{
		sprintf(filename, "%s", "run0065_LEGe_MSD_XtRa_241Am_inChamber_2mmcollimator_152Eu_outChamber_663min");
		sprintf(rawrootname, "%s%s%s", inpathname, filename, ".root");
		sprintf(calrootname, "%s%s%s", inpathname, filename, "_cal.root");
	}
	if (runnumber == 66)
	{
		sprintf(filename, "%s", "run0066_LEGe_MSD_XtRa_241Am_inChamber_2mmcollimator_152Eu_outChamber_440min");
		sprintf(rawrootname, "%s%s%s", inpathname, filename, ".root");
		sprintf(calrootname, "%s%s%s", inpathname, filename, "_cal.root");
	}
	if (runnumber == 74) // 67, 68, 69, 72, 73, 74
	{
		sprintf(filename, "%s", "run0074_LEGe_MSD_XtRa_241Am_inChamber_152Eu_outChamber_triggerfilter_adjusting_55min");
		sprintf(rawrootname, "%s%s%s", inpathname, filename, ".root");
		sprintf(calrootname, "%s%s%s", inpathname, filename, "_cal.root");
	}
	if (runnumber == 70)
	{
		sprintf(filename, "%s", "run0070_LEGe_MSD_XtRa_241Am_inChamber_2mmcollimator_152Eu_outChamber_CFDdelay_adjusted_1035min");
		sprintf(rawrootname, "%s%s%s", inpathname, filename, ".root");
		sprintf(calrootname, "%s%s%s", inpathname, filename, "_cal.root");
	}
	if (runnumber == 75)
	{
		sprintf(filename, "%s", "run0075_LEGe_MSD_XtRa_241Am_inChamber_nocollimator_152Eu_outChamber_CFDdelay_adjusted_1511min");
		sprintf(rawrootname, "%s%s%s", inpathname, filename, ".root");
		sprintf(calrootname, "%s%s%s", inpathname, filename, "_cal.root");
	}
	if (runnumber == 71)
	{
		sprintf(filename, "%s", "run0071_LEGe_MSD_XtRa_241Am_inChamber_2mmcollimator_152Eu_outChamber_CFDdelay_adjusted_1068min");
		sprintf(rawrootname, "%s%s%s", inpathname, filename, ".root");
		sprintf(calrootname, "%s%s%s", inpathname, filename, "_cal.root");
		// timestampshift = 6.21e13;
	}
	if (runnumber == 76)
	{
		sprintf(filename, "%s", "run0076_LEGe_MSD_241Am_inChamber_window1us_742min");
		sprintf(rawrootname, "%s%s%s", inpathname, filename, ".root");
		sprintf(calrootname, "%s%s%s", inpathname, filename, "_cal.root");
	}
	if (runnumber == 77)
	{
		sprintf(filename, "%s", "run0077_LEGe_MSD_241Am_inChamber_window1us_CFDdelay_adjusting_45min"); // name both input and output root files accordingly
		sprintf(rawrootname, "%s%s%s", inpathname, filename, ".root");
		sprintf(calrootname, "%s%s%s", inpathname, filename, "_cal.root");
	}
	if (runnumber == 78)
	{
		sprintf(filename, "%s", "run0078_LEGe_MSD_241Am_inChamber_window1us_CFDdelay_adjusting_66min"); // name both input and output root files accordingly
		sprintf(rawrootname, "%s%s%s", inpathname, filename, ".root");
		sprintf(calrootname, "%s%s%s", inpathname, filename, "_cal.root");
	}
	if (runnumber == 79)
	{
		sprintf(filename, "%s", "run0079_LEGe_MSD_241Am_inChamber_window1.5us_CFDdelay_adjusted_1560min");
	}
	if (runnumber == 80)
	{
		sprintf(filename, "%s", "run0080_LEGe_MSD_241Am_inChamber_window1.5us_CFDdelay_adjusted_2280min");
	}
	if (runnumber == 81)
	{
		sprintf(filename, "%s", "run0081_LEGe_MSD_241Am_inChamber_window1.5us_CFDdelay_adjusted_2460min");
	}
	if (runnumber == 82)
	{
		sprintf(filename, "%s", "run0082_LEGe_MSD_241Am_inChamber_window1.5us_CFDdelay_adjusted_3722min");
	}
	if (runnumber == 83)
	{
		sprintf(filename, "%s", "run0083_LEGe_MSD_241Am_inChamber_window1.5us_CFDdelay_adjusted_6242min");
	}
	if (runnumber == 84)
	{
		sprintf(filename, "%s", "run0084_LEGe_MSD_241Am_inChamber_window1.5us_CFDdelay_adjusted_4556min");
	}
	if (runnumber == 85)
	{
		sprintf(filename, "%s", "run0085_LEGe_MSD_241Am_inChamber_window1.5us_CFDdelay_adjusted_6336min");
	}
	if (runnumber == 86)
	{
		sprintf(filename, "%s", "run0086_LEGe_MSD_241Am_inChamber_window1.5us_CFDdelay_adjusted_2552min");
	}
	if (runnumber == 87)
	{
		sprintf(filename, "%s", "run0087_LEGe_MSD_241Am_inChamber_window1.5us_CFDdelay_adjusted_1517min");
	}
	if (runnumber == 88)
	{
		sprintf(filename, "%s", "run0088_LEGe_MSD_241Am_inChamber_window1.5us_CFDdelay_adjusted_7771min");
	}
	if (runnumber == 89)
	{
		sprintf(filename, "%s", "run0089_LEGe_MSD_241Am_inChamber_window1.5us_CFDdelay_adjusted_11660min");
	}
	if (runnumber == 90)
	{
		sprintf(filename, "%s", "run0090_LEGe_MSD_241Am_inChamber_window1.5us_CFDdelay_adjusted_4515min");
	}
	if (runnumber == 91)
	{
		sprintf(filename, "%s", "run0091_LEGe_MSD_241Am_inChamber_window1.5us_CFDdelay_adjusted_1780min");
	}
	if (runnumber == 92)
	{
		sprintf(filename, "%s", "run0092_LEGe_MSD26_241Am_inChamber_2mmCollimator_window1.5us_CFDdelay_adjusted_901min");
	}
	if (runnumber == 93)
	{
		sprintf(filename, "%s", "run0093_LEGe_MSD26_241Am_inChamber_2mmCollimator_window1.5us_CFDdelay_adjusted_5000min");
	}
	if (runnumber == 94)
	{
		sprintf(filename, "%s", "run0094_LEGe_MSD26_241Am_inChamber_2mmCollimator_window1.5us_CFDdelay_adjusted_9528min");
	}
	if (runnumber == 95)
	{
		sprintf(filename, "%s", "run0095_LEGe_MSD26_241Am_inChamber_2mmCollimator_window1.5us_CFDdelay_adjusted_4099min");
	}
	if (runnumber == 96)
	{
		sprintf(filename, "%s", "run0096_LEGe_241Am_inChamber_window1.5us_CFDdelay_adjusted_60min");
	}
	if (runnumber == 97)
	{
		sprintf(filename, "%s", "run0097_LEGe_XtRa_152Eu_inChamber_window1.5us_CFDdelay_adjusted_180min");
	}
	if (runnumber == 98)
	{
		sprintf(filename, "%s", "run0098_LEGe_XtRa_152Eu_inChamber_window1.5us_CFDdelay_adjusted_750min");
	}
	if (runnumber == 99)
	{
		sprintf(filename, "%s", "run0099_LEGe_XtRa_152Eu_inChamber_flanges_swapped_window1.5us_CFDdelay_adjusted_420min");
	}
	if (runnumber == 100)
	{
		sprintf(filename, "%s", "run0100_LEGe_MSD26_241Am_inChamber_2mmCollimator_window1.5us_CFDdelay_adjusted_997min");
	}
	if (runnumber == 101)
	{
		sprintf(filename, "%s", "run0101_LEGe_MSD26_241Am_inChamber_2mmCollimator_window1.5us_CFDdelay_adjusted_9000min");
	}
	if (runnumber == 102)
	{
		sprintf(filename, "%s", "run0102_LEGe_MSD26_241Am_inChamber_2mmCollimator_window1.5us_CFDdelay_adjusted_360min");
	}
	if (runnumber == 103)
	{
		sprintf(filename, "%s", "run0103_LEGe_MSD26_241Am_inChamber_2mmCollimator_window1.5us_CFDdelay_adjusted_9000min");
	}
	if (runnumber == 200)
	{
		sprintf(filename, "%s", "run0200_LEGe_MSD26_137Cs_M4038_inChamber_vacuum_North_16mm_South_12mm_away_window1.5us_CFDdelay_adjusted");
	}
	if (runnumber == 201)
	{
		sprintf(filename, "%s", "run0201_LEGe_MSD26_137Cs_M4038_inChamber_vacuum_North_16mm_South_12mm_away_window1.5us_CFDdelay_adjusted");
	}
	if (runnumber == 202)
	{
		sprintf(filename, "%s", "run0202_LEGe_MSD26_137Cs_M4038_inChamber_vacuum_North_16mm_South_12mm_away_window1.5us_CFDdelay_adjusted");
	}
	// name both input and output root files accordingly
	sprintf(rawrootname, "%s%s%s", inpathname, filename, ".root");
	sprintf(calrootname, "%s%s%s", outpathname, filename, "_cal.root");

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

	TFile* fout = new TFile(calrootname, "RECREATE");
	cout << "  output file: " << calrootname << endl;
	TTree* tree = new TTree("tree", "tree");

	// If you have closed the channels in CSRA, please comment out the corresponding branch lines. No other changes are needed.
	tree->Branch("lege_e", &lege_e, "lege_e/D");
	tree->Branch("lege_t", &lege_t, "lege_t/D");
	tree->Branch("north_e", &north_e, "north_e/D");
	tree->Branch("north_t", &north_t, "north_t/D");
	tree->Branch("south_e", &south_e, "south_e/D");
	tree->Branch("south_t", &south_t, "south_t/D");
 	//tree->Branch("msd12_e", &msd12_e, "msd12_e/D");
 	//tree->Branch("msd12_t", &msd12_t, "msd12_t/D");
 	//tree->Branch("msd26_e", &msd26_e, "msd26_e/D");
 	//tree->Branch("msd26_t", &msd26_t, "msd26_t/D");
 	//tree->Branch("msdtotal_e", &msdtotal_e, "msdtotal_e/D");

// 	tree->Branch("lege_e_low", &lege_e_low, "lege_e_low/D");
// 	tree->Branch("lege_t_low", &lege_t_low, "lege_t_low/D");
// 	tree->Branch("north_e_low", &north_e_low, "north_e_low/D");
// 	tree->Branch("north_t_low", &north_t_low, "north_t_low/D");
// 	tree->Branch("south_e_low", &south_e_low, "south_e_low/D");
// 	tree->Branch("south_t_low", &south_t_low, "south_t_low/D");

	TH1D* hlege_e;
	TH1D* hnorth_e;
	TH1D* hsouth_e;
	TH1D* hmsd12_e;
	TH1D* hmsd26_e;
	TH1D* hmsd12_e1keVbin;
	TH1D* hmsd26_e1keVbin;
	TH1D* hmsdtotal_e;
	TH1D* hlege_e_low;
	TH1D* hnorth_e_low;
	TH1D* hsouth_e_low;
	TH1D* htiming_lege_msd12;
	TH1D* htiming_lege_msd26;

	// 	if (runnumber == 0)	hlege_e = new TH1D("hlege_e", "hlege_e", 60000, -0.0270924805469, 436.3668916340480); // 60000 channels, 7.3 eV per channel
	// 	if (runnumber == 1)	hlege_e = new TH1D("hlege_e", "hlege_e", 60000, -0.0584419300571, 435.7108166612020); // 60000 channels, 7.3 eV per channel
	// 	if (runnumber == 55)	hlege_e = new TH1D("hlege_e", "hlege_e", 60000, 0.0080746242334, 437.1574593691120); // 60000 channels, 7.3 eV per channel
	// 	if (runnumber == 152)	hlege_e = new TH1D("hlege_e", "hlege_e", 60000, -0.0001420032512, 436.2768890605690); // 60000 channels, 7.3 eV per channel
	//	if (runnumber == 241)	hlege_e = new TH1D("hlege_e", "hlege_e", 60000, 0.0541091878755, 435.4667414524260); // 60000 channels, 7.3 eV per channel
	//hlege_e = new TH1D("hlege_e", "hlege_e", 60000, 0.0541091878755, 435.4667414524260); // 60000 channels, 7.3 eV per channel
	hlege_e = new TH1D("hlege_e", "hlege_e", 60000, -0.0540888491631, 435.8137733893930); // 60000 channels, 7.3 eV per channel run0060
	// hnorth_e = new TH1D("hnorth_e", "hnorth_e", 60000, -0.1774488889544, 2286.04224413527); // 60000 channels, 38 eV per channel from Excel calibration
	hnorth_e = new TH1D("hnorth_e", "hnorth_e", 60000, -0.1943988436200, 2285.29642151811); // 60000 channels, 38 eV per channel from Excel calibration 1/11/2024 run 97-98
	// hsouth_e = new TH1D("hsouth_e", "hsouth_e", 60000, -0.1622624327090, 2278.22727449908); // 60000 channels, 38 eV per channel from Excel calibration
	hsouth_e = new TH1D("hsouth_e", "hsouth_e", 60000, -0.1195774566938, 2274.16156253889); // 60000 channels, 38 eV per channel from Excel calibration 1/11/2024 run 97-98
	hmsd12_e = new TH1D("hmsd12_e", "hmsd12_e", 60000, 0, 7941.93363844394); // 60000 channels, 512 eV per channel
	hmsd26_e = new TH1D("hmsd26_e", "hmsd26_e", 60000, -0.0016133749182, 6931.63915798418); // 60000 channels, 115 eV per channel from Excel calibration 1/18/2024 run 100
	//hmsdtotal_e = new TH1D("hmsdtotal_e", "hmsdtotal_e", 6900, 0, 6900); // 6900 channels, 1000 eV per channel
	//hmsd12_e1keVbin = new TH1D("hmsd12_e1keVbin", "hmsd12_e1keVbin", 7942, 0, 7942); // 7942 channels, 1 keV per channel
	//hmsd26_e1keVbin = new TH1D("hmsd26_e1keVbin", "hmsd26_e1keVbin", 6935, 0, 6935); // 6935 channels, 1 keV per channel

	//hlege_e_low = new TH1D("hlege_e_low", "hlege_e_low", 60000, 0, 60000); // 60000 channels
	//hnorth_e_low = new TH1D("hnorth_e_low", "hnorth_e_low", 60000, 0, 60000); // 60000 channels
	//hsouth_e_low = new TH1D("hsouth_e_low", "hsouth_e_low", 60000, 0, 60000); // 60000 channels
	//htiming_lege_msd12 = new TH1D("htiming_lege_msd12", "htiming_lege_msd12", 3000, -1500, 1500); // 3000 channels, 1 ns per channel
	//htiming_lege_msd26 = new TH1D("htiming_lege_msd26", "htiming_lege_msd26", 3000, -1500, 1500); // 3000 channels, 1 ns per channel

	for (i = 0; i < totalentries; i++)
	{
		lege_e = 0; lege_t = 0; north_e = 0; north_t = 0; south_e = 0; south_t = 0; msd12_e = 0; msd12_t = 0; 	msd26_e = 0; msd26_t = 0;
		lege_e_low = 0; lege_t_low = 0; north_e_low = 0; north_t_low = 0; south_e_low = 0; south_t_low = 0;
		pTree->GetEntry(i);
		pChan = pEvent->GetData();
		Nchannels = (int)pChan.size();  //get active number of channels
		// cout<<"  number of active channels = "<<Nchannels<<endl;
		for (int j = 0; j < Nchannels; j++)
		{
			//if (pChan[j]->GetChannelID() == 0)
			//{
			//	lege_e_low = pChan[j]->GetEnergy();
			//	hlege_e_low->Fill(lege_e_low);
			//	lege_t_low = pChan[j]->GetTime() + timestampshift;
			//}
			if (pChan[j]->GetChannelID() == 1)
			{
				//lege_e = pChan[j]->GetEnergy();
				//lege_e = pChan[j]->GetEnergy() * 0.0072712838511 - 0.0001420032512;
				lege_e = pChan[j]->GetEnergy() * 0.0072644643706 - 0.054088849163; // based on run0060
				hlege_e->Fill(lege_e);
				lege_t = pChan[j]->GetTime() + timestampshift;
			}
			//if (pChan[j]->GetChannelID() == 2)
			//{
			//	north_e_low = pChan[j]->GetEnergy();
			//	hnorth_e_low->Fill(north_e_low);
			//	north_t_low = pChan[j]->GetTime() + timestampshift;
			//}
			if (pChan[j]->GetChannelID() == 3)
			{
				//north_e = pChan[j]->GetEnergy();
				//north_e = pChan[j]->GetEnergy() * 0.0381036615504 - 0.1774488889544; // based on run0060
				north_e = pChan[j]->GetEnergy() * 0.0380915136727 - 0.1943988436200; // based on run0098
				hnorth_e->Fill(north_e);
				north_t = pChan[j]->GetTime() + timestampshift;
			}
			//if (pChan[j]->GetChannelID() == 4)
			//{
			//	south_e_low = pChan[j]->GetEnergy();
			//	hsouth_e_low->Fill(south_e_low);
			//	south_t_low = pChan[j]->GetTime() + timestampshift;
			//}
			if (pChan[j]->GetChannelID() == 5)
			{
				//south_e = pChan[j]->GetEnergy();
				//south_e = pChan[j]->GetEnergy() * 0.0379731589489 - 0.1622624327090; // based on run0060
				south_e = pChan[j]->GetEnergy() * 0.0379046856666 - 0.1195774566938; // based on run0098
				hsouth_e->Fill(south_e);
				south_t = pChan[j]->GetTime() + timestampshift;
			}
			//if (pChan[j]->GetChannelID() == 6)
			//{
			//	//msd12_e = pChan[j]->GetEnergy();
			//	msd12_e = pChan[j]->GetEnergy() * 0.1323655606407 + 0;
			//	hmsd12_e->Fill(msd12_e);
			//	hmsd12_e1keVbin->Fill(msd12_e + gRandom->Uniform(-0.5, 0.5));
			//	msd12_t = pChan[j]->GetTime() + timestampshift;
			//}
			//if (pChan[j]->GetChannelID() == 8)
			//{
			//	//msd26_e = pChan[j]->GetEnergy();
			//	//msd26_e = pChan[j]->GetEnergy() * 0.1152143446474 + 22.352486402900; // based on 148Gd and 241Am two peaks
			//	msd26_e = pChan[j]->GetEnergy() * 0.1155273461893 - 0.001613374918; // based on 241Am two peaks Run 100
			//	hmsd26_e->Fill(msd26_e);
			//	hmsd26_e1keVbin->Fill(msd26_e + gRandom->Uniform(-0.5, 0.5));
			//	msd26_t = pChan[j]->GetTime() + timestampshift;
			//}

			//if (msd12_e > 200 && msd26_e > 200)
			//{
			//	msdtotal_e = msd12_e + msd26_e;
			//	hmsdtotal_e->Fill(msdtotal_e);
			//}

			//if (runnumber >= 92 && lege_e > 59 && lege_e < 60.1 && msd26_e > 5470 && msd26_e < 5510)
			//{
			//	htiming_lege_msd26->Fill(lege_t - msd26_t);
			//}

			//if (runnumber < 92 && lege_e > 58.5 && lege_e < 60.5 && msd12_e>1700 && msd12_e < 2200 && msd26_e>3300 && msd26_e < 3800 && msd12_e + msd26_e > 5400 && msd12_e + msd26_e < 5560)
			//{
			//	htiming_lege_msd12->Fill(lege_t - msd12_t);
			//	htiming_lege_msd26->Fill(lege_t - msd26_t);
			//}


		} // (int j = 0; j < Nchannels; j++)
		tree->Fill();
		if (i % 1000000 == 0 && i != 0) cout << "entry: " << i << endl;
	} // for (int i = 0; i < totalentries; i++)
	cout << "  Last entry: " << i << endl;
	fout->Write();
	fout->Close();
	pFile->Close();
}
