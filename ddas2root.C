/*
 * DDAStoROOR.cxx
 *
 * CALIBRATION FOR PXCT
 *
 * Created : May , 2017.   Author: Moshe Friedman.  Contact: Friedmam@nscl.msu.edu
 * Update: May 5 , 2018.   Author: David Perez Loureiro.  Contact: dperezlo@utk.edu
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

	time_t start, tim;
	struct tm* at;
	char now[80];
	float speed;


	long totalentries, i, i_outtree;
	int Nchannels;
	char rawrootname[400];
	char calrootname[400];
	char inpathname[300], outpathname[300];
	char filename[400];
	int runnumber;
	cout << "input run number: ";
	cin >> runnumber;
	double timestampshift = 0;

	sprintf(inpathname, "%s", "/mnt/daqtesting/pxct/stagearea/");
	//sprintf(inpathname, "%s", "/user/pxct/readout/rootfile/");

	//sprintf(outpathname, "%s", "/mnt/daqtesting/pxct/stagearea/");
	sprintf(outpathname, "%s", "/user/pxct/readout/rootfile/");

	if (runnumber == 7)	sprintf(filename, "%s", "run0007_LEGe_152Eu_inChamber_vacuum_window1.5us_CFDdelay_adjusted_for_LEGe_efficiency_345min");
	if (runnumber == 9)	sprintf(filename, "%s", "run0009_LEGe_152Eu_inChamber_vacuum_window1.5us_CFDdelay_adjusted_for_LEGe_efficiency_720min");
	if (runnumber == 10)	sprintf(filename, "%s", "run0010_LEGe_152Eu_inChamber_vacuum_window0.5us_CFDdelay_adjusted_for_LEGe_efficiency_440min");
	if (runnumber == 16)	sprintf(filename, "%s", "run0016_LEGe_152Eu_inChamber_vacuum_window0.5us_CFDdelay_adjusted_for_LEGe_efficiency_1000min");
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

	if (runnumber == 39)	sprintf(filename, "%s", "run0039-00_LEGe_MSD26_241Am_2mmcollimator_60min_window1us_CFDdelay0.3us");
	if (runnumber == 40)	sprintf(filename, "%s", "run0040-00_LEGe_MSD26_241Am_2mmcollimator_20min_window1us_CFDdelay0.05us");
	if (runnumber == 41)	sprintf(filename, "%s", "run0041-00_LEGe_MSD26_241Am_2mmcollimator_480min_window1us_CFDdelay0.16us");
	if (runnumber == 42)	sprintf(filename, "%s", "run0042-00_LEGe_MSD26_XtRa_241Am_2mmcollimator_152Eu_11min_window1.5us_CFDdelay0.16us");
	if (runnumber == 43)	sprintf(filename, "%s", "run0043-00_LEGe_MSD26_XtRa_241Am_2mmcollimator_152Eu_160min_window1.5us_CFDdelay0.16us");
	if (runnumber == 45)	sprintf(filename, "%s", "run0045-00_LEGe_MSD_XtRa_241Am_3mmcollimator_152Eu_120min_window1.5us");
	if (runnumber == 46)	sprintf(filename, "%s", "run0046-00_LEGe_MSD_XtRa_241Am_3mmcollimator_152Eu_121min_window1.5us");
	if (runnumber == 47)	sprintf(filename, "%s", "run0047-00_LEGe_MSD_XtRa_Pulser_20min_window1.5us");
	if (runnumber == 54)	sprintf(filename, "%s", "run0054-00_MSD26_241Am_2mmcollimator_60min_window1.5us_positivebias");
	if (runnumber == 55)	sprintf(filename, "%s", "run0055-00_MSD26_148Gd_2mmcollimator_80min_window1.5us");
	if (runnumber == 56)	sprintf(filename, "%s", "run0056-00_MSD26_241Am_2mmcollimator_60min_window1.5us_negativebias");
	if (runnumber == 57)	sprintf(filename, "%s", "run0057-00_MSD_148Gd_2mmcollimator_120min_window1.5us");
	if (runnumber == 60)	sprintf(filename, "%s", "run0060_LEGe_XtRa_152Eu_inChamber_400min_efficiency");
	if (runnumber == 63)	sprintf(filename, "%s", "run0063_MSD_LEGe_XtRa_241Am_inChamber_2mmcollimator_600min");
	if (runnumber == 64)	sprintf(filename, "%s", "run0064_LEGe_MSD_XtRa_241Am_inChamber_2mmcollimator_152Eu_outChamber_357min");
	if (runnumber == 65)	sprintf(filename, "%s", "run0065_LEGe_MSD_XtRa_241Am_inChamber_2mmcollimator_152Eu_outChamber_663min");
	if (runnumber == 66)	sprintf(filename, "%s", "run0066_LEGe_MSD_XtRa_241Am_inChamber_2mmcollimator_152Eu_outChamber_440min");
	if (runnumber == 74)	sprintf(filename, "%s", "run0074_LEGe_MSD_XtRa_241Am_inChamber_152Eu_outChamber_triggerfilter_adjusting_55min");
	if (runnumber == 70)	sprintf(filename, "%s", "run0070_LEGe_MSD_XtRa_241Am_inChamber_2mmcollimator_152Eu_outChamber_CFDdelay_adjusted_1035min");
	if (runnumber == 75)	sprintf(filename, "%s", "run0075_LEGe_MSD_XtRa_241Am_inChamber_nocollimator_152Eu_outChamber_CFDdelay_adjusted_1511min");
	if (runnumber == 71)	sprintf(filename, "%s", "run0071_LEGe_MSD_XtRa_241Am_inChamber_2mmcollimator_152Eu_outChamber_CFDdelay_adjusted_1068min");
	if (runnumber == 76)	sprintf(filename, "%s", "run0076_LEGe_MSD_241Am_inChamber_window1us_742min");
	if (runnumber == 77)	sprintf(filename, "%s", "run0077_LEGe_MSD_241Am_inChamber_window1us_CFDdelay_adjusting_45min");
	if (runnumber == 78)	sprintf(filename, "%s", "run0078_LEGe_MSD_241Am_inChamber_window1us_CFDdelay_adjusting_66min");
	if (runnumber == 79)	sprintf(filename, "%s", "run0079_LEGe_MSD_241Am_inChamber_window1.5us_CFDdelay_adjusted_1560min");
	if (runnumber == 80)	sprintf(filename, "%s", "run0080_LEGe_MSD_241Am_inChamber_window1.5us_CFDdelay_adjusted_2280min");
	if (runnumber == 81)	sprintf(filename, "%s", "run0081_LEGe_MSD_241Am_inChamber_window1.5us_CFDdelay_adjusted_2460min");
	if (runnumber == 82)	sprintf(filename, "%s", "run0082_LEGe_MSD_241Am_inChamber_window1.5us_CFDdelay_adjusted_3722min");
	if (runnumber == 83)	sprintf(filename, "%s", "run0083_LEGe_MSD_241Am_inChamber_window1.5us_CFDdelay_adjusted_6242min");
	if (runnumber == 84)	sprintf(filename, "%s", "run0084_LEGe_MSD_241Am_inChamber_window1.5us_CFDdelay_adjusted_4556min");
	if (runnumber == 85)	sprintf(filename, "%s", "run0085_LEGe_MSD_241Am_inChamber_window1.5us_CFDdelay_adjusted_6336min");
	if (runnumber == 86)	sprintf(filename, "%s", "run0086_LEGe_MSD_241Am_inChamber_window1.5us_CFDdelay_adjusted_2552min");
	if (runnumber == 87)	sprintf(filename, "%s", "run0087_LEGe_MSD_241Am_inChamber_window1.5us_CFDdelay_adjusted_1517min");
	if (runnumber == 88)	sprintf(filename, "%s", "run0088_LEGe_MSD_241Am_inChamber_window1.5us_CFDdelay_adjusted_7771min");
	if (runnumber == 89)	sprintf(filename, "%s", "run0089_LEGe_MSD_241Am_inChamber_window1.5us_CFDdelay_adjusted_11660min");
	if (runnumber == 90)	sprintf(filename, "%s", "run0090_LEGe_MSD_241Am_inChamber_window1.5us_CFDdelay_adjusted_4515min");
	if (runnumber == 91)	sprintf(filename, "%s", "run0091_LEGe_MSD_241Am_inChamber_window1.5us_CFDdelay_adjusted_1780min");
	if (runnumber == 92)	sprintf(filename, "%s", "run0092_LEGe_MSD26_241Am_inChamber_2mmCollimator_window1.5us_CFDdelay_adjusted_901min");
	if (runnumber == 93)	sprintf(filename, "%s", "run0093_LEGe_MSD26_241Am_inChamber_2mmCollimator_window1.5us_CFDdelay_adjusted_5000min");
	if (runnumber == 94)	sprintf(filename, "%s", "run0094_LEGe_MSD26_241Am_inChamber_2mmCollimator_window1.5us_CFDdelay_adjusted_9528min");
	if (runnumber == 95)	sprintf(filename, "%s", "run0095_LEGe_MSD26_241Am_inChamber_2mmCollimator_window1.5us_CFDdelay_adjusted_4099min");
	if (runnumber == 96)	sprintf(filename, "%s", "run0096_LEGe_241Am_inChamber_window1.5us_CFDdelay_adjusted_60min");
	if (runnumber == 97)	sprintf(filename, "%s", "run0097_LEGe_XtRa_152Eu_inChamber_window1.5us_CFDdelay_adjusted_180min");
	if (runnumber == 98)	sprintf(filename, "%s", "run0098_LEGe_XtRa_152Eu_inChamber_window1.5us_CFDdelay_adjusted_750min");
	if (runnumber == 99)	sprintf(filename, "%s", "run0099_LEGe_XtRa_152Eu_inChamber_flanges_swapped_window1.5us_CFDdelay_adjusted_420min");
	if (runnumber == 100)	sprintf(filename, "%s", "run0100_LEGe_MSD26_241Am_inChamber_2mmCollimator_window1.5us_CFDdelay_adjusted_997min");
	if (runnumber == 101)	sprintf(filename, "%s", "run0101_LEGe_MSD26_241Am_inChamber_2mmCollimator_window1.5us_CFDdelay_adjusted_9000min");
	if (runnumber == 102)	sprintf(filename, "%s", "run0102_LEGe_MSD26_241Am_inChamber_2mmCollimator_window1.5us_CFDdelay_adjusted_360min");
	if (runnumber == 103)	sprintf(filename, "%s", "run0103_LEGe_MSD26_241Am_inChamber_2mmCollimator_window1.5us_CFDdelay_adjusted_9000min");

	if (runnumber == 200)	sprintf(filename, "%s", "run0200_LEGe_MSD26_137Cs_M4038_inChamber_vacuum_North_16mm_South_12mm_away_window1.5us_CFDdelay_adjusted");
	if (runnumber == 201)	sprintf(filename, "%s", "run0201_LEGe_MSD26_137Cs_M4038_inChamber_vacuum_North_16mm_South_12mm_away_window1.5us_CFDdelay_adjusted");
	if (runnumber == 202)	sprintf(filename, "%s", "run0202_LEGe_MSD26_137Cs_M4038_inChamber_vacuum_North_16mm_South_12mm_away_window1.5us_CFDdelay_adjusted");
	if (runnumber == 203)	sprintf(filename, "%s", "run0203_LEGe_XtRa_MSD26_152Eu_Z2707_inChamber_vacuum_XtRa_12mm_away_window1.5us_CFDdelay_adjusted");
	if (runnumber == 204)	sprintf(filename, "%s", "run0204_XtRa_152Eu_Z2707_inChamber_vacuum_XtRa_12mm_away_window1.5us_CFDdelay_adjusted");
	if (runnumber == 205)	sprintf(filename, "%s", "run0205_XtRa_152Eu_Z2707_inChamber_vacuum_XtRa_12mm_away_window0.3us_CFDdelay_adjusted");
	if (runnumber == 206)	sprintf(filename, "%s", "run0206_XtRa_152Eu_Z2707_inChamber_vacuum_XtRa_220mm_away_window0.3us_CFDdelay_adjusted");
	if (runnumber == 207)	sprintf(filename, "%s", "run0207_XtRa_152Eu_Z2707_inChamber_vacuum_XtRa_220mm_away_window0.3us_CFDdelay_adjusted");
	if (runnumber == 208)	sprintf(filename, "%s", "run0208_LEGe_152Eu_Z2707_inChamber_vacuum_window0.1us_CFDdelay_adjusted");
	if (runnumber == 209)	sprintf(filename, "%s", "run0209_LEGe_152Eu_Z2707_inChamber_vacuum_window1.5us_CFDdelay_adjusted");
	if (runnumber == 210)	sprintf(filename, "%s", "run0210_LEGe_152Eu_Z2707_inChamber_vacuum_window0.02us_CFDdelay_adjusted");
	if (runnumber == 211)	sprintf(filename, "%s", "run0211_LEGe_152Eu_Z2707_inChamber_vacuum_window10us_CFDdelay_adjusted");
	if (runnumber == 212)	sprintf(filename, "%s", "run0212_XtRa_152Eu_Z2707_inChamber_vacuum_XtRa_12mm_away_window0.3us_CFDdelay_adjusted");
	if (runnumber == 213)	sprintf(filename, "%s", "run0213_North_152Eu_Z2707_On_Cap_window0.02us_CFDdelay_adjusted_Extreme_Summing_for_fun");
	if (runnumber == 214)	sprintf(filename, "%s", "run0214_LEGe_XtRa_MSD26_60Co_I7281_inChamber_vacuum_XtRa_12mm_away_window0.02us_CFDdelay_adjusted");
	if (runnumber == 215)	sprintf(filename, "%s", "run0215_LEGe_XtRa_MSD26_60Co_I7281_inChamber_vacuum_XtRa_12mm_away_window1us_CFDdelay_adjusted");
	if (runnumber == 216)	sprintf(filename, "%s", "run0216_LEGe_XtRa_MSD26_60Co_I7281_inChamber_vacuum_XtRa_12mm_away_window1us_CFDdelay_adjusted_AnalogGain1.0");
	if (runnumber == 217)	sprintf(filename, "%s", "run0217_LEGe_XtRa_MSD26_60Co_I7281_inChamber_vacuum_XtRa_12mm_away_window1us_CFDdelay_adjusted_AnalogGain1.0");
	if (runnumber == 218)	sprintf(filename, "%s", "run0218_LEGe_XtRa_MSD26_60Co_I7281_inChamber_vacuum_XtRa_12mm_away_window1us_CFDdelay_adjusted_AnalogGain1.0");
	if (runnumber == 219)	sprintf(filename, "%s", "run0219_LEGe_XtRa_MSD26_60Co_I7281_inChamber_vacuum_XtRa_12mm_away_window1us_CFDdelay0.5usforXtRa_AnalogGain1.0");
	if (runnumber == 220)	sprintf(filename, "%s", "run0220_LEGe_XtRa_MSD26_60Co_I7281_inChamber_vacuum_XtRa_12mm_away_window1us_CFDdelay0.5us_TriggerRiseandGap0.5usforXtRa_AnalogGain1.0");
	if (runnumber == 221)	sprintf(filename, "%s", "run0221_LEGe_XtRa_152Eu_Z2707_inChamber_atmosphere_XtRa_20ishmm_away_window1us_CFDdelay0.16us_AnalogGain4.0");
	if (runnumber == 222)	sprintf(filename, "%s", "run0222_LEGe_XtRa_152Eu_Z2707_inChamber_atmosphere_XtRa_20ishmm_away_window1us_CFDoff_AnalogGain4.0");
	if (runnumber == 223)	sprintf(filename, "%s", "run0223_LEGe_XtRa_152Eu_Z2707_inChamber_atmosphere_XtRa_12ishmm_away_window1us_CFDdelay0.5us_AnalogGain4.0");
	if (runnumber == 224)	sprintf(filename, "%s", "run0224_LEGe_XtRa_152Eu_Z2707_inChamber_atmosphere_XtRa_12ishmm_away_window1us_CFDoff_AnalogGain4.0");
	if (runnumber == 225)	sprintf(filename, "%s", "run0225_LEGe_XtRa_152Eu_Z2707_inChamber_atmosphere_XtRa_12mm_away_window1us_CFDdelay0.3us_AnalogGain4.0");
	if (runnumber == 228)	sprintf(filename, "%s", "run0228_LEGe_XtRa_MSD26_152Eu_Z2707_inChamber_vacuum_XtRa_12mm_away_window1us_XtRaCFDdelay_0.2us_for_efficiency");
	if (runnumber == 229)	sprintf(filename, "%s", "run0229_LEGe_XtRa_MSD26_152Eu_Z2707_inChamber_vacuum_XtRa_12mm_away_window1us_XtRaCFDdelay_0.2us_for_efficiency");
	if (runnumber == 230)	sprintf(filename, "%s", "run0230_LEGe_XtRa_MSD26_152Eu_Z2707_inChamber_vacuum_XtRa_12mm_away_window1us_XtRaCFDdelay_0.2us_for_efficiency");
	if (runnumber == 233)	sprintf(filename, "%s", "run0233_LEGe_XtRa_Room_background_XtRa_12mm_away_window1us_for_background_subtraction");
	if (runnumber == 234)	sprintf(filename, "%s", "run0234_LEGe_XtRa_60Co_I7281_onChamberWall_window1us_TrigRise0.016us_TrigGap1.000us_CFDdisabled");
	if (runnumber == 235)	sprintf(filename, "%s", "run0235_LEGe_XtRa_60Co_I7281_onChamberWall_window1us_TrigRise0.016us_TrigGap1.000us_CFDdelay0.504us");
	if (runnumber == 236)	sprintf(filename, "%s", "run0236_LEGe_XtRa_60Co_I7281_MiddleXtRa_window1us_TrigRise0.016us_TrigGap1.000us_ThN10000_ThS10000_LED");
	if (runnumber == 237)	sprintf(filename, "%s", "run0237_LEGe_XtRa_60Co_I7281_MiddleXtRa_window1us_TrigRise0.016us_TrigGap1.000us_ThN16384_ThS16384_LED");
	if (runnumber == 239)	sprintf(filename, "%s", "run0239_LEGe_XtRa_60Co_I7281_MiddleXtRa_window1us_TrigRise0.016us_TrigGap1.000us_ThN5000_ThS5000_LED");
	if (runnumber == 240)	sprintf(filename, "%s", "run0240_LEGe_XtRa_60Co_I7281_MiddleXtRa_window1us_TrigRise0.016us_TrigGap1.000us_ThN500_ThS950_LED");
	if (runnumber == 241)	sprintf(filename, "%s", "run0241_LEGe_XtRa_60Co_I7281_MiddleXtRa_window1us_TrigRise0.016us_TrigGap1.000us_ThN2000_ThS2000_LED");
	if (runnumber == 242)	sprintf(filename, "%s", "run0242_LEGe_XtRa_60Co_I7281_MiddleXtRa_window1us_TrigRise0.016us_TrigGap1.000us_ThN530_ThS1050_CFDDelay0.504us");
	if (runnumber == 243)	sprintf(filename, "%s", "run0243_LEGe_XtRa_60Co_I7281_MiddleXtRa_window1us_TrigRise0.016us_TrigGap1.000us_ThN530_ThS1050_CFDDelay0.104us");
	if (runnumber == 244)	sprintf(filename, "%s", "run0244_LEGe_XtRa_60Co_I7281_MiddleXtRa_window1us_TrigRise0.064us_TrigGap0.952us_ThN260_ThS610_LED");
	if (runnumber == 245)	sprintf(filename, "%s", "run0245_LEGe_XtRa_60Co_I7281_MiddleXtRa_window1us_TrigRise0.104us_TrigGap0.912us_ThN210_ThS340_LED");
	if (runnumber == 246)	sprintf(filename, "%s", "run0246_LEGe_XtRa_60Co_I7281_MiddleXtRa_window1us_TrigRise0.200us_TrigGap0.816us_ThN190_ThS330_LED");
	if (runnumber == 247)	sprintf(filename, "%s", "run0247_LEGe_XtRa_60Co_I7281_MiddleXtRa_window1us_TrigRise0.104us_TrigGap0.912us_ThN2520_ThS2520_LED");
	if (runnumber == 248)	sprintf(filename, "%s", "run0248_LEGe_XtRa_60Co_I7281_MiddleXtRa_window1us_TrigRise0.064us_TrigGap0.504us_ThN260_ThS610_CFDDelay0.504us");
	if (runnumber == 249)	sprintf(filename, "%s", "run0249_LEGe_XtRa_60Co_I7281_MiddleXtRa_window1us_TrigRise0.064us_TrigGap0.952us_ThN260_ThS610_CFDDelay0.504us");
	if (runnumber == 250)	sprintf(filename, "%s", "run0250_LEGe_XtRa_60Co_I7281_MiddleXtRa_window1us_TrigRise0.064us_TrigGap0.952us_ThN260_ThS610_CFDDelay0.504us_Scale4");
	if (runnumber == 251)	sprintf(filename, "%s", "run0251_LEGe_XtRa_60Co_I7281_MiddleXtRa_window1us_TrigRise0.064us_TrigGap0.952us_ThN260_ThS610_CFDDelay0.504us_Scale7");
	if (runnumber == 252)	sprintf(filename, "%s", "run0252_LEGe_XtRa_60Co_I7281_MiddleXtRa_window1us_TrigRise0.064us_TrigGap0.952us_ThN260_ThS610_CFDDelay0.104us_Scale7");
	if (runnumber == 253)	sprintf(filename, "%s", "run0253_LEGe_XtRa_60Co_I7281_MiddleXtRa_window1us_TrigRise0.048us_TrigGap0.968us_ThN270_ThS720_LED");
	if (runnumber == 254)	sprintf(filename, "%s", "run0254_LEGe_XtRa_60Co_I7281_MiddleXtRa_window1us_TrigRise0.032us_TrigGap0.984us_ThN330_ThS860_LED");
	if (runnumber == 255)	sprintf(filename, "%s", "run0255_LEGe_XtRa_60Co_I7281_MiddleXtRa_window1us_TrigRise0.064us_TrigGap0.952us_ThN300_ThS700_CFDDelay0.008us_Scale7");
	if (runnumber == 256)	sprintf(filename, "%s", "run0256_LEGe_XtRa_60Co_I7281_MiddleXtRa_window1us_TrigRise0.016us_TrigGap1.000us_ThN550_ThS1500_LED");
	if (runnumber == 257)	sprintf(filename, "%s", "run0257_LEGe_XtRa_60Co_I7281_MiddleXtRa_window1us_TrigRise0.064us_TrigGap0.952us_ThN300_ThS700_CFDDelay0.200us_Scale7");
	if (runnumber == 258)	sprintf(filename, "%s", "run0258_LEGe_XtRa_60Co_I7281_MiddleXtRa_window1us_TrigRise0.064us_TrigGap0.952us_ThN300_ThS700_CFDDelay0.064us_Scale7");
	if (runnumber == 259)	sprintf(filename, "%s", "run0259_LEGe_XtRa_60Co_I7281_MiddleXtRa_window1us_TrigRise0.064us_TrigGap0.952us_ThN300_ThS700_CFDDelay0.064us_Scale4");
	if (runnumber == 260)	sprintf(filename, "%s", "run0260_LEGe_XtRa_60Co_I7281_MiddleXtRa_window1us_TrigRise0.064us_TrigGap0.952us_ThN300_ThS700_CFDDelay0.120us_Scale7");
	if (runnumber == 261)	sprintf(filename, "%s", "run0261_LEGe_XtRa_60Co_I7281_MiddleXtRa_window1us_TrigRise0.064us_TrigGap0.952us_ThN300_ThS700_LED");
	if (runnumber == 262)	sprintf(filename, "%s", "run0262_LEGe_XtRa_60Co_I7281_OnChamberWall_window1us_TrigRise0.064us_TrigGap0.952us_ThL150_N300_S700_LED");
	if (runnumber == 263)	sprintf(filename, "%s", "run0263_LEGe_XtRa_60Co_I7281_OnLEGeCap_window1us_TrigRise0.064us_TrigGap0.952us_ThL150_N300_S700_CFDDelay0.120us_Scale7");
	if (runnumber == 264)	sprintf(filename, "%s", "run0264_LEGe_XtRa_152Eu_Z2707_OnLEGeCap_window1us_TrigRise0.064us_TrigGap0.952us_ThL310_N300_S730_LED");
	if (runnumber == 265)	sprintf(filename, "%s", "run0265_LEGe_XtRa_152Eu_Z2707_OnLEGeCap_window1us_TrigRise0.064us_TrigGap0.952us_ThL310_N300_S730_CFDDelay0.120us_Scale7");
	if (runnumber == 266)	sprintf(filename, "%s", "run0266_LEGe_XtRa_152Eu_Z2707_OnLEGeCap_window1us_TrigRise0.064us_TrigGap0.952us_ThL310_N300_S730_CFDDelay0.200us_Scale7");
	if (runnumber == 267)	sprintf(filename, "%s", "run0267_LEGe_XtRa_152Eu_Z2707_OnLEGeCap_window1us_TrigRise0.064us_TrigGap0.952us_ThL310_N300_S730_CFDDelay0.304us_Scale7");
	if (runnumber == 268)	sprintf(filename, "%s", "run0268_LEGe_XtRa_152Eu_Z2707_OnLEGeCap_window1us_TrigRise0.064us_TrigGap0.952us_ThL310_N300_S730_CFDDelay0.400us_Scale7");
	if (runnumber == 269)	sprintf(filename, "%s", "run0269_LEGe_XtRa_152Eu_Z2707_OnLEGeCap_window1us_TrigRise0.200_0.064us_TrigGap0.816_0.952us_ThL310_N300_S730_CFDDelay0.200us_Scale7");
	if (runnumber == 270)	sprintf(filename, "%s", "run0270_LEGe_XtRa_152Eu_Z2707_OnLEGeCap_window1us_TrigRise0.064_0.064us_TrigGap0.600_0.952us_ThL310_N300_S730_CFDDelay0.304us_Scale7");
	if (runnumber == 271)	sprintf(filename, "%s", "run0271_LEGe_XtRa_152Eu_Z2707_OnLEGeCap_window1us_TrigRise0.064us_TrigGap0.952us_ThL310_N300_S730_CFDDelay0.120_0.304us_Scale3_7");
	if (runnumber == 272)	sprintf(filename, "%s", "run0272_LEGe_XtRa_152Eu_Z2707_OnLEGeCap_window1us_TrigRise0.064us_TrigGap0.952us_ThL310_N300_S730_CFDDelay0.304_0.120us_Scale7");
	if (runnumber == 273)	sprintf(filename, "%s", "run0273_LEGe_XtRa_152Eu_Z2707_OnLEGeCap_window1us_TrigRise0.064us_TrigGap0.952us_ThL310_N300_S730_CFDDelay0.120us_Scale0_7");
	if (runnumber == 274)	sprintf(filename, "%s", "run0274_LEGe_XtRa_152Eu_Z2707_OnLEGeCap_window1us_TrigRise0.064us_TrigGap0.952us_ThL310_N300_S730_CFDDelay0.120us_Scale7");
	if (runnumber == 275)	sprintf(filename, "%s", "run0275_LEGe_XtRa_152Eu_Z2707_OnLEGeCap_window1us_TrigRise0.064us_TrigGap0.952us_ThL310_N300_S730_CFDDelay0.120us_Scale2_7");
	if (runnumber == 276)	sprintf(filename, "%s", "run0276_LEGe_XtRa_152Eu_Z2707_OnLEGeCap_window1us_TrigRise0.064us_TrigGap0.952us_ThL310_N300_S730_CFDDelay0.304us_Scale7");
	if (runnumber == 277)	sprintf(filename, "%s", "run0277_LEGe_XtRa_152Eu_Z2707_OnLEGeCap_window1us_TrigRise0.064us_TrigGap0.304us_0.952us_ThL310_N300_S730_CFDDelay0.304us_Scale7");
	if (runnumber == 278)	sprintf(filename, "%s", "run0278_LEGe_XtRa_152Eu_Z2707_OnLEGeCap_window1us_TrigRise0.200us_0.064us_TrigGap0.104us_0.952us_ThL250_N300_S730_CFDDelay0.304us_Scale0_7");
	if (runnumber == 279)	sprintf(filename, "%s", "run0279_LEGe_XtRa_152Eu_Z2707_OnLEGeCap_window1us_TrigRise0.064us_TrigGap0.952us_ThL310_N300_S730_CFDDelay0.104_0.304us_Scale7");
	if (runnumber == 280)	sprintf(filename, "%s", "run0280_Pulser_window1us_TrigRise0.064us_TrigGap0.952us_ThN300_S730_LED");
	if (runnumber == 281)	sprintf(filename, "%s", "run0281_Pulser_window1us_TrigRise0.064us_TrigGap0.952us_ThN300_S730_CFDDelay0.304us_Scale7");
	if (runnumber == 282)	sprintf(filename, "%s", "run0282_Pulser_window1us_TrigRise0.064us_TrigGap0.952us_ThN300_S300_LED");
	if (runnumber == 283)	sprintf(filename, "%s", "run0283_Pulser_window1us_TrigRise0.064us_TrigGap0.952us_ThN300_S300_CFDDelay0.304us_Scale7");
	if (runnumber == 284)	sprintf(filename, "%s", "run0284_LEGe_XtRa_152Eu_Z2707_OnLEGeCap_window1us_TrigRise0.064us_TrigGap0.952us_ThL310_N300_S730_CFDDelay0.304us_Scale7");
	if (runnumber == 285)	sprintf(filename, "%s", "run0285_LEGe_XtRa_152Eu_Z2707_OnLEGeCap_window1us_TrigRise0.064us_TrigGap0.952us_ThL310_N300_S730_CFDDelay0.200_0.304us_Scale7");
	if (runnumber == 286)	sprintf(filename, "%s", "run0286_LEGe_XtRa_152Eu_Z2707_OnLEGeCap_window1us_TrigRise0.064us_TrigGap0.952us_ThL310_N300_S730_CFDDelay0.200_0.304us_Scale4_7");
	if (runnumber == 287)	sprintf(filename, "%s", "run0287_LEGe_XtRa_152Eu_Z2707_OnLEGeCap_window1us_TrigRise0.064us_TrigGap0.952us_ThL310_N300_S730_CFDDelay0.200_0.304us_Scale1_7");
	if (runnumber == 288)	sprintf(filename, "%s", "run0288_LEGe_XtRa_152Eu_Z2707_OnLEGeCap_window1us_TrigRise0.200us_0.064us_TrigGap0.104us_0.952us_ThL310_N300_S730_CFDDelay0.304us_Scale0_7");
	if (runnumber == 289)	sprintf(filename, "%s", "run0289_LEGe_XtRa_152Eu_Z2707_OnLEGeCap_window1us_TrigRise0.064us_TrigGap0.952us_ThL310_N300_S730_LED");
	if (runnumber == 290)	sprintf(filename, "%s", "run0290_LEGe_XtRa_152Eu_Z2707_OnLEGeCap_window1us_TrigRise0.016_0.064us_TrigGap1.000_0.952us_ThL330_N300_S730_LED");
	if (runnumber == 291)	sprintf(filename, "%s", "run0291_LEGe_XtRa_152Eu_Z2707_OnLEGeCap_window1us_TrigRise0.112_0.064us_TrigGap0.904_0.952us_ThL270_N300_S730_LED");
	if (runnumber == 292)	sprintf(filename, "%s", "run0292_LEGe_XtRa_152Eu_Z2707_OnLEGeCap_window1us_TrigRise0.200_0.064us_TrigGap0.816_0.952us_ThL240_N300_S730_LED");
	if (runnumber == 293)	sprintf(filename, "%s", "run0293_LEGe_XtRa_152Eu_Z2707_OnLEGeCap_window1us_TrigRise0.304_0.064us_TrigGap0.712_0.952us_ThL210_N300_S730_LED");
	if (runnumber == 294)	sprintf(filename, "%s", "run0294_LEGe_XtRa_152Eu_Z2707_OnLEGeCap_window1us_TrigRise0.400_0.064us_TrigGap0.616_0.952us_ThL190_N300_S730_LED");
	if (runnumber == 295)	sprintf(filename, "%s", "run0295_LEGe_MSD_241Am_Z7117_ChamberCenter_window1us_TrigRise0.064_0.064_0.112us_TrigGap0.952_0.952_0.904us_Th350_2400_500_LED");
	if (runnumber == 296)	sprintf(filename, "%s", "run0296_LEGe_MSD_241Am_Z7117_ChamberCenter_window1us_TrigRise0.064_0.064_0.112us_TrigGap0.952_0.952_0.904us_Th350_2400_500_CFDDelay0.304us_Scale7");
	if (runnumber == 297)	sprintf(filename, "%s", "run0297_LEGe_MSD_241Am_Z7117_ChamberCenter_window1us_TrigRise0.064_0.016_0.016us_TrigGap0.952_1.000_1.000us_Th350_2700_1000_LED");
	if (runnumber == 298)	sprintf(filename, "%s", "run0298_LEGe_MSD_241Am_Z7117_ChamberCenter_window1us_TrigRise0.064_0.016_0.016us_TrigGap0.952_1.000_1.000us_Th350_2700_1000_MSDCFDDelay0.304us_Scale7");
	if (runnumber == 299)	sprintf(filename, "%s", "run0299_LEGe_MSD_241Am_Z7117_ChamberCenter_window1us_TrigRise0.064_0.016_0.016us_TrigGap0.952_1.000_1.000us_Th350_2700_1000_CFDDelay0.304us_Scale7");
	if (runnumber == 328)	sprintf(filename, "%s", "run0328_LEGe_MSD_241Am_Z7117_ChamberCenter_window1us_TrigRise0.064_0.016_0.016us_TrigGap0.952_1.000_1.000us_Th350_2700_1000_CFDDelay0.104us_Scale7");
	if (runnumber == 329)	sprintf(filename, "%s", "run0329_LEGe_MSD_241Am_Z7117_ChamberCenter_window1us_TrigRise0.064_0.016_0.016us_TrigGap0.952_1.000_1.000us_Th350_2700_1000_CFDDelay0.504us_Scale7");
	if (runnumber == 330)	sprintf(filename, "%s", "run0330_LEGe_MSD_241Am_Z7117_ChamberCenter_window1.5us_TrigRise0.064_0.016_0.016us_TrigGap0.952_1.000_1.000us_Th350_2700_1000_CFDDelay0.304us_Scale7");
	if (runnumber == 331)	sprintf(filename, "%s", "run0331_LEGe_MSD_241Am_Z7117_ChamberCenter_window1.5us_TrigRise0.064_0.016_0.016us_TrigGap0.952_1.000_1.000us_Th350_2700_1000_CFDDelay0.304us_Scale7");
	if (runnumber == 332)	sprintf(filename, "%s", "run0332_LEGe_MSD_241Am_Z7117_ChamberCenter_window1.5us_TrigRise0.064_0.016_0.016us_TrigGap0.952_1.000_1.000us_Th350_2700_1000_CFDDelay0.304us_Scale7");
	if (runnumber == 333)	sprintf(filename, "%s", "run0333_LEGe_241Am_Z7117_ChamberCenter_window1.5us_TrigRise0.064us_TrigGap0.952us_Th350_CFDDelay0.304us_Scale7");
	if (runnumber == 334)	sprintf(filename, "%s", "run0334_LEGe_241Am_Z7117_ChamberCenter_window1.5us_TrigRise0.064us_TrigGap0.952us_Th350_CFDDelay0.304us_Scale7_for_X_efficiency");
	if (runnumber == 335)	sprintf(filename, "%s", "run0335_LEGe_MSD_137Cs_M4038_facingMSD_inChamber_vacuum_window1us");
	if (runnumber == 336)	sprintf(filename, "%s", "run0336_LEGe_MSD_137Cs_M4038_facingMSD_inChamber_vacuum_window1us_TriseTgap_adjusted");
	if (runnumber == 337)	sprintf(filename, "%s", "run0337_LEGe_MSD26_137Cs_M4038_facingMSD26_inChamber_vacuum_window1us");
	if (runnumber == 338)	sprintf(filename, "%s", "run0338_LEGe_MSD26_137Cs_M4038_facingMSD26_inChamber_vacuum_window1us");
	
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

	
	// name both input and output root files accordingly
	sprintf(rawrootname, "%s%s%s", inpathname, filename, ".root");
	sprintf(calrootname, "%s%s%s", outpathname, filename, "_cal.root");

	TFile* pFile = new TFile(rawrootname);
	cout << "input file: " << rawrootname << endl;
	TTree* pTree;
	pFile->GetObject("dchan", pTree);
	totalentries = pTree->GetEntries();
	cout << "Total Entries=" << totalentries << endl;

	DDASEvent* pEvent = new DDASEvent;
	std::vector <ddaschannel*> pChan;
	pTree->SetBranchAddress("ddasevent", &pEvent);

	Double_t lege_e;
	Double_t lege_t;
	Double_t lege_tled;
	Double_t north_e;
	Double_t north_t;
	Double_t north_tled;
	Double_t south_e;
	Double_t south_t;
	Double_t south_tled;
	Double_t msd12_e;
	Double_t msd12_t;
	Double_t msd12_tled;
	Double_t msd26_e;
	Double_t msd26_t;
	Double_t msd26_tled;
	Double_t msdtotal_e;

	Double_t lege_e_low;
	Double_t lege_t_low;
	Double_t north_e_low;
	Double_t north_t_low;
	Double_t south_e_low;
	Double_t south_t_low;
	Double_t pulser_e;
	Double_t pulser_t;

	TFile* fout = new TFile(calrootname, "RECREATE");
	cout << "output file: " << calrootname << endl;
	TTree* tree = new TTree("tree", "tree");

	// If you have closed these channels in CSRA, please comment out the corresponding tree->branch lines. No other changes are needed.
	//tree->Branch("lege_e", &lege_e, "lege_e/D");
	//tree->Branch("lege_t", &lege_t, "lege_t/D");
	//tree->Branch("lege_tled", &lege_tled, "lege_tled/D");
	//tree->Branch("north_e", &north_e, "north_e/D");
	//tree->Branch("north_t", &north_t, "north_t/D");
	//tree->Branch("north_tled", &north_tled, "north_tled/D");
	//tree->Branch("south_e", &south_e, "south_e/D");
	//tree->Branch("south_t", &south_t, "south_t/D");
	//tree->Branch("south_tled", &south_tled, "south_tled/D");
  	//tree->Branch("msd12_e", &msd12_e, "msd12_e/D");
  	//tree->Branch("msd12_t", &msd12_t, "msd12_t/D");
	//tree->Branch("msd12_tled", &msd12_tled, "msd12_tled/D"); // likely useless
  	//tree->Branch("msd26_e", &msd26_e, "msd26_e/D");
  	//tree->Branch("msd26_t", &msd26_t, "msd26_t/D");
	//tree->Branch("msd26_tled", &msd26_tled, "msd26_tled/D"); // likely useless

// 	tree->Branch("lege_e_low", &lege_e_low, "lege_e_low/D");
// 	tree->Branch("lege_t_low", &lege_t_low, "lege_t_low/D");
 	//tree->Branch("north_e_low", &north_e_low, "north_e_low/D");
 	//tree->Branch("north_t_low", &north_t_low, "north_t_low/D");
 	//tree->Branch("south_e_low", &south_e_low, "south_e_low/D");
 	//tree->Branch("south_t_low", &south_t_low, "south_t_low/D");
// 	tree->Branch("pulser_e", &pulser_e, "pulser/D");
// 	tree->Branch("pulser_t", &pulser_t, "pulser_t/D");
//  	TBranch* b_pulser_e = tree->Branch("pulser_e", &pulser_e, "pulser_e/D");
//  	TBranch* b_pulser_t = tree->Branch("pulser_t", &pulser_t, "pulser_t/D");

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
	TH1D* htiming_msd12_msd26;

	// 	if (runnumber == 0)	hlege_e = new TH1D("hlege_e", "hlege_e", 60000, -0.0270924805469, 436.3668916340480); // 60000 channels, 7.3 eV per channel
	// 	if (runnumber == 1)	hlege_e = new TH1D("hlege_e", "hlege_e", 60000, -0.0584419300571, 435.7108166612020); // 60000 channels, 7.3 eV per channel
	// 	if (runnumber == 55)	hlege_e = new TH1D("hlege_e", "hlege_e", 60000, 0.0080746242334, 437.1574593691120); // 60000 channels, 7.3 eV per channel
	// 	if (runnumber == 152)	hlege_e = new TH1D("hlege_e", "hlege_e", 60000, -0.0001420032512, 436.2768890605690); // 60000 channels, 7.3 eV per channel
	//	if (runnumber == 241)	hlege_e = new TH1D("hlege_e", "hlege_e", 60000, 0.0541091878755, 435.4667414524260); // 60000 channels, 7.3 eV per channel
	//hlege_e = new TH1D("hlege_e", "hlege_e", 60000, 0.0541091878755, 435.4667414524260); // 60000 channels, 7.3 eV per channel
	//hlege_e = new TH1D("hlege_e", "hlege_e", 60000, -0.0540888491631, 435.8137733893930); // 60000 channels, 7.3 eV per channel run0060
	hlege_e = new TH1D("hlege_e", "hlege_e", 60000, 0.7583911439124, 1679.872782287820); // 60000 channels, 28 eV per channel run0216 60Co
	// hnorth_e = new TH1D("hnorth_e", "hnorth_e", 60000, -0.1774488889544, 2286.04224413527); // 60000 channels, 38 eV per channel from Excel calibration
	hnorth_e = new TH1D("hnorth_e", "hnorth_e", 60000, -0.1943988436200, 2285.29642151811); // 60000 channels, 38 eV per channel from Excel calibration 1/11/2024 run97-98
	//hnorth_e = new TH1D("hnorth_e", "hnorth_e", 60000, -1.1199018269051, 8815.85740025835); // 60000 channels, 147 eV per channel from Excel calibration run0216 60Co 2/25/2024
	// hsouth_e = new TH1D("hsouth_e", "hsouth_e", 60000, -0.1622624327090, 2278.22727449908); // 60000 channels, 38 eV per channel from Excel calibration
	hsouth_e = new TH1D("hsouth_e", "hsouth_e", 60000, -0.1195774566938, 2274.16156253889); // 60000 channels, 38 eV per channel from Excel calibration 1/11/2024 run97-98
	//hsouth_e = new TH1D("hsouth_e", "hsouth_e", 60000, -0.6039596145015, 8770.24505369435); // 60000 channels, 146 eV per channel from Excel calibration run0216 60Co 2/25/2024
	hmsd12_e = new TH1D("hmsd12_e", "hmsd12_e", 60000, 0, 7941.93363844394); // 60000 channels, 512 eV per channel
	//hmsd26_e = new TH1D("hmsd26_e", "hmsd26_e", 60000, -0.0016133749182, 6931.63915798418); // 60000 channels, 115 eV per channel from Excel calibration 1/18/2024 run 100
	hmsd26_e = new TH1D("hmsd26_e", "hmsd26_e", 60000, 0, 7874.6844); // 60000 channels, 130 eV per channel run 200
	hmsdtotal_e = new TH1D("hmsdtotal_e", "hmsdtotal_e", 6900, 0, 6900); // 6900 channels, 1000 eV per channel
	hmsd12_e1keVbin = new TH1D("hmsd12_e1keVbin", "hmsd12_e1keVbin", 7942, 0, 7942); // 7942 channels, 1 keV per channel
	hmsd26_e1keVbin = new TH1D("hmsd26_e1keVbin", "hmsd26_e1keVbin", 6935, 0, 6935); // 6935 channels, 1 keV per channel

	//hlege_e_low = new TH1D("hlege_e_low", "hlege_e_low", 60000, 0, 60000); // 60000 channels
	//hnorth_e_low = new TH1D("hnorth_e_low", "hnorth_e_low", 60000, 0, 60000); // 60000 channels
	//hsouth_e_low = new TH1D("hsouth_e_low", "hsouth_e_low", 60000, 0, 60000); // 60000 channels
	htiming_lege_msd12 = new TH1D("htiming_lege_msd12", "htiming_lege_msd12", 3000, -1500, 1500); // 3000 channels, 1 ns per channel
	htiming_lege_msd26 = new TH1D("htiming_lege_msd26", "htiming_lege_msd26", 3000, -1500, 1500); // 3000 channels, 1 ns per channel
	htiming_msd12_msd26 = new TH1D("htiming_msd12_msd26", "htiming_msd12_msd26", 3000, -1500, 1500); // 3000 channels, 1 ns per channel
	//i_outtree = 0;

	for (i = 0; i < totalentries; i++)
	{
		lege_e = 0; lege_t = 0; lege_tled = 0; north_e = 0; north_t = 0; north_tled = 0; south_e = 0; south_t = 0; south_tled = 0; msd12_e = 0; msd12_t = 0; msd12_tled = 0; msd26_e = 0; msd26_t = 0; msd26_tled = 0; msdtotal_e = 0;
		lege_e_low = 0; lege_t_low = 0; north_e_low = 0; north_t_low = 0; south_e_low = 0; south_t_low = 0; pulser_e = 0; pulser_t = 0;
		pTree->GetEntry(i);
		pChan = pEvent->GetData();
		Nchannels = (int)pChan.size();  //get active number of channels
		
		for (int j = 1; j < Nchannels; j++)
		{
			if (pChan[j]->GetChannelID() == pChan[j - 1]->GetChannelID())
			{
				cout << "!!! Multi events in one entry: " << i << " ch: " << pChan[j]->GetChannelID() << ", e: " << pChan[j]->GetEnergy() << ", t: " << pChan[j]->GetTime() << endl;
			}
		}
		for (int j = 0; j < Nchannels; j++)
		{
// 			if (pChan[j]->GetChannelID() == 9)
// 			{
// 				pulser_e = pChan[j]->GetEnergy();
// 				pulser_t = pChan[j]->GetTime();
// 				//b_pulser_e->Fill();
// 				//b_pulser_t->Fill();
// 				//tree->SetEntries(i_outtree++);
// 				//if (pulser_e > 0) tree->Fill();
// 				//tree->GetBranch("pulser_e")->Fill();
// 				//tree->GetBranch("pulser_t")->Fill();
// 			}
// 			if (pChan[j]->GetChannelID() == 0)
// 			{
// 				lege_e_low = pChan[j]->GetEnergy();
// 				//hlege_e_low->Fill(lege_e_low);
// 				lege_t_low = pChan[j]->GetTime() + timestampshift;
// 			}
			if (pChan[j]->GetChannelID() == 1)
			{
				//lege_e = pChan[j]->GetEnergy();
				//lege_e = pChan[j]->GetEnergy() * 0.0072712838511 - 0.0001420032512;
				//lege_e = pChan[j]->GetEnergy() * 0.0072644643706 - 0.054088849163; // based on run0060 Gain 4.0
				lege_e = pChan[j]->GetEnergy() * 0.0279852398524 - 0.758391143912; // based on run0216 Gain 1.0
				hlege_e->Fill(lege_e);
				lege_t = pChan[j]->GetTime() + timestampshift;
				//lege_tled = pChan[j]->GetCoarseTime() + timestampshift;
			}
 			//if (pChan[j]->GetChannelID() == 2)
 			//{
 			//	north_e_low = pChan[j]->GetEnergy();
 			//	//hnorth_e_low->Fill(north_e_low);
 			//	north_t_low = pChan[j]->GetTime() + timestampshift;
 			//}
			if (pChan[j]->GetChannelID() == 3)
			{
				//north_e = pChan[j]->GetEnergy();
				//north_e = pChan[j]->GetEnergy() * 0.0381036615504 - 0.1774488889544; // based on run0060
				north_e = pChan[j]->GetEnergy() * 0.0380915136727 - 0.1943988436200; // based on run0098
				//north_e = pChan[j]->GetEnergy() * 0.1469496217014 - 1.1199018269051; // based on run0216
				hnorth_e->Fill(north_e);
				north_t = pChan[j]->GetTime() + timestampshift;
				//north_tled = pChan[j]->GetCoarseTime() + timestampshift;
			}
 			//if (pChan[j]->GetChannelID() == 4)
 			//{
 			//	south_e_low = pChan[j]->GetEnergy();
 			//	//hsouth_e_low->Fill(south_e_low);
 			//	south_t_low = pChan[j]->GetTime() + timestampshift;
 			//}
			if (pChan[j]->GetChannelID() == 5)
			{
				//south_e = pChan[j]->GetEnergy();
				//south_e = pChan[j]->GetEnergy() * 0.0379731589489 - 0.1622624327090; // based on run0060
				south_e = pChan[j]->GetEnergy() * 0.0379046856666 - 0.1195774566938; // based on run0098
				//south_e = pChan[j]->GetEnergy() * 0.1461808168885 - 0.6039596145015; // based on run0216
				hsouth_e->Fill(south_e);
				south_t = pChan[j]->GetTime() + timestampshift;
				//south_tled = pChan[j]->GetCoarseTime() + timestampshift;
			}
 			if (pChan[j]->GetChannelID() == 6)
 			{
 				//msd12_e = pChan[j]->GetEnergy();
 				msd12_e = pChan[j]->GetEnergy() * 0.1323655606407 + 0;
 				hmsd12_e->Fill(msd12_e);
 				//hmsd12_e1keVbin->Fill(msd12_e + gRandom->Uniform(-0.5, 0.5));
 				msd12_t = pChan[j]->GetTime() + timestampshift;
 				//msd12_tled = pChan[j]->GetCoarseTime() + timestampshift;
 			}
			if (pChan[j]->GetChannelID() == 8)
			{
				//msd26_e = pChan[j]->GetEnergy();
				//msd26_e = pChan[j]->GetEnergy() * 0.1152143446474 + 22.352486402900; // based on 148Gd and 241Am two peaks
				//msd26_e = pChan[j]->GetEnergy() * 0.1155273461893 - 0.001613374918; // based on 241Am two peaks Run 100
				msd26_e = pChan[j]->GetEnergy() * 0.13124474; // based on 137Cs internal conversion electronic Run 335
				hmsd26_e->Fill(msd26_e+gRandom->Uniform(-0.26,0.26));
				//hmsd26_e1keVbin->Fill(msd26_e + gRandom->Uniform(-0.5, 0.5));
				msd26_t = pChan[j]->GetTime() + timestampshift;
				//msd26_tled = pChan[j]->GetCoarseTime() + timestampshift;
			}

			if (msd12_e > 100 && msd26_e > 100)
			{
				hmsd12_e1keVbin->Fill(msd12_e + gRandom->Uniform(-0.5, 0.5));
				hmsd26_e1keVbin->Fill(msd26_e + gRandom->Uniform(-0.5, 0.5));
				msdtotal_e = msd12_e + msd26_e;
				hmsdtotal_e->Fill(msdtotal_e);
			}

			//if lege_e > 59.0 && lege_e < 60.1 && msd26_e > 5470 && msd26_e < 5510)
			//{
			//	htiming_lege_msd26->Fill(lege_t - msd26_t);
			//}

			if (lege_e > 59.0 && lege_e < 60.1 && msd12_e>1700 && msd12_e < 2200 && msd26_e>3300 && msd26_e < 3800 && msd12_e + msd26_e > 5400 && msd12_e + msd26_e < 5560)
			{
				htiming_lege_msd12->Fill(lege_t - msd12_t);
				htiming_lege_msd26->Fill(lege_t - msd26_t);
			}

			if (msd12_e>1700 && msd12_e < 2200 && msd26_e>3300 && msd26_e < 3800 && msd12_e + msd26_e > 5400 && msd12_e + msd26_e < 5560)
			{
				htiming_msd12_msd26->Fill(msd12_t - msd26_t);
			}

		} // (int j = 0; j < Nchannels; j++)
		tree->Fill();

		if (i == 0)
		{
			time(&start);
			at = localtime(&start);
			strftime(now, 79, "%Y-%m-%d %H:%M:%S", at);
			cout << now << " started." << endl;
		}
		if (i % 1000000 == 0 && i != 0)
		{
			time(&tim);
			at = localtime(&tim);
			strftime(now, 79, "%Y-%m-%d %H:%M:%S", at);
			cout << now;
			speed = (float)i / (float)difftime(tim, start);
			printf(" %s%ld%s%.1f%s%.0f%s%.0f%s\n", "entry: ", i, " ",100 * (float)i / (float)totalentries, "% done. Speed is ", speed, " e/s. Still need ", float(totalentries - i) / speed, " s.");
		}
		//break at
// 		if (msd26_t >= 21205617870819.0-2e9)
// 		{
// 			cout << "break at: " << i << "	" << msd26_t << endl;
// 			break;
// 		}
	} // for (int i = 0; i < totalentries; i++)
	cout << "Last entry: " << i << endl;
	fout->Write();
	fout->Close();
	pFile->Close();
}
