#include "TFile.h"  
#include "TTree.h"

void tree_scan_print_with_more_digits() {

	TFile fin("F:/e21010/pxct/run0039-00_LEGe_MSD26_241Am_2mmcollimator_60min_window1us_CFDdelay0.3us_cal.root");

	TTree* tree = (TTree*)fin.Get("tree");

	if (!tree) {
		cout << "Tree not found!" << endl;
		return;
	}

	tree->Print("lege_t:%.10f:msd26_t:%.10f", "abs(lege_t-msd26_t)>0&&abs(lege_t-msd26_t)<1000");
	fin.Close();
}