void print_histogram() {
	TFile* fin = new TFile("F:/e21010/pxct/run0200_0201_0202_LEGe_MSD26_137Cs_M4038_inChamber_vacuum_North_16mm_South_12mm_away_window1.5us_CFDdelay_adjusted_cal.root");
	TH1F* histo = (TH1F*)fin->Get("hsouth_e"); //Get spectrum
	histo->Print("all");
}
//usage: root -l -q print_histogram.C > h1.log