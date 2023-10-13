void print_histogram() {
	TFile* fin = new TFile("F:/e21010/pxct/run0012-00_LEGe_55Fe_217min_cal.root");
	TH1F* histo = (TH1F*)fin->Get("hlege_e"); //Get spectrum
	histo->Print("all");
}
//usage: root -l -q print_histogram.C > h1.log