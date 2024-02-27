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

void get_last_entry(const char* fileName)
{
	TFile* file = TFile::Open(fileName);
	if (!file || file->IsZombie()) {
		std::cerr << "Error opening file: " << fileName << std::endl;
		return;
	}

	TTree* pTree;
	file->GetObject("dchan", pTree);
	int totalentries = pTree->GetEntries();

	DDASEvent* pEvent = new DDASEvent;
	std::vector<ddaschannel*> pChan;
	pTree->SetBranchAddress("ddasevent", &pEvent);
	pTree->GetEntry(totalentries - 1); // Get the last entry
	pChan = pEvent->GetData();

	// Open the output file in append mode
	std::ofstream outfile("timestamp_output.txt", std::ios_base::app);
	outfile << "File: " << fileName << ", Total Entries: " << totalentries;
	if (!pChan.empty()) {
		outfile << ", Last Time Stamp: " << std::fixed << std::setprecision(1) << pChan[0]->GetTime() << std::endl;
	}
	else {
		outfile << ", No channels in the last event." << std::endl;
	}

	// Close the file
	outfile.close();
	delete file; // Properly close and delete the TFile object
}
