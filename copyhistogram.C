#include <iostream>
#include <fstream>
#include <iomanip.h>
#include <math.h>
#include <time.h>
#include <map>
#include <TROOT.h>
#include <TSystem.h>
#include <TStyle.h>
#include <TCutG.h>
#include<TChain.h>
#include "TChain.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TEllipse.h"
#include "TString.h"
#include "TGraph.h"
#include "TSpectrum.h"
#include "TVirtualFitter.h"
#include "TRandom.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TMath.h"
#include "TPaveText.h"
#include "TPad.h"
#include "string.h"
#include "stdlib.h"
#include "stdio.h"
#include "TLegend.h"
using namespace std;
void copyhistogram()//copyhistogram from a root to another root
{
	time_t tim;
	struct tm *at;
	char now[80];
	int irawroot;
	char rawrootname[80];
	char simurootname[80];
	
	//sprintf(rawrootname,"X:/T999/P26_0154_0345_pz.root");//modify if change path or nuclide
	//sprintf(simurootname,"%s","X:/T999/P26_pz.root");//modify if change path or nuclide
	sprintf(rawrootname,"X:/T999/other/Al22_0373_0538_pz.root");//modify if change path or nuclide
	sprintf(simurootname,"%s","X:/T999/other/Al22_pz.root");//modify if change path or nuclide
	TFile *fin = new TFile(rawrootname);//after this statement, you can use any ROOT command for this rootfile
	unsigned long nentries=T999->GetEntries();//读事件数
	cout<<rawrootname;
	cout<<"  Entries="<<nentries<<endl;
	TFile *fout = new TFile(simurootname,"RECREATE");//输出文件//It's better to define histograms and then define fout, in case of draw bugs.
 	TH1D *h142px = (TH1D*)fin->Get("himpl142px");//read a histogram from another Root file
 	TH1D *h142py = (TH1D*)fin->Get("himpl142py");//read a histogram from another Root file
 	TH1D *h142pz = (TH1D*)fin->Get("himpl142pz");//read a histogram from another Root file
	TH1D *h40px= (TH1D*)fin->Get("himpl40px");//read a histogram from another Root file
 	TH1D *h40py = (TH1D*)fin->Get("himpl40py");//read a histogram from another Root file
 	TH1D *h40pz = (TH1D*)fin->Get("himpl40pz");//read a histogram from another Root file
	TH1D *h304px = (TH1D*)fin->Get("himpl304px");//read a histogram from another Root file
 	TH1D *h304py = (TH1D*)fin->Get("himpl304py");//read a histogram from another Root file
	TH1D *h304pz = (TH1D*)fin->Get("himpl304pz");//read a histogram from another Root file
	//h304px->SetName("P26");
	h142px->Write();
	h142py->Write();
	h142pz->Write();
	h40px->Write();
	h40py->Write();
	h40pz->Write();
	h304px->Write();
	h304py->Write();
	h304pz->Write();

	fout->Close();
}