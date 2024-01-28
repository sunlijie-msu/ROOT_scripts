#include <iostream>
#include <fstream>
//#include <iomanip.h>
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
void Fill_histogram_from_txt_byFill_DSL2_Bayesian()//read bincontent from a txt file and create a histo and save it as a new root file.
{// input bin by bin, not event by event
	int ii,jj;
	float Channel[500];
	float x,y,xe,ye;
	long Count[500];
	char outrootname[200], datname[200], tname[200];
	TCanvas *canvas;
	int Entries=300;//modify
	sprintf(tname, "%s", "D:/X/out/Bayesian_VS/Bayesian_DSL/DSL_31S4156_3fs/31S4156_samples");
	sprintf(datname, "%s%s", tname, ".dat");
	sprintf(outrootname, "%s%s", tname, ".root");
	ifstream infile(datname,ios::in);
	TFile *fout = new TFile(outrootname,"RECREATE");//输出文件//It's better to define histograms and then define fout, in case of draw bugs.
	canvas=new TCanvas("c1","c1", 900,600);//建立画布
	TH1F *h1 = new TH1F("h1","h1",3000,0,30);//create a histogram
	string line;
	while ( getline(infile, line) )
	{
		stringstream(line) >> x;
		h1->Fill(x);
	}
	canvas->cd();//进入画布
	h1->GetXaxis()->SetTitle("Lifetime (fs)");//轴名
	h1->GetYaxis()->SetTitle("Counts per 0.01 fs");//轴名
	h1->GetXaxis()->CenterTitle();//居中
	h1->GetYaxis()->CenterTitle();//居中
	h1->Draw();
	canvas->Update();
	cout<<"RMS= "<<setiosflags(ios::fixed)<<setprecision(9)<<h1->GetRMS();//you can also find the RMS in the rootfile
	fout->Write();
	fout->Close();
}