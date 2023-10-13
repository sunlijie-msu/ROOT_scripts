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
#include "TChain.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TEllipse.h"
#include "TString.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TSpectrum.h"
#include "TVirtualFitter.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"
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
using namespace std;
void Graphfit()//read x,y from a txt file, and then fit with a exp function
{
	int ii,jj;
	int icalroot;
	float Channel[500];
	long Count[500];
	char rootname[80];
	ifstream infile("C:/Si24/Si22peakcali/Xu/Cali.dat",ios::in);//The data that need to be fitted
	TGraph *g = new TGraph();
	TH1F *h1 = new TH1F("h1","h1", 322,0,9600);
	TCanvas *c = new TCanvas("canvas", "canvas", 800, 500);
	sprintf(rootname,"%s","C:/Si24/Si22peakcali/Xu/fit.root");
	TFile *fout = new TFile(rootname,"RECREATE");//Êä³öÎÄ¼þ//It's better to define histograms and then define fout, in case of draw bugs.
	//gPad->SetLogx();
	double x;
	double y;
	int i = 0;
	string line;
	
	while ( getline(infile, line) )
	{
		stringstream(line) >> x >> y;
		g->SetPoint(i,x,y);
		i++;
	}
	TF1 *eff=new TF1("eff","exp([0]+[1]*x+exp([2]+[3]*x))",1,8000);
	eff->SetParNames("A","B","C","D");
	eff->SetParameters(-4.5,-2.6,0.03,-0.17);
	//TF1 *eff=new TF1("eff","[0]*x^[1]",1,8000);
	//eff->SetParNames("A","B");
	//eff->SetParameters(0.5,-0.5);
	g->Fit("eff","","",450,1800);
	g->SetMarkerStyle(21);
	g->Draw("AP");
	//h1->Fill(g);
	fout->Write();
	fout->Close();
}