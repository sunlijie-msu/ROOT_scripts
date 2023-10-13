#include <iostream>
#include <fstream>
#include <iomanip.h>
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
using namespace std;
void Fill1D()// read data from a txt file and fill a TH1F
{
	int Eg=2234;//DSLmodify  1248, 2234, 3076, 4971, 5156
	int Ea=47;//DSLmodify         49,     47,      46,     42,     42
	int Egfirst,Egbinwidth;
	int Taufirst,Taubinwidth;
	double Nbinsx_Tau, Xlow, Xup, Nbinsy_Eg, Ylow, Yup;
	double input_value[80][80];
	double Sum_input_X[80]={0};
	double Sum_input_Y[80]={0};
	int i,j,k;
	double c;
	char h_name[300],b_name[300],rootname[300],datname[300];
	sprintf(h_name,"%s%d%s%d%s","Prob_Ea",Ea,"Eg",Eg,"_1DX");
	sprintf(b_name,"%s%s","D:/X/out/DSL/",h_name);
	sprintf(rootname,"%s%s",b_name,".root");
	sprintf(datname,"%s%s",b_name,".dat");
	TFile *fout = new TFile(rootname,"RECREATE");
	if (Eg==1248){ Nbinsx_Tau=27, Xlow=100, Xup=5500, Nbinsy_Eg=17, Ylow=1231, Yup=1265; Taufirst=200; Egfirst=1232; Taubinwidth=200; Egbinwidth=2; }//Eg=1232->1264, 32 keV, 17 bins //low and up should be half-bin broader in order to make bin display correctly.
	if (Eg==2234){ Nbinsx_Tau=43, Xlow=75, Xup=505, Nbinsy_Eg=17, Ylow=2216, Yup=2250; Taufirst=80; Egfirst=2217; Taubinwidth=10; Egbinwidth=2; }//Eg=2217->2249, 32 keV, 17 bins
	if (Eg==3076){ Nbinsx_Tau=25, Xlow=-0.5, Xup=24.5, Nbinsy_Eg=17, Ylow=3058, Yup=3092; Taufirst=0; Egfirst=3059; Taubinwidth=1; Egbinwidth=2; }//Eg=3059->3091, 32 keV, 17 bins
	if (Eg==4971){ Nbinsx_Tau=25, Xlow=-0.5, Xup=24.5, Nbinsy_Eg=17, Ylow=4951, Yup=4985; Taufirst=0; Egfirst=4952; Taubinwidth=1; Egbinwidth=2; }//Eg=4952->4984, 32 keV, 17 bins
	if (Eg==5156){ Nbinsx_Tau=25, Xlow=-0.5, Xup=24.5, Nbinsy_Eg=17, Ylow=5137, Yup=5171; Taufirst=0; Egfirst=5138; Taubinwidth=1; Egbinwidth=2; }//Eg=5138->5170, 32 keV, 17 bins
	TH2F *hChisquare2D = new TH2F(h_name,h_name, Nbinsx_Tau, Xlow, Xup, Nbinsy_Eg, Ylow, Yup);//create and name a histogram
	TH1F *hChisquare1DX = new TH1F("Tau","Tau", Nbinsx_Tau, Xlow, Xup);//create and name a histogram
	TH1F *hChisquare1DY = new TH1F("Eg","Eg", Nbinsy_Eg, Ylow, Yup);//create and name a histogram
	//TH2F (const char *name, const char *title, Int_t nbinsx, Double_t xlow, Double_t xup, Int_t nbinsy, Double_t ylow, Double_t yup)
	fstream InputByChar_infile;
	InputByChar_infile.open(datname,ios::in);
	while(!InputByChar_infile.eof())
	{
		for(i=0;i<Nbinsx_Tau;i++)//i is row, j is column
				InputByChar_infile>>Sum_input_X[i];
	}

	for(i=0;i<Nbinsx_Tau;i++)//i is row, j is column
	{
		for (k=0;k<Sum_input_X[i];k++)
		{
			hChisquare1DX->Fill(Taufirst+i*Taubinwidth);//Fill histograms z by event
		}
	}

	for(i=0;i<Nbinsx_Tau;i++)
	{
		cout<<Taufirst+i*Taubinwidth<<"	"<<Sum_input_X[i]<<endl;
	}
	TCanvas* c1   = new TCanvas("c1","c1",1000,700);
	c1->cd();

	hChisquare1DX->Draw("");

	c1->Update();
	//hChisquare2D->Write();
	fout->Write();
	//fout->Close();
}//Fill2D main