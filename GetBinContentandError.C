#include "TH1F.h"
//#include <cmath> //can't use pow() with this header
#include <stdlib.h>
#include "TMinuit.h"
#include "TFumili.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TStyle.h"
#include "TLine.h"
#include <iostream>
#include <fstream>
#include <math.h>
#include <time.h>
#include <map>
#include <TROOT.h>
#include <TSystem.h>
#include <TStyle.h>
#include <TCutG.h>
#include <TChain.h>
#include <TMinuit.h>
#include "TChain.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TEllipse.h"
#include "TString.h"
#include "TGraph.h"
#include "TSpectrum.h"
#include "TVirtualFitter.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "TRandom.h"
#include <TRandom3.h>
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TF1.h"
#include "TMath.h"
#include "TPaveText.h"
#include "TPad.h"
#include "string.h"
#include "stdlib.h"
#include "stdio.h"
#include "TLegend.h"
#include "TPaletteAxis.h"
using namespace std;
void GetBinContentandError()//get BinContent and BinError from DSL gamma-ray spectra data histograms for MADAI observables
{
	int ii,jj,i,j;
	char calrootname[80];
	char hname[80],gname[80];
	char txtname[80];
	int runstart,runstop;
	int binwidth = 2;
	const int factor_rebin=1;
	int Eg[13]={   1248, 2234, 3076, 4970, 5156, 3435, 2186, 2838, 4270, 3541, 5141, 4156, 5293};
int Eg_start[13]={1234, 2222, 3246, 5268, 5450, 3636, 2314, 2998, 4520, 3742, 5464, 4404, 5600};
int Eg_end[13]={1352, 2420, 3326, 5370, 5584, 3716, 2364, 3070, 4610, 3834, 5554, 4494, 5720};
	int down[13]={ 100,   140,   -70,    -154,  -198,   -81,    -40,    -66,   -126,    -71,   -181,  -110, -151};
	TH1F *hFit,*hSega,*hSega_data,*bkgshort,*background;

	sprintf(calrootname,"%s","D:/X/out/DSL/nice.root");//center sega
	TFile *fin_data = new TFile(calrootname);//after this statement, you can use any ROOT command for this rootfile

	for (int Ea=36; Ea<=50; Ea++)
	{
		sprintf(hname,"%s%d%s","D:/X/out/DSL/nice_Ea",Ea,".dat");//alpha gate
		ofstream outfile(hname,ios::out);
		sprintf(hname,"%s%d%s","hGriffinAddback_alpha",Ea,"_t0");//alpha gate
		hSega_data= (TH1F*)fin_data->Get(hname); //Get Gamma ray spectrum
		hSega_data->SetBinErrorOption(TH1::kPoisson);
		for (int i=1; i<=hSega_data->GetNbinsX(); i++)
		{
			outfile<<hSega_data->GetBinCenter(i)<<"	";
			outfile<<hSega_data->GetBinContent(i)<<"	";
			outfile<<hSega_data->GetBinErrorUp(i)<<endl;
		}	
	}


	for (int j=0; j<13; j++)
	{
		sprintf(hname,"%s%d%s","D:/X/out/DSL/S31_Gamma",Eg[j],".dat");//Eg
		ofstream outfile(hname,ios::out);

		sprintf(hname,"%s%d%s","D:/X/out/DSL/S31_Gamma",Eg[j],".root");//Eg
		TFile *fin = new TFile(hname);//after this statement, you can use any ROOT command for this rootfile

		sprintf(hname,"%s","hSega");//gamma spec
		hSega= (TH1F*)fin->Get(hname); //Get Gamma ray spectrum
		hSega->SetBinErrorOption(TH1::kPoisson);
		for (int i=1; i<=hSega->GetNbinsX(); i++)
		{
			outfile<<hSega->GetBinCenter(i)<<"	";
			outfile<<hSega->GetBinContent(i)<<"	";
			outfile<<hSega->GetBinErrorUp(i)<<endl;
		}	
	}

	for (int j=0; j<13; j++)
	{
		sprintf(hname,"%s%d%s","D:/X/out/DSL/experimental_results_MADAI_",Eg[j],".dat");//Eg
		ofstream outfile2(hname,ios::out);

		sprintf(hname,"%s%d%s","D:/X/out/DSL/S31_Gamma",Eg[j],".root");//Eg
		TFile *fin = new TFile(hname);//after this statement, you can use any ROOT command for this rootfile

		sprintf(hname,"%s","hSega");//gamma spec
		hSega= (TH1F*)fin->Get(hname); //Get Gamma ray spectrum
		hSega->Rebin(factor_rebin);
		hSega->SetBinErrorOption(TH1::kPoisson);
		for (int i = (Eg_start[j] - (Eg[j] - down[j])) / binwidth + 1; i<=(Eg_end[j] - (Eg[j] - down[j])) / binwidth; i++)
		{
			outfile2<<hSega->GetBinCenter(i)<<"	";
			outfile2<<hSega->GetBinContent(i)<<"	";
			outfile2<<hSega->GetBinErrorUp(i);
			if (i!=(Eg_end[j] - (Eg[j] - down[j])) / binwidth)
			outfile2<<endl;
		}
	}

		for (int j=0; j<13; j++)
		{
			sprintf(hname,"%s%d%s","D:/X/out/DSL/observable_names_MADAI_",Eg[j],".dat");//Eg
			ofstream outfile3(hname,ios::out);

			sprintf(hname,"%s%d%s","D:/X/out/DSL/S31_Gamma",Eg[j],".root");//Eg
			TFile *fin = new TFile(hname);//after this statement, you can use any ROOT command for this rootfile

			sprintf(hname,"%s","hSega");//gamma spec
			hSega= (TH1F*)fin->Get(hname); //Get Gamma ray spectrum
			hSega->Rebin(factor_rebin);
			hSega->SetBinErrorOption(TH1::kPoisson);
			for (int i = (Eg_start[j] - (Eg[j] - down[j])) / binwidth + 1; i<=(Eg_end[j] - (Eg[j] - down[j])) / binwidth; i++)
			{
				outfile3<<hSega->GetBinCenter(i);
				//outfile3<<hSega->GetBinContent(i)<<"	";
				//outfile3<<hSega->GetBinErrorUp(i);
				if (i!=(Eg_end[j] - (Eg[j] - down[j])) / binwidth)
					outfile3<<endl;
			}
		}

	fin_data->Close();
	fin->Close();
}