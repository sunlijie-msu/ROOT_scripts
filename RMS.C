#include <iostream>
#include <iomanip.h>
#include <fstream>
#include <math.h>
#include <map>
#include <TROOT.h>
#include "TChain.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TEllipse.h"
#include "TString.h"
#include "TGraph.h"
#include "TSpectrum.h"
#include "TVirtualFitter.h"
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
#include <time.h>
using namespace std;
void RMS()
{
	int ii,jj;
	float Channel[500];
	long Count[500];
	char rootname[80];
	TCanvas *canvas;
	int Entries=14;//modify
	ifstream infile("C:/Si24/Si22peakcali/Xu/inputforRMS.dat",ios::in);
	sprintf(rootname,"%s","C:/Si24/Si22peakcali/Xu/outputforRMS.root");
	TFile *fout = new TFile(rootname,"RECREATE");//输出文件
	canvas=new TCanvas("c1","c1", 900,600);//建立画布
	TH1F *h1 = new TH1F("h1","h1", 200000,-0.5,0.5);//bin影响不大
	for(ii=0;ii<Entries;ii++)
	{
		infile>>Channel[ii]>>Count[ii];
		cout<<' '<<Channel[ii]<<' '<<Count[ii]<<endl;
	}//output for check
	for(ii=0;ii<Entries;ii++)
	{
		for(jj=0;jj<Count[ii];jj++)
		{
			h1->Fill(Channel[ii]);//每个channel即每一个bin的计数Fill Count次，算RMS的情况下，一般Count都=1
		}
	}
	//canvas->cd();//进入画布
	//h1->Draw("same");
	cout<<"RMS= "<<setiosflags(ios::fixed)<<setprecision(9)<<h1->GetRMS();//you can also find the RMS in the rootfile
	fout->Write();
	fout->Close();
}