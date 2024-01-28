#include <iostream>
#include <fstream>
#include <iomanip>
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
void Fill_histogram_from_txt_byBinContent()//read bincontent from a txt file and create a histo and save it as a new root file.
{// input bin by bin, not event by event
	int ii,jj;
	float Channel[500];
	float x,y,xe,ye;
	long Count[500];
	char rootname[80];
	char rawrootname[80];
	TCanvas *canvas;
	int Entries=300;//modify
	ifstream infile("C:/Si24/root_scripts/Xu/inputforh1.dat",ios::in);
	sprintf(rootname,"%s","C:/Si24/root_scripts/Xu/outputforh1.root");
//	sprintf(rawrootname,"C:/Si24/root_scripts/Xu/S27_pz_new.root");//modify if change path or nuclide
//	TFile *fin = new TFile(rawrootname);//after this statement, you can use any ROOT command for this rootfile
	//unsigned long nentries=T999->GetEntries();//读事件数
//	cout<<rawrootname<<endl;
	//cout<<"  Entries="<<nentries<<endl;
	TFile *fout = new TFile(rootname,"RECREATE");//输出文件//It's better to define histograms and then define fout, in case of draw bugs.
	canvas=new TCanvas("c1","c1", 900,600);//建立画布
	TH1F *h1 = new TH1F("depth","depth",100,-1.5,298.5);//create a histogram
// 	TH1D *h40pz = (TH1D*)fin->Get("himpl40pz");//read a histogram from another Root file
// 	TH1D *h304pz = (TH1D*)fin->Get("himpl304pz");//read a histogram from another Root file
	string line;
	int ibin = 0;
// 	while ( getline(infile, line) )
// 	{
// 		stringstream(line) >> x >> y;
// 		ibin=h1->FindBin(x);
// 		h1->SetBinContent(ibin,y);//SetBinContent(ibin starts from 1, bincontent);
// //		h1->Fill(x);
// 		ibin++;
// 	}

	while ( getline(infile, line) )
	{
		stringstream(line) >> x >> y;
		ibin=h1->FindBin(x);
		cout<<ibin<<' '<<x<<' '<<y<<endl;
		h1->SetBinContent(ibin,y);//SetBinContent(ibin starts from 1, bincontent);
		//h1->SetBinError(ibin,xe);//SetBinContent(ibin starts from 1, bincontent);
		//		h1->Fill(x);
		ibin++;
	}
	
// 	h40pz->Write();//this histo was copied, not new, this write statement cannot be omitted.
// 	h304pz->Write();//this histo was copied, not new, this write statement cannot be omitted.
	h1->Write();// this write statement should be omitted.
// 	for(ii=0;ii<Entries;ii++)
// 	{
// 		infile>>Channel[ii]>>Count[ii];
// 		cout<<' '<<Channel[ii]<<' '<<Count[ii]<<endl;
// 	}//output for check
// 	for(ii=0;ii<Entries;ii++)
// 	{
// 		for(jj=0;jj<Count[ii];jj++)
// 		{
// 			h1->Fill(Channel[ii]);//每个channel即每一个bin的计数Fill Count次，算RMS的情况下，一般Count都=1
// 		}
// 	}
	canvas->cd();//进入画布
	h1->GetXaxis()->SetTitle("Depth (nm)");//轴名
	h1->GetYaxis()->SetTitle("Counts per 3 nm");//轴名
	h1->GetXaxis()->CenterTitle();//居中
	h1->GetYaxis()->CenterTitle();//居中
	h1->Draw();
	canvas->Update();
	cout<<"RMS= "<<setiosflags(ios::fixed)<<setprecision(9)<<h1->GetRMS();//you can also find the RMS in the rootfile
	fout->Write();
	fout->Close();
}