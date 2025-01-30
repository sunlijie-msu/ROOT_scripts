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
void Fill_histogram_from_txt_by_Frequency()//read x from a txt file and fill a histo event-by-event and save it as a new root file.
{// input bin by bin, not event by event
	int ii,jj;
	float Channel[500];
	float x,y,xe,ye;
	long Count[500];
	char outrootname[300];
	char rawrootname[300];
	TCanvas *canvas;
	int Entries=300;//modify
	ifstream infile("D:/X/Geant4PKU/Geant4_DSL/Build_Geant4_DSL/out.txt",ios::in);
	sprintf(outrootname,"%s","D:/X/Geant4PKU/Geant4_DSL/Build_Geant4_DSL/out.root");
//	sprintf(rawrootname,"C:/Si24/root_scripts/Xu/S27_pz_new.root");//modify if change path or nuclide
//	TFile *fin = new TFile(rawrootname);//after this statement, you can use any ROOT command for this rootfile
	//unsigned long nentries=T999->GetEntries();//读事件数
//	cout<<rawrootname<<endl;
	//cout<<"  Entries="<<nentries<<endl;
	TFile *fout = new TFile(outrootname,"RECREATE");//输出文件//It's better to define histograms and then define fout, in case of draw bugs.
	canvas=new TCanvas("c1","c1", 900,600);//建立画布
	TH1F * angular_distribution = new TH1F("angular_distribution","angular_distribution",1800,-0.05,179.95);//create a histogram
// 	TH1D *h40pz = (TH1D*)fin->Get("himpl40pz");//read a histogram from another Root file
// 	TH1D *h304pz = (TH1D*)fin->Get("himpl304pz");//read a histogram from another Root file
	string line;
	int ibin = 0;
	//while (getline(infile, line))
	//{
	//	stringstream(line) >> x >> y;
	//	ibin = angular_distribution->FindBin(x);
	//	cout << ibin << ' ' << x << ' ' << y << endl;
	//	angular_distribution->SetBinContent(ibin, y);//SetBinContent(ibin starts from 1, bincontent);
	//	//h1->SetBinError(ibin,xe);//SetBinContent(ibin starts from 1, bincontent);
	//	//		h1->Fill(x);
	//	ibin++;
	//}

	while ( getline(infile, line) )
	{
		stringstream(line) >> x;
		angular_distribution->Fill(x);
	}
	
// 	h40pz->Write();//this histo was copied, not new, this write statement cannot be omitted.
// 	h304pz->Write();//this histo was copied, not new, this write statement cannot be omitted.
	angular_distribution->Write();// this write statement should be omitted.
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
	angular_distribution->GetXaxis()->SetTitle("Angle (deg)");//轴名
	angular_distribution->GetYaxis()->SetTitle("Counts");//轴名
	angular_distribution->GetXaxis()->CenterTitle();//居中
	angular_distribution->GetYaxis()->CenterTitle();//居中
	angular_distribution->Draw();
	canvas->Update();
	cout<<"RMS= "<<setiosflags(ios::fixed)<<setprecision(9)<< angular_distribution->GetRMS();//you can also find the RMS in the rootfile to verify
	fout->Write();
	//fout->Close();
}