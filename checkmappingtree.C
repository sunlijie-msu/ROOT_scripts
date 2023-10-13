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
void checkmappingtree()
{
	int iroot;
	int i, j;
	int peakch[10],icchlow[10],icchhigh[10],dssdchlow[10],dssdchhigh[10];
	char command1BS[80],command2AS[80],command3BS[80],command1AL[80];
	char command1BL[80],command2BL[80],command3AL[80];
	char Command[80];
	long counts[128];
	char inrootname[80];
	char resultname[80];
	char filename[80];
	char txtname[80];
	int runstart,runstop;
	cout<<"input runstart:";
	cin>>runstart;
	cout<<"input runstop:";
	cin>>runstop;
	TChain* fChain = new TChain("tree");//root输入文件中的tree名
	//sprintf(resultname,"%s%04d%s%04d%s","/media/60A4A303A42DB34/data/runrootfiles/ana_F17_",runstart,"_",runstop,".txt");
	sprintf(resultname,"%s","X:/T777/checkmappingtree.txt");
	ofstream outfile(resultname,ios::out);
	//sprintf(inrootname,"%s%04d%s","/media/60A4A303A42DB34/data/runrootfiles/ana_F17_",iroot,".root");
	for(iroot=runstart; iroot<=runstop; iroot++)
	{
		sprintf(inrootname,"%s%04d%s","X:/rootfile/data",iroot,".root");
		fChain->Add(inrootname);
	}
	long nentries=fChain->GetEntries();//读事件数
	//sprintf(filename,"%s","/media/60A4A303A42DB34/data/runrootfiles/icch.dat");
// 	if(iroot>=215)
// 	{
// 		for(i=0;i<112;i++)
// 		{
// 			sprintf(Command,"%s%d%s%d%s","gdc[0][",i,"][0]>18000&&gdc[0][",i,"][0]<28000");
// 			puts(Command);
// 			counts[i]=fChain->GetEntries(Command);
// 			outfile<<counts[i]<<endl;
// 		}
// 	}
	if((iroot>=174&&iroot<=208)||(iroot>=213&&iroot<=214))
	{
		for(i=0;i<96;i++)
		{
			sprintf(Command,"%s%d%s%d%s","gdc[0][",i,"][0]>18000&&gdc[0][",i,"][0]<28000");
			puts(Command);
			counts[i]=fChain->GetEntries(Command);
			outfile<<counts[i]<<endl;
		}
	}
// 	if(iroot>=151&&iroot<=173)
// 	{
// 		for(i=0;i<48;i++)
// 		{
// 			sprintf(Command,"%s%d%s%d%s","gdc[0][",i,"][0]>18000&&gdc[0][",i,"][0]<28000");
// 			puts(Command);
// 			counts[i]=fChain->GetEntries(Command);
// 			outfile<<counts[i]<<endl;
// 		}
// 	}
		//fin->Close();
}