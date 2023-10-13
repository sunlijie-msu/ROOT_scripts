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
void SCA()
{
	int icalroot;
	char calrootname[80];
	char resultname[80];
	char txtname[80];
	int runstart,runstop;
	cout<<"input runstart:";
	cin>>runstart;
	cout<<"input runstop:";
	cin>>runstop;
	for(icalroot=runstart; icalroot<=runstop; icalroot++)
	{
		if(icalroot==624)
			continue;
		sprintf(resultname,"%s%04d%s%04d%s","X:/SCAprint_",runstart,"_",runstop,".txt");
		ofstream outfile(resultname,ios::app);
		sprintf(calrootname,"%s%04d%s","X:/rootfile/data",icalroot,".root");
		TFile *fin = new TFile(calrootname);//after this statement, you can use any ROOT command for this rootfile
		TTree *tree = (TTree*)fin->Get("tree");
		long nentries=tree->GetEntries();//读事件数
		//Firstly, redefine the variables to hold the read values.
		//Float_t SCA[3];
		//tree->SetBranchAddress("SCA", SCA);
		ULong64_t sdc[32];
		Long64_t nevt;
		tree->SetBranchAddress("sdc", sdc);//cannot be omitted
		tree->SetBranchAddress("nevt", &nevt);
		//用SetBranchAddress函数将tree的Branch TOFC与重定义好的变量地址&TOFC联系起来
		tree->GetEntry(nentries-1);
		//获得输入root文件的第某个entry相应的Branch变量数据
		outfile<<"Rn"<<icalroot<<'	'<<nentries<<'	'<<nevt<<'	'<<sdc[0]<<'	'<<sdc[1]<<'	'<<sdc[2]<<endl;//time stamp, trigger, event
		//GetEntries and nevt are the same. slightly different from sdc[2]
		fin->Close();
	}
}