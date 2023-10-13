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
void copytree()//Copy one branch of of a Tree in a file to a new Tree in a new file.
{
	time_t start,tim;
	struct tm *at;
	char now[80];
	float speed;
		//Get old file, old tree and set top branch address
		TFile *fin = new TFile("X:/T777/Rn559_152Eu_DSSD3.root");
		TTree *T777 = (TTree*)fin->Get("T777");
		long totalentries=T777->GetEntries();//读事件数
		Float_t Clover[20];
		long i;
		T777->SetBranchAddress("Clover",&Clover);
		TFile *fout = new TFile("X:/T777/copysmall.root","RECREATE");//It's better to define histograms and then define fout, in case of draw bugs.
		TTree *T111 = new TTree("T111","T111");
		T111->Branch("Clover",Clover,"Clover[20]/F");
		for(i=0;i<totalentries;i++)
		{
			memset(Clover,0,sizeof(Clover));
			T777->GetEntry(i);
			if(Clover[1]<50||Clover[1]>500)continue;
			Clover[1]=Clover[1]+gRandom->Uniform(-0.8,0.8);//随机数消除谱上的锯齿
			T111->Fill();
			if(i%100000==0&&i!=0)
			{
				time(&tim);
				at=localtime(&tim);
				strftime(now,79,"%Y-%m-%d %H:%M:%S",at);
				cout<<now;
				speed=(float)i/(float)difftime(tim,start);
				printf(" %.1f%s%.0f%s%.0f%s\n",100*(float)i/(float)totalentries,"%. Speed is ",speed," e/s. Still need ",float(totalentries-i)/speed," s.");
			}
		}
		fout->Write();
		fout->Close();
		fin->Close();
}