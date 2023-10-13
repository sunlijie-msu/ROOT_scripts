#include <iostream>
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
#include<TChain.h>
using namespace std;
void Branch()//almost useless
{
	char D_name[30],D_command[30];
	char anarootname[80],filename[80];
	float peakch[10],chlow[10],chhigh[10];
	int peakarea[10];
	int runstart=636,runstop=650,ianaroot;
	ofstream outfile("C:/Si24/Si22peakcali/Branch.dat",ios::out);
	int iT=1;
	while (iT<=7)
	{
		outfile<<'\n'<<"iT="<<iT;
		for(ianaroot=runstart; ianaroot<=runstop; ianaroot++)
		{
			sprintf(anarootname,"%s%04d%s%04d%s","V:/RIBLL2015/data24/Si24calabcd_",ianaroot,"_",ianaroot,".root");
			//sprintf(anarootname,"%s%04d%s%04d%s","V:/RIBLL2015/data24/Si24calabcd_",runstart,"_",runstop,".root");
			TFile *fin = new TFile(anarootname);//after this statement, you can use any ROOT command for this rootfile
			TTree *T999 = (TTree*)fin->Get("T999");
			unsigned long nentries=T999->GetEntries();
			cout<<anarootname;
			cout<<"  Entries="<<nentries<<endl;
			//cout<<"input DSSD for branching ratio: ";
			//cin>>D_name;
			sprintf(D_name,"%s","D60Bne");//change point
			//cout<<D_name;
			while(strcmp(D_name,"D300Ane")!=0&&strcmp(D_name,"D60Bne")!=0)
			{
				cout<<"Wrong D_name!"<<endl;
				cout<<"input D_name for branching ratio again: ";
				cin>>D_name;
			}

			Float_t D300Ane[16];
			Float_t D60Bne[16];
			Float_t TAB300;
			Float_t TAB60;
			T999->SetBranchAddress("D300Ane", D300Ane);
			T999->SetBranchAddress("D60Bne", D60Bne);
			T999->SetBranchAddress("TAB300", &TAB300);
			T999->SetBranchAddress("TAB60", &TAB60);
			sprintf(filename,"%s%s%s","C:/Si24/Si22peakcali/Branch",D_name,"_peakch.dat");
			ifstream infile(filename,ios::in);
			for(int ipeak=0;ipeak<10;ipeak++)
			{
				infile>>peakch[ipeak]>>chlow[ipeak]>>chhigh[ipeak];
				//cout<<peakch[ipeak]<<' '<<chlow[ipeak]<<' '<<chhigh[ipeak]<<endl;
			}
			outfile<<'\n';
			for(int ipeak=0;ipeak<10;ipeak++)
			{
				sprintf(D_command,"%s%s%f%s%s%s%f%s%d",D_name,">",chlow[ipeak],"&&",D_name,"<",chhigh[ipeak],"&&TAB60<",70*iT);//change point
				puts(D_command);
				peakarea[ipeak]=T999->GetEntries(D_command);
				outfile<<peakarea[ipeak]<<'	';
			}
		}//ianaroot
		iT++;
	}//while (iT<=7)
}//Tf main