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
void checkmappingT777()
{
	int iroot;
	int i, j;
	int peakch[10],icchlow[10],icchhigh[10],dssdchlow[10],dssdchhigh[10];
	char command1BS[80],command2AS[80],command3BS[80],command1AL[80];
	char command1BL[80],command2AL[80],command2BL[80],command3AL[80],command3BL[80];
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
	TChain* fChain = new TChain("T777");//root输入文件中的tree名
	//sprintf(resultname,"%s%04d%s%04d%s","/media/60A4A303A42DB34/data/runrootfiles/ana_F17_",runstart,"_",runstop,".txt");
	sprintf(resultname,"%s","X:/T777/checkmapping.txt");
	ofstream outfile(resultname,ios::out);
	//sprintf(inrootname,"%s%04d%s","/media/60A4A303A42DB34/data/runrootfiles/ana_F17_",iroot,".root");
	for(iroot=runstart; iroot<=runstop; iroot++)
	{
		sprintf(inrootname,"%s%04d%s","X:/T777/S27calabcd",iroot,".root");
		fChain->Add(inrootname);
	}
	long nentries=fChain->GetEntries();//读事件数
	//sprintf(filename,"%s","/media/60A4A303A42DB34/data/runrootfiles/icch.dat");
	if((iroot>=174&&iroot<=208)||(iroot>=213&&iroot<=214))
	{
		for(i=0;i<8;i++)
		{
			sprintf(command1AL,"%s","TDSSD1AL");
			sprintf(command2AL,"%s","TDSSD2AL");
			sprintf(command3AL,"%s","TDSSD3AL");
			sprintf(command1BL,"%s","TDSSD1BL");
			sprintf(command2BL,"%s","TDSSD2BL");
			sprintf(command3BL,"%s","TDSSD3BL");
		}
		for(i=0;i<16;i++)
		{
			sprintf(command1AL,"%s%d%s%d%s","TDSSD1AL[",15-i,"]>18000&&TDSSD1AL[",15-i,"]<28000");
			puts(command1AL);
			counts[i]=fChain->GetEntries(command1AL);
			outfile<<counts[i]<<endl;
		}
		for(i=16;i<32;i++)
		{
			sprintf(command2AL,"%s%d%s%d%s","TDSSD2AL[",31-i,"]>18000&&TDSSD2AL[",31-i,"]<28000");
			puts(command2AL);
			counts[i]=fChain->GetEntries(command2AL);
			outfile<<counts[i]<<endl;
		}
		for(i=32;i<48;i++)
		{
			sprintf(command3AL,"%s%d%s%d%s","TDSSD3AL[",i-32,"]>18000&&TDSSD3AL[",i-32,"]<28000");
			puts(command3AL);
			counts[i]=fChain->GetEntries(command3AL);
			outfile<<counts[i]<<endl;
		}
		for(i=8*6;i<8*7;i++)
		{
			sprintf(command1BL,"%s%d%s%d%s","TDSSD1BL[",63-i,"]>18000&&TDSSD1BL[",63-i,"]<28000");
			puts(command1BL);
			counts[i]=fChain->GetEntries(command1BL);
			outfile<<counts[i]<<endl;
		}
		for(i=8*7;i<8*8;i++)
		{
			sprintf(command1BL,"%s%d%s%d%s","TDSSD1BL[",i-56,"]>18000&&TDSSD1BL[",i-56,"]<28000");
			puts(command1BL);
			counts[i]=fChain->GetEntries(command1BL);
			outfile<<counts[i]<<endl;
		}
		for(i=8*8;i<8*9;i++)
		{
			sprintf(command2BL,"%s%d%s%d%s","TDSSD2BL[",79-i,"]>18000&&TDSSD2BL[",79-i,"]<28000");
			puts(command2BL);
			counts[i]=fChain->GetEntries(command2BL);
			outfile<<counts[i]<<endl;
		}
		for(i=8*9;i<8*10;i++)
		{
			sprintf(command2BL,"%s%d%s%d%s","TDSSD2BL[",i-72,"]>18000&&TDSSD2BL[",i-72,"]<28000");
			puts(command2BL);
			counts[i]=fChain->GetEntries(command2BL);
			outfile<<counts[i]<<endl;
		}
		for(i=8*10;i<8*11;i++)
		{
			sprintf(command3BL,"%s%d%s%d%s","TDSSD3BL[",i-80,"]>18000&&TDSSD3BL[",i-80,"]<28000");
			puts(command3BL);
			counts[i]=fChain->GetEntries(command3BL);
			outfile<<counts[i]<<endl;
		}
		for(i=8*11;i<8*12;i++)
		{
			sprintf(command3BL,"%s%d%s%d%s","TDSSD3BL[",103-i,"]>18000&&TDSSD3BL[",103-i,"]<28000");
			puts(command3BL);
			counts[i]=fChain->GetEntries(command3BL);
			outfile<<counts[i]<<endl;
		}
	}
	if(iroot>=151&&iroot<=173)
	{
		for(i=0;i<8;i++)
		{
			sprintf(command1BS,"%s","TDSSD1BS");
			sprintf(command2AS,"%s","TDSSD2AS");
			sprintf(command3BS,"%s","TDSSD3BS");
		}
		for(i=8*0;i<8*1;i++)
		{
			sprintf(command1BS,"%s%d%s%d%s","TDSSD1BS[",15-i,"]>18000&&TDSSD1BS[",15-i,"]<28000");
			puts(command1BS);
			counts[i]=fChain->GetEntries(command1BS);
			outfile<<counts[i]<<endl;
		}
		for(i=8*1;i<8*2;i++)
		{
			sprintf(command1BS,"%s%d%s%d%s","TDSSD1BS[",i-8,"]>18000&&TDSSD1BS[",i-8,"]<28000");
			puts(command1BS);
			counts[i]=fChain->GetEntries(command1BS);
			outfile<<counts[i]<<endl;
		}
		for(i=8*2;i<8*3;i++)
		{
			sprintf(command2AS,"%s%d%s%d%s","TDSSD2AS[",31-i,"]>18000&&TDSSD2AS[",31-i,"]<28000");
			puts(command2AS);
			counts[i]=fChain->GetEntries(command2AS);
			outfile<<counts[i]<<endl;
		}
		for(i=8*3;i<8*4;i++)
		{
			sprintf(command2AS,"%s%d%s%d%s","TDSSD2AS[",31-i,"]>18000&&TDSSD2AS[",31-i,"]<28000");
			puts(command2AS);
			counts[i]=fChain->GetEntries(command2AS);
			outfile<<counts[i]<<endl;
		}
		for(i=8*4;i<8*5;i++)
		{
			sprintf(command3BS,"%s%d%s%d%s","TDSSD3BS[",i-32,"]>18000&&TDSSD3BS[",i-32,"]<28000");
			puts(command3BS);
			counts[i]=fChain->GetEntries(command3BS);
			outfile<<counts[i]<<endl;
		}
		for(i=8*5;i<8*6;i++)
		{
			sprintf(command3BS,"%s%d%s%d%s","TDSSD3BS[",55-i,"]>18000&&TDSSD3BS[",55-i,"]<28000");
			puts(command3BS);
			counts[i]=fChain->GetEntries(command3BS);
			outfile<<counts[i]<<endl;
		}
	}
	if(iroot>=215)
	{
		for(i=0;i<8;i++)
		{
			sprintf(command1BS,"%s","TDSSD1BS");
			sprintf(command2AS,"%s","TDSSD2AS");
			sprintf(command3BS,"%s","TDSSD3BS");
			sprintf(command1AL,"%s","TDSSD1AL");
			sprintf(command1BL,"%s","TDSSD1BL");
			sprintf(command2BL,"%s","TDSSD2BL");
			sprintf(command3AL,"%s","TDSSD3AL");
		}
		for(i=8*0;i<8*1;i++)
		{
			sprintf(command1BS,"%s%d%s%d%s","TDSSD1BS[",15-i,"]>18000&&TDSSD1BS[",15-i,"]<28000");
			puts(command1BS);
			counts[i]=fChain->GetEntries(command1BS);
			outfile<<counts[i]<<endl;
		}
		for(i=8*1;i<8*2;i++)
		{
			sprintf(command1BS,"%s%d%s%d%s","TDSSD1BS[",i-8,"]>18000&&TDSSD1BS[",i-8,"]<28000");
			puts(command1BS);
			counts[i]=fChain->GetEntries(command1BS);
			outfile<<counts[i]<<endl;
		}
		for(i=8*2;i<8*3;i++)
		{
			sprintf(command2AS,"%s%d%s%d%s","TDSSD2AS[",31-i,"]>18000&&TDSSD2AS[",31-i,"]<28000");
			puts(command2AS);
			counts[i]=fChain->GetEntries(command2AS);
			outfile<<counts[i]<<endl;
		}
		for(i=8*3;i<8*4;i++)
		{
			sprintf(command2AS,"%s%d%s%d%s","TDSSD2AS[",31-i,"]>18000&&TDSSD2AS[",31-i,"]<28000");
			puts(command2AS);
			counts[i]=fChain->GetEntries(command2AS);
			outfile<<counts[i]<<endl;
		}
		for(i=8*4;i<8*5;i++)
		{
			sprintf(command3BS,"%s%d%s%d%s","TDSSD3BS[",i-32,"]>18000&&TDSSD3BS[",i-32,"]<28000");
			puts(command3BS);
			counts[i]=fChain->GetEntries(command3BS);
			outfile<<counts[i]<<endl;
		}
		for(i=8*5;i<8*6;i++)
		{
			sprintf(command3BS,"%s%d%s%d%s","TDSSD3BS[",55-i,"]>18000&&TDSSD3BS[",55-i,"]<28000");
			puts(command3BS);
			counts[i]=fChain->GetEntries(command3BS);
			outfile<<counts[i]<<endl;
		}
		for(i=8*6;i<8*7;i++)
		{
			sprintf(command1AL,"%s%d%s%d%s","TDSSD1AL[",63-i,"]>18000&&TDSSD1AL[",63-i,"]<28000");
			puts(command1AL);
			counts[i]=fChain->GetEntries(command1AL);
			outfile<<counts[i]<<endl;
		}
		for(i=8*7;i<8*8;i++)
		{
			sprintf(command1AL,"%s%d%s%d%s","TDSSD1AL[",63-i,"]>18000&&TDSSD1AL[",63-i,"]<28000");
			puts(command1AL);
			counts[i]=fChain->GetEntries(command1AL);
			outfile<<counts[i]<<endl;
		}
		for(i=8*8;i<8*9;i++)
		{
			sprintf(command1BL,"%s%d%s%d%s","TDSSD1BL[",79-i,"]>18000&&TDSSD1BL[",79-i,"]<28000");
			puts(command1BL);
			counts[i]=fChain->GetEntries(command1BL);
			outfile<<counts[i]<<endl;
		}
		for(i=8*9;i<8*10;i++)
		{
			sprintf(command1BL,"%s%d%s%d%s","TDSSD1BL[",i-72,"]>18000&&TDSSD1BL[",i-72,"]<28000");
			puts(command1BL);
			counts[i]=fChain->GetEntries(command1BL);
			outfile<<counts[i]<<endl;
		}
		for(i=8*10;i<8*11;i++)
		{
			sprintf(command2BL,"%s%d%s%d%s","TDSSD2BL[",95-i,"]>18000&&TDSSD2BL[",95-i,"]<28000");
			puts(command2BL);
			counts[i]=fChain->GetEntries(command2BL);
			outfile<<counts[i]<<endl;
		}
		for(i=8*11;i<8*12;i++)
		{
			sprintf(command2BL,"%s%d%s%d%s","TDSSD2BL[",i-88,"]>18000&&TDSSD2BL[",i-88,"]<28000");
			puts(command2BL);
			counts[i]=fChain->GetEntries(command2BL);
			outfile<<counts[i]<<endl;
		}
		for(i=8*12;i<8*13;i++)
		{
			sprintf(command3AL,"%s%d%s%d%s","TDSSD3AL[",i-96,"]>18000&&TDSSD3AL[",i-96,"]<28000");
			puts(command3AL);
			counts[i]=fChain->GetEntries(command3AL);
			outfile<<counts[i]<<endl;
		}
		for(i=8*13;i<8*14;i++)
		{
			sprintf(command3AL,"%s%d%s%d%s","TDSSD3AL[",i-96,"]>18000&&TDSSD3AL[",i-96,"]<28000");
			puts(command3AL);
			counts[i]=fChain->GetEntries(command3AL);
			outfile<<counts[i]<<endl;
		}
	}
		//fin->Close();
}