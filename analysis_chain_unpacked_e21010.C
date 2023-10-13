// Analysis for e21010
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
#include "TRandom.h"
#include <TRandom3.h>
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
#include "TPaletteAxis.h"
using namespace std;
//在包含C++头文件时一般不用后缀。如果用户自己编写头文件，可以用.h为后缀。
void analysis_chain_unpacked_e21010() // chain unpacked runs, old tree is h101 and new tree is h102.
{
	time_t start, tim;
	struct tm* at;
	char now[80];
	float speed;

	TChain* h101Chain = new TChain("h101");//root文件中的tree名
	// 	T888Chain->Add("Ca407.root");//路径+root文件名（+tree名）
	//  cout<<"Entries="<<T888Chain->GetEntries()<<endl;//总事件数
	//在定义字符串变量时不需指定长度，长度随其中的字符串长度而改变。
	long nentries[652];
	long totalentries;
	memset(nentries, 0, sizeof(nentries));
	int icalroot, ilarge;
	char calrootname[150];
	char anarootname[150];
	char pathname[150];
	char filename[150];
	int runstart, runstop;
	//	Tstring rootname;
	cout << "input runstart: may be 85?  ";
	cin >> runstart;
	cout << "input runstop: may be 117?  ";
	cin >> runstop;
	for (icalroot = runstart; icalroot <= runstop; icalroot++)
	{
		if ((icalroot == 98 && icalroot <= 102) || (icalroot == 115))
		{
			nentries[icalroot] = h101Chain->GetEntries();//if wrong RN was input, nentries[wrong RN]=last nentries
			continue;
		}
		//sprintf(pathname, "%s", "F:/e21010/sorted/");
		sprintf(pathname, "%s", "/mnt/analysis/e21010/sun/");
		sprintf(calrootname, "%s%s%03d%s", pathname, "raw_", icalroot, ".root");
		h101Chain->Add(calrootname);
		cout << calrootname;
		nentries[icalroot] = h101Chain->GetEntries();//unsigned long是最大范围的整数，相当于ULong64_t
		cout << " Total nentries[" << icalroot << "]=" << nentries[icalroot] << endl;
	}
	for (ilarge = runstop + 1; ilarge < 652; ilarge++)
	{
		nentries[ilarge] = 500000000;//in case of no appropriate ellipse loop to get in
	}
	// 	for (int i=0;i<700;i++)
	// 	{
	// 		cout<<i<<" "<<nentries[i]<<endl;
	// 	}

	totalentries = nentries[icalroot - 1];
	long nentriesmax;
	cout << "input the max number of entries (not longer than 2147483647): " << endl;
	//cin>>nentriesmax;
	nentriesmax = 500000000;

	//	T888Chain->Draw("T2");//画出某个leaf,T777Chain当成T777
	//重定义变量cal.root里的变量, Firstly, redefine the variables to hold the read values.
   // Declaration of leaf types
	UInt_t          TRIGGER;
	UInt_t          EVENTNO;
	UInt_t          MEVENTNO;
	UInt_t          BD;
	UInt_t          BDI[8];   //[BD]
	UInt_t          BD_N[8];   //[BD]
	UInt_t          AD;
	UInt_t          ADI[8];   //[AD]
	UInt_t          AD_N[8];   //[AD]
	UInt_t          AR;
	UInt_t          ARI[8];   //[AR]
	UInt_t          AR_N[8];   //[AR]
	UInt_t          CLOCK;
	UInt_t          TPAT;
	UInt_t          TLOW;
	UInt_t          THIGH;
	UInt_t          U1F;
	UInt_t          U1FI[16];   //[U1F]
	UInt_t          U1F_E[16];   //[U1F]
	UInt_t          U1B;
	UInt_t          U1BI[16];   //[U1B]
	UInt_t          U1B_E[16];   //[U1B]
	UInt_t          U1B_T[16];   //[U1B]
	UInt_t          U2F;
	UInt_t          U2FI[16];   //[U2F]
	UInt_t          U2F_E[16];   //[U2F]
	UInt_t          U2B;
	UInt_t          U2BI[16];   //[U2B]
	UInt_t          U2B_E[16];   //[U2B]
	UInt_t          U2B_T[16];   //[U2B]
	UInt_t          U3F;
	UInt_t          U3FI[16];   //[U3F]
	UInt_t          U3F_E[16];   //[U3F]
	UInt_t          U3B;
	UInt_t          U3BI[16];   //[U3B]
	UInt_t          U3B_E[16];   //[U3B]
	UInt_t          U3B_T[16];   //[U3B]
	UInt_t          U4F;
	UInt_t          U4FI[16];   //[U4F]
	UInt_t          U4F_E[16];   //[U4F]
	UInt_t          U4B;
	UInt_t          U4BI[16];   //[U4B]
	UInt_t          U4B_E[16];   //[U4B]
	UInt_t          U4B_T[16];   //[U4B]
	UInt_t          U5F;
	UInt_t          U5FI[16];   //[U5F]
	UInt_t          U5F_E[16];   //[U5F]
	UInt_t          U5B;
	UInt_t          U5BI[16];   //[U5B]
	UInt_t          U5B_E[16];   //[U5B]
	UInt_t          U5B_T[16];   //[U5B]
	UInt_t          U6F;
	UInt_t          U6FI[16];   //[U6F]
	UInt_t          U6F_E[16];   //[U6F]
	UInt_t          U6B;
	UInt_t          U6BI[16];   //[U6B]
	UInt_t          U6B_E[16];   //[U6B]
	UInt_t          U6B_T[16];   //[U6B]
	UInt_t          P1E;
	UInt_t          P1T;
	UInt_t          P2E;
	UInt_t          P2T;
	UInt_t          P3E;
	UInt_t          P3T;
	UInt_t          P4E;
	UInt_t          P4T;
	UInt_t          P5E;
	UInt_t          P5T;
	UInt_t          P6E;
	UInt_t          P6T;
	UInt_t          G1E;
	UInt_t          G1T;
	UInt_t          G2E;
	UInt_t          G2T;

	// List of branches
	TBranch* b_TRIGGER;   //!
	TBranch* b_EVENTNO;   //!
	TBranch* b_MEVENTNO;   //!
	TBranch* b_BD;   //!
	TBranch* b_BDI;   //!
	TBranch* b_BD_N;   //!
	TBranch* b_AD;   //!
	TBranch* b_ADI;   //!
	TBranch* b_AD_N;   //!
	TBranch* b_AR;   //!
	TBranch* b_ARI;   //!
	TBranch* b_AR_N;   //!
	TBranch* b_CLOCK;   //!
	TBranch* b_TPAT;   //!
	TBranch* b_TLOW;   //!
	TBranch* b_THIGH;   //!
	TBranch* b_U1F;   //!
	TBranch* b_U1FI;   //!
	TBranch* b_U1F_E;   //!
	TBranch* b_U1B;   //!
	TBranch* b_U1BI;   //!
	TBranch* b_U1B_E;   //!
	TBranch* b_U1B_T;   //!
	TBranch* b_U2F;   //!
	TBranch* b_U2FI;   //!
	TBranch* b_U2F_E;   //!
	TBranch* b_U2B;   //!
	TBranch* b_U2BI;   //!
	TBranch* b_U2B_E;   //!
	TBranch* b_U2B_T;   //!
	TBranch* b_U3F;   //!
	TBranch* b_U3FI;   //!
	TBranch* b_U3F_E;   //!
	TBranch* b_U3B;   //!
	TBranch* b_U3BI;   //!
	TBranch* b_U3B_E;   //!
	TBranch* b_U3B_T;   //!
	TBranch* b_U4F;   //!
	TBranch* b_U4FI;   //!
	TBranch* b_U4F_E;   //!
	TBranch* b_U4B;   //!
	TBranch* b_U4BI;   //!
	TBranch* b_U4B_E;   //!
	TBranch* b_U4B_T;   //!
	TBranch* b_U5F;   //!
	TBranch* b_U5FI;   //!
	TBranch* b_U5F_E;   //!
	TBranch* b_U5B;   //!
	TBranch* b_U5BI;   //!
	TBranch* b_U5B_E;   //!
	TBranch* b_U5B_T;   //!
	TBranch* b_U6F;   //!
	TBranch* b_U6FI;   //!
	TBranch* b_U6F_E;   //!
	TBranch* b_U6B;   //!
	TBranch* b_U6BI;   //!
	TBranch* b_U6B_E;   //!
	TBranch* b_U6B_T;   //!
	TBranch* b_P1E;   //!
	TBranch* b_P1T;   //!
	TBranch* b_P2E;   //!
	TBranch* b_P2T;   //!
	TBranch* b_P3E;   //!
	TBranch* b_P3T;   //!
	TBranch* b_P4E;   //!
	TBranch* b_P4T;   //!
	TBranch* b_P5E;   //!
	TBranch* b_P5T;   //!
	TBranch* b_P6E;   //!
	TBranch* b_P6T;   //!
	TBranch* b_G1E;   //!
	TBranch* b_G1T;   //!
	TBranch* b_G2E;   //!
	TBranch* b_G2T;   //!

	// Set branch addresses and branch pointers
	h101Chain->SetBranchAddress("TRIGGER", &TRIGGER, &b_TRIGGER);
	h101Chain->SetBranchAddress("EVENTNO", &EVENTNO, &b_EVENTNO);
	h101Chain->SetBranchAddress("MEVENTNO", &MEVENTNO, &b_MEVENTNO);
	h101Chain->SetBranchAddress("BD", &BD, &b_BD);
	h101Chain->SetBranchAddress("BDI", BDI, &b_BDI);
	h101Chain->SetBranchAddress("BD_N", BD_N, &b_BD_N);
	h101Chain->SetBranchAddress("AD", &AD, &b_AD);
	h101Chain->SetBranchAddress("ADI", ADI, &b_ADI);
	h101Chain->SetBranchAddress("AD_N", AD_N, &b_AD_N);
	h101Chain->SetBranchAddress("AR", &AR, &b_AR);
	h101Chain->SetBranchAddress("ARI", ARI, &b_ARI);
	h101Chain->SetBranchAddress("AR_N", AR_N, &b_AR_N);
	h101Chain->SetBranchAddress("CLOCK", &CLOCK, &b_CLOCK);
	h101Chain->SetBranchAddress("TPAT", &TPAT, &b_TPAT);
	h101Chain->SetBranchAddress("TLOW", &TLOW, &b_TLOW);
	h101Chain->SetBranchAddress("THIGH", &THIGH, &b_THIGH);
	h101Chain->SetBranchAddress("U1F", &U1F, &b_U1F);
	h101Chain->SetBranchAddress("U1FI", U1FI, &b_U1FI);
	h101Chain->SetBranchAddress("U1F_E", U1F_E, &b_U1F_E);
	h101Chain->SetBranchAddress("U1B", &U1B, &b_U1B);
	h101Chain->SetBranchAddress("U1BI", U1BI, &b_U1BI);
	h101Chain->SetBranchAddress("U1B_E", U1B_E, &b_U1B_E);
	h101Chain->SetBranchAddress("U1B_T", U1B_T, &b_U1B_T);
	h101Chain->SetBranchAddress("U2F", &U2F, &b_U2F);
	h101Chain->SetBranchAddress("U2FI", U2FI, &b_U2FI);
	h101Chain->SetBranchAddress("U2F_E", U2F_E, &b_U2F_E);
	h101Chain->SetBranchAddress("U2B", &U2B, &b_U2B);
	h101Chain->SetBranchAddress("U2BI", U2BI, &b_U2BI);
	h101Chain->SetBranchAddress("U2B_E", U2B_E, &b_U2B_E);
	h101Chain->SetBranchAddress("U2B_T", U2B_T, &b_U2B_T);
	h101Chain->SetBranchAddress("U3F", &U3F, &b_U3F);
	h101Chain->SetBranchAddress("U3FI", U3FI, &b_U3FI);
	h101Chain->SetBranchAddress("U3F_E", U3F_E, &b_U3F_E);
	h101Chain->SetBranchAddress("U3B", &U3B, &b_U3B);
	h101Chain->SetBranchAddress("U3BI", U3BI, &b_U3BI);
	h101Chain->SetBranchAddress("U3B_E", U3B_E, &b_U3B_E);
	h101Chain->SetBranchAddress("U3B_T", U3B_T, &b_U3B_T);
	h101Chain->SetBranchAddress("U4F", &U4F, &b_U4F);
	h101Chain->SetBranchAddress("U4FI", U4FI, &b_U4FI);
	h101Chain->SetBranchAddress("U4F_E", U4F_E, &b_U4F_E);
	h101Chain->SetBranchAddress("U4B", &U4B, &b_U4B);
	h101Chain->SetBranchAddress("U4BI", U4BI, &b_U4BI);
	h101Chain->SetBranchAddress("U4B_E", U4B_E, &b_U4B_E);
	h101Chain->SetBranchAddress("U4B_T", U4B_T, &b_U4B_T);
	h101Chain->SetBranchAddress("U5F", &U5F, &b_U5F);
	h101Chain->SetBranchAddress("U5FI", U5FI, &b_U5FI);
	h101Chain->SetBranchAddress("U5F_E", U5F_E, &b_U5F_E);
	h101Chain->SetBranchAddress("U5B", &U5B, &b_U5B);
	h101Chain->SetBranchAddress("U5BI", U5BI, &b_U5BI);
	h101Chain->SetBranchAddress("U5B_E", U5B_E, &b_U5B_E);
	h101Chain->SetBranchAddress("U5B_T", U5B_T, &b_U5B_T);
	h101Chain->SetBranchAddress("U6F", &U6F, &b_U6F);
	h101Chain->SetBranchAddress("U6FI", U6FI, &b_U6FI);
	h101Chain->SetBranchAddress("U6F_E", U6F_E, &b_U6F_E);
	h101Chain->SetBranchAddress("U6B", &U6B, &b_U6B);
	h101Chain->SetBranchAddress("U6BI", U6BI, &b_U6BI);
	h101Chain->SetBranchAddress("U6B_E", U6B_E, &b_U6B_E);
	h101Chain->SetBranchAddress("U6B_T", U6B_T, &b_U6B_T);
	h101Chain->SetBranchAddress("P1E", &P1E, &b_P1E);
	h101Chain->SetBranchAddress("P1T", &P1T, &b_P1T);
	h101Chain->SetBranchAddress("P2E", &P2E, &b_P2E);
	h101Chain->SetBranchAddress("P2T", &P2T, &b_P2T);
	h101Chain->SetBranchAddress("P3E", &P3E, &b_P3E);
	h101Chain->SetBranchAddress("P3T", &P3T, &b_P3T);
	h101Chain->SetBranchAddress("P4E", &P4E, &b_P4E);
	h101Chain->SetBranchAddress("P4T", &P4T, &b_P4T);
	h101Chain->SetBranchAddress("P5E", &P5E, &b_P5E);
	h101Chain->SetBranchAddress("P5T", &P5T, &b_P5T);
	h101Chain->SetBranchAddress("P6E", &P6E, &b_P6E);
	h101Chain->SetBranchAddress("P6T", &P6T, &b_P6T);
	h101Chain->SetBranchAddress("G1E", &G1E, &b_G1E);
	h101Chain->SetBranchAddress("G1T", &G1T, &b_G1T);
	h101Chain->SetBranchAddress("G2E", &G2E, &b_G2E);
	h101Chain->SetBranchAddress("G2T", &G2T, &b_G2T);

	//用SetBranchAddress函数将tree的Branch TOFC与重定义好的变量地址&TOFC联系起来
	//The first parameter is the branch name, and the second is the address of the variable where the branch data is to be placed.

	//float TOF也一样，这些变量用来填新branch，与原始cal-tree中的变量可以混杂使用
	Long64_t nevent;
	long i;
	char h_name[50];

	//definition for private use
	sprintf(anarootname, "%s%s%03d%s%03d%s", pathname, "sum_", runstart, "_", runstop, "_unpacked.root");//oftenmodify
	TFile* fout = new TFile(anarootname, "RECREATE");//输出文件。It's better to define histograms and then define fout, in case of draw bugs.
	TTree* h102 = new TTree("h102", "h102");//or TTree *T888 = new TTree("T888","Treetitle");
	//e.g. TTree(const char* name, const char* title, Int_t splitlevel = 99);//输出tree
	TH1F* hG[2];

	for (int ii = 0; ii < 2; ii++)
	{
		sprintf(h_name, "%s%d", "hG", ii);
		hG[ii] = new TH1F(h_name, h_name, 4000, 10, 4010);//name sequence decide the TH1F sequence in rootfile
	}

	//e.g. TH2F *hist_name = new TH2F("hist_name","hist_title",num_bins_x,x_low,x_high,num_bins_y,y_low,y_high);
	//e.g. TH1F *hist_name = new TH1F("hist_name","hist_title",num_bins,x_low,x_high);
	//TCanvas *C = new TCanvas("C");

	//e.g. Branch(const char* name, void* address, const char* leaflist(i.e. variable list));
	//设置Branch带数组leaf, With this Branch method, you can also add a leaf that holds an entire array of variables.

	for (i = 0; i < totalentries; i++)
	{
		nevent = i;
		h101Chain->GetEntry(i);//获得输入root文件的第i个entry相应的Branch变量数据，then the redefined variables could be used.
		//smear the bin counts cuz histograms sometimes have jagged features
// 		if (G1E > 0)	hG[0]->Fill(G1E);// + gRandom->Uniform(-0.5, 0.5));
// 		if (G2E > 0)	hG[1]->Fill(G2E);// + gRandom->Uniform(-0.5, 0.5));
		if (G1E > 0)	hG[0]->Fill(G1E * 0.82741082 - 77.642761 + gRandom->Uniform(-1.5, 1.5));// calibration based on 152Eu, 226Ra, and room background
		if (G2E > 0)	hG[1]->Fill(G2E * 0.79025883 - 76.111357 + gRandom->Uniform(-1.5, 1.5));// calibration based on 152Eu, 226Ra, and room background

		if (i == 0)
		{
			time(&start);
			at = localtime(&start);
			strftime(now, 79, "%Y-%m-%d %H:%M:%S", at);
			cout << now << " started." << endl;
		}
		if (i % 600000 == 0 && i != 0)
		{
			time(&tim);
			at = localtime(&tim);
			strftime(now, 79, "%Y-%m-%d %H:%M:%S", at);
			cout << now;
			speed = (float)i / (float)difftime(tim, start);
			printf(" %.1f%s%.0f%s%.0f%s\n", 100 * (float)i / (float)totalentries, "%. Speed is ", speed, " e/s. Still need ", float(totalentries - i) / speed, " s.");
		}
		if (i >= nentriesmax) break; // can be i>= 300 or i>=nentriesmax
	}//for(i=0;i<nentries;i++)

	fout->Write();//等效于把所有的tree和新一维谱（非copy旧文件的一维谱）都写入文件。file Write，tree Write保留一个就行，因为file里只有一个tree，两个都写生成的root文件会大一点，而且里面有两个b101，画图并无区别。
	fout->Close();//关文件指针最好不要省，不然首次打开root文件时会有warning
}