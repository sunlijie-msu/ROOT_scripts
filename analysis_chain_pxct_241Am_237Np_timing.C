#include <iostream>
#include <fstream>
#include <math.h>
#include <time.h>
#include <map>
#include <TROOT.h>
#include <TSystem.h>
#include <TStyle.h>
#include <TCutG.h>
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
void analysis_chain_pxct_241Am_237Np_timing()// chain pxct 237Np 59-keV lifetime measurement runs, old tree is tree and new tree is tree2
{
	time_t start, tim;
	struct tm* at;
	char now[80];
	float speed;

	TChain* chain = new TChain("tree");//root文件中的tree名
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
	char txtfilename[150];
	int runstart, runstop;
	//	Tstring rootname;
	cout << "input runstart: may be 79 or 92?  ";
	cin >> runstart;
	cout << "input runstop: may be 91 or 100?  ";
	cin >> runstop;
	for (icalroot = runstart; icalroot <= runstop; icalroot++)
	{
		if ((icalroot >= 96 && icalroot <= 99) || (icalroot == 1215))
		{
			nentries[icalroot] = chain->GetEntries();//if wrong RN was input, nentries[wrong RN]=last nentries
			continue;
		}
		sprintf(pathname, "%s", "F:/e21010/pxct/");
		//sprintf(pathname, "%s", "/mnt/analysis/e21010/sun/");
		if (icalroot <= 91)
			sprintf(calrootname, "%s%s%04d%s", pathname, "run", icalroot, "_LEGe_MSD_241Am_inChamber_window1.5us_CFDdelay_adjusted_cal.root"); // MSD12 and MSD26 both
		if (icalroot >= 92)
			sprintf(calrootname, "%s%s%04d%s", pathname, "run", icalroot, "_LEGe_MSD26_241Am_inChamber_2mmCollimator_window1.5us_CFDdelay_adjusted_cal.root"); // MSD26 only
		chain->Add(calrootname);
		cout << calrootname;
		nentries[icalroot] = chain->GetEntries();//unsigned long是最大范围的整数，相当于ULong64_t
		cout << " Total nentries[" << icalroot << "]=" << nentries[icalroot] << endl;
	}
	for (ilarge = runstop + 1; ilarge < 652; ilarge++)
	{
		nentries[ilarge] = 2100000000;//in case of no appropriate ellipse loop to get in
	}

	icalroot = icalroot - 1;
	totalentries = nentries[icalroot];
	long nentriesmax;
	cout << "input the max number of entries (not longer than 2147483647): " << endl;
	//cin>>nentriesmax;
	nentriesmax = 2100000000;

	//	T888Chain->Draw("T2");//画出某个leaf,T777Chain当成T777
	//重定义变量cal.root里的变量, Firstly, redefine the variables to hold the read values.
   // Declaration of leaf types
	Double_t lege_e;
	Double_t lege_t;
	Double_t msd12_e;
	Double_t msd12_t;
	Double_t msd26_e;
	Double_t msd26_t;

	// List of branches
	TBranch* b_lege_e;   //!
	TBranch* b_lege_t;   //!
	TBranch* b_msd12_e;   //!
	TBranch* b_msd12_t;   //!
	TBranch* b_msd26_e;   //!
	TBranch* b_msd26_t;   //!

	// Set branch addresses and branch pointers
	chain->SetBranchAddress("lege_e", &lege_e, &b_lege_e);
	chain->SetBranchAddress("lege_t", &lege_t, &b_lege_t);
	if (icalroot <= 91) chain->SetBranchAddress("msd12_e", &msd12_e, &b_msd12_e);
	if (icalroot <= 91) chain->SetBranchAddress("msd12_t", &msd12_t, &b_msd12_t);
	chain->SetBranchAddress("msd26_e", &msd26_e, &b_msd26_e);
	chain->SetBranchAddress("msd26_t", &msd26_t, &b_msd26_t);

	//用SetBranchAddress函数将tree的Branch TOFC与重定义好的变量地址&TOFC联系起来
	//The first parameter is the branch name, and the second is the address of the variable where the branch data is to be placed.

	int Ea_central = 0, msd_e_cut_low = 0, msd_e_cut_high = 0, Ea_gate_start = 0, Ea_gate_end = 0;
	if (icalroot <= 91)
	{
		Ea_central = 5469;
		Ea_gate_start = 3;
		Ea_gate_end = 60; // 58 Ea gate choices
	}
	if (icalroot >= 92)
	{
		Ea_central = 5486;
		Ea_gate_start = 3;
		Ea_gate_end = 30; // 28 Ea gate choices
	}
	
	for (int ianaroot = Ea_gate_start; ianaroot <= Ea_gate_end; ianaroot++)
	{
		msd_e_cut_low = Ea_central - ianaroot;
		msd_e_cut_high = Ea_central + ianaroot;
		cout << "msd_e_cut_low = " << msd_e_cut_low << "	msd_e_cut_high = " << msd_e_cut_high << endl;
		sprintf(anarootname, "%s%s%04d%s%04d%s%d%s%d%s", pathname, "sum_", runstart, "_", runstop, "_msd_e_", msd_e_cut_low, "_", msd_e_cut_high, ".root");//modify
		TFile* fout = new TFile(anarootname, "RECREATE");//输出文件。It's better to define histograms and then define fout, in case of draw bugs.
		TTree* tree2 = new TTree("tree2", "tree2");//or TTree *T888 = new TTree("T888","Treetitle");
		//e.g. TTree(const char* name, const char* title, Int_t splitlevel = 99);//输出tree
		TH1D* htiming_lege_msd12 = new TH1D("htiming_lege_msd12", "htiming_lege_msd12", 3000, -1500, 1500); // 3000 channels, 1 ns per channel
		TH1D* htiming_lege_msd26 = new TH1D("htiming_lege_msd26", "htiming_lege_msd26", 3000, -1500, 1500); // 3000 channels, 1 ns per channel

		//e.g. TH2F *hist_name = new TH2F("hist_name","hist_title",num_bins_x,x_low,x_high,num_bins_y,y_low,y_high);
		//e.g. TH1F *hist_name = new TH1F("hist_name","hist_title",num_bins,x_low,x_high);
		//TCanvas *C = new TCanvas("C");

		//e.g. Branch(const char* name, void* address, const char* leaflist(i.e. variable list));
		long i;
		sprintf(txtfilename, "%s%s", pathname, "lege_t-msd12_t-msd26_t.dat");
		ofstream outfile(txtfilename, ios::out);
		for (i = 0; i < totalentries; i++)
		{
			chain->GetEntry(i);//获得输入root文件的第i个entry相应的Branch变量数据，then the redefined variables could be used.
			if (icalroot <= 91)
			{
				if (lege_e > 59 && lege_e < 60.1 && msd12_e > 1700 && msd12_e < 2200 && msd26_e > 3300 && msd26_e < 3800 && msd12_e + msd26_e > msd_e_cut_low && msd12_e + msd26_e < msd_e_cut_high)
				{
					htiming_lege_msd12->Fill(lege_t - msd12_t);
					htiming_lege_msd26->Fill(lege_t - msd26_t);
					//outfile << msd12_e << "	" << lege_t - msd12_t << "	" << msd26_e << "	" << lege_t - msd26_t << endl;
				}
			}
			
			if (icalroot >= 92)
			{
				if (lege_e > 59 && lege_e < 60.1 && msd26_e > msd_e_cut_low && msd26_e < msd_e_cut_high) // for MSD26 only runs
				{
					htiming_lege_msd26->Fill(lege_t - msd26_t);
					//outfile << msd12_e << "	" << lege_t - msd12_t << "	" << msd26_e << "	" << lege_t - msd26_t << endl;
				}
			}
			
			if (i == 0)
			{
				time(&start);
				at = localtime(&start);
				strftime(now, 79, "%Y-%m-%d %H:%M:%S", at);
				cout << now << "	" << ianaroot <<" started." << endl;
			}
			if (i % 10000000 == 0 && i != 0)
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
	} // ianaroot
}