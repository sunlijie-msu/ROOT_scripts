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
void analysis_chain_pxct_241Am_237Np_timing_Run2()
// chain pxct 237Np 59-keV lifetime measurement runs. Later I hadded the cal.root files so chain is only for 1 0000_cal.root file now.
// old tree is tree and new tree is tree2. tree2 is empty though.
// To generate many Run1_timing_msdtotal_e_5361_5481_msdtotal_t.root files corresponding to α gates, each contains a htiming_lege_msd12 and 1 htiming_lege_msd26 gated by the α energy gate. Bin counts are also output to timing_msdtotal_e_5361_5481_msd26_t_bin1ns.csv files for later emcee fit.
// Upstream code: /user/pxct/readout/rootfile/run0079_0091_ddas2root.C
// Downstream code: peakfit_expdecay_band_lifetime_pxct_241Am_237Np.C
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
	char calrootname[500];
	char anarootname[500];
	char pathname[500];
	char filename[500];
	char txtfilename[500];
	int runstart = 0, runstop = 0;
	//	Tstring rootname;
// 	cout << "input runstart: Type 0:  ";
// 	cin >> runstart;
// 	cout << "input runstop: Type 191 for run0092-0095, 0100-0105; Type 0 for run0079-0091 or run0330-0332:  ";
// 	cin >> runstop;
	int Which_Dataset = 2; // 1 for MSDtotal; 2 for MSD26; 3 for MSDtotal data in May 2024
	if (Which_Dataset == 1) sprintf(filename, "%s", "run0079_0091_LEGe_MSD_241Am_inChamber_window1.5us_CFDdelay_adjusted"); // input root
	if (Which_Dataset == 2) sprintf(filename, "%s", "run0092_0095_0100_0109_LEGe_MSD26_241Am_window1.5us"); // input root
	if (Which_Dataset == 3) sprintf(filename, "%s", "run0330_0332_LEGe_MSD_241Am_Z7117_ChamberCenter_window1.5us_TrigRise0.064_0.016_0.016us_TrigGap0.952_1.000_1.000us_Th350_2700_1000_CFDDelay0.304us_Scale7"); // input root
	//sprintf(filename, "%s", "run0092_0095_0100_0105_LEGe_MSD26_241Am_inChamber_2mmCollimator_window1.5us_CFDdelay_adjusted"); // input root
	
	for (icalroot = runstart; icalroot <= runstop; icalroot++)
	{
		if ((icalroot >= 10000 && icalroot <= 10000) || (icalroot == 10000))
		{
			nentries[icalroot] = chain->GetEntries();//if wrong RN was input, nentries[wrong RN]=last nentries
			continue;
		}
		
		sprintf(pathname, "%s", "F:/e21010/pxct/");
		sprintf(calrootname, "%s%s%s%04d%s", pathname, filename, "_", icalroot, "_cal.root"); // input root files. Because of gain matching by 4 hours for MSD26 runs, there are hundreds of small files. MSDtotal run is just 1 file with icalroot=0000.
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
	cout << "input the max number of entries (not longer than 2147483647, 32-bit integer): " << endl;
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
// 	TBranch* b_msd12_e;   //!
// 	TBranch* b_msd12_t;   //!
	TBranch* b_msd26_e;   //!
	TBranch* b_msd26_t;   //!

	// Set branch addresses and branch pointers
	chain->SetBranchAddress("lege_e", &lege_e, &b_lege_e);
	chain->SetBranchAddress("lege_t", &lege_t, &b_lege_t);
// 	chain->SetBranchAddress("msd12_e", &msd12_e, &b_msd12_e);
// 	chain->SetBranchAddress("msd12_t", &msd12_t, &b_msd12_t);
	chain->SetBranchAddress("msd26_e", &msd26_e, &b_msd26_e);
	chain->SetBranchAddress("msd26_t", &msd26_t, &b_msd26_t);

	//用SetBranchAddress函数将tree的Branch TOFC与重定义好的变量地址&TOFC联系起来
	//The first parameter is the branch name, and the second is the address of the variable where the branch data is to be placed.

	int Ea_central = 0, msd_e_cut_low = 0, msd_e_cut_high = 0, Ea_gate_start = 0, Ea_gate_end = 0;
	
// 	Ea_central = 5418;
// 	Ea_gate_start = 4;
// 	Ea_gate_end = 60; // Ea gate choices for MSD12+MSD26 runs

	Ea_central = 5479;
	Ea_gate_start = 60;
	Ea_gate_end = 60; // Ea gate choices for MSD26 only runs
	
	for (int ianaroot = Ea_gate_start; ianaroot <= Ea_gate_end; ianaroot+=1)
	{
		msd_e_cut_low = Ea_central - ianaroot;
		msd_e_cut_high = Ea_central + ianaroot;
		msd_e_cut_low = 5460; // for single peak fit
		msd_e_cut_high = 5500; // for single peak fit
		cout << "msd_e_cut_low = " << msd_e_cut_low << "	msd_e_cut_high = " << msd_e_cut_high << endl;
		sprintf(anarootname, "%s%s%d%s%d%s%d%s", pathname, "Run", Which_Dataset, "_timing_msd26_e_", msd_e_cut_low, "_", msd_e_cut_high, "_msd26_t.root");//output root modify

		//sprintf(txtfilename, "%s%s%d%s%d%s%d%s", pathname, "Run", Which_Dataset, "_timing_msdtotal_e_", msd_e_cut_low, "_", msd_e_cut_high, "_msd12_t_bin01ns.csv");//output csv
		//ofstream outfile1(txtfilename, ios::out);
		//sprintf(txtfilename, "%s%s%d%s%d%s%d%s", pathname, "Run", Which_Dataset, "_timing_msdtotal_e_", msd_e_cut_low, "_", msd_e_cut_high, "_msd12_t_bin1ns.csv");//output csv
		//ofstream outfile2(txtfilename, ios::out);
		//sprintf(txtfilename, "%s%s%d%s%d%s%d%s", pathname, "Run", Which_Dataset, "_timing_msdtotal_e_", msd_e_cut_low, "_", msd_e_cut_high, "_msd12_t_bin2ns.csv");//output csv
		//ofstream outfile3(txtfilename, ios::out);
		//sprintf(txtfilename, "%s%s%d%s%d%s%d%s", pathname, "Run", Which_Dataset, "_timing_msdtotal_e_", msd_e_cut_low, "_", msd_e_cut_high, "_msd12_t_bin5ns.csv");//output csv
		//ofstream outfile4(txtfilename, ios::out);

		//sprintf(txtfilename, "%s%s%d%s%d%s%d%s", pathname, "Run", Which_Dataset, "_timing_msdtotal_e_", msd_e_cut_low, "_", msd_e_cut_high, "_msd26_t_bin01ns.csv");//output csv
		//ofstream outfile5(txtfilename, ios::out);
		sprintf(txtfilename, "%s%s%d%s%d%s%d%s", pathname, "Run", Which_Dataset, "_timing_msd26_e_", msd_e_cut_low, "_", msd_e_cut_high, "_msd26_t_bin1ns.csv");//output csv
		ofstream outfile6(txtfilename, ios::out);
		//sprintf(txtfilename, "%s%s%d%s%d%s%d%s", pathname, "Run", Which_Dataset, "_timing_msdtotal_e_", msd_e_cut_low, "_", msd_e_cut_high, "_msd26_t_bin2ns.csv");//output csv
		//ofstream outfile7(txtfilename, ios::out);
		//sprintf(txtfilename, "%s%s%d%s%d%s%d%s", pathname, "Run", Which_Dataset, "_timing_msdtotal_e_", msd_e_cut_low, "_", msd_e_cut_high, "_msd26_t_bin5ns.csv");//output csv
		//ofstream outfile8(txtfilename, ios::out);

		TFile* fout = new TFile(anarootname, "RECREATE");//输出文件。It's better to define histograms and then define fout, in case of draw bugs.
		TTree* tree2 = new TTree("tree2", "tree2");//or TTree *T888 = new TTree("T888","Treetitle");
		//e.g. TTree(const char* name, const char* title, Int_t splitlevel = 99);//输出tree
		TH1D* hlege_e = new TH1D("hlege_e", "hlege_e", 60000, -0.0540888491631, 435.8137733893930); // 60000 channels, 7.3 eV per channel
		//TH1D* hmsd12_e = new TH1D("hmsd12_e", "hmsd12_e", 14000, 0, 7000); // 7000 channels, 0.5 keV per channel
		TH1D* hmsd26_e = new TH1D("hmsd26_e", "hmsd26_e", 14000, 0, 7000); // 7000 channels, 0.5 keV per channel
		//TH1D* hmsdtotal_e = new TH1D("hmsdtotal_e", "hmsdtotal_e", 14000, 0, 7000); // 7000 channels, 0.5 keV per channel

		//TH1D* htiming_lege_msd12_bin01ns = new TH1D("htiming_lege_msd12_bin01ns", "htiming_lege_msd12_bin01ns", 30000, -1500, 1500); // 6000 channels, 0.1 ns per channel
		//TH1D* htiming_lege_msd12_bin1ns = new TH1D("htiming_lege_msd12_bin1ns", "htiming_lege_msd12_bin1ns", 3000, -1500, 1500); // 3000 channels, 1 ns per channel
		//TH1D* htiming_lege_msd12_bin2ns = new TH1D("htiming_lege_msd12_bin2ns", "htiming_lege_msd12_bin2ns", 1500, -1500, 1500); // 1500 channels, 2 ns per channel
		//TH1D* htiming_lege_msd12_bin5ns = new TH1D("htiming_lege_msd12_bin5ns", "htiming_lege_msd12_bin5ns", 600, -1500, 1500); // 600 channels, 5 ns per channel

		//TH1D* htiming_lege_msd26_bin01ns = new TH1D("htiming_lege_msd26_bin01ns", "htiming_lege_msd26_bin01ns", 30000, -1500, 1500); // 6000 channels, 0.1 ns per channel
		TH1D* htiming_lege_msd26_bin1ns = new TH1D("htiming_lege_msd26_bin1ns", "htiming_lege_msd26_bin1ns", 3000, -1500, 1500); // 3000 channels, 1 ns per channel
		//TH1D* htiming_lege_msd26_bin2ns = new TH1D("htiming_lege_msd26_bin2ns", "htiming_lege_msd26_bin2ns", 1500, -1500, 1500); // 1500 channels, 2 ns per channel
		//TH1D* htiming_lege_msd26_bin5ns = new TH1D("htiming_lege_msd26_bin5ns", "htiming_lege_msd26_bin5ns", 600, -1500, 1500); // 600 channels, 5 ns per channel

		//e.g. TH2F *hist_name = new TH2F("hist_name","hist_title",num_bins_x,x_low,x_high,num_bins_y,y_low,y_high);
		//e.g. TH1F *hist_name = new TH1F("hist_name","hist_title",num_bins,x_low,x_high);
		//TCanvas *C = new TCanvas("C");

		//e.g. Branch(const char* name, void* address, const char* leaflist(i.e. variable list));
		long i;
		//sprintf(txtfilename, "%s%s%s", pathname, filename, "_lege_t-msd12_t-msd26_t.dat"); // 
		//ofstream outfile(txtfilename, ios::out);
		for (i = 0; i < totalentries; i++)
		{
			chain->GetEntry(i);//获得输入root文件的第i个entry相应的Branch变量数据，then the redefined variables could be used.

			if (lege_e > 0)	hlege_e->Fill(lege_e);
			//if (msd12_e > 0)	hmsd12_e->Fill(msd12_e + gRandom->Uniform(-0.2, 0.2));
			if (msd26_e > 0)	hmsd26_e->Fill(msd26_e + gRandom->Uniform(-0.2, 0.2));
			//if (msd12_e > 0 && msd26_e > 0)	hmsdtotal_e->Fill(msd12_e + msd26_e);

			//if (lege_e > 59.0 && lege_e < 60.1 && msd12_e > 1600 && msd12_e < 2050 && msd26_e > 3300 && msd26_e < 3800 && msd12_e + msd26_e > msd_e_cut_low && msd12_e + msd26_e < msd_e_cut_high)
			//{
			//	//htiming_lege_msd12_bin01ns->Fill(lege_t - msd12_t);
			//	htiming_lege_msd12_bin1ns->Fill(lege_t - msd12_t);
			//	//htiming_lege_msd12_bin2ns->Fill(lege_t - msd12_t);
			//	//htiming_lege_msd12_bin5ns->Fill(lege_t - msd12_t);

			//	//htiming_lege_msd26_bin01ns->Fill(lege_t - msd26_t);
			//	htiming_lege_msd26_bin1ns->Fill(lege_t - msd26_t);
			//	//htiming_lege_msd26_bin2ns->Fill(lege_t - msd26_t);
			//	//htiming_lege_msd26_bin5ns->Fill(lege_t - msd26_t);

			//	//outfile << msd12_e << "	" << lege_t - msd12_t << "	" << msd26_e << "	" << lege_t - msd26_t << endl;
			//}

			if (lege_e > 59.0 && lege_e < 60.1 && msd26_e > msd_e_cut_low && msd26_e < msd_e_cut_high) // for MSD26 only runs
			{
				//htiming_lege_msd26_bin01ns->Fill(lege_t - msd26_t);
				htiming_lege_msd26_bin1ns->Fill(lege_t - msd26_t);
				//htiming_lege_msd26_bin2ns->Fill(lege_t - msd26_t);
				//htiming_lege_msd26_bin5ns->Fill(lege_t - msd26_t);

				//outfile << msd12_e << "	" << lege_t - msd12_t << "	" << msd26_e << "	" << lege_t - msd26_t << endl;
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

		//htiming_lege_msd12_bin01ns->SetBinErrorOption(TH1::kPoisson);
		//htiming_lege_msd12_bin1ns->SetBinErrorOption(TH1::kPoisson);
		//htiming_lege_msd12_bin2ns->SetBinErrorOption(TH1::kPoisson);
		//htiming_lege_msd12_bin5ns->SetBinErrorOption(TH1::kPoisson);
		//htiming_lege_msd26_bin01ns->SetBinErrorOption(TH1::kPoisson);
		htiming_lege_msd26_bin1ns->SetBinErrorOption(TH1::kPoisson);
		//htiming_lege_msd26_bin2ns->SetBinErrorOption(TH1::kPoisson);
		//htiming_lege_msd26_bin5ns->SetBinErrorOption(TH1::kPoisson);

		//for (int i = 1; i <= htiming_lege_msd12_bin01ns->GetNbinsX(); i++)
		//{
		//	outfile1 << htiming_lege_msd12_bin01ns->GetBinCenter(i) << ",";
		//	outfile1 << htiming_lege_msd12_bin01ns->GetBinContent(i) << ",";
		//	outfile1 << htiming_lege_msd12_bin01ns->GetBinErrorLow(i) << ",";
		//	outfile1 << htiming_lege_msd12_bin01ns->GetBinErrorUp(i) << endl;
		//}
		//for (int i = 1; i <= htiming_lege_msd12_bin1ns->GetNbinsX(); i++)
		//{
		//	outfile2 << htiming_lege_msd12_bin1ns->GetBinCenter(i) << ",";
		//	outfile2 << htiming_lege_msd12_bin1ns->GetBinContent(i) << ",";
		//	outfile2 << htiming_lege_msd12_bin1ns->GetBinErrorLow(i) << ",";
		//	outfile2 << htiming_lege_msd12_bin1ns->GetBinErrorUp(i) << endl;
		//}
		//for (int i = 1; i <= htiming_lege_msd12_bin2ns->GetNbinsX(); i++)
		//{
		//	outfile3 << htiming_lege_msd12_bin2ns->GetBinCenter(i) << ",";
		//	outfile3 << htiming_lege_msd12_bin2ns->GetBinContent(i) << ",";
		//	outfile3 << htiming_lege_msd12_bin2ns->GetBinErrorLow(i) << ",";
		//	outfile3 << htiming_lege_msd12_bin2ns->GetBinErrorUp(i) << endl;
		//}
		//for (int i = 1; i <= htiming_lege_msd12_bin5ns->GetNbinsX(); i++)
		//{
		//	outfile4 << htiming_lege_msd12_bin5ns->GetBinCenter(i) << ",";
		//	outfile4 << htiming_lege_msd12_bin5ns->GetBinContent(i) << ",";
		//	outfile4 << htiming_lege_msd12_bin5ns->GetBinErrorLow(i) << ",";
		//	outfile4 << htiming_lege_msd12_bin5ns->GetBinErrorUp(i) << endl;
		//}

		//for (int i = 1; i <= htiming_lege_msd26_bin01ns->GetNbinsX(); i++)
		//{
		//	outfile5 << htiming_lege_msd26_bin01ns->GetBinCenter(i) << ",";
		//	outfile5 << htiming_lege_msd26_bin01ns->GetBinContent(i) << ",";
		//	outfile5 << htiming_lege_msd26_bin01ns->GetBinErrorLow(i) << ",";
		//	outfile5 << htiming_lege_msd26_bin01ns->GetBinErrorUp(i) << endl;
		//}
		for (int i = 1; i <= htiming_lege_msd26_bin1ns->GetNbinsX(); i++)
		{
			outfile6 << htiming_lege_msd26_bin1ns->GetBinCenter(i) << ",";
			outfile6 << htiming_lege_msd26_bin1ns->GetBinContent(i) << ",";
			outfile6 << htiming_lege_msd26_bin1ns->GetBinErrorLow(i) << ",";
			outfile6 << htiming_lege_msd26_bin1ns->GetBinErrorUp(i) << endl;
		}
		//for (int i = 1; i <= htiming_lege_msd26_bin2ns->GetNbinsX(); i++)
		//{
		//	outfile7 << htiming_lege_msd26_bin2ns->GetBinCenter(i) << ",";
		//	outfile7 << htiming_lege_msd26_bin2ns->GetBinContent(i) << ",";
		//	outfile7 << htiming_lege_msd26_bin2ns->GetBinErrorLow(i) << ",";
		//	outfile7 << htiming_lege_msd26_bin2ns->GetBinErrorUp(i) << endl;
		//}
		//for (int i = 1; i <= htiming_lege_msd26_bin5ns->GetNbinsX(); i++)
		//{
		//	outfile8 << htiming_lege_msd26_bin5ns->GetBinCenter(i) << ",";
		//	outfile8 << htiming_lege_msd26_bin5ns->GetBinContent(i) << ",";
		//	outfile8 << htiming_lege_msd26_bin5ns->GetBinErrorLow(i) << ",";
		//	outfile8 << htiming_lege_msd26_bin5ns->GetBinErrorUp(i) << endl;
		//}
		//outfile1.close();
		//outfile2.close();
		//outfile3.close();
		//outfile4.close();
		//outfile5.close();
		outfile6.close();
		//outfile7.close();
		//outfile8.close();
		fout->Write();//等效于把所有的tree和新一维谱（非copy旧文件的一维谱）都写入文件。file Write，tree Write保留一个就行，因为file里只有一个tree，两个都写生成的root文件会大一点，而且里面有两个b101，画图并无区别。
		fout->Close();//关文件指针最好不要省，不然首次打开root文件时会有warning
	} // ianaroot
}