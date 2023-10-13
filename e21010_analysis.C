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
void e21010_analysis()
{
	time_t start, tim;
	struct tm* at;
	char now[80];
	float speed;

	TChain* a101Chain = new TChain("a101");//root文件中的tree名
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
			nentries[icalroot] = a101Chain->GetEntries();//if wrong RN was input, nentries[wrong RN]=last nentries
			continue;
		}
		//sprintf(pathname, "%s", "F:/e21010/sorted/");
		sprintf(pathname, "%s", "/mnt/analysis/e21010/sun/");
		sprintf(calrootname, "%s%s%03d%s", pathname, "run_", icalroot, ".root");
		a101Chain->Add(calrootname);
		cout << calrootname;
		nentries[icalroot] = a101Chain->GetEntries();//unsigned long是最大范围的整数，相当于ULong64_t
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
	UInt_t          mul;
	UInt_t          FI[66];   //[mul]
	UInt_t          FT[66];   //[mul]
	Double_t        FE[66];   //[mul]
	UInt_t          BI[66];   //[mul]
	UInt_t          BT[66];   //[mul]
	Double_t        BE[66];   //[mul]
	Double_t        theta[66];   //[mul]
	Double_t        phi[66];   //[mul]
	UInt_t          id[66];   //[mul]
	UInt_t          CLOCK;
	UInt_t          TRIGGER;

	// List of branches
	TBranch* b_mul;   //!
	TBranch* b_FI;   //!
	TBranch* b_FT;   //!
	TBranch* b_FE;   //!
	TBranch* b_BI;   //!
	TBranch* b_BT;   //!
	TBranch* b_BE;   //!
	TBranch* b_theta;   //!
	TBranch* b_phi;   //!
	TBranch* b_id;   //!
	TBranch* b_CLOCK;   //!
	TBranch* b_TRIGGER;   //!

	a101Chain->SetBranchAddress("mul", &mul, &b_mul);
	a101Chain->SetBranchAddress("FI", FI, &b_FI);
	a101Chain->SetBranchAddress("FT", FT, &b_FT);
	a101Chain->SetBranchAddress("FE", FE, &b_FE);
	a101Chain->SetBranchAddress("BI", BI, &b_BI);
	a101Chain->SetBranchAddress("BT", BT, &b_BT);
	a101Chain->SetBranchAddress("BE", BE, &b_BE);
	a101Chain->SetBranchAddress("theta", theta, &b_theta);
	a101Chain->SetBranchAddress("phi", phi, &b_phi);
	a101Chain->SetBranchAddress("id", id, &b_id);
	a101Chain->SetBranchAddress("CLOCK", &CLOCK, &b_CLOCK);
	a101Chain->SetBranchAddress("TRIGGER", &TRIGGER, &b_TRIGGER);

	//用SetBranchAddress函数将tree的Branch TOFC与重定义好的变量地址&TOFC联系起来
	//The first parameter is the branch name, and the second is the address of the variable where the branch data is to be placed.

	//float TOF也一样，这些变量用来填新branch，与原始cal-tree中的变量可以混杂使用
	Long64_t nevent;
	long i;
	char h_name[50];

	//definition for private use
	sprintf(anarootname, "%s%s%03d%s%03d%s", pathname, "sum_", runstart, "_", runstop, "_test.root");//oftenmodify
	TFile* fout = new TFile(anarootname, "RECREATE");//输出文件。It's better to define histograms and then define fout, in case of draw bugs.
	TTree* b101 = new TTree("b101", "b101");//or TTree *T888 = new TTree("T888","Treetitle");
	//e.g. TTree(const char* name, const char* title, Int_t splitlevel = 99);//输出tree
	TH1F* hG[2];
	TH1F* hbG[2];
	TH1F* hpG[2];

	for (int ii = 0; ii < 2; ii++)
	{
		sprintf(h_name, "%s%d", "hG", ii);
		hG[ii] = new TH1F(h_name, h_name, 8800, 100, 4500);//name sequence decide the TH1F sequence in rootfile
	}
	for (int ii = 0; ii < 2; ii++)
	{
		sprintf(h_name, "%s%d", "hbG", ii);
		hbG[ii] = new TH1F(h_name, h_name, 8800, 100, 4500);//name sequence decide the TH1F sequence in rootfile
	}
	for (int ii = 0; ii < 2; ii++)
	{
		sprintf(h_name, "%s%d", "hpG", ii);
		hpG[ii] = new TH1F(h_name, h_name, 8800, 100, 4500);//name sequence decide the TH1F sequence in rootfile
	}

	//e.g. TH2F *hist_name = new TH2F("hist_name","hist_title",num_bins_x,x_low,x_high,num_bins_y,y_low,y_high);
	//e.g. TH1F *hist_name = new TH1F("hist_name","hist_title",num_bins,x_low,x_high);
	//TCanvas *C = new TCanvas("C");

	//e.g. Branch(const char* name, void* address, const char* leaflist(i.e. variable list));
	//设置Branch带数组leaf, With this Branch method, you can also add a leaf that holds an entire array of variables.

	for (i = 0; i < totalentries; i++)
	{
		nevent = i;
		a101Chain->GetEntry(i);//获得输入root文件的第i个entry相应的Branch变量数据，then the redefined variables could be used.

		if (mul > 0 && mul < 66)
			for (int ii = 0; ii < mul; ii++) //mul is the number of hits (recorded channels) in each event
			{
				if (FE[ii] > 0 && id[ii] == 12)	hG[0]->Fill(FE[ii]);
				if (FE[ii] > 0 && id[ii] == 13)	hG[1]->Fill(FE[ii]);
				if (i % 10000 == 0)
				{
					cout << "Event: " << i << " Hit: " << ii << " FE: " << FE[ii] << " id: " << id[ii] << endl;
				}
			}

		if (mul >= 2 && mul < 66)
			for (int ii = 0; ii < mul; ii++)
			{
				if (FE[ii] > 0 && id[ii] == 12)	hbG[0]->Fill(FE[ii]);
				if (FE[ii] > 0 && id[ii] == 13)	hbG[1]->Fill(FE[ii]);
			}

		if (mul >= 3 && mul < 66)
			for (int ii = 0; ii < mul; ii++)
			{
				if (FE[ii] > 0 && id[ii] == 12)	hpG[0]->Fill(FE[ii]);
				if (FE[ii] > 0 && id[ii] == 13)	hpG[1]->Fill(FE[ii]);
			}

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