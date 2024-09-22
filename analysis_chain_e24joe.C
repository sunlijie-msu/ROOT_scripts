// Analysis for e24joe
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
void analysis_chain_e24joe()// chain sorted runs, old tree is a101 and new tree is tree_newchain.
{
	time_t start, tim;
	struct tm* at;
	char now[80];
	float speed;

	TChain* tree_chain = new TChain("tree");//root文件中的tree名
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
	cout << "input runstart:  ";
	cin >> runstart;
	cout << "input runstop:  ";
	cin >> runstop;
	for (icalroot = runstart; icalroot <= runstop; icalroot++)
	{
		if (icalroot<=100)
		{
			nentries[icalroot] = tree_chain->GetEntries();//if wrong RN was input, nentries[wrong RN]=last nentries
			continue;
		}
		sprintf(pathname, "%s", "/mnt/daqtesting/protondetector/stagearea/");
		sprintf(calrootname, "%s%s%04d%s", pathname, "run-", icalroot, "_cal.root");
		tree_chain->Add(calrootname);
		cout << calrootname;
		nentries[icalroot] = tree_chain->GetEntries();//unsigned long是最大范围的整数，相当于ULong64_t
		cout << " Total nentries[" << icalroot << "]=" << nentries[icalroot] << endl;
	}

	totalentries = nentries[icalroot - 1];
	long nentriesmax;
	//cout << "input the max number of entries (not longer than 2147483647): " << endl;
	//cin>>nentriesmax;
	nentriesmax = 2147483646;

	//	T888Chain->Draw("T2");//画出某个leaf,T777Chain当成T777
	//重定义变量cal.root里的变量, Firstly, redefine the variables to hold the read values.
	// Declaration of leaf types
	Double_t north_e;
	Double_t north_t;
	Double_t south_e;
	Double_t south_t;
	Double_t mesh_e;
	Double_t mesh_t;
	Double_t mesh_logic_e;
	Double_t mesh_logic_t;

	tree_chain->SetBranchAddress("north_e", &north_e);
	tree_chain->SetBranchAddress("north_t", &north_t);
	tree_chain->SetBranchAddress("south_e", &south_e);
	tree_chain->SetBranchAddress("south_t", &south_t);
	tree_chain->SetBranchAddress("mesh_e", &mesh_e);
	tree_chain->SetBranchAddress("mesh_t", &mesh_t);
	tree_chain->SetBranchAddress("mesh_logic_e", &mesh_logic_e);
	tree_chain->SetBranchAddress("mesh_logic_t", &mesh_logic_t);

	//用SetBranchAddress函数将tree的Branch TOFC与重定义好的变量地址&TOFC联系起来
	//The first parameter is the branch name, and the second is the address of the variable where the branch data is to be placed.

	//float TOF也一样，这些变量用来填新branch，与原始cal-tree中的变量可以混杂使用
	Long64_t nevent;
	long i;

	//definition for private use
	sprintf(anarootname, "%s%s%04d%s%04d%s", pathname, "sum_", runstart, "_", runstop, "_cal.root");//oftenmodify
	TFile* fout = new TFile(anarootname, "RECREATE");//输出文件。It's better to define histograms and then define fout, in case of draw bugs.
	//TTree* tree_newchain = new TTree("tree_newchain", "tree_newchain");//or TTree *T888 = new TTree("T888","Treetitle");
	//e.g. TTree(const char* name, const char* title, Int_t splitlevel = 99);//输出tree
	TH1F* hnorthmesht = new TH1F("hnorthmesht", "Mesh Timestamp - North Germanium Timestamp", 1100, -11000, 11000);
	hnorthmesht->GetXaxis()->SetTitle("North Germanium Mesh Time Difference (ns)");

	TH1F* hsouthmesht = new TH1F("hsouthmesht", "Mesh Timestamp - South Germanium Timestamp", 1100, -11000, 11000);
	hsouthmesht->GetXaxis()->SetTitle("South Germanium Mesh Time Difference (ns)");

	TH1F* hnorthsoutht = new TH1F("hnorthsoutht", "South Germanium Timestamp - North Germanium Timestamp", 1100, -11000, 11000);
	hnorthsoutht->GetXaxis()->SetTitle("North Germanium South Time Difference (ns)");

	TH1F* hnorthe = new TH1F("hnorthe", "North Germanium Energy", 1400, 0, 2800);
	hnorthe->GetXaxis()->SetTitle("North Germanium Energy (keV)");

	TH1F* hsouthe = new TH1F("hsouthe", "South Germanium Energy", 1400, 0, 2800);
	hsouthe->GetXaxis()->SetTitle("South Germanium Energy (keV)");

	TH1F* hmeshe = new TH1F("hmeshe", "Raw Mesh Energy", 3200, 0, 32000);
	hmeshe->GetXaxis()->SetTitle("Raw Mesh Energy (channels)");

	TH2F* hnorthemeshe = new TH2F("hnorthemeshe", "North Ge v. Mesh Signal 10 us Coincidence", 280, 0, 2800, 320, 0, 32000);
	hnorthemeshe->GetYaxis()->SetTitle("Mesh Energy");
	hnorthemeshe->GetXaxis()->SetTitle("North Germanium Energy");

	TH2F* hsouthemeshe = new TH2F("hsouthemeshe", "South Ge v. Mesh Signal 10 us Coincidence", 280, 0, 2800, 320, 0, 32000);
	hsouthemeshe->GetYaxis()->SetTitle("Mesh Energy");
	hsouthemeshe->GetXaxis()->SetTitle("South Germanium Energy");

	for (int i = 0; i < totalentries; i++)
	{
		tree_chain->GetEntry(i);//获得输入root文件的第i个entry相应的Branch变量数据，then the redefined variables could be used.
		if (mesh_e>20 && north_e>20)
		{
			hnorthmesht->Fill(mesh_t - north_t);
		}
		if (mesh_e > 20 && south_e > 20)
		{
			hsouthmesht->Fill(mesh_t - south_t);
		}
		if (north_e > 20 && south_e > 20)
		{
			hnorthsoutht->Fill(south_t - north_t);
		}
		if (mesh_e>0)
		{
			hmeshe->Fill(mesh_e);
		}
		if (north_e>0)
		{
			hnorthe->Fill(north_e);
		}
		if (south_e > 0)
		{
			hsouthe->Fill(south_e);
		}
		if (north_e > 20 && mesh_e > 20)
		{
			hnorthemeshe->Fill(north_e, mesh_e);
		}
		if (south_e > 20 && mesh_e > 20)
		{
			hsouthemeshe->Fill(south_e, mesh_e);
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