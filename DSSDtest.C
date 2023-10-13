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
using namespace std;
void DSSDtest()
{
	time_t tim;
	struct tm *at;
	char now[80];
	char rawrootname[80];
	char calrootname[80];
	for(int irawroot=19;irawroot<=19;irawroot++)
	{
		sprintf(rawrootname,"%s%d%s","C:/Si24/DSSDtest2017/N1068_test_Rn",irawroot,".root");
		sprintf(calrootname,"%s%d%s","C:/Si24/DSSDtest2017/TEST2017_Rn",irawroot,".root");
		TFile *fin = new TFile(rawrootname);
		TTree *T777 = (TTree*)fin->Get("T777");
		int nentries=T777->GetEntries();//读事件数
		cout<<"Entries="<<nentries<<endl;

		//	T777->Draw("T2");//画出某个leaf,fChain当成T777
		Float_t         A16;
		Float_t         A15;
		Float_t         A14;
		Float_t         A13;
		Float_t         A12;
		Float_t         A11;
		Float_t         A10;
		Float_t         A09;
		Float_t         A08;
		Float_t         A07;
		Float_t         A06;
		Float_t         A05;
		Float_t         A04;
		Float_t         A03;
		Float_t         A02;
		Float_t         A01;
		Float_t         B16;
		Float_t         B15;
		Float_t         B14;
		Float_t         B13;
		Float_t         B12;
		Float_t         B11;
		Float_t         B10;
		Float_t         B09;
		Float_t         B01;
		Float_t         B02;
		Float_t         B03;
		Float_t         B04;
		Float_t         B05;
		Float_t         B06;
		Float_t         B07;
		Float_t         B08;
		Float_t         adccnt1;
		Float_t         SCA0;
		Float_t         SCA1;
		TBranch        *b_A16;   //!
		TBranch        *b_A15;   //!
		TBranch        *b_A14;   //!
		TBranch        *b_A13;   //!
		TBranch        *b_A12;   //!
		TBranch        *b_A11;   //!
		TBranch        *b_A10;   //!
		TBranch        *b_A09;   //!
		TBranch        *b_A08;   //!
		TBranch        *b_A07;   //!
		TBranch        *b_A06;   //!
		TBranch        *b_A05;   //!
		TBranch        *b_A04;   //!
		TBranch        *b_A03;   //!
		TBranch        *b_A02;   //!
		TBranch        *b_A01;   //!
		TBranch        *b_B16;   //!
		TBranch        *b_B15;   //!
		TBranch        *b_B14;   //!
		TBranch        *b_B13;   //!
		TBranch        *b_B12;   //!
		TBranch        *b_B11;   //!
		TBranch        *b_B10;   //!
		TBranch        *b_B09;   //!
		TBranch        *b_B01;   //!
		TBranch        *b_B02;   //!
		TBranch        *b_B03;   //!
		TBranch        *b_B04;   //!
		TBranch        *b_B05;   //!
		TBranch        *b_B06;   //!
		TBranch        *b_B07;   //!
		TBranch        *b_B08;   //!
		TBranch        *b_adccnt1;   //!
		TBranch        *b_SCA0;   //!
		TBranch        *b_SCA1;   //!

		T777->SetBranchAddress("A16", &A16, &b_A16);
		T777->SetBranchAddress("A15", &A15, &b_A15);
		T777->SetBranchAddress("A14", &A14, &b_A14);
		T777->SetBranchAddress("A13", &A13, &b_A13);
		T777->SetBranchAddress("A12", &A12, &b_A12);
		T777->SetBranchAddress("A11", &A11, &b_A11);
		T777->SetBranchAddress("A10", &A10, &b_A10);
		T777->SetBranchAddress("A09", &A09, &b_A09);
		T777->SetBranchAddress("A08", &A08, &b_A08);
		T777->SetBranchAddress("A07", &A07, &b_A07);
		T777->SetBranchAddress("A06", &A06, &b_A06);
		T777->SetBranchAddress("A05", &A05, &b_A05);
		T777->SetBranchAddress("A04", &A04, &b_A04);
		T777->SetBranchAddress("A03", &A03, &b_A03);
		T777->SetBranchAddress("A02", &A02, &b_A02);
		T777->SetBranchAddress("A01", &A01, &b_A01);
		T777->SetBranchAddress("B16", &B16, &b_B16);
		T777->SetBranchAddress("B15", &B15, &b_B15);
		T777->SetBranchAddress("B14", &B14, &b_B14);
		T777->SetBranchAddress("B13", &B13, &b_B13);
		T777->SetBranchAddress("B12", &B12, &b_B12);
		T777->SetBranchAddress("B11", &B11, &b_B11);
		T777->SetBranchAddress("B10", &B10, &b_B10);
		T777->SetBranchAddress("B09", &B09, &b_B09);
		T777->SetBranchAddress("B01", &B01, &b_B01);
		T777->SetBranchAddress("B02", &B02, &b_B02);
		T777->SetBranchAddress("B03", &B03, &b_B03);
		T777->SetBranchAddress("B04", &B04, &b_B04);
		T777->SetBranchAddress("B05", &B05, &b_B05);
		T777->SetBranchAddress("B06", &B06, &b_B06);
		T777->SetBranchAddress("B07", &B07, &b_B07);
		T777->SetBranchAddress("B08", &B08, &b_B08);
		T777->SetBranchAddress("adccnt1", &adccnt1, &b_adccnt1);
		T777->SetBranchAddress("SCA0", &SCA0, &b_SCA0);
		T777->SetBranchAddress("SCA1", &SCA1, &b_SCA1);
		
		//用SetBranchAddress函数将tree的Branch T2与重定义好的变量地址&T2联系起来

		//	Float_t TOF;//float TOF也一样，这些变量与原始tree中的变量T2可以混杂使用
		//	Float_t DE;//存刻度后的变量数值
		Float_t DSSD2A[16];//存刻度后的变量数值
		Float_t DSSD2B[16];//存刻度后的变量数值
		Float_t cnt;

		int i,j,ii,jj;//自用变量
		TFile *fout = new TFile(calrootname,"RECREATE");//输出文件
		TTree *tree = new TTree("tree","tree");//or TTree *T777 = new TTree("T777","Treetitle");
		//e.g. TTree(const char* name, const char* title, Int_t splitlevel = 99);//输出tree
		//	TH2F *detof = new TH2F("detof","DE-TOF",120,203,233,540,335,470);//定义输出一二维谱
		//	TH1F *ES1L = new TH1F("ES1L","Energy",500,0,10000);//ES1L
		//e.g. TH2F *hist_name = new TH2F("hist_name","hist_title",num_bins_x,x_low,x_high,num_bins_y,y_low,y_high);
		//e.g. TH1F *hist_name = new TH1F("hist_name","hist_title",num_bins,x_low,x_high);
		TH2I *hposiDSSD = new TH2I("hposiDSSD","hposiDSSD",16,0,16,16,0,16);//create a histogram
		tree->Branch("DSSD2A",DSSD2A,"DSSD2A[16]/F");
		tree->Branch("DSSD2B",DSSD2B,"DSSD2B[16]/F");
		tree->Branch("cnt",&cnt,"cnt/F");

		//e.g. Branch(const char* name, void* address, const char* leaflist(i.e. variable list));
		for(i=0;i<nentries;i++)
		{
			memset(DSSD2A,0,sizeof(DSSD2A));
			memset(DSSD2B,0,sizeof(DSSD2B));
			//变量初始化

			T777->GetEntry(i);//获得输入root文件的第i个entry相应的Branch变量数据
			DSSD2A[0]=A01;
			DSSD2A[1]=A02;
			DSSD2A[2]=A03;
			DSSD2A[3]=A04;
			DSSD2A[4]=A05;
			DSSD2A[5]=A06;
			DSSD2A[6]=A07;
			DSSD2A[7]=A08;
			DSSD2A[8]=A09;
			DSSD2A[9]=A10;
			DSSD2A[10]=A11;
			DSSD2A[11]=A12;
			DSSD2A[12]=A13;
			DSSD2A[13]=A14;
			DSSD2A[14]=A15;
			DSSD2A[15]=A16;
			DSSD2B[0]=B01;
			DSSD2B[1]=B02;
			DSSD2B[2]=B03;
			DSSD2B[3]=B04;
			DSSD2B[4]=B05;
			DSSD2B[5]=B06;
			DSSD2B[6]=B07;
			DSSD2B[7]=B08;
			DSSD2B[8]=B09;
			DSSD2B[9]=B10;
			DSSD2B[10]=B11;
			DSSD2B[11]=B12;
			DSSD2B[12]=B13;
			DSSD2B[13]=B14;
			DSSD2B[14]=B15;
			DSSD2B[15]=B16;
			for(ii=0;ii<16;ii++)
			{
				for(jj=0;jj<16;jj++)
				{
					if(DSSD2A[ii]>800&&DSSD2B[jj]>600&&DSSD2A[ii]<3000&&DSSD2B[jj]<3000)
					{
						hposiDSSD->Fill(jj,ii);
					}
				}
			}//for(ii=0;ii<16;ii++)
			
			tree->Fill();//逐事件填充tree，各Branch统一填充
			if(i%100000==0)
			{
				time(&tim);
				at=localtime(&tim);
				strftime(now,79,"%Y-%m-%d %H:%M:%S",at);
				cout<<now;
				printf(" complete %.1f%s\n",100*(float)i/(float)nentries,"%");
			}
			//		if(i+1>=1000)break;
		}//for(i=0;i<nentries;i++)
		hposiDSSD->GetXaxis()->SetTitle("DSSDB");
		hposiDSSD->GetYaxis()->SetTitle("DSSDA");
		hposiDSSD->GetXaxis()->CenterTitle();
		hposiDSSD->GetYaxis()->CenterTitle();
		hposiDSSD->GetYaxis()->SetTitleOffset(1);
		hposiDSSD->Draw("colz");
		//	fout.cd();//写入文件，只开了一个fout，该cd命令可以省略
		//tree->Write();//将tree统一写入root文件，root文件中有多个T888时，以第一个为准
		fout->Write();//等效于把所有的tree都写入文件。file Write，tree Write保留一个就行，因为file里只有一个tree，两个都写生成的root文件会大一点，而且里面有两个T888。
		fout->Close();//关文件指针最好不要省，不然首次打开root文件时会有warning
		fin->Close();
	}//for(int irawroot=7;irawroot<=15;irawroot++)
}//T777->tree