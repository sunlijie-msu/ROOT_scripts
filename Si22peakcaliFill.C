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
void Si22peakcaliFill()
{
	time_t tim;
	struct tm *at;
	char now[80];
	int irawroot;
	char rawrootname[80];
	double slope[300];
	double intercept[300];
	int ID=16;
	int binmin=500,binmax=3900;
	double sigma=12,thresh=0.4;
	float gaplow=40.,gaphigh=40.;//随分辨不同调整
	unsigned long i;
	int jj,ii;
	char paraprint[30],h_name[50],histo_name[50],b_name[50],hfit_name[50],hcali_name[50];
	//string paraprint[30],h_name[80],histo_name[80],b_name[80],hfit_name[80],hcali_name[80];
	TCanvas *canvaspeak[300];
	TCanvas *canvascali[300];
	TH1D *histo[16];//TH1F peak search+gauss fit,creat histograms, 
	//TH1D *h300AL[16] cannot be TH1D *h300AL[ii]= new TH1D(h_name,h_name,4000,0,4000);in loop
	TF1 *g[10];//creat function
	TGraph *graph[300];//TGraph刻度
	float *peakch;
	float peakmin;
	int k;
	float chtemp;
	int peaknum;
	const int nummax=3;//search peaknumber can be many, fitnumber can be limited to 3
	const long iprint=100000;//useless for TH1Dget, only for treeFill
	const long imax=800000;//useless for TH1Dget, only for treeFill
	double par[10][3];//随峰数不同改动
	float energy[3]={1320,1862,2036};
	//float energy[3]={5157,5486,5805};
	//float energy[3]={797,1670,2700};
	//int runstart,runstop;
	//cout<<"input runstart:";
	//cin>>runstart;
	//cout<<"input runstop:";
	//cin>>runstop;
	//for(irawroot=runstart; irawroot<=runstop; irawroot++)
	//{
		//if()continue;
		sprintf(rawrootname,"V:/RIBLL2015/Si22_0205_0285_20.root");
		//sprintf(calirootname,"%s%04d%s","X:/RIBLL2015/cali",irawroot,".root");
		TFile *fin = new TFile(rawrootname);
		TTree *T999 = (TTree*)fin->Get("T999");
		unsigned long nentries=T999->GetEntries();//读事件数
		cout<<rawrootname;
		cout<<"  Entries="<<nentries<<endl;
		Float_t D300ALch[16];//define the variables to hold the read values //△treeFill
		Float_t D300BLch[16];//define the variables to hold the read values //△treeFill
		Float_t D60ALch[16];//define the variables to hold the read values //△treeFill
		Float_t D60BLch[16];//define the variables to hold the read values //△treeFill
		T999->SetBranchAddress("D300ALch", D300ALch);//other channels can be fitted by rename D300ALch //△treeFill
		T999->SetBranchAddress("D300BLch", D300BLch);//other channels can be fitted by rename D300ALch //△treeFill
		T999->SetBranchAddress("D60ALch", D60ALch);//other channels can be fitted by rename D300ALch //△treeFill
		T999->SetBranchAddress("D60BLch", D60BLch);//other channels can be fitted by rename D300ALch //△treeFill
		cout<<"input h_name for cali: ";
		cin>>h_name;
		while(strcmp(h_name,"D300AL")!=0&&strcmp(h_name,"D300BL")!=0&&strcmp(h_name,"D60AL")!=0&&strcmp(h_name,"D60BL")!=0)
		{
			cout<<"Wrong h_name!"<<endl;
			cout<<"input h_name for cali again: ";
			cin>>h_name;
		}
		for(ii=0;ii<16;ii++)
		{
			sprintf(histo_name,"%s%d",h_name,ii);
			histo[ii] = new TH1D(histo_name,histo_name,3900,100,4000);//name a histogram
		}//name sequence decide the TH1F sequence in rootfile, name histograms if use Fill(branch)
		
		for(i=0;i<nentries;i++)
		{
			memset(D300ALch,0,sizeof(D300ALch));
			memset(D300BLch,0,sizeof(D300BLch));
			memset(D60ALch,0,sizeof(D60ALch));
			memset(D60ALch,0,sizeof(D60BLch));
			for(ii=0;ii<16;ii++)
			{
				T999->GetEntry(i);//Read the entry i in the tree (the i th event in the tree)
				//sprintf(b_name,"%s%s%d%s",h_name,"ch[",ii,"]");
				if(strcmp(h_name,"D300AL")==0)
					histo[ii]->Fill(D300ALch[ii]);//Fill Branch for cali
				else if(strcmp(h_name,"D300BL")==0)
					histo[ii]->Fill(D300BLch[ii]);//Fill Branch for cali
				else if(strcmp(h_name,"D60AL")==0)
					histo[ii]->Fill(D60ALch[ii]);//Fill Branch for cali
				else if(strcmp(h_name,"D60BL")==0)
					histo[ii]->Fill(D60BLch[ii]);//Fill Branch for cali
			}//Fill(Float_t), not Fill(char *)
			if(i%iprint==0)
			{
				time(&tim);
				at=localtime(&tim);
				strftime(now,79,"%Y-%m-%d %H:%M:%S",at);
				cout<<now;
				printf(" complete %d ...\n",i);
			}
			if(i+1>=imax)break;
		}//for(i=0;i<nentries;i++) △treeFill
		//从tree Fill TH1F或者从root文件Get TH1F形成新的TH1F用于寻峰拟合
		//for(ii=0;ii<16;ii++)
		//{
		//sprintf(h_name,"%s%d","h300ALch",ii);
		//h300AL[ii] = (TH1D*)fin->Get(h_name);
		//}//read histograms from a file, no need to name histograms, also cannot rebin
		// 	h300AL[0]=(TH1D*)fin->Get("D300ALch0");

		ofstream outfile("C:/Si24/Si23peakcali/Sipeakcali.dat",ios::out);
		for(i=0;i<ID;i++)
		{
			sprintf(hfit_name,"peak_%s%d",h_name,i);
			histo[i]->GetXaxis()->SetRangeUser(binmin,binmax);//zoom the axis
			canvaspeak[i]=new TCanvas(hfit_name,hfit_name,600,480);//建立画布
			canvaspeak[i]->cd();//进入画布
			histo[i]->SetTitle(hfit_name);//图名
			histo[i]->GetXaxis()->SetTitle("Channel");//轴名
			histo[i]->GetYaxis()->SetTitle("Counts");//轴名
			histo[i]->GetXaxis()->CenterTitle();//居中
			histo[i]->GetYaxis()->CenterTitle();//居中
			histo[i]->GetXaxis()->SetLabelFont(132);//坐标字体
			histo[i]->GetYaxis()->SetLabelFont(132);//坐标字体
			histo[i]->GetXaxis()->SetTitleFont(132);//轴名字体
			histo[i]->GetYaxis()->SetTitleFont(132);//轴名字体
			histo[i]->GetXaxis()->SetTitleOffset(1.2);//轴名偏移
			histo[i]->GetYaxis()->SetTitleOffset(1.3);//轴名偏移
			histo[i]->Draw();
			TSpectrum *s = new TSpectrum();
			//TSpectrum contains advanced spectra processing functions for 1- and 2-dimensional background estimation, smoothing, deconvolution, peak search and fitting, and orthogonal transformations.
			peaknum = s->Search(histo[i],sigma,"",thresh);//△寻峰个数
			// Int_t Search(const TH1* hist, Double_t sigma = 2, Option_t* option = "", Double_t threshold = 0.05)
			//hin: pointer to the histogram of source spectrum
			//sigma: sigma of searched peaks, for details we refer to manual
			//threshold: (default=0.05) peaks with amplitude less than threshold*highest_peak are discarded. 0 By default, the background is removed before deconvolution. Specify the option "nobackground" to not remove the background.
			//printf("Found %d candidate peaks to fit in DSSD%d.\n",peaknum,i);
			peakch = s->GetPositionX();//△寻峰道数
			for(ii=0;ii<peaknum;ii++)
			{
				k=ii;
				peakmin=peakch[k];
				for(jj=ii+1;jj<peaknum;jj++)
				{
					if(peakch[jj]<peakmin)
					{
						k=jj;
						peakmin=peakch[k];
						chtemp=peakch[jj];
						peakch[jj]=peakch[ii];
						peakch[ii]=chtemp;
					}
				}	
			}
			for(jj=0;jj<peaknum;jj++)
			{
				outfile<<hfit_name<<" canvaspeak"<<jj<<" "<<peakch[jj]<<endl;//输出文本查看
			}
			if(peaknum>nummax) peaknum=nummax;//not more than 10 peaks
			for(ii=0;ii<peaknum;ii++)
			{
				g[ii]=new TF1("g","gaus",peakch[ii]-gaplow,peakch[ii]+gaphigh);//高斯拟合，调节合适的拟合range
				histo[i]->Fit("g","R+");//restrict the fit to the range specified in the TF1 constructor
				g[ii]->GetParameters(par[ii]);//g[ii],pointer to the TF1, GetParameters的数组得是double类型
			}
			//TF1 *g2=new TF1("g2","gaus",peakch[1]-gaplow/2,peakch[1]+gaphigh/2);//不同大小的峰有不同的拟合range，不适合写自动拟合循环
			//gStyle->SetFitColor(3);
			//h300AL[i]->Fit("g2","R+");//
			//TF1 *gt=new TF1("gt","gaus(0)+gaus(3)",peakch[0]-gaplow,peakch[1]+gaphigh);
			//g2->GetParameters(&par[3]);//Get到par数组，同时也做gt的Set数组，更自动
			//gt->SetParameters(par);//自动，两个高斯的参数做gt的初始化，In the more complicated case of the sum of 3 Gaussian functions, the initial values of parameters must be set. In this particular case, the initial values are taken from the result of the individual fits.
			//h300AL[i]->SetLineColor(1);
			//h300AL[i]->Fit("gt","R+");//+ means adding this new fitted function to the list of fitted functions (by default, the previous function is deleted and only the last one is kept)
			//double Sigma1=gt->GetParameter(2);//Obtaining the value of the 3rd parameter (Sigma)
			//double Mean1=gt->GetParameter(1);//Obtaining the value of the 2nd parameter (Mean)
			//double Sigma2=gt->GetParameter(5);//Obtaining the value of the 3rd parameter (Sigma)
			//double Mean2=gt->GetParameter(4);//Obtaining the value of the 2nd parameter (Mean)
			TPaveText *textgaus = new TPaveText(0.13,0.50,0.31,0.85,"brNDC");//加标注
			textgaus->SetBorderSize(1);//边框宽度
			textgaus->SetFillColor(0);//填充颜色
			textgaus->SetTextAlign(12);//align = 10*HorizontalAlign + VerticalAlign, 12 means水平左对齐、垂直居中对齐
			textgaus->SetTextFont(132);//font = 10 * fontID + precision, 12+2,12 means Symbol; 13+2, 13 means News Times Roman
			//text->SetTextColor(2);//文本颜色
			for(ii=0;ii<peaknum;ii++)
			{
				sprintf(paraprint,"Mean%d=%.1f",ii,par[ii][1]);//par数组还保持着g1、g2的参数
				textgaus->AddText(paraprint);
				sprintf(paraprint,"Sigma%d=%.1f",ii,par[ii][2]);
				textgaus->AddText(paraprint);
				sprintf(paraprint,"Res%d=%.2f%%",ii,par[ii][2]/par[ii][1]*2.355*100);
				textgaus->AddText(paraprint);
			}
			//sprintf(paraprint,"Mean1=%4.2f",Mean1);//gt拟合参数，似乎寻峰不太准，拟合的Mean作峰位更准
			textgaus->Draw();
			sprintf(hfit_name,"%s.png",hfit_name);
			canvaspeak[i]->SaveAs(hfit_name);//存图
			for(ii=0;ii<peaknum;ii++)
			{
				outfile<<"h300AL"<<i<<" Mean"<<ii<<"="<<par[ii][1]<<endl;//输出文本查看
			}
			sprintf(hcali_name,"cali_%s%d",h_name,i);
			canvascali[i]=new TCanvas(hcali_name,hcali_name,600,480);//建立画布
			canvascali[i]->cd();//进入画布
			//canvascali[i]->SetGrid();//显示网格
			graph[i]=new TGraph(peaknum,peakch,energy);//TGraph *gr1=new TGraph(n,x,y);
			//float eenergy[2]={200,280};
			//float epeakch[2]={0,0};
			//graph[i]= new TGraphErrors(peaknum,peakch,energy,epeakch,eenergy);//画error bars
			graph[i]->SetTitle(hcali_name);
			graph[i]->GetXaxis()->SetTitle("Channel");//轴名
			graph[i]->GetYaxis()->SetTitle("Energy");//轴名
			graph[i]->GetXaxis()->CenterTitle();//居中
			graph[i]->GetYaxis()->CenterTitle();//居中
			graph[i]->GetXaxis()->SetLabelFont(132);//坐标字体
			graph[i]->GetYaxis()->SetLabelFont(132);//坐标字体
			graph[i]->GetXaxis()->SetTitleFont(132);//轴名字体
			graph[i]->GetYaxis()->SetTitleFont(132);//轴名字体
			//graph[i]->GetYaxis()->SetLabelSize(0.05);//坐标字号
			//graph[i]->GetYaxis()->SetTitleSize(0.05);//轴名字号
			graph[i]->GetXaxis()->SetTitleOffset(1.2);//轴名偏移
			graph[i]->GetYaxis()->SetTitleOffset(1.3);//轴名偏移
			graph[i]->SetMarkerStyle(21);
			graph[i]->SetMarkerColor(1);
			TF1 *pol1=new TF1("pol1","pol1",200,4000);//多项式拟合，调节合适的拟合range
			graph[i]->Fit("pol1");//pol1 can be used directly without TF1 constructor in CINT
			intercept[i]=pol1->GetParameter(0);
			slope[i]=pol1->GetParameter(1);
			graph[i]->Draw("AP"); //A-Axis around the graph,AP is suitable
			TPaveText *textpol1 = new TPaveText(0.13,0.70,0.33,0.85,"brNDC");
			textpol1->SetBorderSize(1);
			textpol1->SetFillColor(0);
			textpol1->SetTextAlign(12);//align = 10*HorizontalAlign + VerticalAlign, 12 means水平左对齐、垂直居中对齐
			textpol1->SetTextFont(132);//font = 10 * fontID + precision, 12+2,12 means Symbol; 13+2, 13 means News Times Roman
			sprintf(paraprint,"slope=%.2f",slope[i]);
			textpol1->AddText(paraprint);
			sprintf(paraprint,"intercept=%.2f",intercept[i]);
			textpol1->AddText(paraprint);
			textpol1->Draw();
			sprintf(hcali_name,"%s.png",hcali_name);
			canvascali[i]->SaveAs(hcali_name);//存图
			outfile<<"h300AL"<<i<<" slope="<<slope[i]<<endl;//输出文本查看
			outfile<<"h300AL"<<i<<" intercept="<<intercept[i]<<endl;//输出文本查看
		}//for (i=0;i<ID;i++)
	//}//for(irawroot=runstart; irawroot<=runstop; irawroot++)
}//peakcali main