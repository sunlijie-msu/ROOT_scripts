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
void Peakcali_treedraw()
{
	time_t tim;
	struct tm *at;
	char now[80];
	int irawroot;
	char rawrootname[80];
	double slope[300];
	double intercept[300];
	const int ID1=7;//adjust
	const int ID2=8;//adjust
	const int binmin=300,binmax=6000;//adjust//search for peaks in this range
	const int binwidth=10;
	int binnum=(binmax-binmin)/binwidth;
	double sigma=7,thresh=0.3,sigmab=20;//adjust
	float gaplow=100.,gaphigh=100.;//fitting range随分辨不同调整//adjust
	unsigned long i;
	int jj,ii;
	char paraprint[30],b_name[50],histo_name[50],h_name[50],hfit_name[50],hcali_name[50];
	char command[60];
	//string paraprint[30],h_name[80],histo_name[80],b_name[80],hfit_name[80],hcali_name[80];
	TCanvas *canvaspeak[300];
	TCanvas *canvascali[300];
	TH1F *histo[ID2];//TH1F peak search+gauss fit,creat histograms
	//TH1D *h300AL[16] cannot be TH1D *h300AL[ii]= new TH1D(h_name,h_name,4000,0,4000);in loop
	TF1 *g[12];//creat function
	TF1 *p[12];
	TGraph *graph[300];//TGraph刻度
	float *peakx;
	float *peakxn;
	float *peaky;
	float peakmin;
	int k;
	float chtemp;
	int peaknum;
	//float *peakmean;
	const int nummax=12;//search peaknumber can be many, fitnumber can be limited to 3//adjust
	double par[12][5];//随峰数不同改动par[peaknum][3]//adjust
	//float energy[6]={778.90,867.39,964.06,1085.84,1112.09,1408.02};//adjust
	//float energy[8]={121.78,244.70,344.28,778.90,964.06,1085.84,1112.09,1408.02};
	float energy[8]={401,943,1917,2162,2307,5805,3463,4252};
	//float energy[3]={797,1670,2700};
	//int runstart,runstop;
	//cout<<"input runstart:";
	//cin>>runstart;
	//cout<<"input runstop:";
	//cin>>runstop
	//for(irawroot=runstart; irawroot<=runstop; irawroot++)
	//{
		//if()continue;
		sprintf(rawrootname,"X:/T999/Si25_0154_0345_p230_c250_t5T_Erep.root");//modify if change path
		//sprintf(calirootname,"%s%04d%s","X:/RIBLL2015/cali",irawroot,".root");
		TFile *fin = new TFile(rawrootname);//after this statement, you can use any ROOT command for this rootfile
		//TTree *tree = (TTree*)fin->Get("tree");
		//TTree *T999 = (TTree*)fin->Get("T999");//此句可有可无
		unsigned long nentries=T999->GetEntries();//读事件数
		cout<<rawrootname;
		cout<<"  Entries="<<nentries<<endl;
		//Int_t adc[5][ID];//define the variables to hold the read values //△treeFill
		//tree->SetBranchAddress("adc", adc);//△treeFill
		cout<<"input branch name for cali: ";
		sprintf(b_name,"%s","D40BLch");
		//cin>>b_name;//input the branch name you want to search its peaks
		while(strcmp(b_name,"D142ALch")!=0&&strcmp(b_name,"D142BLch")!=0&&strcmp(b_name,"D40ALch")!=0&&strcmp(b_name,"D40BLch")!=0&&strcmp(b_name,"D304ALch")!=0&&strcmp(b_name,"D304BLch")!=0)
		{
			cout<<"Wrong b_name!"<<endl;
			cout<<"input b_name for cali again: ";
			cin>>b_name;
		}
		for(ii=ID1;ii<ID2;ii++)
		{
			sprintf(histo_name,"h%s%s%d%s",b_name,"[",ii,"]");
			histo[ii] = new TH1F(histo_name,histo_name,binnum,binmin,binmax);//name a histogram,hDSSD1AL[0]
		}//name sequence decide the TH1F sequence in rootfile, name histograms if use tree->draw
		for(ii=ID1;ii<ID2;ii++)
		{
				sprintf(command,"%s%s%d%s%s%s%d%s",b_name,"[",ii,"]>>h",b_name,"[",ii,"]");
				puts(command);
				T999->Draw(command);//T999-Draw("DSSD1AL[0]>>hDSSD1AL[0]")
		}
// 		for(i=0;i<nentries;i++)
// 		{
// 			memset(adc,0,sizeof(adc));
// 			for(ii=0;ii<32;ii++)
// 			{
// 				tree->GetEntry(i);//Read the entry i in the tree (the i th event in the tree)
// 				//sprintf(b_name,"%s%s%d%s",h_name,"ch[",ii,"]");
// 				if(strcmp(h_name,"adc[2]")==0)
// 					histo[ii]->Fill(adc[2][ii]);//Fill Branch for cali
// 				else if(strcmp(h_name,"adc[3]")==0)
// 					histo[ii]->Fill(adc[3][ii]);//Fill Branch for cali
// 			}//Fill(Float_t), not Fill(char *)
// 			if(i%iprint==0)
// 			{
// 				time(&tim);
// 				at=localtime(&tim);
// 				strftime(now,79,"%Y-%m-%d %H:%M:%S",at);
// 				cout<<now;
// 				printf(" complete %.1f%s\n",100*(float)i/(float)nentries,"%");
// 			}
// 			if(i+1>=imax)break;
// 		}//for(i=0;i<nentries;i++) △treedraw
		//从tree Draw TH1F或者从root文件Get TH1F形成新的TH1F用于寻峰拟合
		//for(ii=0;ii<16;ii++)
		//{
		//sprintf(h_name,"%s%d","h300ALch",ii);
		//h300AL[ii] = (TH1D*)fin->Get(h_name);
		//}//read histograms from a file, no need to name histograms, also cannot rebin
		// 	h300AL[0]=(TH1D*)fin->Get("D300ALch0");

		ofstream outfile("X:/T999/cali/peakcali.dat",ios::out);//modify if change path
		for(i=ID1;i<ID2;i++)
		{
			sprintf(hfit_name,"peak_%s%s%d%s",b_name,"[",i,"]");
			histo[i]->GetXaxis()->SetRangeUser(binmin,binmax);//zoom the axis
			canvaspeak[i]=new TCanvas(hfit_name,hfit_name,900,600);//建立画布
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
			TH1F *histob = s->Background(histo[i], sigmab, "same");
			histo[i]->Add(histob, -1);
			//TSpectrum contains advanced spectra processing functions for 1- and 2-dimensional background estimation, smoothing, deconvolution, peak search and fitting, and orthogonal transformations.
			peaknum = s->Search(histo[i],sigma,"",thresh);//△寻峰个数
			// Int_t Search(const TH1* hist, Double_t sigma = 2, Option_t* option = "", Double_t threshold = 0.05)
			//hin: pointer to the histogram of source spectrum
			//sigma: sigma of searched peaks, for details we refer to manual
			//threshold: (default=0.05) peaks with amplitude less than threshold*highest_peak are discarded. 0 By default, the background is removed before deconvolution. Specify the option "nobackground" to not remove the background.
			//printf("Found %d candidate peaks to fit in DSSD%d.\n",peaknum,i);
			peakx = s->GetPositionX();//△寻峰道数
			peaky = s->GetPositionY();//△寻峰道数
			for(ii=0;ii<peaknum;ii++)
			{
				k=ii;
				peakmin=peakx[k];
				for(jj=ii+1;jj<peaknum;jj++)
				{
					if(peakx[jj]<peakmin)
					{
						//cout<<peakx[jj]<<" "<<peakmin<<endl;
						k=jj;
						peakmin=peakx[k];
						chtemp=peakx[jj];
						peakx[jj]=peakx[ii];
						peakx[ii]=chtemp;

						chtemp=peaky[jj];
						peaky[jj]=peaky[ii];
						peaky[ii]=chtemp;
					}
				}	
			}
			//TMath::Sort(peaknum, peakxn, peakx);
			for(jj=0;jj<peaknum;jj++)
			{
				outfile<<hfit_name<<" canvaspeak"<<jj<<"	xpos="<<peakx[jj]<<"	ypos="<<peaky[jj]<<endl;//输出文本查看
			}
			if(peaknum>nummax) peaknum=nummax;//not more than 10 peaks
			for(ii=0;ii<peaknum;ii++)
			{
				g[ii]=new TF1("g","gaus(0)+pol1(3)",peakx[ii]-gaplow,peakx[ii]+gaphigh);//高斯拟合，调节合适的拟合range
				p[ii]=new TF1("p","pol1",peakx[ii]-gaplow,peakx[ii]+gaphigh);
				g[ii]->SetParameters(peaky[ii],peakx[ii],10,0,0);
				g[ii]->SetParLimits(0,5,300000);//Constant,min,max
				//g[ii]->SetParLimits(1,5,300000);//Mean
				//g[ii]->SetParLimits(2,5,300000);//Sigma
				//g[ii]->SetParLimits(3,0,0);//a
				//g[ii]->SetParLimits(4,0,0);//b
				histo[i]->Fit("g","R+");//restrict the fit to the range specified in the TF1 constructor: peakx[ii]-gaplow,peakx[ii]+gaphigh
				g[ii]->GetParameters(par[ii]);//g[ii],pointer to the TF1, GetParameters的数组得是double类型Obtaining the value of parameters and saving them to par[]; 
				p[ii]->SetParameters(par[ii][3],par[ii][4]);
				p[ii]->SetLineColor(8);
				p[ii]->Draw("same");
			}
			//TF1 *g2=new TF1("g2","gaus",peakx[1]-gaplow/2,peakx[1]+gaphigh/2);//不同大小的峰有不同的拟合range，不适合写自动拟合循环
			//gStyle->SetFitColor(3);
			//h300AL[i]->Fit("g2","R+");//
			//TF1 *gt=new TF1("gt","gaus(0)+gaus(3)",peakx[0]-gaplow,peakx[1]+gaphigh);
			//g2->GetParameters(&par[3]);//Get到par数组，同时也做gt的Set数组，更自动
			//gt->SetParameters(par);//自动，两个高斯的参数做gt的初始化，In the more complicated case of the sum of 3 Gaussian functions, the initial values of parameters must be set. In this particular case, the initial values are taken from the result of the individual fits.
			//h300AL[i]->SetLineColor(1);
			//h300AL[i]->Fit("gt","R+");//+ means adding this new fitted function to the list of fitted functions (by default, the previous function is deleted and only the last one is kept)
			//double Sigma1=gt->GetParameter(2);//Obtaining the value of the 3rd parameter (Sigma)
			//double Mean1=gt->GetParameter(1);//Obtaining the value of the 2nd parameter (Mean)
			//double Sigma2=gt->GetParameter(5);//Obtaining the value of the 3rd parameter (Sigma)
			//double Mean2=gt->GetParameter(4);//Obtaining the value of the 2nd parameter (Mean)
			TPaveText *textgaus = new TPaveText(0.83,0.15,0.99,0.77,"brNDC");//加标注left, down, right, up
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
			//sprintf(hfit_name,"%s.png",hfit_name);
			sprintf(hfit_name,"%s%s%s","X:/T999/cali/",hfit_name,".png");//modify if change path
			canvaspeak[i]->SaveAs(hfit_name);//存图
			for(ii=0;ii<peaknum;ii++)
			{
				outfile<<b_name<<i<<"	Constant"<<ii<<"="<<par[ii][0]<<"	Mean"<<ii<<"="<<par[ii][1]<<"	Sigma"<<ii<<"="<<par[ii][2]<<"	Area"<<ii<<"="<<sqrt(3.14159*2)*par[ii][2]*par[ii][0]/binwidth<<endl;//输出文本查看
				peakx[ii]=par[ii][1];// For graph pol fit, mean is more accurate than peakx
			}
			sprintf(hcali_name,"cali_%s%s%d%s",b_name,"[",i,"]");
			canvascali[i]=new TCanvas(hcali_name,hcali_name,900,600);//建立画布
			canvascali[i]->cd();//进入画布
			//canvascali[i]->SetGrid();//显示网格
			graph[i]=new TGraph(peaknum,peakx,energy);//TGraph *gr1=new TGraph(n,x,y);
			//float eenergy[2]={200,280};
			//float epeakch[2]={0,0};
			//graph[i]= new TGraphErrors(peaknum,peakx,energy,epeakch,eenergy);//画error bars
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
			TF1 *pol1=new TF1("pol1","pol1",20,5000);//多项式拟合，调节合适的拟合range
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
			//sprintf(hcali_name,"%s.png",hcali_name);
			sprintf(hcali_name,"%s%s%s","X:/T999/cali/",hcali_name,".png");//modify if change path
			canvascali[i]->SaveAs(hcali_name);//存图
			outfile<<b_name<<i<<" slope="<<slope[i]<<endl;//输出文本查看
			outfile<<b_name<<i<<" intercept="<<intercept[i]<<endl;//输出文本查看
		}//for (i=0;i<ID;i++)
	//}//for(irawroot=runstart; irawroot<=runstop; irawroot++)
}//peakcali main