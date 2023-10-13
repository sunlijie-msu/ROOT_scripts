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
#include "TChain.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TEllipse.h"
#include "TString.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TSpectrum.h"
#include "TVirtualFitter.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "TRandom.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH1D.h"
#include "TH2.h"
#include "TF1.h"
#include "TMath.h"
#include "TPaveText.h"
#include "TPad.h"
#include "TAxis.h"
#include "TClass.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "string.h"
#include "stdlib.h"
#include "stdio.h"
#include "TLegend.h"
using namespace std;
void gausnerfcpol1_peakfit()// old way, peak search then fit, not good. peak draw-fit (gausn+pol1) used for 25Si SeGA 16 channels calibration (after gain matching)
{
	time_t tim;
	struct tm *at;
	char now[80];
	int irawroot;
	char rawrootname[300];
	double slope[300];
	double intercept[300];
	const int ID1=0;//adjust i=ID1
	const int ID2=0;//adjust i<=ID2
	const int binmin=480,binmax=520;//adjust//search for peaks in this range
	const int binwidth=1;
	int binnum=(binmax-binmin)/binwidth;
	double sigma=1,thresh=0.7,sigmab=20;//adjust
	float gaplow=7.,gaphigh=7;//fitting range随分辨不同调整//adjust
	unsigned long i;
	int jj,ii;
	char paraprint[30],b_name[50],histo_name[50],h_name[50],hfit_name[50],hcali_name[50];
	char command1[100],command2[100];
	//string paraprint[30],h_name[80],histo_name[80],b_name[80],hfit_name[80],hcali_name[80];
	TCanvas *canvaspeak[300];
	TCanvas *canvascali[300];
	TH1F *histo[ID2+1];//TH1F peak search+gauss fit,creat histograms
	//TH1D *h300AL[16] cannot be TH1D *h300AL[ii]= new TH1D(h_name,h_name,4000,0,4000);in loop
	TF1 *total[12];//creat function
	TF1 *p[12], *g[12], *b[12];
	TGraph *graph[300];//TGraph刻度
	float *peakx;
	float *peakxn;
	float *peaky;
	float peakmin;
	int k;
	float chtemp;
	int peaknum;
	//float *peakmean;
	const int nummax=1;//search peaknumber can be many, fitnumber can be limited to 3//adjust
	double par[1][5],parerr[5];//随峰数不同改动par[peaknum][3]//adjust
	double parchi,parNDF;

	//float energy[6]={778.90,867.39,964.06,1085.84,1112.09,1408.02};//adjust
	//float energy[8]={121.78,244.70,344.28,778.90,964.06,1085.84,1112.09,1408.02};
	//float energy[11]={121.8, 244.7, 344.3, 411.1, 444.0, 778.9, 867.4, 964.1, 1085.8, 1112.0, 1408.0}//wangjianguo 152Eu
	//float intensity[11]={28.37,7.53,26.57,2.23,3.12,12.97,4.214,14.63,10.13,13.54,20.85};//wangjianguo 152Eu
	//float energy[5]={80.9, 276.4,302.8,356.0,383.8};//wangjianguo 133Ba
	//float intensity[5]={34.11,7.147,18.30,61.94,8.905};//wangjianguo 133Ba
	//float energy[2]={1173.2,1332.5};//wangjianguo 60Co
	//float intensity[2]={99.857,99.983};//wangjianguo 60Co
	//float energy[8]={401,943,1917,2162,2307,5805,3463,4252};
	float energy[7]={407,944,1922,2163,2306,3473,4257};
//	int runstart=123,runstop=140;
	//cout<<"input runstart:";
	//cin>>runstart;
	//cout<<"input runstop:";
	//cin>>runstop;
	ofstream outfile("D:/X/peakcali.dat",ios::out);//modify if change path
	//for(irawroot=runstart; irawroot<=runstop; irawroot++)
	//{
	//if()continue;
	sprintf(rawrootname,"%s","D:/X/out/Si25/100eV-bin-ungated_sun_25Si.root");//modify if change path
	TFile *fin = new TFile(rawrootname);//after this statement, you can use any ROOT command1 for this rootfile
	//TTree *tree = (TTree*)fin->Get("tree");
	//TTree *tree = (TTree*)fin->Get("tree");//此句可有可无
	//	unsigned long nentries=tree->GetEntries();//读事件数
	cout<<rawrootname<<endl;
	//	cout<<"  Entries="<<nentries<<endl;
	//Int_t adc[5][ID];//define the variables to hold the read values //△treeFill
	//tree->SetBranchAddress("adc", adc);//△treeFill
	//cout<<"input branch name for cali: ";
	//sprintf(b_name,"%s","Clover");
//	sprintf(b_name,"%s","SeGA_ungated");
	//cin>>b_name;//input the branch name you want to search its peaks
	// 	while(strcmp(b_name,"DSSD2e")!=0&&strcmp(b_name,"DSSD3e")!=0&&strcmp(b_name,"QSD1e")!=0)
	// 	{
	// 		cout<<"Wrong b_name!"<<endl;
	// 		cout<<"input b_name for cali again: ";
	// 		cin>>b_name;
	// 	}
	// tree draw *****************
// 	sprintf(b_name,"%s","SeGAenergy");
// 	for(ii=ID1;ii<=ID2;ii++)
// 	{
// 		sprintf(histo_name,"h%s%s%d%s",b_name,"[",ii,"]");
// 		histo[ii] = new TH1F(histo_name,histo_name,binnum,binmin,binmax);//name a histogram,hDSSD1AL[0]
// 	}//name sequence decide the TH1F sequence in rootfile, name histograms if use tree->draw
// 	for(ii=ID1;ii<=ID2;ii++)
// 	{
// 		sprintf(command1,"%s%s%d%s%s%s%d%s",b_name,"[",ii,"]>>h",b_name,"[",ii,"]");
// 		puts(command1);
// 		sprintf(command2,"%s%s%d%s",b_name,"[",ii,"]>0");
// 		puts(command2);
// 		outTree->Draw(command1,command2);//tree-Draw("DSSD1AL[0]>>hDSSD1AL[0]")
// 	}
	// tree draw *******************
	// get histogram *******************
	for(ii=ID1;ii<=ID2;ii++)
	{
		sprintf(b_name,"%s%02d","SeGA_",ii);
		histo[ii] = (TH1F*)fin->Get(b_name);
		histo[ii]->Rebin(1);
		histo[ii]->GetXaxis()->SetRangeUser(binmin,binmax);
	}
	// get histogram *******************
	// 	{
	// 		sprintf(command1,"%s%s%s%s%d%s",b_name,">>h",b_name,"[",ii,"]");
	// 		puts(command1);
	// 		sprintf(command2,"%s","DSSD2e<100&&QSD1e<100");
	// 		puts(command2);
	// 		tree->Draw(command1,command2);//tree-Draw("DSSD1AL[0]>>hDSSD1AL[0]")
	// 	}
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

	for(i=ID1;i<=ID2;i++)
	{
		if(i==8) continue;
		sprintf(hfit_name,"%s%d","Rn_all_SeGA_",i);
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
		//TH1F *histob = s->Background(histo[i], sigmab, "same");
		//histo[i]->Add(histob, -1);
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
// 		for(jj=0;jj<peaknum;jj++)
// 		{
// 			outfile<<hfit_name<<" canvaspeak"<<jj<<"	xpos="<<peakx[jj]<<"	ypos="<<peaky[jj]<<endl;//输出文本查看
// 		}
		if(peaknum>nummax) peaknum=nummax;//not more than 10 peaks
		for(ii=0;ii<peaknum;ii++)
		{
			// 			total[ii]=new TF1("total","[0]*x+[1]+[2]/2/[3]*exp(0.5*([4]*[4]/([3]*[3]))+(x-[5])/[3])*ROOT::Math::erfc(1/sqrt(2)*([4]/[3]+(x-[5])/[4]))",peakx[ii]-gaplow,peakx[ii]+gaphigh);//高斯拟合，调节合适的拟合range / Glassman PLB2018 low-energy tail
			// 			g[ii]=new TF1("g","gausn",peakx[ii]-gaplow,peakx[ii]+gaphigh);// The [2]-N parameter in total is equivalent to the Constant in gausn
			// 			p[ii]=new TF1("p","[0]*x+[1]-[0]*x-[1]+[2]/2/[3]*exp(0.5*([4]*[4]/([3]*[3]))+(x-[5])/[3])*ROOT::Math::erfc(1/sqrt(2)*([4]/[3]+(x-[5])/[4]))",peakx[ii]-gaplow,peakx[ii]+gaphigh);//pure peak
			// 			b[ii]=new TF1("b","[0]*x+[1]",peakx[ii]-gaplow,peakx[ii]+gaphigh);//pure bkg

			total[ii]=new TF1("total","[0]*x+[1]+sqrt(3.141592654/2)*[2]/[3]*[4]*exp(0.5*([4]*[4]/([3]*[3]))+(x-[5])/[3])*ROOT::Math::erfc(1/sqrt(2)*([4]/[3]+(x-[5])/[4]))",peakx[ii]-gaplow,peakx[ii]+gaphigh);//高斯拟合，调节合适的拟合range / Glassman PRC2019 low-energy tail
			g[ii]=new TF1("g","gausn",peakx[ii]-gaplow,peakx[ii]+gaphigh);// The [2]-N parameter in total is equivalent to the Constant in gausn
			p[ii]=new TF1("p","[0]*x+[1]-[0]*x-[1]+sqrt(3.141592654/2)*[2]/[3]*[4]*exp(0.5*([4]*[4]/([3]*[3]))+(x-[5])/[3])*ROOT::Math::erfc(1/sqrt(2)*([4]/[3]+(x-[5])/[4]))",peakx[ii]-gaplow,peakx[ii]+gaphigh);//pure peak
			b[ii]=new TF1("b","[0]*x+[1]",peakx[ii]-gaplow,peakx[ii]+gaphigh);//pure bkg

			// 			total[ii]=new TF1("total","[0]*x+[1]+[2]/2/[3]*exp([4]*[4]/2/[3]/[3]+(x-[5])/[3])*(1-ROOT::Math::erf([4]*[4]+[3]*(x-[5])/sqrt(2)/[4]/[3]))",peakx[ii]-gaplow,peakx[ii]+gaphigh);//高斯拟合，调节合适的拟合range / Bennett
			// 			g[ii]=new TF1("g","gausn",peakx[ii]-gaplow,peakx[ii]+gaphigh);
			// 			p[ii]=new TF1("p","[0]*x+[1]-[0]*x-[1]+[2]/2/[3]*exp([4]*[4]/2/[3]/[3]+(x-[5])/[3])*(1-ROOT::Math::erf([4]*[4]+[3]*(x-[5])/sqrt(2)/[4]/[3]))",peakx[ii]-gaplow,peakx[ii]+gaphigh);
			//  			b[ii]=new TF1("b","[0]*x+[1]",peakx[ii]-gaplow,peakx[ii]+gaphigh);//pure bkg

			//total[ii]->SetParameters(0.1,15,peaky[ii],10,10,peakx[ii]);//initial value [0]-A, [1]-B, [2]-N, [3]-τ, [4]-σ, [5]-μ
			// 			total[ii]->SetParameters(-0.6,200,40000,1000,3,390);//initial value [0]-A, [1]-B, [2]-N, [3]-τ, [4]-σ, [5]-μ  //Pérez-Loureiro
			// 			total[ii]->SetParLimits(0,-50,50);//Bkg A
			// 			total[ii]->SetParLimits(1,-5000,230);//Bkg B
			// 			total[ii]->SetParLimits(2,1000,1000000);//Constant,min,max
			// 			total[ii]->SetParLimits(3,0,5000);//Lambda
			// 			total[ii]->SetParLimits(4,0,20);//Sigma
			// 			total[ii]->SetParLimits(5,350,420);//Mean

			//total[ii]->SetParameters(0.1,15,peaky[ii],10,10,peakx[ii]);//initial value [0]-A, [1]-B, [2]-N, [3]-λ, [4]-σ, [5]-μ
			total[ii]->SetParameters(0.1,15,peaky[ii],10,10,peakx[ii]);//initial value [0]-A, [1]-B, [2]-N, [3]-τ, [4]-σ, [5]-μ
			total[ii]->SetParameters(0,500,1500000,100,1.5,511);//initial value [0]-A, [1]-B, [2]-N, [3]-λ, [4]-σ, [5]-μ  //Glassman
			total[ii]->SetParLimits(0,-500,500);//Bkg A  //Glassman
			total[ii]->SetParLimits(1,-50000,300000);//Bkg B  //Glassman
			total[ii]->SetParLimits(2,5,1000000);//Constant,min,max  //Glassman
			total[ii]->SetParLimits(3,0.01,1000);//Lambda  //Glassman
			total[ii]->SetParLimits(4,1,10);//Sigma  //Glassman
			total[ii]->SetParLimits(5,400,1730);//Mean  //Glassman

			total[ii]->SetParNames("BkgA","BkgB","Const*bin","Lamda","Sigma","Mean");
			histo[i]->Fit("total","MLE","",peakx[ii]-gaplow,peakx[ii]+gaphigh);///specify a range in the Fit
			//histo[i]->Fit("total","R+");//or restrict the fit to the range specified in the TF1 constructor: peakx[ii]-gaplow,peakx[ii]+gaphigh
			total[ii]->GetParameters(par[ii]);//二维数组的par[ii]是地址,pointer to the TF1, GetParameters的数组得是double类型Obtaining the value of parameters and saving them to par[]; 
			parerr[0]=total[ii]->GetParError(0);//Obtaining the error of the 1st parameter
			parerr[1]=total[ii]->GetParError(1);//Obtaining the error of the 2nd parameter
			parerr[2]=total[ii]->GetParError(2);//Obtaining the error of the 3rd parameter
			parerr[3]=total[ii]->GetParError(3);//Obtaining the error of the 4th parameter
			parerr[4]=total[ii]->GetParError(4);//Obtaining the error of the 5th parameter
			parerr[5]=total[ii]->GetParError(5);//Obtaining the error of the 6th parameter
			parchi=total[ii]->GetChisquare();
			parNDF=total[ii]->GetNDF();
			g[ii]->SetParameters(par[ii][2],par[ii][5],par[ii][4]);//set parameters for drawing gausn
			g[ii]->SetLineColor(4);
			p[ii]->SetParameters(par[ii][0],par[ii][1],par[ii][2],par[ii][3],par[ii][4],par[ii][5]);//set parameters for drawing peak
			p[ii]->SetLineColor(6);
			b[ii]->SetParameters(par[ii][0],par[ii][1]);//set parameters for drawing bkg
			b[ii]->SetLineColor(8);
			g[ii]->Draw("same");
			p[ii]->Draw("same");
			b[ii]->Draw("same");
			histo[i]->GetXaxis()->SetRangeUser(peakx[ii]-1*gaplow,peakx[ii]+1*gaphigh);//zoom the axis
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
		TPaveText *textgaus = new TPaveText(0.7,0.25,0.99,0.87,"brNDC");//加标注left, down, right, up
		textgaus->SetBorderSize(1);//边框宽度
		textgaus->SetFillColor(0);//填充颜色
		textgaus->SetTextAlign(12);//align = 10*HorizontalAlign + VerticalAlign, 12 means水平左对齐、垂直居中对齐
		textgaus->SetTextFont(132);//font = 10 * fontID + precision, 12+2,12 means Symbol; 13+2, 13 means News Times Roman
		//text->SetTextColor(2);//文本颜色

		for(ii=0;ii<peaknum;ii++)
		{
			sprintf(paraprint,"Constant%d=%.1f%s%.3f",ii,par[ii][2]/binwidth,"+/-",parerr[2]/binwidth);//par数组还保持着刚才的参数
			textgaus->AddText(paraprint);
			sprintf(paraprint,"Mean%d=%.1f%s%.3f",ii,par[ii][5],"+/-",parerr[5]);
			textgaus->AddText(paraprint);
			sprintf(paraprint,"Sigma%d=%.3f%s%.3f",ii,par[ii][4],"+/-",parerr[4]);
			textgaus->AddText(paraprint);
			sprintf(paraprint,"Res%d=%.2f%%",ii,par[ii][4]/par[ii][5]*2.355*100);
			textgaus->AddText(paraprint);
			sprintf(paraprint,"Lambda%d=%.2f%s%.3f",ii,par[ii][3],"+/-",parerr[3]);
			textgaus->AddText(paraprint);
			sprintf(paraprint,"A%d(%.2f%s%.3f)*x+B%d(%.1f%s%.3f)",ii,par[ii][0],"+/-",parerr[0],ii,par[ii][1],"+/-",parerr[1]);
			textgaus->AddText(paraprint);
			sprintf(paraprint,"Chisquare%d=%.2f",ii,parchi);
			textgaus->AddText(paraprint);
			sprintf(paraprint,"NDF%d=%.2f",ii,parNDF);
			textgaus->AddText(paraprint);
			sprintf(paraprint,"Centroid%d=%.1f%s",ii,peakx[ii]);//par数组还保持着刚才的参数
			textgaus->AddText(paraprint);
		}
		//sprintf(paraprint,"Mean1=%4.2f",Mean1);//gt拟合参数，似乎寻峰不太准，拟合的Mean作峰位更准
		textgaus->Draw();
		sprintf(h_name,"D:/X/out/cali/%s_%.1f.png",hfit_name,peakx[0]);//modify if change path
		canvaspeak[i]->SaveAs(h_name);//存图
		// 		for(ii=0;ii<peaknum;ii++)
		// 		{
		// 			outfile<<b_name<<i<<"	Constant"<<ii<<"="<<par[ii][0]<<"	Mean"<<ii<<"="<<par[ii][1]<<"	Sigma"<<ii<<"="<<par[ii][2]<<"	Area"<<ii<<"="<<sqrt(3.14159*2)*par[ii][2]*par[ii][0]/binwidth<<endl;//输出文本查看
		// 			peakx[ii]=par[ii][1];// For graph pol fit, mean is more accurate than peakx
		// 		}
		// 		sprintf(hcali_name,"cali_%s%s%d%s",b_name,"[",i,"]");
		// 		canvascali[i]=new TCanvas(hcali_name,hcali_name,900,600);//建立画布
		// 		canvascali[i]->cd();//进入画布
		// 		//canvascali[i]->SetGrid();//显示网格
		// 		graph[i]=new TGraph(peaknum,peakx,energy);//TGraph *gr1=new TGraph(n,x,y);
		// 		//float eenergy[2]={200,280};
		// 		//float epeakch[2]={0,0};
		// 		//graph[i]= new TGraphErrors(peaknum,peakx,energy,epeakch,eenergy);//画error bars
		// 		graph[i]->SetTitle(hcali_name);
		// 		graph[i]->GetXaxis()->SetTitle("Channel");//轴名
		// 		graph[i]->GetYaxis()->SetTitle("Energy");//轴名
		// 		graph[i]->GetXaxis()->CenterTitle();//居中
		// 		graph[i]->GetYaxis()->CenterTitle();//居中
		// 		graph[i]->GetXaxis()->SetLabelFont(132);//坐标字体
		// 		graph[i]->GetYaxis()->SetLabelFont(132);//坐标字体
		// 		graph[i]->GetXaxis()->SetTitleFont(132);//轴名字体
		// 		graph[i]->GetYaxis()->SetTitleFont(132);//轴名字体
		// 		//graph[i]->GetYaxis()->SetLabelSize(0.05);//坐标字号
		// 		//graph[i]->GetYaxis()->SetTitleSize(0.05);//轴名字号
		// 		graph[i]->GetXaxis()->SetTitleOffset(1.2);//轴名偏移
		// 		graph[i]->GetYaxis()->SetTitleOffset(1.3);//轴名偏移
		// 		graph[i]->SetMarkerStyle(21);
		// 		graph[i]->SetMarkerColor(1);
		// 		TF1 *pol1=new TF1("pol1","pol1",20,5000);//多项式拟合，调节合适的拟合range
		// 		graph[i]->Fit("pol1");//pol1 can be used directly without TF1 constructor in CINT
		// 		intercept[i]=pol1->GetParameter(0);
		// 		slope[i]=pol1->GetParameter(1);
		// 		graph[i]->Draw("AP"); //A-Axis around the graph,AP is suitable
		// 		TPaveText *textpol1 = new TPaveText(0.13,0.70,0.33,0.85,"brNDC");
		// 		textpol1->SetBorderSize(1);
		// 		textpol1->SetFillColor(0);
		// 		textpol1->SetTextAlign(12);//align = 10*HorizontalAlign + VerticalAlign, 12 means水平左对齐、垂直居中对齐
		// 		textpol1->SetTextFont(132);//font = 10 * fontID + precision, 12+2,12 means Symbol; 13+2, 13 means News Times Roman
		// 		sprintf(paraprint,"slope=%.2f",slope[i]);
		// 		textpol1->AddText(paraprint);
		// 		sprintf(paraprint,"intercept=%.2f",intercept[i]);
		// 		textpol1->AddText(paraprint);
		// 		textpol1->Draw();
		// 		//sprintf(hcali_name,"%s.png",hcali_name);
		// 		sprintf(hcali_name,"%s%s%s","D:/X/",hcali_name,".png");//modify if change path
		// 		canvascali[i]->SaveAs(hcali_name);//存图
		// 		outfile<<b_name<<i<<" slope="<<slope[i]<<endl;//输出文本查看
		// 		outfile<<b_name<<i<<" intercept="<<intercept[i]<<endl;//输出文本查看
	}//for (i=0;i<ID;i++)
	//}//for(irawroot=runstart; irawroot<=runstop; irawroot++)
}//peakcali main