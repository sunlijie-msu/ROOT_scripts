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
void respfunc_derivative()// peak draw-fit (gausn+pol1) used for 25Si SeGA 16 channels calibration (after gain matching+calibration with 6 peaks)
{
	time_t tim;
	struct tm *at;
	char now[80];
	int irawroot;
	char rawrootname[300];
	double slope[300];
	double intercept[300];
	const int ID1=0;// no need to modify i=ID1//which detector
	const int ID2=15;// no need to modify i<=ID2//which detector
	
	int binwidth=1;
//	int binnum=(upperlimit-lowerlimit)/binwidth;
	int minrange=0,maxrange=0,minbin,maxbin;
	double sigma=1,thresh=0.8,sigmab=20;//adjust
	float gaplow=70.,gaphigh=70.;//fitting range随分辨不同调整//no need to modify
	unsigned long i;
	int jj,ii,ibin,k,ipeak;
	char paraprint[30],b_name[50],histo_name[50],h_name[50],hfit_name[50],hcali_name[50];
	char command1[100],command2[100];
	//string paraprint[30],h_name[80],histo_name[80],b_name[80],hfit_name[80],hcali_name[80];
	TCanvas *canvaspeak[300];
	TCanvas *canvascali[300];
	TH1F *histo[ID2+1];//TH1F peak search+gauss fit,creat histograms
	const int nummax=6;//search peaknumber can be many, fitnumber can be limited to 3//no need to modify
	TF1 *total[nummax];//creat function
	TF1 *derivative[nummax];
	TF1 *p[nummax], *g[nummax], *b[nummax];
	TGraph *graph[300];//TGraph刻度
	int peaknum=nummax;
	Double_t peakx[nummax];
	Double_t peaky[nummax];
// 	float *peakx;
// 	float *peaky;
	float peakmin;
	float chtemp;
	double par[nummax][5],parerr[5];//随峰数不同改动par[peaknum][3]//adjust
	double parchi[nummax],parNDF[nummax];
	ifstream infilelowerlimit("X:/aaron/lowerlimitforresponsefunction.dat",ios::in);
	ifstream infileupperlimit("X:/aaron/upperlimitforresponsefunction.dat",ios::in);
	int lowerlimit[ID2+1][nummax],upperlimit[ID2+1][nummax];
	for(i=ID1;i<=ID2;i++)//which detector
	{
		if(i==8) continue;
		for(ii=0;ii<nummax;ii++)//which peak
		{
			infilelowerlimit>>lowerlimit[i][ii];
			infileupperlimit>>upperlimit[i][ii];
		}
	}
	for(i=ID1;i<=ID2;i++)//which detector
	{
		if(i==8) continue;
		for(ii=0;ii<nummax;ii++)//which peak
		{
			cout<<upperlimit[i][ii]<<"	";
		}
		cout<<endl;
	}
	// 	int lowerlimit[6]={3400,3700,7200,11100,12300,20000};//portal//search for peaks in this range
	// 	int upperlimit[6]={3550,3850,7330,11350,12500,20220};//search for peaks in this range
	//float energy[6]={778.90,867.39,964.06,1085.84,1112.09,1408.02};//adjust
	//float energy[8]={121.78,244.70,344.28,778.90,964.06,1085.84,1112.09,1408.02};
	//float energy[11]={121.8, 244.7, 344.3, 411.1, 444.0, 778.9, 867.4, 964.1, 1085.8, 1112.0, 1408.0}//wangjianguo 152Eu
	//float intensity[11]={28.37,7.53,26.57,2.23,3.12,12.97,4.214,14.63,10.13,13.54,20.85};//wangjianguo 152Eu
	//float energy[5]={80.9, 276.4,302.8,356.0,383.8};//wangjianguo 133Ba
	//float intensity[5]={34.11,7.147,18.30,61.94,8.905};//wangjianguo 133Ba
	//float energy[2]={1173.2,1332.5};//wangjianguo 60Co
	//float intensity[2]={99.857,99.983};//wangjianguo 60Co
	//float energy[8]={401,943,1917,2162,2307,5805,3463,4252};
	//float energy[7]={407,944,1922,2163,2306,3473,4257};
//	int runstart=123,runstop=140;
	//cout<<"input runstart:";
	//cin>>runstart;
	//cout<<"input runstop:";
	//cin>>runstop;
	ofstream outfile("X:/aaron/peakcali.dat",ios::out);
	//FILE *outfile=fopen ("D:/X/peakcali.dat","a");
	//for(irawroot=runstart; irawroot<=runstop; irawroot++)
	//{
	//if()continue;
	sprintf(rawrootname,"%s","X:/aaron/rootfiles/run-all_after_cali.root");
//	sprintf(rawrootname,"D:/X/out/Si25/good-dets-ungated.root");//
	//	sprintf(rawrootname,"X:/T777/copysmall.root");//
	//sprintf(calirootname,"%s%04d%s","X:/RIBLL2015/cali",irawroot,".root");
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
// 		histo[ii] = new TH1F(histo_name,histo_name,binnum,lowerlimit,upperlimit);//name a histogram,hDSSD1AL[0]
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
	for(i=ID1;i<=ID2;i++)//which detector
	{
		sprintf(b_name,"%s%02d","SeGA_",i);
		histo[i] = (TH1F*)fin->Get(b_name);
		//histo[i]->Rebin(10);
		//histo[i]->GetYaxis()->SetRangeUser(0,4000);//zoom the axis
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

	for(i=12;i<=13;i++)//which SeGA detector modify
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
//		TSpectrum *s = new TSpectrum();
		//TH1F *histob = s->Background(histo[i], sigmab, "same");
		//histo[i]->Add(histob, -1);
		//TSpectrum contains advanced spectra processing functions for 1- and 2-dimensional background estimation, smoothing, deconvolution, peak search and fitting, and orthogonal transformations.
//		peaknum = s->Search(histo[i],sigma,"",thresh);//△寻峰个数
		// Int_t Search(const TH1* hist, Double_t sigma = 2, Option_t* option = "", Double_t threshold = 0.05)
		//hin: pointer to the histogram of source spectrum
		//sigma: sigma of searched peaks, for details we refer to manual
		//threshold: (default=0.05) peaks with amplitude less than threshold*highest_peak are discarded. 0 By default, the background is removed before deconvolution. Specify the option "nobackground" to not remove the background.
		//printf("Found %d candidate peaks to fit in DSSD%d.\n",peaknum,i);
// 		peakx = s->GetPositionX();//△寻峰道数
// 		peaky = s->GetPositionY();//△寻峰道数
// 		for(ii=0;ii<peaknum;ii++)
// 		{
// 			k=ii;
// 			peakmin=peakx[k];
// 			for(jj=ii+1;jj<peaknum;jj++)
// 			{
// 				if(peakx[jj]<peakmin)
// 				{
// 					//cout<<peakx[jj]<<" "<<peakmin<<endl;
// 					k=jj;
// 					peakmin=peakx[k];
// 					chtemp=peakx[jj];
// 					peakx[jj]=peakx[ii];
// 					peakx[ii]=chtemp;
// 
// 					chtemp=peaky[jj];
// 					peaky[jj]=peaky[ii];
// 					peaky[ii]=chtemp;
// 				}
// 			}	
// 		}
		//TMath::Sort(peaknum, peakxn, peakx);
// 		for(jj=0;jj<peaknum;jj++)
// 		{
// 			outfile<<hfit_name<<" canvaspeak"<<jj<<"	xpos="<<peakx[jj]<<"	ypos="<<peaky[jj]<<endl;//输出文本查看
// 		}
//		if(peaknum>nummax) peaknum=nummax;//not more than 10 peaks
		for(ii=0;ii<6;ii++)//which peak in one SeGA detector modify
		{
			peaky[ii]=0;
			histo[i]->GetXaxis()->SetRangeUser(lowerlimit[i][ii],upperlimit[i][ii]);
			peaky[ii]=histo[i]->GetMaximum();
			peakx[ii]=histo[i]->GetBinCenter(histo[i]->GetMaximumBin());
			if(ii==0){	gaplow=7; gaphigh=7;}
			if(ii==1){	gaplow=7; gaphigh=7;}
			if(ii==2){	gaplow=8; gaphigh=8;}
			if(ii==3){	gaplow=10; gaphigh=10;}
			if(ii==4){	gaplow=10; gaphigh=10;}
			if(ii==5){	gaplow=14; gaphigh=14;}
			histo[i]->GetXaxis()->SetRangeUser(peakx[ii]-gaplow,peakx[ii]+gaphigh);//zoom the axis
			//cout<<"************"<<peakx[ii]<<"	"<<peaky[ii]<<endl;
			minrange=peakx[ii]-gaplow;
			maxrange=peakx[ii]+gaphigh;
			minbin=histo[i]->FindBin(minrange);
			maxbin=histo[i]->FindBin(maxrange);
			ibin=minbin;
			while(histo[i]->GetBinContent(ibin)<(peaky[ii]/2))
			{
				ibin++;
				if(ibin>=maxbin)break;
			}
			double sigmaguess=2*(peakx[ii]-histo[i]->GetBinCenter(ibin))/2.355;
			float highcounts=0,lowcounts=0;
			for (jj=0;jj<10;jj++)
			{
				highcounts+=histo[i]->GetBinContent(maxbin-jj);
				lowcounts+=histo[i]->GetBinContent(minbin+jj);
			}
			highcounts=highcounts/10; lowcounts=lowcounts/10;
			double aguess=(highcounts-lowcounts)/(maxrange-minrange);
			double bguess=lowcounts-minrange*aguess;
			cout<<aguess<<"	"<<bguess<<"	"<<sigmaguess<<"	"<<peakx[ii]<<"	"<<peaky[ii]<<endl;
// 			total[ii]=new TF1("total","[0]*x+[1]+[2]/2/[3]*exp(0.5*([4]*[4]/([3]*[3]))+(x-[5])/[3])*ROOT::Math::erfc(1/sqrt(2)*([4]/[3]+(x-[5])/[4]))",peakx[ii]-gaplow,peakx[ii]+gaphigh);//高斯拟合，调节合适的拟合range / Glassman PLB2018 low-energy tail
// 			g[ii]=new TF1("g","gausn",peakx[ii]-gaplow,peakx[ii]+gaphigh);// The [2]-N parameter in total is equivalent to the Constant in gausn
// 			p[ii]=new TF1("p","[0]*x+[1]-[0]*x-[1]+[2]/2/[3]*exp(0.5*([4]*[4]/([3]*[3]))+(x-[5])/[3])*ROOT::Math::erfc(1/sqrt(2)*([4]/[3]+(x-[5])/[4]))",peakx[ii]-gaplow,peakx[ii]+gaphigh);//pure peak
// 			b[ii]=new TF1("b","[0]*x+[1]",peakx[ii]-gaplow,peakx[ii]+gaphigh);//pure bkg
			[0]*x+[1]+sqrt(3.141592654/2)*[2]/[3]*[4]*exp(0.5*[4]*[4]/([3]*[3])+(x-[5])/[3])*ROOT::Math::erfc(1/sqrt(2)*([4]/[3]+(x-[5])/[4]))
			total[ii]=new TF1("total","[0]*x+[1]+sqrt(3.141592654/2)*[2]/[3]*[4]*exp(0.5*[4]*[4]/([3]*[3])+(x-[5])/[3])*ROOT::Math::erfc(1/sqrt(2)*([4]/[3]+(x-[5])/[4]))",peakx[ii]-gaplow,peakx[ii]+gaphigh);//高斯拟合，调节合适的拟合range / Glassman PRC2019 low-energy tail
			derivative[ii]=new TF1("derivative","[0]+exp(0.5*[4]*[4]/([3]*[3])+(x-[5])/[3])/[3]*ROOT::Math::erfc(1/sqrt(2)*([4]/[3]+(x-[5])/[4]))+exp(0.5*[4]*[4]/([3]*[3])+(x-[5])/[3])*(-sqrt(3.141592654/2)/[4])*exp(-(1/sqrt(2)*([4]/[3]+(x-[5])/[4]))*(1/sqrt(2)*([4]/[3]+(x-[5])/[4])))",peakx[ii]-gaplow,peakx[ii]+gaphigh);
			
			g[ii]=new TF1("g","gausn",peakx[ii]-gaplow,peakx[ii]+gaphigh);// The [2]-N parameter in total is equivalent to the Constant in gausn
			p[ii]=new TF1("p","[0]*x+[1]-[0]*x-[1]+sqrt(3.141592654/2)*[2]/[3]*[4]*exp(0.5*[4]*[4]/([3]*[3])+(x-[5])/[3])*ROOT::Math::erfc(1/sqrt(2)*([4]/[3]+(x-[5])/[4]))",peakx[ii]-gaplow,peakx[ii]+gaphigh);//pure peak
			b[ii]=new TF1("b","[0]*x+[1]",peakx[ii]-gaplow,peakx[ii]+gaphigh);//pure bkg
			total[ii]->SetNpx((gaphigh+gaplow)*10);
			g[ii]->SetNpx((gaphigh+gaplow)*10);
			p[ii]->SetNpx((gaphigh+gaplow)*10);
			b[ii]->SetNpx((gaphigh+gaplow)*10);
			// 			total[ii]=new TF1("total","[0]*x+[1]+[2]/2/[3]*exp([4]*[4]/2/[3]/[3]+(x-[5])/[3])*(1-ROOT::Math::erf([4]*[4]+[3]*(x-[5])/sqrt(2)/[4]/[3]))",peakx[ii]-gaplow,peakx[ii]+gaphigh);//高斯拟合，调节合适的拟合range / Bennett
			// 			g[ii]=new TF1("g","gausn",peakx[ii]-gaplow,peakx[ii]+gaphigh);
			// 			p[ii]=new TF1("p","[0]*x+[1]-[0]*x-[1]+[2]/2/[3]*exp([4]*[4]/2/[3]/[3]+(x-[5])/[3])*(1-ROOT::Math::erf([4]*[4]+[3]*(x-[5])/sqrt(2)/[4]/[3]))",peakx[ii]-gaplow,peakx[ii]+gaphigh);
			//  			b[ii]=new TF1("b","[0]*x+[1]",peakx[ii]-gaplow,peakx[ii]+gaphigh);//pure bkg

			//total[ii]->SetParameters(0.1,15,peaky[ii],10,10,peakx[ii]);//initial value [0]-A, [1]-B, [2]-N, [3]-λ, [4]-σ, [5]-μ
			total[ii]->SetParameters(aguess,bguess,peaky[ii],0.7,sigmaguess,peakx[ii]);//initial value [0]-A, [1]-B, [2]-N, [3]-τ, [4]-σ, [5]-μ
//			total[ii]->SetParameters(0,500,150000,0.7,1.5,3480);//initial value [0]-A, [1]-B, [2]-N, [3]-λ, [4]-σ, [5]-μ  //Glassman
// 			total[ii]->SetParLimits(0,-500,500);//Bkg A  //Glassman
// 			total[ii]->SetParLimits(1,-50000,300000);//Bkg B  //Glassman
 			total[ii]->SetParLimits(2,50,500000);//Constant,min,max  //Glassman
			total[ii]->SetParLimits(3,0.1,2);//Tau  //Glassman portal
//			total[ii]->SetParLimits(3,0.7,0.7);//Tau  //Glassman portal
 			total[ii]->SetParLimits(4,0.4,20);//Sigma  //Glassman
 			total[ii]->SetParLimits(5,400,4500);//Mean  //Glassman
			total[ii]->SetParNames("BkgA","BkgB","Const*bin","Tau","Sigma","Mean");
			
			histo[i]->Fit("total","","",peakx[ii]-gaplow,peakx[ii]+gaphigh);///specify a range in the Fit
			//histo[i]->Fit("total","R+");//or restrict the fit to the range specified in the TF1 constructor: peakx[ii]-gaplow,peakx[ii]+gaphigh
			total[ii]->GetParameters(par[ii]);//二维数组的par[ii]是地址,pointer to the TF1, GetParameters的数组得是double类型Obtaining the value of parameters and saving them to par[]; 
			parerr[0]=total[ii]->GetParError(0);//Obtaining the error of the 1st parameter
			parerr[1]=total[ii]->GetParError(1);//Obtaining the error of the 2nd parameter
			parerr[2]=total[ii]->GetParError(2);//Obtaining the error of the 3rd parameter
			parerr[3]=total[ii]->GetParError(3);//Obtaining the error of the 4th parameter
			parerr[4]=total[ii]->GetParError(4);//Obtaining the error of the 5th parameter
			parerr[5]=total[ii]->GetParError(5);//Obtaining the error of the 6th parameter
			parchi[ii]=total[ii]->GetChisquare();
			parNDF[ii]=total[ii]->GetNDF();
			g[ii]->SetParameters(par[ii][2]*par[ii][4]*sqrt(3.141592654*2),par[ii][5],par[ii][4]);//set parameters for drawing gausn
			g[ii]->SetLineColor(4);
			p[ii]->SetParameters(par[ii][0],par[ii][1],par[ii][2],par[ii][3],par[ii][4],par[ii][5]);//set parameters for drawing peak
			p[ii]->SetLineColor(6);
			b[ii]->SetParameters(par[ii][0],par[ii][1]);//set parameters for drawing bkg
			b[ii]->SetLineColor(8);
			total[ii]->Draw("same");
			g[ii]->Draw("same");
			p[ii]->Draw("same");
			b[ii]->Draw("same");

// 			fprintf(outfile,"Constant%d=	%.1f%s%.3f\n",ii,par[ii][2]/binwidth,"	",parerr[2]/binwidth);//par数组还保持着刚才的参数
// 			fprintf(outfile,"Mean%d=	%.1f%s%.3f\n",ii,par[ii][5],"	",parerr[5]);
// 			fprintf(outfile,"Sigma%d=	%.3f%s%.3f\n",ii,par[ii][4],"	",parerr[4]);
// 			//fprintf(outfile,"Res%d=	%.2f%%\n",ii,par[ii][4]/par[ii][5]*2.355*100);
// 			fprintf(outfile,"Tau%d=	%.2f%s%.3f\n",ii,par[ii][3],"	",parerr[3]);
// 			fprintf(outfile,"A%d=	%.2f%s%.3f\nB%d=	%.1f%s%.3f\n",ii,par[ii][0],"	",parerr[0],ii,par[ii][1],"	",parerr[1]);
// 			fprintf(outfile,"Chisquare%d=	%.2f\n",ii,parchi);
// 			fprintf(outfile,"NDF%d=	%.2f\n",ii,parNDF);
// 			fprintf(outfile,"Maximum%d=	%.1f%s\n",ii,peakx[ii]);//par数组还保持着刚才的参数
			outfile<<"SeGA_"<<i<<"	Constant"<<ii<<"=	"<<par[ii][2]<<"	+/-	"<<parerr[2]<<"	Mean"<<ii<<"=	"<<par[ii][5]<<"	+/-	"<<parerr[5]<<"	Maximum"<<ii<<"=	"<<peakx[ii]<<"	Sigma"<<ii<<"=	"<<par[ii][4]<<"	+/-	"<<parerr[4]<<"	Tau"<<ii<<"=	"<<par[ii][3]<<"	+/-	"<<parerr[3]<<"	A"<<ii<<"=	"<<par[ii][0]<<"	+/-	"<<parerr[0]<<"	B"<<ii<<"=	"<<par[ii][1]<<"	+/-	"<<parerr[1]<<"	Chi2"<<ii<<"=	"<<parchi[ii]<<"	NDF"<<ii<<"=	"<<parNDF[ii]<<"	Area"<<ii<<"=	"<<sqrt(3.1415926536*2)*par[ii][2]*par[ii][4]<<endl;//输出文本查看
			//peakx[ii]=par[ii][1];// For graph pol fit, mean is more accurate than peakx
		
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
			sprintf(paraprint,"Constant%d=%.1f%s%.3f",ii,par[ii][2]/binwidth,"+/-",parerr[2]/binwidth);//par数组还保持着刚才的参数
			textgaus->AddText(paraprint);
			sprintf(paraprint,"Mean%d=%.1f%s%.3f",ii,par[ii][5],"+/-",parerr[5]);
			textgaus->AddText(paraprint);
			sprintf(paraprint,"Sigma%d=%.3f%s%.3f",ii,par[ii][4],"+/-",parerr[4]);
			textgaus->AddText(paraprint);
			sprintf(paraprint,"Res%d=%.2f%%",ii,par[ii][4]/par[ii][5]*2.355*100);
			textgaus->AddText(paraprint);
			sprintf(paraprint,"Tau%d=%.2f%s%.3f",ii,par[ii][3],"+/-",parerr[3]);
			textgaus->AddText(paraprint);
			sprintf(paraprint,"A%d(%.2f%s%.3f)*x+B%d(%.1f%s%.3f)",ii,par[ii][0],"+/-",parerr[0],ii,par[ii][1],"+/-",parerr[1]);
			textgaus->AddText(paraprint);
			sprintf(paraprint,"Chisquare%d=%.2f",ii,parchi[ii]);
			textgaus->AddText(paraprint);
			sprintf(paraprint,"NDF%d=%.2f",ii,parNDF[ii]);
			textgaus->AddText(paraprint);
			sprintf(paraprint,"Maximum%d=%.1f%s",ii,peakx[ii]);//par数组还保持着刚才的参数
			textgaus->AddText(paraprint);
			textgaus->Draw();
			sprintf(h_name,"D:/X/out/cali/%s_%.1f.png",hfit_name,par[ii][5]);
			canvaspeak[i]->SaveAs(h_name);//存图
		}//for(ii=0;ii<peaknum;ii++)
		outfile<<"\n\n"<<endl;
		//sprintf(paraprint,"Mean1=%4.2f",Mean1);//gt拟合参数，似乎寻峰不太准，拟合的Mean作峰位更准
		
		
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
		// 		sprintf(hcali_name,"%s%s%s","D:/X/",hcali_name,".png");
		// 		canvascali[i]->SaveAs(hcali_name);//存图
		// 		outfile<<b_name<<i<<" slope="<<slope[i]<<endl;//输出文本查看
		// 		outfile<<b_name<<i<<" intercept="<<intercept[i]<<endl;//输出文本查看
	}//for (i=0;i<ID;i++)
	//}//for(irawroot=runstart; irawroot<=runstop; irawroot++)
}//peakcali main