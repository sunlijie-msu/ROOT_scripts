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
void gausnerfcpol1_peakfit_bp2164_2310()// peak draw-fit (gausn+erfc+pol1) used for fit one peak, used to formally fit β-delayed proton to get the energy and intensity.
{//for 25Si known bp
	//save each peak fit, and tau, sigma as figures
	time_t tim;
	struct tm *at;
	char now[80];
	int irawroot;
	char rawrootname[300];
	double slope[300],slopeerr[300];
	double intercept[300],intercepterr[300];
	const int ID1=0;// i=ID1//which detector
	const int ID2=0;// i<=ID2//which detector

	int binwidth=1;
	//	int binnum=(upperlimit-lowerlimit)/binwidth;
	int minrange=0,maxrange=0,minbin,maxbin;
	double sigma=1,thresh=0.8,sigmab=20;//adjust
	float gaplow=70.,gaphigh=70., peakdiff1=0., peakdiff2=0.;//fitting range随分辨不同调整
	unsigned long i;
	int jj,ii,ibin,k,ipeak;
	char paraprint[30],b_name[200],histo_name[200],h_name[200],hfit_name[200],hcali_name[200],tauflag[200];
	char command1[100],command2[100];
	//string paraprint[30],h_name[80],histo_name[80],b_name[80],hfit_name[80],hcali_name[80];
	TCanvas *canvaspeak[300];
	TCanvas *canvascali[300];
	TH1F *histo[ID2+1];//TH1F peak search+gauss fit,creat histograms
	const int nummax=6;//search peaknumber can be many, fitnumber can be limited to 3//no need to modify
	int peaknum=nummax;
	TF1 *total[nummax];//creat function
	TF1 *p1[nummax], *g1[nummax], *p2[nummax], *g2[nummax], *b[nummax];
	TGraph *graph1[300], *graph2[300];//TGraph刻度
	Double_t peakx[nummax],peakxerr[nummax];
	Double_t peaky[nummax],peakyerr[nummax];
	Double_t sig[nummax],sigerr[nummax];
	Double_t tau[nummax],tauerr[nummax];
	// 	float *peakx;//if you don't know how many peaks will be found, use this
	// 	float *peaky;//if you don't know how many peaks will be found, use this
	double energylit[6]={451.7820, 493.1847, 944.6500, 1460.8200, 1612.6869, 2614.5110};//25Si
	double energyliterr[6]={0.384111, 0.455554, 	0.353553, 0.005000, 0.384111, 0.010000};//25Si
//	double energylit[6]={450.6924, 1460.8200, 1599.5500, 2614.5110, 2904.7794, 7801.4706};//23Al
//	double energyliterr[6]={0.137758, 0.005000, 0.365148, 0.010000, 0.603657, 0.485071};//23Al
	float peakmin;
	float chtemp;
	double par[nummax][5],parerr[nummax][5];//随峰数不同改动par[peaknum][3]//adjust
	double parchi[nummax],parNDF[nummax],p_value[nummax];
 	ifstream infilelowerlimit("D:/X/out/Si25/lowerlimitforresponsefunction25Si.dat",ios::in);//12 SeGA detectors
 	ifstream infileupperlimit("D:/X/out/Si25/upperlimitforresponsefunction25Si.dat",ios::in);//12 SeGA detectors
// 	ifstream infilelowerlimit("D:/X/out/Si25/lowerlimitforresponsefunction25Si_newbg.dat",ios::in);//12 SeGA detectors 
// 	ifstream infileupperlimit("D:/X/out/Si25/upperlimitforresponsefunction25Si_newbg.dat",ios::in);//12 SeGA detectors
	//ifstream infilelowerlimit("D:/X/out/Si25/lowerlimitforresponsefunction23Al.dat",ios::in);//12 SeGA detectors
	//ifstream infileupperlimit("D:/X/out/Si25/upperlimitforresponsefunction23Al.dat",ios::in);//12 SeGA detectors
	int lowerlimit[ID2+1][nummax],upperlimit[ID2+1][nummax];
	for(i=ID1;i<=ID2;i++)//which detector
	{
		if(i==8||i==11||i==14|i==15) continue;
		for(ii=0;ii<nummax;ii++)//which peak
		{
			infilelowerlimit>>lowerlimit[i][ii];
			infileupperlimit>>upperlimit[i][ii];
			lowerlimit[i][ii]=2100;
			upperlimit[i][ii]=2370;
		}
	}
	for(i=ID1;i<=ID2;i++)//which detector
	{
		if(i==8||i==11||i==14|i==15) continue;
		for(ii=0;ii<nummax;ii++)//which peak
		{
			cout<<lowerlimit[i][ii]<<"	";
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
	ofstream outfile("D:/X/out/Si25/outfile/peakcali.dat",ios::out);
	ofstream outfile2("D:/X/out/Si25/outfile/peakcalipara.dat",ios::out);
	//FILE *outfile=fopen ("D:/X/peakcali.dat","a");
	//for(irawroot=runstart; irawroot<=runstop; irawroot++)
	//{
	//if()continue;
//	sprintf(rawrootname,"%s","D:/X/out/Si25/100eV-bin-ungated_sun_25Si.root");//modify formal Eg
	sprintf(rawrootname,"%s","D:/X/out/Moshe/secondAnalysis/protons-si25_final.root");
//	sprintf(rawrootname,"%s","D:/X/out/Si25/100eV-bin-gated_sun_25Si.root");//modify
//	sprintf(rawrootname,"%s","D:/X/out/Si25/100eV-bin-ungated_sun_25Si_10peakscali.root");
//	sprintf(rawrootname,"%s","D:/X/out/Si25/100eV-bin-gated_sun_25Si_10peakscali.root");
	//sprintf(rawrootname,"%s","X:/aaron/rootfiles/run-all_after_cali_25Si.root");
	//sprintf(rawrootname,"%s","X:/aaron/rootfiles/run-all_after_cali_23Al.root");

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
		histo[i]= (TH1F*)fin->Get("h8"); //Get Gamma ray spectrum ungated modify
//		histo[i]= (TH1F*)fin->Get("gbCoincidences"); //Get Gamma ray spectrum gated
		histo[i]->Rebin(5);//rebin 10 will get 1-keV bin
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

	for(i=0;i<=0;i++)//don't change which SeGA detector =0<=15. 0-10 are usually pretty good to fit, 12,13 are troublesome
	{//don't change i=0 to other values
		if(i==8||i==11||i==14||i==15) continue;
		sprintf(hfit_name,"%s%d","fit_proton_3pads_",i);
		canvaspeak[i]=new TCanvas(hfit_name,hfit_name,1100,600);//建立画布
		canvaspeak[i]->cd();//进入画布
		histo[i]->SetTitle(hfit_name);//图名
		histo[i]->GetXaxis()->SetTitle("Energy");//轴名
		histo[i]->GetYaxis()->SetTitle("Counts");//轴名
		histo[i]->GetXaxis()->CenterTitle();//居中
		histo[i]->GetYaxis()->CenterTitle();//居中
		histo[i]->GetXaxis()->SetLabelFont(132);//坐标字体
		histo[i]->GetYaxis()->SetLabelFont(132);//坐标字体
		histo[i]->GetXaxis()->SetTitleFont(132);//轴名字体
		histo[i]->GetYaxis()->SetTitleFont(132);//轴名字体
		histo[i]->GetXaxis()->SetTitleOffset(0.95);//轴名偏移
		histo[i]->GetYaxis()->SetTitleOffset(1.1);//轴名偏移		
		histo[i]->Draw("e");
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
		for(ii=0; ii<=0; ii++)//which peak in one SeGA detector modify =0<=5
		{
			//if(ii==4)continue;
			peaky[ii]=0; peakx[ii]=0; sig[ii]=0; tau[ii]=0;
			histo[i]->GetXaxis()->SetRangeUser(lowerlimit[i][ii],upperlimit[i][ii]);
			peaky[ii]=histo[i]->GetMaximum();
			peakx[ii]=histo[i]->GetBinCenter(histo[i]->GetMaximumBin());
			if(ii==0){	gaplow=100; gaphigh=100; peakdiff1=0; peakdiff2=140;}//25Si 452 old bg modify peakdiff1=0 if peak1 is bigger; peakdiff2=0 if peak2 is bigger
			histo[i]->GetXaxis()->SetRangeUser(peakx[ii]-peakdiff1-gaplow,peakx[ii]+peakdiff2+gaphigh);//zoom the axis
			//cout<<"************"<<peakx[ii]<<"	"<<peaky[ii]<<endl;
			minrange=peakx[ii]-peakdiff1-gaplow;
			maxrange=peakx[ii]+peakdiff2+gaphigh;
			cout<<minrange<<maxrange<<endl;
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
			total[ii]=new TF1("total","[0]*x+[1]+[2]/2/[3]*exp(0.5*([4]*[4]/([3]*[3]))-(x-[5])/[3])*ROOT::Math::erfc(1/sqrt(2)*([4]/[3]-(x-[5])/[4]))+[6]/2/[7]*exp(0.5*([8]*[8]/([7]*[7]))-(x-[9])/[7])*ROOT::Math::erfc(1/sqrt(2)*([8]/[7]-(x-[9])/[8]))",peakx[ii]-peakdiff1-gaplow,peakx[ii]+peakdiff2+gaphigh);// Glassman PLB2018 N=area high-energy tail
			g1[ii]=new TF1("g1","gausn",peakx[ii]-peakdiff1-gaplow,peakx[ii]+peakdiff2+gaphigh);// The [2]-N parameter in total is equivalent to the Constant in gausn
			p1[ii]=new TF1("p1","[0]*x+[1]-[0]*x-[1]+[2]/2/[3]*exp(0.5*([4]*[4]/([3]*[3]))-(x-[5])/[3])*ROOT::Math::erfc(1/sqrt(2)*([4]/[3]-(x-[5])/[4]))",peakx[ii]-peakdiff1-gaplow,peakx[ii]+peakdiff2+gaphigh);//pure peak
			g2[ii]=new TF1("g2","gausn",peakx[ii]-peakdiff1-gaplow,peakx[ii]+peakdiff2+gaphigh);// The [2]-N parameter in total is equivalent to the Constant in gausn
			p2[ii]=new TF1("p2","[0]*x+[1]-[0]*x-[1]+[2]/2/[3]*exp(0.5*([4]*[4]/([3]*[3]))-(x-[5])/[3])*ROOT::Math::erfc(1/sqrt(2)*([4]/[3]-(x-[5])/[4]))",peakx[ii]-peakdiff1-gaplow,peakx[ii]+peakdiff2+gaphigh);//pure peak
			b[ii]=new TF1("b","[0]*x+[1]",peakx[ii]-peakdiff1-gaplow,peakx[ii]+peakdiff2+gaphigh);//pure bkg

// 			total[ii]=new TF1("total","[0]*x+[1]+sqrt(3.141592654/2)*[2]/[3]*[4]*exp(0.5*([4]*[4]/([3]*[3]))+(x-[5])/[3])*ROOT::Math::erfc(1/sqrt(2)*([4]/[3]+(x-[5])/[4]))",peakx[ii]-gaplow,peakx[ii]+gaphigh);// Glassman PRC2019 low-energy tail
// 			g1[ii]=new TF1("g1","gausn",peakx[ii]-gaplow,peakx[ii]+gaphigh);// The [2]-N parameter in total is equivalent to the Constant in gausn
// 			p1[ii]=new TF1("p1","[0]*x+[1]-[0]*x-[1]+sqrt(3.141592654/2)*[2]/[3]*[4]*exp(0.5*([4]*[4]/([3]*[3]))+(x-[5])/[3])*ROOT::Math::erfc(1/sqrt(2)*([4]/[3]+(x-[5])/[4]))",peakx[ii]-gaplow,peakx[ii]+gaphigh);//pure peak
// 			b[ii]=new TF1("b","[0]*x+[1]",peakx[ii]-gaplow,peakx[ii]+gaphigh);//pure bkg
			total[ii]->SetNpx((gaphigh+peakdiff2+gaplow+peakdiff1)*10);
			g1[ii]->SetNpx((gaphigh+peakdiff2+gaplow+peakdiff1)*10);
			p1[ii]->SetNpx((gaphigh+peakdiff2+gaplow+peakdiff1)*10);
			g2[ii]->SetNpx((gaphigh+peakdiff2+gaplow+peakdiff1)*10);
			p2[ii]->SetNpx((gaphigh+peakdiff2+gaplow+peakdiff1)*10);
			b[ii]->SetNpx((gaphigh+peakdiff2+gaplow+peakdiff1)*10);
			// 			total[ii]=new TF1("total","[0]*x+[1]+[2]/2/[3]*exp([4]*[4]/2/[3]/[3]+(x-[5])/[3])*(1-ROOT::Math::erf([4]*[4]+[3]*(x-[5])/sqrt(2)/[4]/[3]))",peakx[ii]-gaplow,peakx[ii]+gaphigh);//高斯拟合，调节合适的拟合range / Bennett
			// 			g1[ii]=new TF1("g1","gausn",peakx[ii]-gaplow,peakx[ii]+gaphigh);
			// 			p1[ii]=new TF1("p1","[0]*x+[1]-[0]*x-[1]+[2]/2/[3]*exp([4]*[4]/2/[3]/[3]+(x-[5])/[3])*(1-ROOT::Math::erf([4]*[4]+[3]*(x-[5])/sqrt(2)/[4]/[3]))",peakx[ii]-gaplow,peakx[ii]+gaphigh);
			//  			b[ii]=new TF1("b","[0]*x+[1]",peakx[ii]-gaplow,peakx[ii]+gaphigh);//pure bkg

			//total[ii]->SetParameters(0.1,15,peaky[ii],10,10,peakx[ii]);//initial value [0]-A, [1]-B, [2]-N, [3]-λ, [4]-σ, [5]-μ
			total[ii]->SetParameters(aguess/2,bguess/2,peaky[ii],2,sigmaguess,peakx[ii]-peakdiff1, peaky[ii],2,sigmaguess,peakx[ii]+peakdiff2);//initial value [0]-A, [1]-B, [2]-N1, [3]-τ1, [4]-σ1, [5]-μ1, [6]-N2, [7]-τ2, [8]-σ2, [9]-μ2
			//			total[ii]->SetParameters(0,500,150000,0.7,1.5,3480);//initial value [0]-A, [1]-B, [2]-N, [3]-λ, [4]-σ, [5]-μ  //Glassman
			// 			total[ii]->SetParLimits(0,-500,500);//Bkg A  //Glassman
			// 			total[ii]->SetParLimits(1,-50000,300000);//Bkg B  //Glassman
			total[ii]->SetParLimits(2,1,5000000);//Constant,min,max  //Glassman
			total[ii]->SetParLimits(3,0.3,55);//Tau  //Glassman portal
			total[ii]->SetParLimits(4,2,29);//Sigma  //Glassman
			total[ii]->SetParLimits(5,2100,2220);//Mean  //Glassman
			total[ii]->SetParLimits(6,1,5000000);//Constant,min,max  //Glassman
			total[ii]->SetParLimits(7,0.3,55);//Tau  //Glassman portal
			total[ii]->SetParLimits(8,2,29);//Sigma  //Glassman
			total[ii]->SetParLimits(9,2250,2350);//Mean  //Glassman
			total[ii]->SetParNames("BkgA","BkgB","Const1*bin","Tau1","Sigma1","Mean1","Const2*bin","Tau2","Sigma2","Mean2");
			//			sprintf(tauflag,"%s","tau0.7");//
			sprintf(tauflag,"%s","taufree");
//			TFitResultPtr r_fit = histo[i]->Fit("total","MLES","",peakx[ii]-gaplow,peakx[ii]+gaphigh);//"S" means the result of the fit is returned in the TFitResultPtr
			//If the option "S" is instead used, TFitResultPtr contains the TFitResult and behaves as a smart pointer to it.
			//The fit parameters, error and chi2 (but not covariance matrix) can be retrieved also from the fitted function.
			histo[i]->Fit("total","MLE","",peakx[ii]-peakdiff1-gaplow,peakx[ii]+peakdiff2+gaphigh);///specify a range in the Fit
			histo[i]->Fit("total","MLE","",peakx[ii]-peakdiff1-gaplow,peakx[ii]+peakdiff2+gaphigh);///specify a range in the Fit
			//histo[i]->Fit("total","R+");//or restrict the fit to the range specified in the TF1 constructor: peakx[ii]-gaplow,peakx[ii]+gaphigh
			total[ii]->GetParameters(par[ii]);//二维数组的par[ii]是地址,pointer to the TF1, GetParameters的数组得是double类型Obtaining the value of parameters and saving them to par[]; 
			parerr[ii][0]=total[ii]->GetParError(0);//Obtaining the error of the 1st parameter
			parerr[ii][1]=total[ii]->GetParError(1);//Obtaining the error of the 2nd parameter
			parerr[ii][2]=total[ii]->GetParError(2);//Obtaining the error of the 3rd parameter
			parerr[ii][3]=total[ii]->GetParError(3);//Obtaining the error of the 4th parameter
			parerr[ii][4]=total[ii]->GetParError(4);//Obtaining the error of the 5th parameter
			parerr[ii][5]=total[ii]->GetParError(5);//Obtaining the error of the 6th parameter
			parerr[ii][6]=total[ii]->GetParError(6);//Obtaining the error of the 7th parameter
			parerr[ii][7]=total[ii]->GetParError(7);//Obtaining the error of the 8th parameter
			parerr[ii][8]=total[ii]->GetParError(8);//Obtaining the error of the 9th parameter
			parerr[ii][9]=total[ii]->GetParError(9);//Obtaining the error of the 10th parameter
			parchi[ii]=total[ii]->GetChisquare();
			parNDF[ii]=total[ii]->GetNDF();
			p_value[ii]=total[ii]->GetProb();
			g1[ii]->SetParameters(par[ii][2]*par[ii][4]*sqrt(3.141592654*2),par[ii][5],par[ii][4]);//set parameters for drawing gausn
			g1[ii]->SetLineColor(4);
			p1[ii]->SetParameters(par[ii][0],par[ii][1],par[ii][2],par[ii][3],par[ii][4],par[ii][5]);//set parameters for drawing peak
			p1[ii]->SetLineColor(4);
			g2[ii]->SetParameters(par[ii][6]*par[ii][8]*sqrt(3.141592654*2),par[ii][9],par[ii][8]);//set parameters for drawing gausn
			g2[ii]->SetLineColor(4);
			p2[ii]->SetParameters(par[ii][0],par[ii][1],par[ii][6],par[ii][7],par[ii][8],par[ii][9]);//set parameters for drawing peak
			p2[ii]->SetLineColor(4);
			b[ii]->SetParameters(par[ii][0],par[ii][1]);//set parameters for drawing bkg
			b[ii]->SetLineColor(8);
			total[ii]->Draw("same");
//			g1[ii]->Draw("same");
			p1[ii]->Draw("same");
//			g2[ii]->Draw("same");
			p2[ii]->Draw("same");
			b[ii]->Draw("same");

		TVirtualFitter * fitter = TVirtualFitter::GetFitter();//The method TVirtualFitter::GetFitter())->Get the parameters of your fitting function after having it fitted to an histogram.
//		TMatrixD cor = r_fit->GetCorrelationMatrix();//get rou
//		cor.Print();
		double Utau_Utau = fitter->GetCovarianceMatrixElement(3,3);//(Utau)^2
		cout<<sqrt(Utau_Utau)<<endl;
		double rou_Utau_Umean = fitter->GetCovarianceMatrixElement(3,5);//(ρ*Utau*Umean), (3,5) (5,3) doesn't matter.
		cout<<rou_Utau_Umean<<endl;
		double Umean_Umean = fitter->GetCovarianceMatrixElement(5,5);//(Umean)^2
		cout<<sqrt(Umean_Umean)<<endl;
		cout<<sqrt(Utau_Utau+Umean_Umean+2*rou_Utau_Umean)<<endl;
			//cout<<"lalalala="<<E_mu_mu<<endl;
			sig[ii]=par[ii][4]; sigerr[ii]=parerr[ii][4];
			tau[ii]=par[ii][3]; tauerr[ii]=parerr[ii][3];
			peaky[ii]=total[ii]->GetMaximum();
			peakx[ii]=total[ii]->GetMaximumX();
			peakxerr[ii]=sqrt(Utau_Utau+Umean_Umean+2*rou_Utau_Umean);
			// 			fprintf(outfile,"Constant%d=	%.1f%s%.3f\n",ii,par[ii][2]/binwidth,"	",parerr[2]/binwidth);//par数组还保持着刚才的参数
			// 			fprintf(outfile,"Mean%d=	%.1f%s%.3f\n",ii,par[ii][5],"	",parerr[5]);
			// 			fprintf(outfile,"Sigma%d=	%.3f%s%.3f\n",ii,par[ii][4],"	",parerr[4]);
			// 			//fprintf(outfile,"Res%d=	%.2f%%\n",ii,par[ii][4]/par[ii][5]*2.355*100);
			// 			fprintf(outfile,"Tau%d=	%.2f%s%.3f\n",ii,par[ii][3],"	",parerr[3]);
			// 			fprintf(outfile,"A%d=	%.2f%s%.3f\nB%d=	%.1f%s%.3f\n",ii,par[ii][0],"	",parerr[0],ii,par[ii][1],"	",parerr[1]);
			// 			fprintf(outfile,"Chisquare%d=	%.2f\n",ii,parchi);
			// 			fprintf(outfile,"NDF%d=	%.2f\n",ii,parNDF);
			// 			fprintf(outfile,"Maximum%d=	%.1f%s\n",ii,peakx[ii]);//par数组还保持着刚才的参数
			outfile<<"SeGA_"<<i<<"	Constant*bin"<<ii<<"=	"<<par[ii][2]<<"	+/-	"<<parerr[ii][2]<<"	Mean"<<ii<<"=	"<<par[ii][5]<<"	+/-	"<<parerr[ii][5]<<"	Maximum"<<ii<<"=	"<<peakx[ii]<<"	+/-	"<<peakxerr[ii]<<"	Sigma"<<ii<<"=	"<<par[ii][4]<<"	+/-	"<<parerr[ii][4]<<"	Tau"<<ii<<"=	"<<par[ii][3]<<"	+/-	"<<parerr[ii][3]<<"	A"<<ii<<"=	"<<par[ii][0]<<"	+/-	"<<parerr[ii][0]<<"	B"<<ii<<"=	"<<par[ii][1]<<"	+/-	"<<parerr[ii][1]<<"	Chi2"<<ii<<"=	"<<parchi[ii]<<"	NDF"<<ii<<"=	"<<parNDF[ii]<<"	Area"<<ii<<"=	"<<sqrt(3.1415926536*2)*par[ii][2]*par[ii][4]<<endl;//输出文本查看
			//peakx[ii]=par[ii][1];// For graph1 pol fit, mean is more accurate than peakx

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
			sprintf(paraprint,"1Constant*bin%d=%.2f%s%.3f",ii,par[ii][2]/binwidth,"+/-",parerr[ii][2]/binwidth);//par数组还保持着刚才的参数
			textgaus->AddText(paraprint);
			sprintf(paraprint,"1Mean%d=%.4f%s%.4f",ii,par[ii][5],"+/-",parerr[ii][5]);
			textgaus->AddText(paraprint);
			sprintf(paraprint,"1Sigma%d=%.3f%s%.3f",ii,par[ii][4],"+/-",parerr[ii][4]);
			textgaus->AddText(paraprint);
			sprintf(paraprint,"1Res%d=%.2f%%",ii,par[ii][4]/par[ii][5]*2.355*100);
			textgaus->AddText(paraprint);
			sprintf(paraprint,"1Tau%d=%.2f%s%.3f",ii,par[ii][3],"+/-",parerr[ii][3]);
			textgaus->AddText(paraprint);

			sprintf(paraprint,"2Constant*bin%d=%.2f%s%.3f",ii,par[ii][6]/binwidth,"+/-",parerr[ii][6]/binwidth);//par数组还保持着刚才的参数
			textgaus->AddText(paraprint);
			sprintf(paraprint,"2Mean%d=%.4f%s%.4f",ii,par[ii][9],"+/-",parerr[ii][9]);
			textgaus->AddText(paraprint);
			sprintf(paraprint,"2Sigma%d=%.3f%s%.3f",ii,par[ii][8],"+/-",parerr[ii][8]);
			textgaus->AddText(paraprint);
			sprintf(paraprint,"2Res%d=%.2f%%",ii,par[ii][8]/par[ii][9]*2.355*100);
			textgaus->AddText(paraprint);
			sprintf(paraprint,"2Tau%d=%.2f%s%.3f",ii,par[ii][7],"+/-",parerr[ii][7]);
			textgaus->AddText(paraprint);

			sprintf(paraprint,"A%d(%.4f%s%.4f)*x+B%d(%.3f%s%.3f)",ii,par[ii][0],"+/-",parerr[ii][0],ii,par[ii][1],"+/-",parerr[ii][1]);
			textgaus->AddText(paraprint);
			sprintf(paraprint,"Chisquare%d=%.2f",ii,parchi[ii]);
			textgaus->AddText(paraprint);
			sprintf(paraprint,"NDF%d=%.0f",ii,parNDF[ii]);
			textgaus->AddText(paraprint);
			sprintf(paraprint,"p1-val%d=%.5f",ii,p_value[ii]);
			textgaus->AddText(paraprint);
			sprintf(paraprint,"Maximum%d=%.4f%s%.4f",ii,peakx[ii],"+/-",peakxerr[ii]);//par数组还保持着刚才的参数
			textgaus->AddText(paraprint);
			textgaus->Draw();
			//sprintf(h_name,"D:/X/out/response23Al/%s%s%d_%s.png",hfit_name,"_peak",ii,tauflag);//modify
			sprintf(h_name,"D:/X/out/bg/%s%s%d_%s.png",hfit_name,"_peak",ii,tauflag);
			canvaspeak[i]->SaveAs(h_name);//存图
		}//for(ii=0;ii<peaknum;ii++)
		outfile<<"\n\n"<<endl;
		//sprintf(paraprint,"Mean1=%4.2f",Mean1);//gt拟合参数，似乎寻峰不太准，拟合的Mean作峰位更准

		// 		sprintf(hcali_name,"%s%d","sigma_Rn_all_SeGA_",i);
		// 		canvascali[i]=new TCanvas(hcali_name,hcali_name,900,600);//建立画布
		// 		canvascali[i]->cd();//进入画布
		// 		//canvascali[i]->SetGrid();//显示网格
		// //		graph1[i]=new TGraph(peaknum,peakx,sig);//TGraph *gr1=new TGraph(n,x,y);
		// 		//float eenergy[2]={200,280};
		// 		//float epeakch[2]={0,0};
		// 		graph1[i]= new TGraphErrors(peaknum,peakx,sig,peakxerr,sigerr);//画error bars TGraph(n,x,y,ex,ey);
		// 		graph1[i]->SetTitle(hcali_name);
		// 		graph1[i]->GetXaxis()->SetTitle("Energy");//轴名
		// 		graph1[i]->GetYaxis()->SetTitle("Sigma");//轴名
		// 		graph1[i]->GetXaxis()->CenterTitle();//居中
		// 		graph1[i]->GetYaxis()->CenterTitle();//居中
		// 		graph1[i]->GetXaxis()->SetLabelFont(132);//坐标字体
		// 		graph1[i]->GetYaxis()->SetLabelFont(132);//坐标字体
		// 		graph1[i]->GetXaxis()->SetTitleFont(132);//轴名字体
		// 		graph1[i]->GetYaxis()->SetTitleFont(132);//轴名字体
		// 		//graph1[i]->GetYaxis()->SetLabelSize(0.05);//坐标字号
		// 		//graph1[i]->GetYaxis()->SetTitleSize(0.05);//轴名字号
		// 		graph1[i]->GetXaxis()->SetTitleOffset(1.2);//轴名偏移
		// 		graph1[i]->GetYaxis()->SetTitleOffset(1.3);//轴名偏移
		// 		graph1[i]->SetMarkerStyle(21);
		// 		graph1[i]->SetMarkerColor(1);
		// 		TF1 *pol1=new TF1("pol1","pol1",0,50000);//多项式拟合，调节合适的拟合range
		// 		graph1[i]->Fit("pol1");//pol1 can be used directly without TF1 constructor in CINT
		// 		intercept[i]=pol1->GetParameter(0);
		// 		slope[i]=pol1->GetParameter(1);
		// 		intercepterr[i]=pol1->GetParError(0);
		// 		slopeerr[i]=pol1->GetParError(1);
		// 		graph1[i]->Draw("AP"); //A-Axis around the graph1,AP is suitable
		// 		TPaveText *textpol1 = new TPaveText(0.13,0.74,0.55,0.89,"brNDC");//left, down, right, up
		// 		textpol1->SetBorderSize(1);
		// 		textpol1->SetFillColor(0);
		// 		textpol1->SetTextAlign(12);//align = 10*HorizontalAlign + VerticalAlign, 12 means水平左对齐、垂直居中对齐
		// 		textpol1->SetTextFont(132);//font = 10 * fontID + precision, 12+2,12 means Symbol; 13+2, 13 means News Times Roman
		// 		sprintf(paraprint,"slope=%.9f+/-%.9f",slope[i],slopeerr[i]);
		// 		textpol1->AddText(paraprint);
		// 		sprintf(paraprint,"intercept=%.5f+/-%.5f",intercept[i],intercepterr[i]);
		// 		textpol1->AddText(paraprint);
		// 		textpol1->Draw();
		// 		//sprintf(hcali_name,"%s.png",hcali_name);
		// 		sprintf(hcali_name,"%s%s%s%s%s","D:/X/out/response23Al/",hcali_name,"_",tauflag,".png");//modify
		// //		sprintf(hcali_name,"%s%s%s%s%s","D:/X/out/response/",hcali_name,"_",tauflag,".png");
		// 		canvascali[i]->SaveAs(hcali_name);//存图
		// 		outfile2<<"SeGA_"<<i<<"	sigslope=	"<<slope[i]<<"	+/-	"<<slopeerr[i]<<"	sigintercept=	"<<intercept[i]<<"	+/-	"<<intercepterr[i]<<endl;
		// 
		// 		sprintf(hcali_name,"%s%d","tau_Rn_all_SeGA_",i);
		// 		canvascali[i]=new TCanvas(hcali_name,hcali_name,900,600);//建立画布
		// 		canvascali[i]->cd();//进入画布
		// 		//canvascali[i]->SetGrid();//显示网格
		// //		graph1[i]=new TGraph(peaknum,peakx,tau);//TGraph *gr1=new TGraph(n,x,y);
		// 		//float eenergy[2]={200,280};
		// 		//float epeakch[2]={0,0};
		// 		graph2[i]= new TGraphErrors(peaknum,peakx,tau,peakxerr,tauerr);//画error bars TGraph(n,x,y,ex,ey);
		// 		graph2[i]->SetTitle(hcali_name);
		// 		graph2[i]->GetXaxis()->SetTitle("Energy");//轴名
		// 		graph2[i]->GetYaxis()->SetTitle("Tau");//轴名
		// 		graph2[i]->GetXaxis()->CenterTitle();//居中
		// 		graph2[i]->GetYaxis()->CenterTitle();//居中
		// 		graph2[i]->GetXaxis()->SetLabelFont(132);//坐标字体
		// 		graph2[i]->GetYaxis()->SetLabelFont(132);//坐标字体
		// 		graph2[i]->GetXaxis()->SetTitleFont(132);//轴名字体
		// 		graph2[i]->GetYaxis()->SetTitleFont(132);//轴名字体
		// 		//graph1[i]->GetYaxis()->SetLabelSize(0.05);//坐标字号
		// 		//graph1[i]->GetYaxis()->SetTitleSize(0.05);//轴名字号
		// 		graph2[i]->GetXaxis()->SetTitleOffset(1.2);//轴名偏移
		// 		graph2[i]->GetYaxis()->SetTitleOffset(1.3);//轴名偏移
		// 		graph2[i]->SetMarkerStyle(21);
		// 		graph2[i]->SetMarkerColor(1);
		// 		TF1 *pol1=new TF1("pol1","pol1",0,50000);//多项式拟合，调节合适的拟合range
		// 		graph2[i]->Fit("pol1");//pol1 can be used directly without TF1 constructor in CINT
		// 		intercept[i]=pol1->GetParameter(0);
		// 		slope[i]=pol1->GetParameter(1);
		// 		intercepterr[i]=pol1->GetParError(0);
		// 		slopeerr[i]=pol1->GetParError(1);
		// 		graph2[i]->Draw("AP"); //A-Axis around the graph1,AP is suitable
		// 		TPaveText *textpol1 = new TPaveText(0.13,0.74,0.55,0.89,"brNDC");//left, down, right, up
		// 		textpol1->SetBorderSize(1);
		// 		textpol1->SetFillColor(0);
		// 		textpol1->SetTextAlign(12);//align = 10*HorizontalAlign + VerticalAlign, 12 means水平左对齐、垂直居中对齐
		// 		textpol1->SetTextFont(132);//font = 10 * fontID + precision, 12+2,12 means Symbol; 13+2, 13 means News Times Roman
		// 		sprintf(paraprint,"slope=%.9f+/-%.9f",slope[i],slopeerr[i]);
		// 		textpol1->AddText(paraprint);
		// 		sprintf(paraprint,"intercept=%.5f+/-%.5f",intercept[i],intercepterr[i]);
		// 		textpol1->AddText(paraprint);
		// 		textpol1->Draw();
		// 		//sprintf(hcali_name,"%s.png",hcali_name);
		// 		sprintf(hcali_name,"%s%s%s%s%s","D:/X/out/response23Al/",hcali_name,"_",tauflag,".png");//modify
		// //		sprintf(hcali_name,"%s%s%s%s%s","D:/X/out/response/",hcali_name,"_",tauflag,".png");
		// 		canvascali[i]->SaveAs(hcali_name);//存图
		// 		outfile2<<"SeGA_"<<i<<"	tauslope=	"<<slope[i]<<"	+/-	"<<slopeerr[i]<<"	tauintercept=	"<<intercept[i]<<"	+/-	"<<intercepterr[i]<<endl;
		//save each peak fit, and tau, sigma
	}//for (i=0;i<ID;i++)
	//}//for(irawroot=runstart; irawroot<=runstop; irawroot++)
}//peakcali main