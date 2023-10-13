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
#include "TMatrixDSymfwd.h"
#include "TRandom.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH1D.h"
#include "TH2.h"
#include "TF1.h"
#include "TMath.h"
#include "TPaveText.h"
#include "TSpline.h"
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
void graphpol6fit_LBayesian_DSL_varyTau()// graph-pol6 fit Likelihood of Bayesian
{
	time_t tim;
	struct tm *at;
	char now[80];
	int irawroot;
	char rawrootname[300];
	double slope[300],slopeerr[300];
	double intercept[300],intercepterr[300];
	double p7[300],p6[300],p5[300],p4[300],p3[300],p2[300],p1[300],p0[300];//for output
	double p7err[300],p6err[300],p5err[300],p4err[300],p3err[300],p2err[300],p1err[300],p0err[300];
	const int ID1=0;// no need to change i=ID1//which detector
	const int ID2=0;// no need to change i<=ID2//which detector

	int binwidth=1;
	//	int binnum=(upperlimit-lowerlimit)/binwidth;
	int minrange=0,maxrange=0,minbin,maxbin;
	double sigma=1,thresh=0.8,sigmab=20;//adjust
	double gaplow=70.,gaphigh=70.;//fitting range随分辨不同调整//no need to change
	unsigned long i;
	int jj,ii,k,ipeak;
	char paraprint[100],o_name[200],b_name[200],histo_name[200],h_name[200],hfit_name[200],hcali_name[200],tauflag[200];
	char command1[100],command2[100];
	//string paraprint[30],h_name[80],histo_name[80],b_name[80],hfit_name[80],hcali_name[80];
	TCanvas *canvaspeak[300];
	TCanvas *canvascali[300];
	TH1F *histo[ID2+1];//TH1F peak search+gauss fit,creat histograms
	const int nummax=100;//search peaknumber can be many, fitnumber can be limited to 3//no need to change
	int peaknum=nummax;
	TF1 *total[nummax];//creat function
	TF1 *p[nummax], *g[nummax], *b[nummax];
	TGraph *graph1[300], *graph1bandhigh[300], *graph1bandlow[300];//TGraph
	TGraph *graph2[300], *graph2bandhigh[300], *graph2bandlow[300];//TGraph
	TGraph *graph3[300], *graph3bandhigh[300], *graph3bandlow[300];//TGraph
	TGraph *graph4[300], *graph4bandhigh[300], *graph4bandlow[300];//TGraph
	double Xvalue[ID2+1][nummax],meanerr[ID2+1][nummax];
	double peaky[nummax],peakyerr[nummax];
	double Yvalue[ID2+1][nummax],Chierr[ID2+1][nummax];
	double XvalueSum[ID2+1],YvalueSum[ID2+1];
	double tau[ID2+1][nummax],tauerr[ID2+1][nummax];
	double range_sig1[ID2+1], range_sig2[ID2+1], range_tau1[ID2+1], range_tau2[ID2+1];
	// 	double *Xvalue;//if you don't know how many peaks will be found, use this
	// 	double *peaky;//if you don't know how many peaks will be found, use this
//	double energylit[6]={451.7, 493.3, 944.9, 1460.820, 1612.4, 2614.511};//25Si
//	double energyliterr[6]={0.5, 0.7, 0.5, 0.005, 0.5, 0.010};//25Si
// 	double energylit[6]={450.7, 1460.820, 1599, 2614.511, 2908, 7801};//23Al
// 	double energyliterr[6]={0.15, 0.005, 2, 0.010, 3, 2};//23Al

// 	const int Npoint=4;
// 	double energy_point[Npoint] = {1369,2754,2870,4238};// set location of point for single value
// 	double err_point[Npoint];  // error on the function at point x0 for single value
	
// 	const int Npoint=230;
// 	double energy_point[Npoint];// set location of point for drawing
// 	double err_point[Npoint];  // error on the function at point x0 for drawing
// 	for(jj=0;jj<Npoint;jj++){	energy_point[jj]=jj*40;}

	double peakmin;
	double chtemp;
	double par[nummax][5],parerr[nummax][5];//随峰数不同改动par[peaknum][3]
	double parchi[ID2+1],parNDF[ID2+1],p_value[ID2+1];
	int Eg=1248; //DSLmodify always finite value
	//1248, 2234, 3076, 4971, 5156, 5141, 4156, 4045, 3435, 2186, 2838, 4270, 3541, 5294
	//  49,     47,      46,     42,     42,     39,      39,     39,      45,      45,     44,     41,     40,      39
	int Ea=0; // no need to change
	double Egv=0; // no need to change
	char flag[20]; // no need to change
	sprintf(flag,"base"); // no need to change
	double ibin=0,Lowbin=0,Centralbin=0,Highbin=0;
	double totalarea=0, Larea=0;
	double istart=0, istop=0, istep=0.01;

	ipeak=2; //choose a peak DSLmodify
	if (ipeak==0) { Eg=1248; Ea=49; Egv=1248; minrange=700; maxrange=5700; istart=0; istop=5700; istep=1; }
	if (ipeak==1) { Eg=2234; Ea=47; Egv=2234; minrange=120; maxrange=500; istart=80; istop=500; istep=1; }
	if (ipeak==2) { Eg=3076; Ea=46; Egv=3076; minrange=0; maxrange=20; istart=0; istop=50; istep=0.1; }
	if (ipeak==3) { Eg=4971; Ea=42; Egv=4970; minrange=0; maxrange=30; istart=0; istop=50; istep=0.1; }
	if (ipeak==4) { Eg=5156; Ea=42; Egv=5156; minrange=0; maxrange=30; istart=0; istop=50; istep=0.1; }
	if (ipeak==5) { Eg=5141; Ea=39; Egv=5141; minrange=0; maxrange=40; istart=0; istop=50; istep=0.1; }
	if (ipeak==6) { Eg=4156; Ea=39; Egv=4156; minrange=0; maxrange=30; istart=0; istop=50; istep=0.1; }
	if (ipeak==7) { Eg=3435; Ea=45; Egv=3435; minrange=1; maxrange=50; istart=0; istop=50; istep=0.1; }
	if (ipeak==8) { Eg=2186; Ea=45; Egv=2186; minrange=0; maxrange=30; istart=0; istop=50; istep=0.1; }
	if (ipeak==9) { Eg=2838; Ea=44; Egv=2838; minrange=0; maxrange=30; istart=0; istop=50; istep=0.1; }
	if (ipeak==10) { Eg=4270; Ea=41; Egv=4270; minrange=0; maxrange=50; istart=0; istop=50; istep=0.1; }
	if (ipeak==11) { Eg=3541; Ea=40; Egv=3541; minrange=0; maxrange=30; istart=0; istop=50; istep=0.1; }
	if (ipeak==12) { Eg=5294; Ea=39; Egv=5293; minrange=0; maxrange=50; istart=0; istop=50; istep=0.1; }

	string line;
//	ifstream infile("C:/Si24/Si22peakcali/Clover.dat",ios::in);//The data that need to be fitted
//	int lowerlimit[ID2+1][nummax],upperlimit[ID2+1][nummax];
	sprintf(o_name,"%s%d%s%d%s%.0f%s","D:/X/out/DSL/UL/LBayesian_Gamma",Eg,"_Ea",Ea,"_Egv",Egv,"slice_vTau.dat");
	cout<<o_name<<endl;
	ofstream outfile2(o_name,ios::out);
	for(int irun=15;irun<=15;irun++)//loop each Gamma's Tau to do LBayesian Integration
	{// DSLmodify choose a flag
		if (irun==0) sprintf(flag,"base");
		if (irun==1) sprintf(flag,"bkglow");
		if (irun==2) sprintf(flag,"bkghigh");
		if (irun==3) sprintf(flag,"SP0.9");
		if (irun==4) sprintf(flag,"SP1.1");
		if (irun==5) sprintf(flag,"AC2-0.8");
		if (irun==6) sprintf(flag,"respsigma1");
		if (irun==7) sprintf(flag,"respsigma-1");
		if (irun==8) sprintf(flag,"resptau1");
		if (irun==9) sprintf(flag,"resptau-1");
		if (irun==10) sprintf(flag,"SP0.95");
		if (irun==11) sprintf(flag,"SP1.05");
		if (irun==12) sprintf(flag,"bkg");
		if (irun==13) sprintf(flag,"SP");
		if (irun==14) sprintf(flag,"AC2");
		if (irun==15) sprintf(flag,"coll14");//modify

	for(i=ID1;i<=ID2;i++)//which detector //no need to change
	{
// 		if(i==8||i==11||i==14|i==15) continue;
// 		if(i==5)peaknum=nummax-1;//
// 		if(i==0)peaknum=nummax-4;

//		if(Eg==1369) sprintf(b_name,"%s%s%s%s%s%d%s%d%s%d%s%d%s%d%s","D:/X/out/",histo_name,"/Chi2_gated_",histo_name,"_",Eg,"_st",Flag_st,"_sig",Flag_sig,"_tau",Flag_tau,"_order",Flag_or,".dat");//temporary

//	  	if(Eg==1369) sprintf(b_name,"%s%s%s%s%s%d%s%d%s%d%s%d%s","D:/X/out/",histo_name,"/Chi2_gated_",histo_name,"_",Eg,"_st",Flag_st,"_sig",Flag_sig,"_tau",Flag_tau,"_August.dat");//formal 4bgcali
//	  	if(Eg!=1369) sprintf(b_name,"%s%s%s%s%s%d%s%d%s%d%s%d%s","D:/X/out/",histo_name,"/Chi2_ungated_",histo_name,"_",Eg,"_st",Flag_st,"_sig",Flag_sig,"_tau",Flag_tau,".dat");//formal 4bgcali
		
//		sprintf(b_name,"%s%d%s%d%s","D:/X/out/DSL/Chi2_Ea",Ea,"Eg",Eg,"projection_vTau.dat");//DSLmodify
		sprintf(b_name,"%s%d%s%d%s%.0f%s%s%s","D:/X/out/DSL/LBayesian/LBayesian_Gamma",Eg,"_Ea",Ea,"_Egv",Egv,"slice_vTau_",flag,".dat");
		cout<<b_name<<endl;
		ifstream infile(b_name,ios::in);//The data that need to be fitted
		ii=0; XvalueSum[i]=0; YvalueSum[i]=0;
		while ( getline(infile, line) )
		{
			stringstream(line)>>Xvalue[i][ii]>>Yvalue[i][ii];
			cout<<Xvalue[i][ii]<<'	'<<Yvalue[i][ii]<<endl;
			YvalueSum[i]+=Yvalue[i][ii];
			ii++;
		}
		cout<<ii<<endl;
		peaknum=ii-1;// the last line is empty and should be discarded, or it will add a (0,0) point to the graph
	}
	// 	const int Npoint=4;
	// 	double energy_point[Npoint] = {1369,2754,2870,4238};// set location of point for single value
	// 	double err_point[Npoint];  // error on the function at point x0 for single value
	//sprintf(tauflag,"%s","taufree");
	//ofstream outfile("D:/X/out/Si25/outfile/peakcali.dat",ios::out);
	//ofstream outfile2("D:/X/out/Si25/outfile/peakcalipara.dat",ios::app);





	for(i=0;i<=0;i++)//which SeGA detector =0<=15
	{
// 		if(i==8||i==11||i==14|i==15) continue;
// 		if(i==5)peaknum=nummax-1;//
// 		if(i==0)peaknum=nummax-4;


		//************ sigma fit by linear func *****************************************

		sprintf(hcali_name,"%s%d%s%d%s%.0f%s%s","LBayesian_Gamma",Eg,"_Ea",Ea,"_Egv",Egv,"slice_vTau_",flag);
		canvascali[i]=new TCanvas(hcali_name,hcali_name,900,600);//建立画布
		canvascali[i]->cd();//进入画布
		//canvascali[i]->SetGrid();//显示网格
		for (ii=0; ii<peaknum; ii++)
		{
			Yvalue[i][ii]=Yvalue[i][ii]/YvalueSum[i];//normalized p.d.f. Without normalization, the Yvalues are too small to be correctly smoothed
		}
		graph1[i]=new TGraph(peaknum,Xvalue[i],Yvalue[i]);//TGraph *gr1=new TGraph(n,x,y);
		//graph1[i]= new TGraphErrors(peaknum,Xvalue[i],Yvalue[i],meanerr[i],Chierr[i]);//画error bars TGraph(n,x,y,ex,ey);
		graph1[i]->SetTitle(hcali_name);
		graph1[i]->SetName(hcali_name);
		//graph1[i]->GetXaxis()->SetTitle("Lifetime (fs)");//轴名
		graph1[i]->GetXaxis()->SetTitle("#tau (fs)");//轴名
		graph1[i]->GetYaxis()->SetTitle("Likelihood");//轴名
		graph1[i]->GetXaxis()->CenterTitle();//居中
		graph1[i]->GetYaxis()->CenterTitle();//居中
		graph1[i]->GetXaxis()->SetLabelFont(132);//坐标字体
		graph1[i]->GetYaxis()->SetLabelFont(132);//坐标字体
		graph1[i]->GetXaxis()->SetTitleFont(132);//轴名字体
		graph1[i]->GetYaxis()->SetTitleFont(132);//轴名字体
		//graph1[i]->GetYaxis()->SetLabelSize(0.05);//坐标字号
		//graph1[i]->GetYaxis()->SetTitleSize(0.05);//轴名字号
		graph1[i]->GetXaxis()->SetTitleOffset(1.2);//轴名偏移
		graph1[i]->GetYaxis()->SetTitleOffset(1.3);//轴名偏移
		graph1[i]->SetMarkerStyle(21);
		graph1[i]->SetMarkerColor(1);
		gPad->SetRightMargin(0.03);
		graph1[i]->GetXaxis()->SetRangeUser(minrange,maxrange);
		//graph1[i]->Draw("AP");
		
		TSpline3 *s3 = new TSpline3("s3",Xvalue[i],Yvalue[i],peaknum);
		s3->SetLineColor(kRed);
		graph1[i]->SetLineWidth(2);
 		graph1[i]->Draw("AL");
		s3->SetLineWidth(2);
 		s3->Draw("l same");

		totalarea=0; Larea=0;
		if (ipeak==0||ipeak==1||ipeak==10||ipeak==12) // get finite lifetime values
		{
			for (ibin=istart; ibin<=istop; ibin+=istep)
			{
				totalarea+=graph1[i]->Eval(ibin)*istep; //Get total area under pdf curve
			}
			for (ibin=istart; ibin<=istop; ibin+=istep)
			{
				Larea+=graph1[i]->Eval(ibin)*istep;
				if (Larea>=totalarea*0.17) break; // Get 10% CL
			}
			Lowbin=ibin; outfile2<<ibin;
			Larea=0;
			for (ibin=istart; ibin<=istop; ibin+=istep)
			{
				Larea+=graph1[i]->Eval(ibin)*istep;
				if (Larea>=totalarea*0.5) break; // Get 10% CL
			}
			Centralbin=ibin; outfile2<<"	"<<ibin;
			Larea=0;
			for (ibin=istart; ibin<=istop; ibin+=istep)
			{
				Larea+=graph1[i]->Eval(ibin)*istep;
				if (Larea>=totalarea*0.83) break; // Get 10% CL
			}
			Highbin=ibin; outfile2<<"	"<<ibin<<endl;
		}


		totalarea=0; Larea=0;
		if (ipeak!=0&&ipeak!=1&&ipeak!=10&&ipeak!=12) // get upper limits
		{
			for (ibin=istart; ibin<=istop; ibin+=istep)
			{
				totalarea+=graph1[i]->Eval(ibin)*istep; //Get total area under pdf curve
			}

			for (ibin=istart; ibin<=istop; ibin+=istep)
			{
				Larea+=graph1[i]->Eval(ibin)*istep;
				if (Larea>=totalarea*0.1) break; // Get 10% CL
			}
			Lowbin=ibin; outfile2<<ibin;
			Larea=0;
			for (ibin=istart; ibin<=istop; ibin+=istep)
			{
				Larea+=graph1[i]->Eval(ibin)*istep;
				if (Larea>=totalarea*0.5) break; // Get 10% CL
			}
			Centralbin=ibin; outfile2<<"	"<<ibin;
			Larea=0;
			for (ibin=istart; ibin<=istop; ibin+=istep)
			{
				Larea+=graph1[i]->Eval(ibin)*istep;
				if (Larea>=totalarea*0.9) break; // Get 10% CL
			}
			Highbin=ibin; outfile2<<"	"<<ibin<<endl;
		}




		TPaveText *textpol6 = new TPaveText(0.50,0.75,0.97,0.90,"brNDC");//left, down, right, up
		textpol6->SetBorderSize(1);
		textpol6->SetFillColor(0);
		textpol6->SetTextAlign(12);//align = 10*HorizontalAlign + VerticalAlign, 12 means水平左对齐、垂直居中对齐
		textpol6->SetTextFont(132);//font = 10 * fontID + precision, 12+2,12 means Symbol; 13+2, 13 means News Times Roman
		if (ipeak==0||ipeak==1)
			sprintf(paraprint,"#tau=%.0f%s%.0f%s%.0f%s",Centralbin,"+",Highbin-Centralbin,"-",Centralbin-Lowbin," fs");
		if (ipeak==10||ipeak==12)
			sprintf(paraprint,"#tau=%.1f%s%.1f%s%.1f%s",Centralbin,"+",Highbin-Centralbin,"-",Centralbin-Lowbin," fs");
		if (ipeak!=0&&ipeak!=1&&ipeak!=10&&ipeak!=12)
			sprintf(paraprint,"90%% UL=%.1f%s",Highbin," fs"); // modify
		textpol6->AddText(paraprint);
		textpol6->Draw();
		//sprintf(hcali_name,"%s.png",hcali_name);
		sprintf(hcali_name,"%s%s%s","D:/X/out/Chi2/",hcali_name,".png");
		canvascali[i]->SaveAs(hcali_name);//存图

	}//for (i=0;i<ID;i++)
	}//irun
	//}//for(irawroot=runstart; irawroot<=runstop; irawroot++)
}//peakcali main