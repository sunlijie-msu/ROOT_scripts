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
void Fillgraph()// read data from a txt file and draw a plot
{
	time_t tim;
	struct tm *at;
	char now[80];
	int irawroot;
	char rawrootname[300];
	double slope[300],slopeerr[300];
	double intercept[300],intercepterr[300];
	double p2[300],p1[300],p0[300];
	double p2err[300],p1err[300],p0err[300];
	const int ID1=0;// no need to change i=ID1//which detector
	const int ID2=0;// no need to change i<=ID2//which detector

	int binwidth=1;
	//	int binnum=(upperlimit-lowerlimit)/binwidth;
	int minrange=0,maxrange=0,minbin,maxbin;
	double sigma=1,thresh=0.8,sigmab=20;//adjust
	double gaplow=70.,gaphigh=70.;//fitting range随分辨不同调整//no need to change
	unsigned long i;
	int jj,ii,ibin,k,ipeak;
	char paraprint[100],b_name[200],histo_name[200],h_name[200],hfit_name[200],hcali_name[200],tauflag[200];
	char command1[100],command2[100];
	//string paraprint[30],h_name[80],histo_name[80],b_name[80],hfit_name[80],hcali_name[80];
	TCanvas *canvaspeak[300];
	TCanvas *canvascali[300];
	TH1F *histo[ID2+1];//TH1F peak search+gauss fit,creat histograms
	const int nummax=100;//search peaknumber can be many, fitnumber can be limited to 3//no need to change
	int peaknum=0;
	TF1 *total[nummax];//creat function
	TF1 *p[nummax], *g[nummax], *b[nummax];
	TGraph *graph1[300], *graph1bandhigh[300], *graph1bandlow[300];//TGraph
	TGraph *graph2[300], *graph2bandhigh[300], *graph2bandlow[300];//TGraph
	TGraph *graph3[300], *graph3bandhigh[300], *graph3bandlow[300];//TGraph
	TGraph *graph4[300], *graph4bandhigh[300], *graph4bandlow[300];//TGraph
	double Number[ID2+1][nummax],Numbererr[ID2+1][nummax];
	double peaky[nummax],peakyerr[nummax];
	double input_value[ID2+1][nummax],input_error[ID2+1][nummax];
	double tau[ID2+1][nummax],tauerr[ID2+1][nummax];
	double range_sig1[ID2+1], range_sig2[ID2+1], range_tau1[ID2+1], range_tau2[ID2+1];
	// 	double *Number;//if you don't know how many peaks will be found, use this
	// 	double *peaky;//if you don't know how many peaks will be found, use this
//	double energylit[6]={451.7, 493.3, 944.9, 1460.820, 1612.4, 2614.511};//25Si
//	double energyliterr[6]={0.5, 0.7, 0.5, 0.005, 0.5, 0.010};//25Si
// 	double energylit[6]={450.7, 1460.820, 1599, 2614.511, 2908, 7801};//23Al
// 	double energyliterr[6]={0.15, 0.005, 2, 0.010, 3, 2};//23Al

// 	const int Npoint=4;
// 	double channel_point[Npoint] = {1369,2754,2870,4238};// set location of point for single value
// 	double err_point[Npoint];  // error on the function at point x0 for single value
	
	const int Npoint=230;
	double channel_point[Npoint];// set location of point for drawing
	double err_point[Npoint];  // error on the function at point x0 for drawing
	for(jj=0;jj<Npoint;jj++){	channel_point[jj]=jj*40;}

	double peakmin;
	double chtemp;
	double par[nummax][5],parerr[nummax][5];//随峰数不同改动par[peaknum][3]//adjust
	double parchi[ID2+1],parNDF[ID2+1],p_value[ID2+1];
	int pe=452;//modify 1369, 2754, 2870, 4238
	int Flag_sig = 0;//default =0
	int Flag_tau = 0;//default =0
	int Flag_st = 0;//default =0
	int Flag_or = 0;
	sprintf(histo_name,"%s","Robertson");//modify
//	sprintf(histo_name,"%s","Thomas");
	string line;
//	ifstream infile("C:/Si24/Si22peakcali/Clover.dat",ios::in);//The data that need to be fitted
//	int lowerlimit[ID2+1][nummax],upperlimit[ID2+1][nummax];
	for(int ior=0;ior<1;ior++)//<1 only show baseline modify, <7 show baseline+-σ, τ, sp
	{
// 		if(ior==0){	Flag_sig=0; Flag_tau=0; Flag_st=0; }//baseline
// 		if(ior==1){	Flag_sig=1; Flag_tau=0; Flag_st=0; }
// 		if(ior==2){	Flag_sig=-1; Flag_tau=0; Flag_st=0; }
// 		if(ior==3){	Flag_sig=0; Flag_tau=1; Flag_st=0; }
// 		if(ior==4){	Flag_sig=0; Flag_tau=-1; Flag_st=0; }
// 		if(ior==5){	Flag_sig=0; Flag_tau=0; Flag_st=1; }
// 		if(ior==6){	Flag_sig=0; Flag_tau=0; Flag_st=-1; }
	for(i=ID1;i<=ID2;i++)//which detector //don't change
	{
// 		if(i==8||i==11||i==14|i==15) continue;
// 		if(i==5)peaknum=nummax-1;//
// 		if(i==0)peaknum=nummax-4;

		//if(pe==1369) sprintf(b_name,"%s%s%s%s%s%d%s%d%s%d%s%d%s%d%s","D:/X/out/",histo_name,"/Chi2_gated_",histo_name,"_",pe,"_st",Flag_st,"_sig",Flag_sig,"_tau",Flag_tau,"_order",Flag_or,".dat");//temporary

// 		if(pe==1369) sprintf(b_name,"%s%s%s%s%s%d%s%d%s%d%s%d%s","D:/X/out/",histo_name,"/Chi2_gated_",histo_name,"_",pe,"_st",Flag_st,"_sig",Flag_sig,"_tau",Flag_tau,".dat");//formal bgcali
// 		if(pe!=1369) sprintf(b_name,"%s%s%s%s%s%d%s%d%s%d%s%d%s","D:/X/out/",histo_name,"/Chi2_ungated_",histo_name,"_",pe,"_st",Flag_st,"_sig",Flag_sig,"_tau",Flag_tau,".dat");//formal bgcali

// 		if(pe==1369) sprintf(b_name,"%s%s%s%s%s%d%s%d%s%d%s%d%s","D:/X/out/",histo_name,"/bkgcali_Chi2_gated_",histo_name,"_",pe,"_st",Flag_st,"_sig",Flag_sig,"_tau",Flag_tau,".dat");//for bkgcali modify
// 		if(pe!=1369) sprintf(b_name,"%s%s%s%s%s%d%s%d%s%d%s%d%s","D:/X/out/",histo_name,"/bkgcali_Chi2_ungated_",histo_name,"_",pe,"_st",Flag_st,"_sig",Flag_sig,"_tau",Flag_tau,".dat");//for bkgcali
		sprintf(b_name,"%s","C:/Si24/Si22peakcali/Xu/inputforgraph1.dat");//formal bgcali
		cout<<b_name<<endl;
		ifstream infile(b_name,ios::in);//The data that need to be fitted
		ii=0;
		while ( getline(infile, line) )
		{
			stringstream(line)>>input_value[i][ii]>>input_error[i][ii];
			cout<<input_value[i][ii]<<'	'<<input_error[i][ii]<<endl;
			ii++;
		}
		peaknum=ii;
		for(int jj=0;jj<peaknum;jj++)//set x values 0,1,2,3..., doesn't affect the result
		{
			Number[i][jj]=jj+1;
			Numbererr[i][jj]=0;
		}
	}
	
	// 	const int Npoint=4;
	// 	double channel_point[Npoint] = {1369,2754,2870,4238};// set location of point for single value
	// 	double err_point[Npoint];  // error on the function at point x0 for single value

	//ofstream outfile("D:/X/out/Si25/outfile/peakcali.dat",ios::out);
	ofstream outfile2("D:/X/out/Si25/outfile/peakcalipara.dat",ios::app);




	for(i=ID1;i<=ID2;i++)//which SeGA detector =0<=15
	{
// 		if(i==8||i==11||i==14|i==15) continue;
// 		if(i==5)peaknum=nummax-1;//
// 		if(i==0)peaknum=nummax-4;


		//************ sigma fit by linear func *****************************************

		sprintf(hcali_name,"%s","Graph");
		canvascali[i]=new TCanvas(hcali_name,hcali_name,900,700);//建立画布
		canvascali[i]->cd();//进入画布
		//canvascali[i]->SetGrid();//显示网格
		//graph1[i]=new TGraph(peaknum,Number[i],input_value[i]);//TGraph *gr1=new TGraph(n,x,y);
		graph1[i]= new TGraphErrors(peaknum,Number[i],input_value[i],Numbererr[i],input_error[i]);//画error bars TGraph(n,x,y,ex,ey);
		graph1[i]->SetTitle(hcali_name);
		graph1[i]->GetXaxis()->SetTitle("Ex (keV)");//轴名
		graph1[i]->GetYaxis()->SetTitle("{B(GT)");//轴名
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
//		graph1[i]->GetYaxis()->SetRangeUser(1,8);
		graph1[i]->GetXaxis()->SetNdivisions(100);//n = n1 + 100*n2 + 10000*n3
		graph1[i]->SetMarkerStyle(21);
		graph1[i]->SetMarkerColor(1);
		graph1[i]->SetLineWidth(2);
		TF1 *pol0=new TF1("pol0","pol0",0,5000);//多项式拟合，调节合适的拟合range
		pol0->SetNpx(5000);
		pol0->SetParNames("p0");//y=p2*x^2+p1*x+p0
		pol0->SetLineWidth(3);
		//Robertson initial guess from Excel bgcali 8peaks cali modify
// 		if(pe==1369){	pol0->SetParameter(2,5640.8557);	pol0->SetParameter(1,-15444380.6684);	pol0->SetParameter(0,10571486307.3844);	}//y = 5640.8557x2 - 15444380.6684x + 10571486307.3844
// 		if(pe==2754){	pol0->SetParameter(2,159.5469);		pol0->SetParameter(1,-878949.9169);		pol0->SetParameter(0,1210543062.7130);		}//y= 159.5469x2 - 878949.9169x + 1210543062.7130
// 		if(pe==2870){	pol0->SetParameter(2,7.1306);			pol0->SetParameter(1,-40936.7678);			pol0->SetParameter(0,58754924.3600);			}//y= 7.1306x2 - 40936.7678x + 58754924.3600
// 		if(pe==4238){	pol0->SetParameter(2,19.5575);		pol0->SetParameter(1,-165780.8106);		pol0->SetParameter(0,351314059.5030);		}//y=  19.5575x2 - 165780.8106x + 351314059.5030
		//bkgcali 10peaks cali modify
// 		if(pe==1369){	pol0->SetParameter(2,5661.2324);	pol0->SetParameter(1,-15495626.6938);	pol0->SetParameter(0,10603453845.6451);	}//y = 5661.2324x2 - 15495626.6938x + 10603453845.6451
// 		if(pe==2754){	pol0->SetParameter(2,161.8199);		pol0->SetParameter(1,-891253.9988);		pol0->SetParameter(0,1227188537.1820);		}//y= 161.8199x2 - 891253.9988x + 1227188537.1820
// 		if(pe==2870){	pol0->SetParameter(2,6.8576);			pol0->SetParameter(1,-39360.5329);			pol0->SetParameter(0,56479622.4639);			}//y= 6.8576x2 - 39360.5329x + 56479622.4639
// 		if(pe==4238){	pol0->SetParameter(2,19.8617);		pol0->SetParameter(1,-168322.3824);		pol0->SetParameter(0,356622896.6320);		}//y=  19.8617x2 - 168322.3824x + 356622896.6320
		//Thomas initial guess from Excel modify
// 		if(pe==2754){	pol0->SetParameter(2,127.62);		pol0->SetParameter(1,-702956.54);		pol0->SetParameter(0,968041242.88);	}//y= 127.62x2 - 702956.54x + 968041242.88
// 		if(pe==1369){	pol0->SetParameter(2,4971.43);		pol0->SetParameter(1,-13606663.06);		pol0->SetParameter(0,9310267562.46);	}//y= 4971.43x2 - 13606663.06x + 9310267562.46
// 		if(pe==4238){	pol0->SetParameter(2,13.14);		pol0->SetParameter(1,-111379.13);		pol0->SetParameter(0,236028923.12);	}//y= 13.14x2 - 111379.13x + 236028923.12
// 		if(pe==2870){	pol0->SetParameter(2,5.66);		pol0->SetParameter(1,-32483.20);		pol0->SetParameter(0,46616203.77);	}//y= 5.66x2 - 32483.20x + 46616203.77

		//graph1[i]->Fit("pol1");//pol1 can be used directly without TF1 constructor in CINT
		//graph1[i]->Fit("pol0","ME");		graph1[i]->Fit("pol0","ME");		graph1[i]->Fit("pol0","ME");		graph1[i]->Fit("pol0","ME");
//		TFitResultPtr r_chi2 = graph1[i]->Fit("pol0","MES");//"S" means the result of the fit is returned in the TFitResultPtr
		//“E” Perform better errors estimation using the Minos technique.
		//“M” Improve fit results, by using the IMPROVE algorithm of TMinuit.
		//If the option "S" is instead used, TFitResultPtr contains the TFitResult and behaves as a smart pointer to it.
		//The fit parameters, error and chi2 (but not covariance matrix) can be retrieved also from the fitted function.
//		TH1D *hint_chi2 = new TH1D("hint_chi2", "Fitted exponential with conf.band", 40000, -10, 30);//Create a histogram to hold the confidence intervals
//		(TVirtualFitter::GetFitter())->GetConfidenceIntervals(hint_chi2, 0.683);//By default the intervals are corrected using the chi2/ndf value of the fit if a chi2 fit is performed
		//hint_chi2 will contain the CL result that you can draw on top of your fitted graph.
		//where hint_chi2 will hold the errors and could superimpose it on the same canvas where you plot central values.
		//The method TVirtualFitter::GetFitter())->GetConfidenceIntervals(hint_chi2, 0.95) computes the confidence level of your model (fitting) function after having it fitted to an histogram.
//		hint_chi2->SetStats(kFALSE);
//		hint_chi2->SetFillColor(kRed-9);
//		TMatrixDSym cov = r_chi2->GetCovarianceMatrix();//useless for now
//		r_chi2->GetConfidenceIntervals(Npoint, 1, 1, channel_point, err_point, 0.683, false);//(Number of x points, 1, 1, x, err, confidence level, false); norm is a flag to control if the intervals need to be normalized to the chi2/ndf value.
// 		for(ii=0;ii<Npoint;ii++)
// 		{
// 			outfile<<"SeGA_"<<i<<"	Eg=	"<<channel_point[ii]<<"	input_value=	"<<pol1->Eval(channel_point[ii])<<"	err_sig=	"<<err_point[ii]<<endl;
// 		}
		graph1[i]->Draw("AP");//"A": Axis are drawn around the graph, "P": The current marker is plotted at each point
//		hint_chi2->Draw("e3 same");

		double errylow[Npoint],erryhigh[Npoint];

// 		for(ii=0;ii<Npoint;ii++)
// 		{
// 			//x[ii]=channel_point[ii];
// 			//y[ii]=pol1->Eval(channel_point[ii]);//is equal to y[ii]=p1[i]*x[ii]+p0[i];
// 			//dy[ii]=sqrt(x[ii]*x[ii]*p1err[i]*p1err[i]+p0err[i]*p0err[i]+p1[i]*p1[i]*dx[ii]*dx[ii]);
// 			//dy[ii]=err_point[ii];
// 			//cout<<"	"<<channel_point[ii]<<"	"<<err_point[ii]<<endl;
// 			errylow[ii]=pol0->Eval(channel_point[ii])-err_point[ii];
// 			erryhigh[ii]=pol0->Eval(channel_point[ii])+err_point[ii];
// 		}
// 		graph1bandlow[i]= new TGraph(Npoint,channel_point,errylow);//画error bars TGraph(n,x,y,ex,ey);
// 		graph1bandlow[i]->SetLineColor(2);
// 		graph1bandlow[i]->SetLineWidth(1);
// 		graph1bandlow[i]->Draw("C");//"C" A smooth Curve is drawn
// 		graph1bandhigh[i]= new TGraph(Npoint,channel_point,erryhigh);//画error bars TGraph(n,x,y,ex,ey);
// 		graph1bandhigh[i]->SetLineColor(2);
// 		graph1bandhigh[i]->SetLineWidth(1);
// 		graph1bandhigh[i]->Draw("C");//"C" A smooth Curve is drawn
//		graph1[i]->Draw("P same");//draw the points again above error band

		p0[i]=pol0->GetParameter(0);
		// 		p1[i]=pol0->GetParameter(1);
		// 		p2[i]=pol0->GetParameter(2);
		p0err[i]=pol0->GetParError(0);
		// 		p1err[i]=pol0->GetParError(1);
		// 		p2err[i]=pol0->GetParError(2);
		parchi[i]=pol0->GetChisquare();
		parNDF[i]=pol0->GetNDF();
		p_value[i]=pol0->GetProb();//This probability is not the “probability that your fit is good.” If you did many fake experiments (draw many random samples of data points from the assumed distribution (your fit function)), this is the percentage of experiments that would give χ2 values ≥ to the one you got in this experiment.

// 		double x,y,z,dx,xlow,xhigh;
// 		x=-p1[i]/(2*p2[i]);
// 		y=p2[i]*x[ii]*x[ii]+p1[i]*x[ii]+p0[i];
// 		z=p0[i]-y-1;
// 		xlow=(-p1[i]-sqrt(p1[i]*p1[i]-4*p2[i]*z))/(2*p2[i]);
// 		xhigh=(-p1[i]+sqrt(p1[i]*p1[i]-4*p2[i]*z))/(2*p2[i]);
// 		dx=x-xlow;

		TPaveText *textpol2 = new TPaveText(0.32,0.10,0.90,0.35,"brNDC");//left, down, right, up
		textpol2->SetBorderSize(1);
		textpol2->SetFillColor(0);
		textpol2->SetTextAlign(12);//align = 10*HorizontalAlign + VerticalAlign, 12 means水平左对齐、垂直居中对齐
		textpol2->SetTextFont(132);//font = 10 * fontID + precision, 12+2,12 means Symbol; 13+2, 13 means News Times Roman
// 		sprintf(paraprint,"y=p2*x^2+p1*x+p0");
// 		textpol2->AddText(paraprint);
// 		sprintf(paraprint,"p2=%.9f+/-%.9f",p2[i],p2err[i]);
// 		textpol2->AddText(paraprint);
// 		sprintf(paraprint,"p1=%.5f+/-%.5f",p1[i],p1err[i]);
// 		textpol2->AddText(paraprint);
		sprintf(paraprint,"Chi2 = %.2f",parchi[i]);
		textpol2->AddText(paraprint);
		sprintf(paraprint,"NDF = %.0f",parNDF[i]);
		textpol2->AddText(paraprint);
		sprintf(paraprint,"p-val = %.2f",p_value[i]);
		textpol2->AddText(paraprint);
		sprintf(paraprint,"Average = %.4f+/-%.4f",p0[i],p0err[i]);
		textpol2->AddText(paraprint);
		sprintf(paraprint,"Inflation factor = sqrt(%.2f/%.0f)=%.2f",parchi[i],parNDF[i],sqrt(parchi[i]/parNDF[i]));
		textpol2->AddText(paraprint);
		sprintf(paraprint,"Average with inflation = %.4f+/-%.4f",p0[i],p0err[i]*sqrt(parchi[i]/parNDF[i]));
		textpol2->AddText(paraprint);
// 		sprintf(paraprint,"ymin=%.5f",y);
// 		textpol2->AddText(paraprint);
// 		sprintf(paraprint,"Eg=%.5f+/-%.5f",x,dx);
// 		textpol2->AddText(paraprint);
//		textpol2->Draw();
		//sprintf(hcali_name,"%s.png",hcali_name);
		sprintf(hcali_name,"%s%s%s","D:/X/out/Chi2/",hcali_name,".png");
		canvascali[i]->SaveAs(hcali_name);//存图
		//outfile2<<"y=p2*x^2+p1*x+p0	chip2=	"<<setprecision(11)<<p2[i]<<"	+/-	"<<p2err[i]<<"	chip1=	"<<setprecision(8)<<p1[i]<<"	+/-	"<<p1err[i]<<"	chip0=	"<<p0[i]<<"	+/-	"<<p0err[i]<<"	Chi2=	"<<parchi[i]<<"	ndf=	"<<parNDF[i]<<"	p-value=	"<<p_value[i]<<"	Chi2min=	"<<y<<setprecision(9)<<"	Eg=	"<<x<<"	+/-	"<<dx<<endl;

/*




		//************ tau fit by sqrt + const func *****************************************

		sprintf(hcali_name,"%s%d","Tau of SeGA",i);
		canvascali[i]=new TCanvas(hcali_name,hcali_name,900,600);//建立画布
		canvascali[i]->cd();//进入画布
		//canvascali[i]->SetGrid();//显示网格
		//		graph1[i]=new TGraph(peaknum,Number,tau);//TGraph *gr1=new TGraph(n,x,y);
		graph3[i]= new TGraphErrors(peaknum,Number[i],tau[i],Numbererr[i],tauerr[i]);//画error bars TGraph(n,x,y,ex,ey);
		graph3[i]->SetTitle(hcali_name);
		graph3[i]->GetXaxis()->SetTitle("Energy");//轴名
		graph3[i]->GetYaxis()->SetTitle("Tau");//轴名
		graph3[i]->GetXaxis()->CenterTitle();//居中
		graph3[i]->GetYaxis()->CenterTitle();//居中
		graph3[i]->GetXaxis()->SetLabelFont(132);//坐标字体
		graph3[i]->GetYaxis()->SetLabelFont(132);//坐标字体
		graph3[i]->GetXaxis()->SetTitleFont(132);//轴名字体
		graph3[i]->GetYaxis()->SetTitleFont(132);//轴名字体
		//graph1[i]->GetYaxis()->SetLabelSize(0.05);//坐标字号
		//graph1[i]->GetYaxis()->SetTitleSize(0.05);//轴名字号
		graph3[i]->GetXaxis()->SetTitleOffset(1.2);//轴名偏移
		graph3[i]->GetYaxis()->SetTitleOffset(1.3);//轴名偏移
		graph3[i]->GetYaxis()->SetRangeUser(range_tau1[i],range_tau2[i]);
		graph3[i]->SetMarkerStyle(21);
		graph3[i]->SetMarkerColor(1);

		TF1 *func_sqrt_const=new TF1("func_sqrt_const","[0]+[1]*sqrt(x)",0,60000);//多项式拟合，调节合适的拟合range
		func_sqrt_const->SetParNames("p0","p1");//y=p1*sqrt(x)+p0
//		graph3[i]->Fit("func_sqrt_const");//pol1 can be used directly without TF1 constructor in CINT

		TFitResultPtr r_tau_sqrt_const = graph3[i]->Fit("func_sqrt_const","S");//"S" means the result of the fit is returned in the TFitResultPtr
		//If the option "S" is instead used, TFitResultPtr contains the TFitResult and behaves as a smart pointer to it.
		//The fit parameters, error and chi2 (but not covariance matrix) can be retrieved also from the fitted function.
		TH1D *hint_tau_sqrt_const = new TH1D("hint_tau_sqrt_const", "Fitted exponential with conf.band", 9000, 0, 9000);//Create a histogram to hold the confidence intervals
		(TVirtualFitter::GetFitter())->GetConfidenceIntervals(hint_tau_sqrt_const, 0.683);//By default the intervals are corrected using the chi2/ndf value of the fit if a chi2 fit is performed
		//hint_tau_sqrt_const will contain the CL result that you can draw on top of your fitted graph.
		//where hint_chi2 will hold the errors and could superimpose it on the same canvas where you plot central values.
		//The method TVirtualFitter::GetFitter())->GetConfidenceIntervals(hint_chi2, 0.95) computes the confidence level of your model (fitting) function after having it fitted to an histogram.
		hint_tau_sqrt_const->SetStats(kFALSE);
		hint_tau_sqrt_const->SetFillColor(kRed-9);
		TMatrixDSym cov = r_tau_sqrt_const->GetCovarianceMatrix();//useless for now
		r_tau_sqrt_const->GetConfidenceIntervals(Npoint, 1, 1, channel_point, err_point, 0.683, false);//(Number of x points, 1, 1, x, err, confidence level, norm); norm is a flag to control if the intervals need to be normalized to the chi2/ndf value.
		for(ii=0;ii<Npoint;ii++)
		{
			outfile<<"SeGA_"<<i<<"	Eg=	"<<channel_point[ii]<<"	tau=	"<<func_sqrt_const->Eval(channel_point[ii])<<"	err_tau=	"<<err_point[ii]<<endl;
		}
		graph3[i]->Draw("AP");//"A": Axis are drawn around the graph, "P": The current marker is plotted at each point
		hint_tau_sqrt_const->Draw("e3 same");
		graph3[i]->Draw("P same");//draw the points again above error band

		p0[i]=func_sqrt_const->GetParameter(0);
		p0err[i]=func_sqrt_const->GetParError(0);
		p1[i]=func_sqrt_const->GetParameter(1);
		p1err[i]=func_sqrt_const->GetParError(1);
		parchi[i]=func_sqrt_const->GetChisquare();
		parNDF[i]=func_sqrt_const->GetNDF();
		p_value[i]=func_sqrt_const->GetProb();//This probability is not the “probability that your fit is good.” If you did many fake experiments (draw many random samples of data points from the assumed distribution (your fit function)), this is the percentage of experiments that would give χ2 values ≥ to the one you got in this experiment.
		double errylow[Npoint],erryhigh[Npoint];
		for(ii=0;ii<Npoint;ii++)
		{
			//x[ii]=channel_point[ii];
			//y[ii]=func_sqrt_const->Eval(channel_point[ii]);//is equal to y[ii]=p1[i]*sqrt(x[ii])+p0[i];
			//dy[ii]=sqrt(sqrt(x[ii])*sqrt(x[ii])*p1err[i]*p1err[i]+p0err[i]*p0err[i]);
			//dy[ii]=err_point[ii];
			//cout<<"	"<<channel_point[ii]<<"	"<<err_point[ii]<<endl;
			errylow[ii]=func_sqrt_const->Eval(channel_point[ii])-err_point[ii];
			erryhigh[ii]=func_sqrt_const->Eval(channel_point[ii])+err_point[ii];
		}
// 		graph3bandlow[i]= new TGraph(Npoint,channel_point,errylow);//画error bars TGraph(n,x,y,ex,ey);
// 		graph3bandlow[i]->SetLineColor(2);
// 		graph3bandlow[i]->SetLineWidth(1);
// 		graph3bandlow[i]->Draw("C");
// 		graph3bandhigh[i]= new TGraph(Npoint,channel_point,erryhigh);//画error bars TGraph(n,x,y,ex,ey);
// 		graph3bandhigh[i]->SetLineColor(2);
// 		graph3bandhigh[i]->SetLineWidth(1);
// 		graph3bandhigh[i]->Draw("C");

		TPaveText *textpol3 = new TPaveText(0.10,0.64,0.34,0.89,"brNDC");//left, down, right, up
		textpol3->SetBorderSize(1);
		textpol3->SetFillColor(0);
		textpol3->SetTextAlign(12);//align = 10*HorizontalAlign + VerticalAlign, 12 means水平左对齐、垂直居中对齐
		textpol3->SetTextFont(132);//font = 10 * fontID + precision, 12+2,12 means Symbol; 13+2, 13 means News Times Roman
		sprintf(paraprint,"y=p1*sqrt(x)+p0");
		textpol3->AddText(paraprint);
		sprintf(paraprint,"p1=%.5f+/-%.5f",p1[i],p1err[i]);
		textpol3->AddText(paraprint);
		sprintf(paraprint,"p0=%.5f+/-%.5f",p0[i],p0err[i]);
		textpol3->AddText(paraprint);
		sprintf(paraprint,"Chisquare=%.2f",parchi[i]);
		textpol3->AddText(paraprint);
		sprintf(paraprint,"NDF=%.0f",parNDF[i]);
		textpol3->AddText(paraprint);
		sprintf(paraprint,"p-value=%e",p_value[i]);
		textpol3->AddText(paraprint);
		textpol3->Draw();
		//sprintf(hcali_name,"%s.png",hcali_name);
		sprintf(hcali_name,"%s%s%s%s%s","D:/X/out/responseforDoppler/",hcali_name,"_",tauflag,"_band.png");
		canvascali[i]->SaveAs(hcali_name);//存图
		outfile2<<"	SeGA_"<<i<<"	y=p1*sqrt(x)+p0	taup1=	"<<setprecision(8)<<p1[i]<<"	+/-	"<<p1err[i]<<"	taup0=	"<<p0[i]<<"	+/-	"<<p0err[i]<<"	Chi2=	"<<parchi[i]<<"	ndf=	"<<parNDF[i]<<"	p-value=	"<<p_value[i]<<endl;





		//************ tau fit by quadratic func *****************************************

		sprintf(hcali_name,"%s%d","AlSi_tau_pol2_SeGA_",i);
		canvascali[i]=new TCanvas(hcali_name,hcali_name,900,600);//建立画布
		canvascali[i]->cd();//进入画布
		//canvascali[i]->SetGrid();//显示网格
		//		graph1[i]=new TGraph(peaknum,Number,tau);//TGraph *gr1=new TGraph(n,x,y);
		graph2[i]= new TGraphErrors(peaknum,Number[i],tau[i],Numbererr[i],tauerr[i]);//画error bars TGraph(n,x,y,ex,ey);
		graph2[i]->SetTitle(hcali_name);
		graph2[i]->GetXaxis()->SetTitle("Energy");//轴名
		graph2[i]->GetYaxis()->SetTitle("Tau");//轴名
		graph2[i]->GetXaxis()->CenterTitle();//居中
		graph2[i]->GetYaxis()->CenterTitle();//居中
		graph2[i]->GetXaxis()->SetLabelFont(132);//坐标字体
		graph2[i]->GetYaxis()->SetLabelFont(132);//坐标字体
		graph2[i]->GetXaxis()->SetTitleFont(132);//轴名字体
		graph2[i]->GetYaxis()->SetTitleFont(132);//轴名字体
		//graph1[i]->GetYaxis()->SetLabelSize(0.05);//坐标字号
		//graph1[i]->GetYaxis()->SetTitleSize(0.05);//轴名字号
		graph2[i]->GetXaxis()->SetTitleOffset(1.2);//轴名偏移
		graph2[i]->GetYaxis()->SetTitleOffset(1.3);//轴名偏移
		graph2[i]->GetYaxis()->SetRangeUser(range_tau1[i],range_tau2[i]);
		graph2[i]->SetMarkerStyle(21);
		graph2[i]->SetMarkerColor(1);
		TF1 *pol0=new TF1("pol0","pol0",0,60000);//多项式拟合，调节合适的拟合range
		pol0->SetParNames("p0","p1","p2");//y=p2*x^2+p1*x+p0
//		graph2[i]->Fit("pol0");//pol1 can be used directly without TF1 constructor in CINT

		TFitResultPtr r_tau_pol2 = graph2[i]->Fit("pol0","S");//"S" means the result of the fit is returned in the TFitResultPtr
		//If the option "S" is instead used, TFitResultPtr contains the TFitResult and behaves as a smart pointer to it.
		//The fit parameters, error and chi2 (but not covariance matrix) can be retrieved also from the fitted function.
		TH1D *hint_tau_pol2 = new TH1D("hint_tau_pol2", "Fitted exponential with conf.band", 9000, 0, 9000);//Create a histogram to hold the confidence intervals
		(TVirtualFitter::GetFitter())->GetConfidenceIntervals(hint_tau_pol2, 0.683);//By default the intervals are corrected using the chi2/ndf value of the fit if a chi2 fit is performed
		//hint_tau_pol2 will contain the CL result that you can draw on top of your fitted graph.
		//where hint_chi2 will hold the errors and could superimpose it on the same canvas where you plot central values.
		//The method TVirtualFitter::GetFitter())->GetConfidenceIntervals(hint_chi2, 0.95) computes the confidence level of your model (fitting) function after having it fitted to an histogram.
		hint_tau_pol2->SetStats(kFALSE);
		hint_tau_pol2->SetFillColor(kRed-9);
		TMatrixDSym cov = r_tau_pol2->GetCovarianceMatrix();//useless for now
		r_tau_pol2->GetConfidenceIntervals(Npoint, 1, 1, channel_point, err_point, 0.683, false);//(Number of x points, 1, 1, x, err, confidence level, false); norm is a flag to control if the intervals need to be normalized to the chi2/ndf value.
		outfile<<"Eg=	"<<channel_point[3]<<"	tau=	"<<pol0->Eval(channel_point[3])<<"	err_tau=	"<<err_point[3]<<endl;
		graph2[i]->Draw("AP");//"A": Axis are drawn around the graph, "P": The current marker is plotted at each point
		hint_tau_pol2->Draw("e3 same");
		graph2[i]->Draw("P same");//draw the points again above error band

		p0[i]=pol0->GetParameter(0);
		p1[i]=pol0->GetParameter(1);
		p2[i]=pol0->GetParameter(2);
		p0err[i]=pol0->GetParError(0);
		p1err[i]=pol0->GetParError(1);
		p2err[i]=pol0->GetParError(2);
		parchi[i]=pol0->GetChisquare();
		parNDF[i]=pol0->GetNDF();
		p_value[i]=pol0->GetProb();//This probability is not the “probability that your fit is good.” If you did many fake experiments (draw many random samples of data points from the assumed distribution (your fit function)), this is the percentage of experiments that would give χ2 values ≥ to the one you got in this experiment.
		double errylow[Npoint],erryhigh[Npoint];
		for(ii=0;ii<Npoint;ii++)
		{
			//x[ii]=channel_point[ii];
			//y[ii]=pol0->Eval(channel_point[ii]);//is equal to y[ii]=p2[i]*x[ii]*x[ii]+p1[i]*x[ii]+p0[i];
			//dy[ii]=sqrt(x[ii]*x[ii]*x[ii]*x[ii]*p2err[i]*p2err[i]+x[ii]*x[ii]*p1err[i]*p1err[i]+p0err[i]*p0err[i]+(2*p2[i]+p1[i])*(2*p2[i]+p1[i])*dx[ii]);
			//dy[ii]=err_point[ii];
			//cout<<"	"<<channel_point[ii]<<"	"<<err_point[ii]<<endl;
			errylow[ii]=pol0->Eval(channel_point[ii])-err_point[ii];
			erryhigh[ii]=pol0->Eval(channel_point[ii])+err_point[ii];
		}
		graph2bandlow[i]= new TGraph(Npoint,channel_point,errylow);//画error bars TGraph(n,x,y,ex,ey);
		graph2bandlow[i]->SetLineColor(2);
		graph2bandlow[i]->SetLineWidth(1);
		graph2bandlow[i]->Draw("C");
		graph2bandhigh[i]= new TGraph(Npoint,channel_point,erryhigh);//画error bars TGraph(n,x,y,ex,ey);
		graph2bandhigh[i]->SetLineColor(2);
		graph2bandhigh[i]->SetLineWidth(1);
		graph2bandhigh[i]->Draw("C");

		TPaveText *textpol2 = new TPaveText(0.10,0.64,0.34,0.89,"brNDC");//left, down, right, up
		textpol2->SetBorderSize(1);
		textpol2->SetFillColor(0);
		textpol2->SetTextAlign(12);//align = 10*HorizontalAlign + VerticalAlign, 12 means水平左对齐、垂直居中对齐
		textpol2->SetTextFont(132);//font = 10 * fontID + precision, 12+2,12 means Symbol; 13+2, 13 means News Times Roman
		sprintf(paraprint,"y=p2*x^2+p1*x+p0");
		textpol2->AddText(paraprint);
		sprintf(paraprint,"p2=%.9f+/-%.9f",p2[i],p2err[i]);
		textpol2->AddText(paraprint);
		sprintf(paraprint,"p1=%.5f+/-%.5f",p1[i],p1err[i]);
		textpol2->AddText(paraprint);
		sprintf(paraprint,"p0=%.5f+/-%.5f",p0[i],p0err[i]);
		textpol2->AddText(paraprint);
		sprintf(paraprint,"Chi2=%.2f",parchi[i]);
		textpol2->AddText(paraprint);
		sprintf(paraprint,"NDF=%.0f",parNDF[i]);
		textpol2->AddText(paraprint);
		sprintf(paraprint,"p-val=%e",p_value[i]);
		textpol2->AddText(paraprint);
		textpol2->Draw();
		//sprintf(hcali_name,"%s.png",hcali_name);
		sprintf(hcali_name,"%s%s%s%s%s","D:/X/out/responseforDoppler/",hcali_name,"_",tauflag,"_band.png");
		canvascali[i]->SaveAs(hcali_name);//存图
		outfile2<<"	SeGA_"<<i<<"	y=p2*x^2+p1*x+p0	taup2=	"<<setprecision(11)<<p2[i]<<"	+/-	"<<p2err[i]<<"	taup1=	"<<setprecision(8)<<p1[i]<<"	+/-	"<<p1err[i]<<"	taup0=	"<<p0[i]<<"	+/-	"<<p0err[i]<<"	Chi2=	"<<parchi[i]<<"	ndf=	"<<parNDF[i]<<"	p-value=	"<<p_value[i];





		//************ tau fit by sqrt func *****************************************

		sprintf(hcali_name,"%s%d","AlSi_tau_sqrt_SeGA_",i);
		canvascali[i]=new TCanvas(hcali_name,hcali_name,900,600);//建立画布
		canvascali[i]->cd();//进入画布
		//canvascali[i]->SetGrid();//显示网格
		//		graph1[i]=new TGraph(peaknum,Number,tau);//TGraph *gr1=new TGraph(n,x,y);
		graph4[i]= new TGraphErrors(peaknum,Number[i],tau[i],Numbererr[i],tauerr[i]);//画error bars TGraph(n,x,y,ex,ey);
		graph4[i]->SetTitle(hcali_name);
		graph4[i]->GetXaxis()->SetTitle("Energy");//轴名
		graph4[i]->GetYaxis()->SetTitle("Tau");//轴名
		graph4[i]->GetXaxis()->CenterTitle();//居中
		graph4[i]->GetYaxis()->CenterTitle();//居中
		graph4[i]->GetXaxis()->SetLabelFont(132);//坐标字体
		graph4[i]->GetYaxis()->SetLabelFont(132);//坐标字体
		graph4[i]->GetXaxis()->SetTitleFont(132);//轴名字体
		graph4[i]->GetYaxis()->SetTitleFont(132);//轴名字体
		//graph1[i]->GetYaxis()->SetLabelSize(0.05);//坐标字号
		//graph1[i]->GetYaxis()->SetTitleSize(0.05);//轴名字号
		graph4[i]->GetXaxis()->SetTitleOffset(1.2);//轴名偏移
		graph4[i]->GetYaxis()->SetTitleOffset(1.3);//轴名偏移
		graph4[i]->GetYaxis()->SetRangeUser(range_tau1[i],range_tau2[i]);
		graph4[i]->SetMarkerStyle(21);
		graph4[i]->SetMarkerColor(1);

		TF1 *func_sqrt=new TF1("func_sqrt","[0]*sqrt(x)",0,60000);//多项式拟合，调节合适的拟合range
		func_sqrt->SetParNames("p0");//y=p0*sqrt(x)
//		graph4[i]->Fit("func_sqrt");//pol1 can be used directly without TF1 constructor in CINT

		TFitResultPtr r_tau_sqrt = graph4[i]->Fit("func_sqrt","S");//"S" means the result of the fit is returned in the TFitResultPtr
		//If the option "S" is instead used, TFitResultPtr contains the TFitResult and behaves as a smart pointer to it.
		//The fit parameters, error and chi2 (but not covariance matrix) can be retrieved also from the fitted function.
		TH1D *hint_tau_sqrt = new TH1D("hint_tau_sqrt", "Fitted exponential with conf.band", 9000, 0, 9000);//Create a histogram to hold the confidence intervals
		(TVirtualFitter::GetFitter())->GetConfidenceIntervals(hint_tau_sqrt, 0.683);//By default the intervals are corrected using the chi2/ndf value of the fit if a chi2 fit is performed
		//hint_tau_sqrt will contain the CL result that you can draw on top of your fitted graph.
		//where hint_tau_sqrt will hold the errors and could superimpose it on the same canvas where you plot central values.
		//The method TVirtualFitter::GetFitter())->GetConfidenceIntervals(hint_chi2, 0.95) computes the confidence level of your model (fitting) function after having it fitted to an histogram.
		hint_tau_sqrt->SetStats(kFALSE);
		hint_tau_sqrt->SetFillColor(kRed-9);
		TMatrixDSym cov = r_tau_sqrt->GetCovarianceMatrix();
		r_tau_sqrt->GetConfidenceIntervals(Npoint, 1, 1, channel_point, err_point, 0.683, false);//(Number of x points, 1, 1, x, err, confidence level, norm); norm is a flag to control if the intervals need to be normalized to the chi2/ndf value.
		outfile<<"Eg=	"<<channel_point[3]<<"	tau=	"<<func_sqrt->Eval(channel_point[3])<<"	err_tau=	"<<err_point[3]<<endl;
		graph4[i]->Draw("AP");//"A": Axis are drawn around the graph, "P": The current marker is plotted at each point
		hint_tau_sqrt->Draw("e3 same");
		graph4[i]->Draw("P same");//draw the points again above error band

		p0[i]=func_sqrt->GetParameter(0);
		p0err[i]=func_sqrt->GetParError(0);
		parchi[i]=func_sqrt->GetChisquare();
		parNDF[i]=func_sqrt->GetNDF();
		p_value[i]=func_sqrt->GetProb();//This probability is not the “probability that your fit is good.” If you did many fake experiments (draw many random samples of data points from the assumed distribution (your fit function)), this is the percentage of experiments that would give χ2 values ≥ to the one you got in this experiment.
		double errylow[Npoint],erryhigh[Npoint];
		for(ii=0;ii<Npoint;ii++)
		{
			//x[ii]=channel_point[ii];
			//y[ii]=func_sqrt->Eval(channel_point[ii]);//is equal to y[ii]=p1[i]*sqrt(x[ii])+p0[i];
			//dy[ii]=sqrt(sqrt(x[ii])*sqrt(x[ii])*p1err[i]*p1err[i]+p0err[i]*p0err[i]);
			//dy[ii]=err_point[ii];
			//cout<<"	"<<channel_point[ii]<<"	"<<err_point[ii]<<endl;
			errylow[ii]=func_sqrt->Eval(channel_point[ii])-err_point[ii];
			erryhigh[ii]=func_sqrt->Eval(channel_point[ii])+err_point[ii];
		}
		graph4bandlow[i]= new TGraph(Npoint,channel_point,errylow);//画error bars TGraph(n,x,y,ex,ey);
		graph4bandlow[i]->SetLineColor(2);
		graph4bandlow[i]->SetLineWidth(1);
		graph4bandlow[i]->Draw("C");
		graph4bandhigh[i]= new TGraph(Npoint,channel_point,erryhigh);//画error bars TGraph(n,x,y,ex,ey);
		graph4bandhigh[i]->SetLineColor(2);
		graph4bandhigh[i]->SetLineWidth(1);
		graph4bandhigh[i]->Draw("C");

		TPaveText *textpol4 = new TPaveText(0.10,0.64,0.34,0.89,"brNDC");//left, down, right, up
		textpol4->SetBorderSize(1);
		textpol4->SetFillColor(0);
		textpol4->SetTextAlign(12);//align = 10*HorizontalAlign + VerticalAlign, 12 means水平左对齐、垂直居中对齐
		textpol4->SetTextFont(132);//font = 10 * fontID + precision, 12+2,12 means Symbol; 13+2, 13 means News Times Roman
		sprintf(paraprint,"y=p0*sqrt(x)");
		textpol4->AddText(paraprint);
		sprintf(paraprint,"p0=%.5f+/-%.5f",p0[i],p0err[i]);
		textpol4->AddText(paraprint);
		sprintf(paraprint,"Chisquare=%.2f",parchi[i]);
		textpol4->AddText(paraprint);
		sprintf(paraprint,"NDF=%.0f",parNDF[i]);
		textpol4->AddText(paraprint);
		sprintf(paraprint,"p-value=%e",p_value[i]);
		textpol4->AddText(paraprint);
		textpol4->Draw();
		//sprintf(hcali_name,"%s.png",hcali_name);
		sprintf(hcali_name,"%s%s%s%s%s","D:/X/out/responseforDoppler/",hcali_name,"_",tauflag,"_band.png");
		canvascali[i]->SaveAs(hcali_name);//存图
		outfile2<<"	SeGA_"<<i<<"	y=p0*sqrt(x)	taup0=	"<<p0[i]<<"	+/-	"<<p0err[i]<<"	Chi2=	"<<parchi[i]<<"	ndf=	"<<parNDF[i]<<"	p-value=	"<<p_value[i]<<endl;
*/

	}//for (i=0;i<ID;i++)
	}//ior
	//}//for(irawroot=runstart; irawroot<=runstop; irawroot++)
}//peakcali main