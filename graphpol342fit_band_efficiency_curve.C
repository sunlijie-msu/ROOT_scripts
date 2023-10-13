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
void graphpol342fit_band_efficiency_curve()// graph-pol342 fit efficiency curve (with colored bands) of Proton Detector
{
	time_t tim;
	struct tm *at;
	char now[80];
	int irawroot;
	char rawrootname[300];
	double slope[300],slopeerr[300];
	double intercept[300],intercepterr[300];
	double p5[300],p4[300],p3[300],p2[300],p1[300],p0[300];
	double p5err[300],p4err[300],p3err[300],p2err[300],p1err[300],p0err[300];
	const int ID1=0;// no need to change i=ID1//which detector
	const int ID2=15;// no need to change i<=ID2//which detector

	int binwidth=1;
	//	int binnum=(upperlimit-lowerlimit)/binwidth;
	int minrange=0,maxrange=0,minbin,maxbin;
	double sigma=1,thresh=0.8,sigmab=20;//adjust
	double gaplow=70.,gaphigh=70.;//fitting range随分辨不同调整//
	unsigned long i;
	int jj,ii,ibin,k,ipeak;
	char paraprint[100],b_name[200],histo_name[200],h_name[200],hfit_name[200],hcali_name[200],tauflag[200];
	char command1[100],command2[100];
	//string paraprint[30],h_name[80],histo_name[80],b_name[80],hfit_name[80],hcali_name[80];
	TCanvas *canvaspeak[300];
	TCanvas *canvascali[300];
	TH1F *histo[ID2+1];//TH1F peak search+gauss fit,creat histograms
	const int nummax=10;//modify
	int peaknum=nummax;
	TF1 *total[nummax];//creat function
	TF1 *p[nummax], *g[nummax], *b[nummax];
	TGraph *graph1[300], *graph1bandhigh[300], *graph1bandlow[300];//TGraph
	TGraph *graph2[300], *graph2bandhigh[300], *graph2bandlow[300];//TGraph
	TGraph *graph3[300], *graph3bandhigh[300], *graph3bandlow[300];//TGraph
	TGraph *graph4[300], *graph4bandhigh[300], *graph4bandlow[300];//TGraph
	double Energy[ID2+1][nummax],Energyerr[ID2+1][nummax];
	double peaky[nummax],peakyerr[nummax];
	double sig[ID2+1][nummax],sigerr[ID2+1][nummax];
	double Efficiency[ID2+1][nummax],Efficiencyerr[ID2+1][nummax];
	double range_x1[ID2+1], range_x2[ID2+1], range_y1[ID2+1], range_y2[ID2+1];
	// 	double *Energy;//if you don't know how many peaks will be found, use this
	// 	double *peaky;//if you don't know how many peaks will be found, use this
//	double energylit[6]={451.7, 493.3, 944.9, 1460.820, 1612.4, 2614.511};//25Si
//	double energyliterr[6]={0.5, 0.7, 0.5, 0.005, 0.5, 0.010};//25Si
// 	double energylit[6]={450.7, 1460.820, 1599, 2614.511, 2908, 7801};//23Al
// 	double energyliterr[6]={0.15, 0.005, 2, 0.010, 3, 2};//23Al

	const int Npoint=11;//modify
	double energy_point[Npoint] = {401, 555, 719, 943, 1040, 1268, 1377, 1804, 1917, 2162, 2307};// set location of point for single value
	double err_point[Npoint];  // error on the function at point x0 for single value
	
// 	const int Npoint=230;
// 	double energy_point[Npoint];// set location of point for drawing
// 	double err_point[Npoint];  // error on the function at point x0 for drawing
// 	for(jj=0;jj<Npoint;jj++){	energy_point[jj]=jj*40;}

	double peakmin;
	double chtemp;
	double par[nummax][5],parerr[nummax][5];//随峰数不同改动par[peaknum][3]//adjust
	double parchi[ID2+1],parNDF[ID2+1],p_value[ID2+1];
//	ifstream infile("C:/Si24/Si22peakcali/Clover.dat",ios::in);//The data that need to be fitted
//	int lowerlimit[ID2+1][nummax],upperlimit[ID2+1][nummax];
	for(i=ID1;i<=ID2;i++)//which detector
	{
		peaknum=nummax;
		if(i==8||i==11||i==14||i==15) continue;
		if(i==1)peaknum=nummax-0;//modify
//		sprintf(b_name,"%s%d%s","D:/X/out/Si25/SeGA_",i,"_sigma_freetau.dat");//modify 12 peaks
		sprintf(b_name,"%s%d%s","D:/X/out/Si25/SeGA_",i,"_sigma_tau_onlybg.dat");//only 8 beta-delayed gamma peaks
		ifstream infile(b_name,ios::in);//The data that need to be fitted
// 		for(ii=0;ii<peaknum;ii++)//which peak in one SeGA detector modify =0<12
// 		{
// 			infile>>Energy[i][ii]>>Energyerr[i][ii]>>sig[i][ii]>>sigerr[i][ii]>>Efficiency[i][ii]>>Efficiencyerr[i][ii];
// 			//cout<<"SeGA_"<<i<<'	'<<Energy[i][ii]<<'	'<<Energyerr[i][ii]<<'	'<<sig[i][ii]<<'	'<<sigerr[i][ii]<<'	'<<Efficiency[i][ii]<<'	'<<Efficiencyerr[i][ii]<<endl;
// 		}
 		if(i==0) {range_x1[i]=0; range_x2[i]=2500; range_y1[i]=0; range_y2[i]=0.70;}
// 		if(i==1) {range_x1[i]=0; range_x2[i]=4; range_y1[i]=2000; range_y2[i]=19000;}
// 		if(i==2) {range_x1[i]=0; range_x2[i]=4; range_y1[i]=0; range_y2[i]=4;}
// 		if(i==3) {range_x1[i]=0; range_x2[i]=4; range_y1[i]=0; range_y2[i]=3;}
// 		if(i==4) {range_x1[i]=0; range_x2[i]=4; range_y1[i]=0; range_y2[i]=3;}
// 		if(i==5) {range_x1[i]=0; range_x2[i]=2; range_y1[i]=0; range_y2[i]=1.6;}
// 		if(i==6) {range_x1[i]=0; range_x2[i]=4; range_y1[i]=0; range_y2[i]=3;}
// 		if(i==7) {range_x1[i]=0; range_x2[i]=4; range_y1[i]=0; range_y2[i]=2.5;}
// 		if(i==9) {range_x1[i]=0; range_x2[i]=4; range_y1[i]=-2.5; range_y2[i]=4;}
// 		if(i==10) {range_x1[i]=0; range_x2[i]=4.5; range_y1[i]=0; range_y2[i]=4;}
// 		if(i==12) {range_x1[i]=-6; range_x2[i]=10; range_y1[i]=-3; range_y2[i]=5;}
// 		if(i==13) {range_x1[i]=-4; range_x2[i]=7; range_y1[i]=-2; range_y2[i]=6;}
		Energy[0][0]=401; Energyerr[0][0]=10;
		Energy[0][1]=555; Energyerr[0][1]=10;
		Energy[0][2]=943; Energyerr[0][2]=10;
		Energy[0][3]=1040; Energyerr[0][3]=10;
		Energy[0][4]=1268; Energyerr[0][4]=10;
		Energy[0][5]=1377; Energyerr[0][5]=10;
		Energy[0][6]=1804; Energyerr[0][6]=10;	
		Energy[0][7]=1917; Energyerr[0][7]=10;
		Energy[0][8]=2162; Energyerr[0][8]=10;
		Energy[0][9]=2307; Energyerr[0][9]=10;
		
		Efficiency[0][0]=53.70/100; Efficiencyerr[0][0]=13.10/100.0;
		Efficiency[0][1]=49.53/100; Efficiencyerr[0][1]=5.21/100.0;
		Efficiency[0][2]=29.91/100; Efficiencyerr[0][2]=8.05/100.0;
		Efficiency[0][3]=15.96/100; Efficiencyerr[0][3]=6.53/100.0;
		Efficiency[0][4]=  9.14/100; Efficiencyerr[0][4]=5.05/100.0;
		Efficiency[0][5]=  6.64/100; Efficiencyerr[0][5]=2.89/100.0;
		Efficiency[0][6]=  4.33/100; Efficiencyerr[0][6]=1.61/100.0;
		Efficiency[0][7]=  3.11/100; Efficiencyerr[0][7]=0.85/100.0;
		Efficiency[0][8]=  1.90/100; Efficiencyerr[0][8]=0.55/100.0;
		Efficiency[0][9]=  1.24/100; Efficiencyerr[0][9]=0.53/100.0;

	}
	sprintf(tauflag,"%s","taufree");
	ofstream outfile("D:/X/out/Si25/outfile/peakcali.dat",ios::out);
	ofstream outfile2("D:/X/out/Si25/outfile/peakcalipara.dat",ios::out);

	for(i=0;i<=0;i++)
	{
		peaknum=nummax;
		if(i==8||i==11||i==14||i==15) continue;
		if(i==1)peaknum=nummax-0;//modify




		//************ Efficiency fit by pol3 *****************************************

		if(i==0)sprintf(hcali_name,"%s","pol3");
		canvascali[i]=new TCanvas(hcali_name,hcali_name,900,600);//建立画布
		canvascali[i]->cd();//进入画布
		//canvascali[i]->SetGrid();//显示网格
		//		graph1[i]=new TGraph(peaknum,Energy,Efficiency);//TGraph *gr1=new TGraph(n,x,y);
		graph1[i]= new TGraphErrors(peaknum,Energy[i],Efficiency[i],Energyerr[i],Efficiencyerr[i]);//画error bars TGraph(n,x,y,ex,ey);
		graph1[i]->SetTitle(hcali_name);
		graph1[i]->GetXaxis()->SetTitle("Energy");//轴名
		graph1[i]->GetYaxis()->SetTitle("Efficiency");//轴名
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
		graph1[i]->GetYaxis()->SetRangeUser(range_y1[i],range_y2[i]);
		graph1[i]->SetMarkerStyle(21);
		graph1[i]->SetMarkerColor(1);
		
		TF1 *pol3=new TF1("pol3","pol3",-60000,60000);//多项式拟合，调节合适的拟合range
		pol3->SetParNames("p0","p1","p2","p3");//y=p3*x^2+p2*x^2+p1*x+p0
		pol3->SetParameter(3,-1.52E-11);
		pol3->SetParameter(2,2.08E-07);
		pol3->SetParameter(1,-6.63E-04);
		pol3->SetParameter(0,6.27E-01);

//		graph2[i]->Fit("pol3");//pol1 can be used directly without TF1 constructor in CINT

		TFitResultPtr r_tau_pol3 = graph1[i]->Fit("pol3","MS");//"S" means the result of the fit is returned in the TFitResultPtr
		//“M” Improve fit results, by using the IMPROVE algorithm of TMinuit.
		//If the option "S" is instead used, TFitResultPtr contains the TFitResult and behaves as a smart pointer to it.
		//The fit parameters, error and chi2 (but not covariance matrix) can be retrieved also from the fitted function.
		TH1D *hint_tau_pol3 = new TH1D("hint_tau_pol3", "Fitted exponential with conf.band", 9000, 0, 9000);//Create a histogram to hold the confidence intervals
		(TVirtualFitter::GetFitter())->GetConfidenceIntervals(hint_tau_pol3, 0.683);//By default the intervals are corrected using the chi2/ndf value of the fit if a chi2 fit is performed
		//hint_tau_pol4 will contain the CL result that you can draw on top of your fitted graph.
		//where hint_sig will hold the errors and could superimpose it on the same canvas where you plot central values.
		//The method TVirtualFitter::GetFitter())->GetConfidenceIntervals(hint_sig, 0.95) computes the confidence level of your model (fitting) function after having it fitted to an histogram.
		hint_tau_pol3->SetStats(kFALSE);
		hint_tau_pol3->SetFillColor(kRed-9);
		//TMatrixDSym cov = r_tau_pol4->GetCovarianceMatrix();//useless for now
		TMatrixD cov = r_tau_pol3->GetCorrelationMatrix();
		TMatrixD cor = r_tau_pol3->GetCovarianceMatrix();
		cov.Print();//get cov=rou*sigma1*sigma2...
		cor.Print();//get rou
		r_tau_pol3->GetConfidenceIntervals(Npoint, 1, 1, energy_point, err_point, 0.683, false);//(Number of x points, 1, 1, x, err, confidence level, norm); norm is a flag to control if the intervals need to be normalized to the chi2/ndf value.
		for(ii=0;ii<Npoint;ii++)
		{
			outfile<<"pol3_"<<i<<"	Energy=	"<<energy_point[ii]<<"	Efficiency=	"<<pol3->Eval(energy_point[ii])<<"	err_Y=	"<<err_point[ii]<<endl;
		}
		graph1[i]->Draw("AP");//"A": Axis are drawn around the graph, "P": The current marker is plotted at each point
		hint_tau_pol3->Draw("e3 same");
		graph1[i]->Draw("P same");//draw the points again above error band

		p0[i]=pol3->GetParameter(0);
		p0err[i]=pol3->GetParError(0);
		p1[i]=pol3->GetParameter(1);
		p1err[i]=pol3->GetParError(1);
		p2[i]=pol3->GetParameter(2);
		p2err[i]=pol3->GetParError(2);
		p3[i]=pol3->GetParameter(3);
		p3err[i]=pol3->GetParError(3);
		parchi[i]=pol3->GetChisquare();
		parNDF[i]=pol3->GetNDF();
		p_value[i]=pol3->GetProb();//This probability is not the “probability that your fit is good.” If you did many fake experiments (draw many random samples of data points from the assumed distribution (your fit function)), this is the percentage of experiments that would give χ2 values ≥ to the one you got in this experiment.
		double errylow[Npoint],erryhigh[Npoint];
// 		for(ii=0;ii<Npoint;ii++)
// 		{
// 			//x[ii]=energy_point[ii];
// 			//y[ii]=func_sqrt_const->Eval(energy_point[ii]);//is equal to y[ii]=p1[i]*sqrt(x[ii])+p0[i];
// 			//dy[ii]=sqrt(sqrt(x[ii])*sqrt(x[ii])*p1err[i]*p1err[i]+p0err[i]*p0err[i]);
// 			//dy[ii]=err_point[ii];
// 			//cout<<"	"<<energy_point[ii]<<"	"<<err_point[ii]<<endl;
// 			errylow[ii]=pol4->Eval(energy_point[ii])-err_point[ii];
// 			erryhigh[ii]=pol4->Eval(energy_point[ii])+err_point[ii];
// 		}
// 		graph1bandlow[i]= new TGraph(Npoint,energy_point,errylow);//画error bars TGraph(n,x,y,ex,ey);
// 		graph1bandlow[i]->SetLineColor(4);
// 		graph1bandlow[i]->SetLineWidth(1);
// 		graph1bandlow[i]->Draw("C");
// 		graph1bandhigh[i]= new TGraph(Npoint,energy_point,erryhigh);//画error bars TGraph(n,x,y,ex,ey);
// 		graph1bandhigh[i]->SetLineColor(4);
// 		graph1bandhigh[i]->SetLineWidth(1);
// 		graph1bandhigh[i]->Draw("C");

		TPaveText *textpol3 = new TPaveText(0.54,0.58,0.90,0.89,"brNDC");//left, down, right, up
		textpol3->SetBorderSize(1);
		textpol3->SetFillColor(0);
		textpol3->SetTextAlign(12);//align = 10*HorizontalAlign + VerticalAlign, 12 means水平左对齐、垂直居中对齐
		textpol3->SetTextFont(132);//font = 10 * fontID + precision, 12+2,12 means Symbol; 13+2, 13 means News Times Roman
		sprintf(paraprint,"y=p3*x^3+p2*x^2+p1*x+p0");
		textpol3->AddText(paraprint);
		sprintf(paraprint,"p3=%e+/-%e",p3[i],p3err[i]);
		textpol3->AddText(paraprint);
		sprintf(paraprint,"p2=%e+/-%e",p2[i],p2err[i]);
		textpol3->AddText(paraprint);
		sprintf(paraprint,"p1=%e+/-%e",p1[i],p1err[i]);
		textpol3->AddText(paraprint);
		sprintf(paraprint,"p0=%e+/-%e",p0[i],p0err[i]);
		textpol3->AddText(paraprint);
		sprintf(paraprint,"Chisquare=%.3f",parchi[i]);
		textpol3->AddText(paraprint);
		sprintf(paraprint,"NDF=%.0f",parNDF[i]);
		textpol3->AddText(paraprint);
		sprintf(paraprint,"p-value=%f",p_value[i]);
		textpol3->AddText(paraprint);
		textpol3->Draw();
		//sprintf(hcali_name,"%s.png",hcali_name);
		sprintf(hcali_name,"%s%s%s%s%s","D:/X/out/",hcali_name,"_",tauflag,"_band.png");
		canvascali[i]->SaveAs(hcali_name);//存图
		outfile2<<"	pol3_"<<i<<"	y=p3*x^3+p2*x^2+p1*x+p0	p3=	"<<setprecision(8)<<p3[i]<<"	+/-	"<<p3err[i]<<"	p2=	"<<p2[i]<<"	+/-	"<<p2err[i]<<"	p1=	"<<p1[i]<<"	+/-	"<<p1err[i]<<"	p0=	"<<p0[i]<<"	+/-	"<<p0err[i]<<"	Chi2=	"<<parchi[i]<<"	ndf=	"<<parNDF[i]<<"	p-value=	"<<p_value[i]<<endl;










		//************ Efficiency fit by pol4 *****************************************

		if(i==0)sprintf(hcali_name,"%s","pol4");
		canvascali[i]=new TCanvas(hcali_name,hcali_name,900,600);//建立画布
		canvascali[i]->cd();//进入画布
		//canvascali[i]->SetGrid();//显示网格
		//		graph1[i]=new TGraph(peaknum,Energy,Efficiency);//TGraph *gr1=new TGraph(n,x,y);
		graph2[i]= new TGraphErrors(peaknum,Energy[i],Efficiency[i],Energyerr[i],Efficiencyerr[i]);//画error bars TGraph(n,x,y,ex,ey);
		graph2[i]->SetTitle(hcali_name);
		graph2[i]->GetXaxis()->SetTitle("Energy");//轴名
		graph2[i]->GetYaxis()->SetTitle("Efficiency");//轴名
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
		graph2[i]->GetYaxis()->SetRangeUser(range_y1[i],range_y2[i]);
		graph2[i]->SetMarkerStyle(21);
		graph2[i]->SetMarkerColor(1);
		
		TF1 *pol4=new TF1("pol4","pol4",-60000,60000);//多项式拟合，调节合适的拟合range
		pol4->SetParNames("p0","p1","p2","p3","p4");//y=p2*x^2+p1*x+p0
		pol4->SetParameter(4,-2.00E-13);
		pol4->SetParameter(3,1.09E-09);
		pol4->SetParameter(2,-1.90E-06);
		pol4->SetParameter(1,9.29E-04);
		pol4->SetParameter(0,2.44E-01);
//		graph2[i]->Fit("pol4");//pol1 can be used directly without TF1 constructor in CINT

		TFitResultPtr r_tau_pol4 = graph2[i]->Fit("pol4","MS");//"S" means the result of the fit is returned in the TFitResultPtr
		//“M” Improve fit results, by using the IMPROVE algorithm of TMinuit.
		//If the option "S" is instead used, TFitResultPtr contains the TFitResult and behaves as a smart pointer to it.
		//The fit parameters, error and chi2 (but not covariance matrix) can be retrieved also from the fitted function.
		TH1D *hint_tau_pol4 = new TH1D("hint_tau_pol4", "Fitted exponential with conf.band", 9000, 0, 9000);//Create a histogram to hold the confidence intervals
		(TVirtualFitter::GetFitter())->GetConfidenceIntervals(hint_tau_pol4, 0.683);//By default the intervals are corrected using the chi2/ndf value of the fit if a chi2 fit is performed
		//hint_tau_pol4 will contain the CL result that you can draw on top of your fitted graph.
		//where hint_sig will hold the errors and could superimpose it on the same canvas where you plot central values.
		//The method TVirtualFitter::GetFitter())->GetConfidenceIntervals(hint_sig, 0.95) computes the confidence level of your model (fitting) function after having it fitted to an histogram.
		hint_tau_pol4->SetStats(kFALSE);
		hint_tau_pol4->SetFillColor(kRed-9);
		//TMatrixDSym cov = r_tau_pol4->GetCovarianceMatrix();//useless for now
		TMatrixD cov = r_tau_pol4->GetCorrelationMatrix();
		TMatrixD cor = r_tau_pol4->GetCovarianceMatrix();
		cov.Print();//get cov=rou*sigma1*sigma2...
		cor.Print();//get rou
		r_tau_pol4->GetConfidenceIntervals(Npoint, 1, 1, energy_point, err_point, 0.683, false);//(Number of x points, 1, 1, x, err, confidence level, norm); norm is a flag to control if the intervals need to be normalized to the chi2/ndf value.
		for(ii=0;ii<Npoint;ii++)
		{
			outfile<<"pol4_"<<i<<"	Energy=	"<<energy_point[ii]<<"	Efficiency=	"<<pol4->Eval(energy_point[ii])<<"	err_y=	"<<err_point[ii]<<endl;
		}
		graph2[i]->Draw("AP");//"A": Axis are drawn around the graph, "P": The current marker is plotted at each point
		hint_tau_pol4->Draw("e3 same");
		graph2[i]->Draw("P same");//draw the points again above error band

		p0[i]=pol4->GetParameter(0);
		p0err[i]=pol4->GetParError(0);
		p1[i]=pol4->GetParameter(1);
		p1err[i]=pol4->GetParError(1);
		p2[i]=pol4->GetParameter(2);
		p2err[i]=pol4->GetParError(2);
		p3[i]=pol4->GetParameter(3);
		p3err[i]=pol4->GetParError(3);
		p4[i]=pol4->GetParameter(4);
		p4err[i]=pol4->GetParError(4);

		parchi[i]=pol4->GetChisquare();
		parNDF[i]=pol4->GetNDF();
		p_value[i]=pol4->GetProb();//This probability is not the “probability that your fit is good.” If you did many fake experiments (draw many random samples of data points from the assumed distribution (your fit function)), this is the percentage of experiments that would give χ2 values ≥ to the one you got in this experiment.
		double errylow[Npoint],erryhigh[Npoint];
// 		for(ii=0;ii<Npoint;ii++)
// 		{
// 			//x[ii]=energy_point[ii];
// 			//y[ii]=func_sqrt_const->Eval(energy_point[ii]);//is equal to y[ii]=p1[i]*sqrt(x[ii])+p0[i];
// 			//dy[ii]=sqrt(sqrt(x[ii])*sqrt(x[ii])*p1err[i]*p1err[i]+p0err[i]*p0err[i]);
// 			//dy[ii]=err_point[ii];
// 			//cout<<"	"<<energy_point[ii]<<"	"<<err_point[ii]<<endl;
// 			errylow[ii]=pol4->Eval(energy_point[ii])-err_point[ii];
// 			erryhigh[ii]=pol4->Eval(energy_point[ii])+err_point[ii];
// 		}
// 		graph3bandlow[i]= new TGraph(Npoint,energy_point,errylow);//画error bars TGraph(n,x,y,ex,ey);
// 		graph3bandlow[i]->SetLineColor(4);
// 		graph3bandlow[i]->SetLineWidth(1);
// 		graph3bandlow[i]->Draw("C");
// 		graph3bandhigh[i]= new TGraph(Npoint,energy_point,erryhigh);//画error bars TGraph(n,x,y,ex,ey);
// 		graph3bandhigh[i]->SetLineColor(4);
// 		graph3bandhigh[i]->SetLineWidth(1);
// 		graph3bandhigh[i]->Draw("C");

		TPaveText *textpol4 = new TPaveText(0.54,0.58,0.90,0.89,"brNDC");//left, down, right, up
		textpol4->SetBorderSize(1);
		textpol4->SetFillColor(0);
		textpol4->SetTextAlign(12);//align = 10*HorizontalAlign + VerticalAlign, 12 means水平左对齐、垂直居中对齐
		textpol4->SetTextFont(132);//font = 10 * fontID + precision, 12+2,12 means Symbol; 13+2, 13 means News Times Roman
		sprintf(paraprint,"y=p4*x^4+p3*x^3+p2*x^2+p1*x+p0");
		textpol4->AddText(paraprint);
		sprintf(paraprint,"p4=%e+/-%e",p4[i],p4err[i]);
		textpol4->AddText(paraprint);
		sprintf(paraprint,"p3=%e+/-%e",p3[i],p3err[i]);
		textpol4->AddText(paraprint);
		sprintf(paraprint,"p2=%e+/-%e",p2[i],p2err[i]);
		textpol4->AddText(paraprint);
		sprintf(paraprint,"p1=%e+/-%e",p1[i],p1err[i]);
		textpol4->AddText(paraprint);
		sprintf(paraprint,"p0=%e+/-%e",p0[i],p0err[i]);
		textpol4->AddText(paraprint);
		sprintf(paraprint,"Chisquare=%.3f",parchi[i]);
		textpol4->AddText(paraprint);
		sprintf(paraprint,"NDF=%.0f",parNDF[i]);
		textpol4->AddText(paraprint);
		sprintf(paraprint,"p-value=%f",p_value[i]);
		textpol4->AddText(paraprint);
		textpol4->Draw();
		//sprintf(hcali_name,"%s.png",hcali_name);
		sprintf(hcali_name,"%s%s%s%s%s","D:/X/out/",hcali_name,"_",tauflag,"_band.png");
		canvascali[i]->SaveAs(hcali_name);//存图
		outfile2<<"	pol4_"<<i<<"	y=p2*x^2+p1*x+p0	p2=	"<<setprecision(8)<<p2[i]<<"	+/-	"<<p2err[i]<<"	p1=	"<<p1[i]<<"	+/-	"<<p1err[i]<<"	p0=	"<<p0[i]<<"	+/-	"<<p0err[i]<<"	Chi2=	"<<parchi[i]<<"	ndf=	"<<parNDF[i]<<"	p-value=	"<<p_value[i]<<endl;
















				//************ Efficiency fit by pol2 *****************************************

		if(i==0)sprintf(hcali_name,"%s","pol2");
		canvascali[i]=new TCanvas(hcali_name,hcali_name,900,600);//建立画布
		canvascali[i]->cd();//进入画布
		//canvascali[i]->SetGrid();//显示网格
		//		graph1[i]=new TGraph(peaknum,Energy,Efficiency);//TGraph *gr1=new TGraph(n,x,y);
		graph3[i]= new TGraphErrors(peaknum,Energy[i],Efficiency[i],Energyerr[i],Efficiencyerr[i]);//画error bars TGraph(n,x,y,ex,ey);
		graph3[i]->SetTitle(hcali_name);
		graph3[i]->GetXaxis()->SetTitle("Energy");//轴名
		graph3[i]->GetYaxis()->SetTitle("Efficiency");//轴名
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
		graph3[i]->GetYaxis()->SetRangeUser(range_y1[i],range_y2[i]);
		graph3[i]->SetMarkerStyle(21);
		graph3[i]->SetMarkerColor(1);
		
		TF1 *pol2=new TF1("pol2","pol2",-60000,60000);//多项式拟合，调节合适的拟合range
		pol2->SetParNames("p0","p1","p2");//y=p2*x^2+p1*x+p0
		pol2->SetParameter(2,1.47E-07);
		pol2->SetParameter(1,-5.92E-04);
		pol2->SetParameter(0,6.04E-01);

//		graph2[i]->Fit("pol4");//pol1 can be used directly without TF1 constructor in CINT

		TFitResultPtr r_tau_pol2 = graph3[i]->Fit("pol2","MS");//"S" means the result of the fit is returned in the TFitResultPtr
		//“M” Improve fit results, by using the IMPROVE algorithm of TMinuit.
		//If the option "S" is instead used, TFitResultPtr contains the TFitResult and behaves as a smart pointer to it.
		//The fit parameters, error and chi2 (but not covariance matrix) can be retrieved also from the fitted function.
		TH1D *hint_tau_pol2 = new TH1D("hint_tau_pol2", "Fitted exponential with conf.band", 9000, 0, 9000);//Create a histogram to hold the confidence intervals
		(TVirtualFitter::GetFitter())->GetConfidenceIntervals(hint_tau_pol2, 0.683);//By default the intervals are corrected using the chi2/ndf value of the fit if a chi2 fit is performed
		//hint_tau_pol4 will contain the CL result that you can draw on top of your fitted graph.
		//where hint_sig will hold the errors and could superimpose it on the same canvas where you plot central values.
		//The method TVirtualFitter::GetFitter())->GetConfidenceIntervals(hint_sig, 0.95) computes the confidence level of your model (fitting) function after having it fitted to an histogram.
		hint_tau_pol2->SetStats(kFALSE);
		hint_tau_pol2->SetFillColor(kRed-9);
		//TMatrixDSym cov = r_tau_pol4->GetCovarianceMatrix();//useless for now
		TMatrixD cov = r_tau_pol2->GetCorrelationMatrix();
		TMatrixD cor = r_tau_pol2->GetCovarianceMatrix();
		cov.Print();//get cov=rou*sigma1*sigma2...
		cor.Print();//get rou
		r_tau_pol2->GetConfidenceIntervals(Npoint, 1, 1, energy_point, err_point, 0.683, false);//(Number of x points, 1, 1, x, err, confidence level, norm); norm is a flag to control if the intervals need to be normalized to the chi2/ndf value.
		for(ii=0;ii<Npoint;ii++)
		{
			outfile<<"pol2_"<<i<<"	Energy=	"<<energy_point[ii]<<"	Efficiency=	"<<pol2->Eval(energy_point[ii])<<"	err_Y=	"<<err_point[ii]<<endl;
		}
		graph3[i]->Draw("AP");//"A": Axis are drawn around the graph, "P": The current marker is plotted at each point
		hint_tau_pol2->Draw("e3 same");
		graph3[i]->Draw("P same");//draw the points again above error band

		p0[i]=pol2->GetParameter(0);
		p0err[i]=pol2->GetParError(0);
		p1[i]=pol2->GetParameter(1);
		p1err[i]=pol2->GetParError(1);
		p2[i]=pol2->GetParameter(2);
		p2err[i]=pol2->GetParError(2);
		parchi[i]=pol2->GetChisquare();
		parNDF[i]=pol2->GetNDF();
		p_value[i]=pol2->GetProb();//This probability is not the “probability that your fit is good.” If you did many fake experiments (draw many random samples of data points from the assumed distribution (your fit function)), this is the percentage of experiments that would give χ2 values ≥ to the one you got in this experiment.
		double errylow[Npoint],erryhigh[Npoint];
// 		for(ii=0;ii<Npoint;ii++)
// 		{
// 			//x[ii]=energy_point[ii];
// 			//y[ii]=func_sqrt_const->Eval(energy_point[ii]);//is equal to y[ii]=p1[i]*sqrt(x[ii])+p0[i];
// 			//dy[ii]=sqrt(sqrt(x[ii])*sqrt(x[ii])*p1err[i]*p1err[i]+p0err[i]*p0err[i]);
// 			//dy[ii]=err_point[ii];
// 			//cout<<"	"<<energy_point[ii]<<"	"<<err_point[ii]<<endl;
// 			errylow[ii]=pol4->Eval(energy_point[ii])-err_point[ii];
// 			erryhigh[ii]=pol4->Eval(energy_point[ii])+err_point[ii];
// 		}
// 		graph3bandlow[i]= new TGraph(Npoint,energy_point,errylow);//画error bars TGraph(n,x,y,ex,ey);
// 		graph3bandlow[i]->SetLineColor(4);
// 		graph3bandlow[i]->SetLineWidth(1);
// 		graph3bandlow[i]->Draw("C");
// 		graph3bandhigh[i]= new TGraph(Npoint,energy_point,erryhigh);//画error bars TGraph(n,x,y,ex,ey);
// 		graph3bandhigh[i]->SetLineColor(4);
// 		graph3bandhigh[i]->SetLineWidth(1);
// 		graph3bandhigh[i]->Draw("C");

		TPaveText *textpol2 = new TPaveText(0.54,0.58,0.90,0.89,"brNDC");//left, down, right, up
		textpol2->SetBorderSize(1);
		textpol2->SetFillColor(0);
		textpol2->SetTextAlign(12);//align = 10*HorizontalAlign + VerticalAlign, 12 means水平左对齐、垂直居中对齐
		textpol2->SetTextFont(132);//font = 10 * fontID + precision, 12+2,12 means Symbol; 13+2, 13 means News Times Roman
		sprintf(paraprint,"y=p2*x^2+p1*x+p0");
		textpol2->AddText(paraprint);
		sprintf(paraprint,"p2=%e+/-%e",p2[i],p2err[i]);
		textpol2->AddText(paraprint);
		sprintf(paraprint,"p1=%e+/-%e",p1[i],p1err[i]);
		textpol2->AddText(paraprint);
		sprintf(paraprint,"p0=%e+/-%e",p0[i],p0err[i]);
		textpol2->AddText(paraprint);
		sprintf(paraprint,"Chisquare=%.3f",parchi[i]);
		textpol2->AddText(paraprint);
		sprintf(paraprint,"NDF=%.0f",parNDF[i]);
		textpol2->AddText(paraprint);
		sprintf(paraprint,"p-value=%f",p_value[i]);
		textpol2->AddText(paraprint);
		textpol2->Draw();
		//sprintf(hcali_name,"%s.png",hcali_name);
		sprintf(hcali_name,"%s%s%s%s%s","D:/X/out/",hcali_name,"_",tauflag,"_band.png");
		canvascali[i]->SaveAs(hcali_name);//存图
		outfile2<<"	pol2_"<<i<<"	y=p3*x^3+p2*x^2+p1*x+p0	p3=	"<<setprecision(8)<<p3[i]<<"	+/-	"<<p3err[i]<<"	p2=	"<<p2[i]<<"	+/-	"<<p2err[i]<<"	p1=	"<<p1[i]<<"	+/-	"<<p1err[i]<<"	p0=	"<<p0[i]<<"	+/-	"<<p0err[i]<<"	Chi2=	"<<parchi[i]<<"	ndf=	"<<parNDF[i]<<"	p-value=	"<<p_value[i]<<endl;










	}//for (i=0;i<ID;i++)
	//}//for(irawroot=runstart; irawroot<=runstop; irawroot++)
}//peakcali main