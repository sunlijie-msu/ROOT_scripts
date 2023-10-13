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
#include "string.h"
#include "stdlib.h"
#include "stdio.h"
#include "TLegend.h"
using namespace std;
void graphpol1fit_band_energy_cali_3rd()// read in 12 .dat files obtained from Si25gausnerfcpol1_peakfit_cali.C and graph-pol1 fit maximum (with colored bands) of SeGA 16 channels for 12 25Si,23Al peaks to extract formal ax+b parameters from fitting.
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
	const int ID2=16;// no need to change i<=ID2//which detector

	int binwidth=1;
	//	int binnum=(upperlimit-lowerlimit)/binwidth;
	int minrange=0,maxrange=0,minbin,maxbin;
	double sigma=1,thresh=0.8,sigmab=20;//adjust
	double gaplow=70.,gaphigh=70.;//fitting range随分辨不同调整
	unsigned long i;
	int jj,ii,ibin,k,ipeak;
	char paraprint[100],b_name[200],histo_name[200],h_name[200],hfit_name[200],hcali_name[200],tauflag[200];
	char command1[100],command2[100];
	//string paraprint[30],h_name[80],histo_name[80],b_name[80],hfit_name[80],hcali_name[80];
	TCanvas *canvaspeak[300];
	TCanvas *canvascali[300];
	TH1F *histo[ID2+1];//TH1F peak search+gauss fit,creat histograms
	const int nummax=8;//search peaknumber can be many, fitnumber can be limited to 3//modify
	int peaknum;
	TF1 *total[nummax];//creat function
	TF1 *p[nummax], *g[nummax], *b[nummax];
	TGraph *graph1[300], *graph2[300], *graph1bandhigh[300], *graph1bandlow[300];//TGraph刻度
	Double_t energy_lit[ID2+1][nummax],energy_literr[ID2+1][nummax];
	Double_t peaky[nummax],peakyerr[nummax];
	Double_t maximum[ID2+1][nummax],meanerr[ID2+1][nummax];
	Double_t tau[ID2+1][nummax],tauerr[ID2+1][nummax];
	// 	double *mean;//if you don't know how many peaks will be found, use this
	// 	double *peaky;//if you don't know how many peaks will be found, use this
	//	double energylit[6]={451.7, 493.3, 944.9, 1460.820, 1612.4, 2614.511};//25Si
	//	double energyliterr[6]={0.5, 0.7, 0.5, 0.005, 0.5, 0.010};//25Si
// 	double energylit[6]={450.7, 1460.820, 1599, 2614.511, 2908, 7801};//23Al
// 	double energyliterr[6]={0.15, 0.005, 2, 0.010, 3, 2};//23Al
	double peakmin;
	double chtemp;
	double par[nummax][5],parerr[nummax][5];//随峰数不同改动par[peaknum][3]//adjust
	double parchi[ID2+1],parNDF[ID2+1],p_value[ID2+1];
	const int Npoint=19;
	double channel_point[Npoint];// set location of point for single value
	double err_point[Npoint];  // error on the function at point x0 for single value
//	double energy_point[Npoint] = {451.7820, 493.1847, 944.6500, 1460.8200, 1612.6869, 2614.5110, 450.6924, 1599.5500, 2904.7794, 7801.4706, 1369, 2754, 2870, 4238};// set location of point for single value
	double energy_point[Npoint] = {451.7, 493.3, 944.9, 1612.4, 1368.626, 2754.007, 2869.5, 4237.96, 844.6, 1337.8, 1789.4, 883.8, 6955, 2221.5, 2673.1, 6288, 7900, 7389, 6878};// all NuDat values, set location of point for single value modify whose uncertainty do you want?
	
//	double energy_point[Npoint] = {1369, 2754, 2870, 4238};// set location of point for single value

//  	const int Npoint=2000;
//  	double channel_point[Npoint];// set location of point for drawing
//  	double err_point[Npoint];  // error on the function at point x0 for drawing
//  	for(jj=0;jj<Npoint;jj++){	channel_point[jj]=jj*35;}


//	ifstream infile("C:/Si24/Si22peakcali/Clover.dat",ios::in);//The data that need to be fitted
//	int lowerlimit[ID2+1][nummax],upperlimit[ID2+1][nummax];
	for(i=ID2;i<=ID2;i++)//which detector
	{
		peaknum=nummax;
		if(i==8) continue;
		//if(i==5)peaknum=nummax-1;
		//sprintf(b_name,"%s%d%s","D:/X/out/Si25/SeGA_",i,"_cali_para_AlSi8peaks.dat");//modify
		sprintf(b_name,"%s","D:/X/out/Si25/SeGA_all_cali_para_Si8peaks.dat");//modify 3rd step
		ifstream infile(b_name,ios::in);//The data that need to be fitted
		for(ii=0;ii<peaknum;ii++)//which peak in one SeGA detector =0<12
		{
			infile>>energy_lit[i][ii]>>energy_literr[i][ii]>>maximum[i][ii]>>meanerr[i][ii];
			cout<<"SeGA_"<<i<<'	'<<energy_lit[i][ii]<<'	'<<energy_literr[i][ii]<<'	'<<maximum[i][ii]<<'	'<<meanerr[i][ii]<<endl;
		}
	}
	sprintf(tauflag,"%s","taufree");
	ofstream outfile("D:/X/out/Si25/outfile/peakcalipara.dat",ios::out);





	for(i=16;i<=16;i++)//which SeGA detector modify =0<=15
	{
		peaknum=nummax;
		if(i==8) continue;
		//if(i==5)peaknum=nummax-1;
		sprintf(hcali_name,"%s%d","Si_Maximum_SeGA_",i);
		canvascali[i]=new TCanvas(hcali_name,hcali_name,900,600);//建立画布
		canvascali[i]->cd();//进入画布
		//canvascali[i]->SetGrid();//显示网格
		//		graph1[i]=new TGraph(peaknum,mean,sig);//TGraph *gr1=new TGraph(n,x,y);
		graph1[i]= new TGraphErrors(peaknum,maximum[i],energy_lit[i],meanerr[i],energy_literr[i]);//画error bars TGraph(n,x,y,ex,ey);
		graph1[i]->SetTitle(hcali_name);
		graph1[i]->GetXaxis()->SetTitle("Channel");//轴名
		graph1[i]->GetYaxis()->SetTitle("Energy");//轴名
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
		TF1 *pol1=new TF1("pol1","pol1",0,70000);//多项式拟合，调节合适的拟合range
		pol1->SetParNames("p0","p1");//y=p1*x+p0
		//graph1[i]->Fit("pol1","V");//pol1 can be used directly without TF1 constructor in CINT
		TFitResultPtr r_cali = graph1[i]->Fit("pol1","MES");//"S" means the result of the fit is returned in the TFitResultPtr
		//If the option "S" is instead used, TFitResultPtr contains the TFitResult and behaves as a smart pointer to it.
		//The fit parameters, error and chi2 (but not covariance matrix) can be retrieved also from the fitted function.
		TH1D *hint_cali = new TH1D("hint_cali", "Fitted exponential with conf.band",70000,0,70000);//Create a histogram to hold the confidence intervals
		(TVirtualFitter::GetFitter())->GetConfidenceIntervals(hint_cali, 0.683);//By default the intervals are corrected using the chi2/ndf value of the fit if a chi2 fit is performed
		//hint_cali will contain the CL result that you can draw on top of your fitted graph.
		//where hint_cali will hold the errors and could superimpose it on the same canvas where you plot central values.
		//The method TVirtualFitter::GetFitter())->GetConfidenceIntervals(hint_cali, 0.95) computes the confidence level of your model (fitting) function after having it fitted to an histogram.
		hint_cali->SetStats(kFALSE);
		hint_cali->SetFillColor(kRed-9);
		//TMatrixDSym cov = r_cali->GetCovarianceMatrix();//useless for now
		TMatrixD cov = r_cali->GetCorrelationMatrix();//parameter correlation coefficients
		TMatrixD cor = r_cali->GetCovarianceMatrix();//error matrix
		cov.Print();
		cor.Print();
		
		p0[i]=pol1->GetParameter(0);//intercept
		p1[i]=pol1->GetParameter(1);//slope
		p0err[i]=pol1->GetParError(0);
		p1err[i]=pol1->GetParError(1);
		parchi[i]=pol1->GetChisquare();
		parNDF[i]=pol1->GetNDF();
		p_value[i]=pol1->GetProb();//This probability is not the “probability that your fit is good.” If you did many fake experiments (draw many random samples of data points from the assumed distribution (your fit function)), this is the percentage of experiments that would give χ2 values ≥ to the one you got in this experiment.
// 		TVirtualFitter * fitter = TVirtualFitter::GetFitter();
// 		p0[i]=fitter->GetParameter(0);
// 		p1[i]=fitter->GetParameter(1);
// 		double E_mu_mu = fitter->GetCovarianceMatrixElement(0,0);
// 		cout<<sqrt(E_mu_mu)<<endl;
// 		double E_mu_mu = fitter->GetCovarianceMatrixElement(1,1);
// 		cout<<sqrt(E_mu_mu)<<endl;
// 		double E_mu_mu = fitter->GetCovarianceMatrixElement(0,1);
// 		cout<<sqrt(abs(E_mu_mu))<<endl;

		for(ii=0;ii<Npoint;ii++)
		{
			channel_point[ii]=(energy_point[ii]-p0[i])/p1[i];
		}

		r_cali->GetConfidenceIntervals(Npoint, 1, 1, channel_point, err_point, 0.683, false);//(Number of x points, 1, 1, x, err, confidence level, false); norm is a flag to control if the intervals need to be normalized to the chi2/ndf value, true=inflate, false=orginal. modify
		outfile<<"SeGA_"<<i<<"	p1=	"<<setprecision(8)<<p1[i]<<"	+/-	"<<p1err[i]<<"	p0=	"<<p0[i]<<"	+/-	"<<p0err[i]<<"	Chi2=	"<<parchi[i]<<"	ndf=	"<<parNDF[i]<<"	p-value=	"<<p_value[i]<<endl;
		for(ii=0;ii<Npoint;ii++)
		{
			outfile<<"SeGA_"<<i<<"	Ch=	"<<channel_point[ii]<<"	Eg=	"<<pol1->Eval(channel_point[ii])<<"	err_Eg=	"<<err_point[ii]<<endl;
		}
		outfile<<endl;
		graph1[i]->Draw("AP");//"A": Axis are drawn around the graph, "P": The current marker is plotted at each point
		hint_cali->Draw("e3 same");
		graph1[i]->Draw("P same");//draw the points again above error band

		double errylow[Npoint],erryhigh[Npoint];
		for(ii=0;ii<Npoint;ii++)
		{
			//x[ii]=energy_point[ii];
			//y[ii]=pol1->Eval(energy_point[ii]);//is equal to y[ii]=p1[i]*x[ii]+p0[i];
			//dy[ii]=sqrt(x[ii]*x[ii]*p1err[i]*p1err[i]+p0err[i]*p0err[i]+p1[i]*p1[i]*dx[ii]*dx[ii]);
			//dy[ii]=err_point[ii];
			//cout<<"	"<<energy_point[ii]<<"	"<<err_point[ii]<<endl;
			errylow[ii]=pol1->Eval(channel_point[ii])-err_point[ii];
			erryhigh[ii]=pol1->Eval(channel_point[ii])+err_point[ii];
		}
// 		graph1bandlow[i]= new TGraph(Npoint,channel_point,errylow);//画error bars TGraph(n,x,y,ex,ey);
// 		graph1bandlow[i]->SetLineColor(2);
// 		graph1bandlow[i]->SetLineWidth(1);
// 		graph1bandlow[i]->Draw("C");//"C" A smooth Curve is drawn
// 		graph1bandhigh[i]= new TGraph(Npoint,channel_point,erryhigh);//画error bars TGraph(n,x,y,ex,ey);
// 		graph1bandhigh[i]->SetLineColor(2);
// 		graph1bandhigh[i]->SetLineWidth(1);
// 		graph1bandhigh[i]->Draw("C");//"C" A smooth Curve is drawn

		TPaveText *textpol1 = new TPaveText(0.13,0.64,0.55,0.89,"brNDC");//left, down, right, up
		textpol1->SetBorderSize(1);
		textpol1->SetFillColor(0);
		textpol1->SetTextAlign(12);//align = 10*HorizontalAlign + VerticalAlign, 12 means水平左对齐、垂直居中对齐
		textpol1->SetTextFont(132);//font = 10 * fontID + precision, 12+2,12 means Symbol; 13+2, 13 means News Times Roman
		sprintf(paraprint,"slope(p1)=%.9f+/-%.9f",p1[i],p1err[i]);
		textpol1->AddText(paraprint);
		sprintf(paraprint,"intercept(p0)=%.5f+/-%.5f",p0[i],p0err[i]);
		textpol1->AddText(paraprint);
		sprintf(paraprint,"Chisquare=%.2f",parchi[i]);
		textpol1->AddText(paraprint);
		sprintf(paraprint,"NDF=%.2f",parNDF[i]);
		textpol1->AddText(paraprint);
		sprintf(paraprint,"p-value=%e",p_value[i]);
		textpol1->AddText(paraprint);
		textpol1->Draw();
		//sprintf(hcali_name,"%s.png",hcali_name);
		sprintf(hcali_name,"%s%s%s%s%s","D:/X/out/cali/",hcali_name,"_",tauflag,"_formal_cali.png");
		canvascali[i]->SaveAs(hcali_name);//存图
		
		//goes to excel match
	}//for (i=0;i<ID;i++)
	//}//for(irawroot=runstart; irawroot<=runstop; irawroot++)
}//peakcali main