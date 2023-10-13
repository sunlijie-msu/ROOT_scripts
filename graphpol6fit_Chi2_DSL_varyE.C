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
void graphpol6fit_Chi2_DSL_varyE()// graph-pol6 fit Chi2
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
	int jj,ii,ibin,k,ipeak;
	char paraprint[100],b_name[200],histo_name[200],h_name[200],hfit_name[200],hcali_name[200],tauflag[200];
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
	double Energy[ID2+1][nummax],meanerr[ID2+1][nummax];
	double peaky[nummax],peakyerr[nummax];
	double Chi[ID2+1][nummax],Chierr[ID2+1][nummax];
	double tau[ID2+1][nummax],tauerr[ID2+1][nummax];
	double range_sig1[ID2+1], range_sig2[ID2+1], range_tau1[ID2+1], range_tau2[ID2+1];
	// 	double *Energy;//if you don't know how many peaks will be found, use this
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
	int Eg=3076;//DSLmodify  1248, 2234, 3076, 4971, 5156
	int Ea=46;//DSLmodify         49,     47,      46,     42,     42
	double Tau=0;//DSLmodify 1200, 270, 
	char flag[20];
	sprintf(flag,"base");//DSLmodify
// 	sprintf(flag,"SP0.9");
// 	sprintf(flag,"SP1.1");

	int Flag_sig = 0;//default =0
	int Flag_tau = 0;//default =0
	int Flag_st = 0;//default =0
	int Flag_or = 0;
//	sprintf(histo_name,"%s","Robertson");
//	sprintf(histo_name,"%s","Thomas");
	string line;
//	ifstream infile("C:/Si24/Si22peakcali/Clover.dat",ios::in);//The data that need to be fitted
//	int lowerlimit[ID2+1][nummax],upperlimit[ID2+1][nummax];
	for(int ior=0;ior<1;ior++)//<1 only show baseline modify, <7 show baseline+-σ, τ, sp
	{
		if(ior==0){	Flag_sig=0; Flag_tau=0; Flag_st=0; }//baseline
		if(ior==1){	Flag_sig=1; Flag_tau=0; Flag_st=0; }
		if(ior==2){	Flag_sig=-1; Flag_tau=0; Flag_st=0; }
		if(ior==3){	Flag_sig=0; Flag_tau=1; Flag_st=0; }
		if(ior==4){	Flag_sig=0; Flag_tau=-1; Flag_st=0; }
		if(ior==5){	Flag_sig=0; Flag_tau=0; Flag_st=1; }
		if(ior==6){	Flag_sig=0; Flag_tau=0; Flag_st=-1; }
//		{	Flag_or=ior+1; }//baseline
//  		if(ior==0){	Flag_or=15; }//baseline
//  		if(ior==1){	Flag_or=18; }
// 		if(ior==2){	Flag_or=21; }
// 		if(ior==3){	Flag_or=24; }
// 		if(ior==4){	Flag_or=30; }
// 		if(ior==5){	Flag_or=40; }
// 		if(ior==6){	Flag_or=50; }
	for(i=ID1;i<=ID2;i++)//which detector //don't change
	{
// 		if(i==8||i==11||i==14|i==15) continue;
// 		if(i==5)peaknum=nummax-1;//
// 		if(i==0)peaknum=nummax-4;

//		if(Eg==1369) sprintf(b_name,"%s%s%s%s%s%d%s%d%s%d%s%d%s%d%s","D:/X/out/",histo_name,"/Chi2_gated_",histo_name,"_",Eg,"_st",Flag_st,"_sig",Flag_sig,"_tau",Flag_tau,"_order",Flag_or,".dat");//temporary

//	  	if(Eg==1369) sprintf(b_name,"%s%s%s%s%s%d%s%d%s%d%s%d%s","D:/X/out/",histo_name,"/Chi2_gated_",histo_name,"_",Eg,"_st",Flag_st,"_sig",Flag_sig,"_tau",Flag_tau,"_August.dat");//formal 4bgcali
//	  	if(Eg!=1369) sprintf(b_name,"%s%s%s%s%s%d%s%d%s%d%s%d%s","D:/X/out/",histo_name,"/Chi2_ungated_",histo_name,"_",Eg,"_st",Flag_st,"_sig",Flag_sig,"_tau",Flag_tau,".dat");//formal 4bgcali
		
//		sprintf(b_name,"%s%d%s%d%s","D:/X/out/DSL/Chi2_Ea",Ea,"Eg",Eg,"projection_vEg.dat");//DSLmodify
		sprintf(b_name,"%s%d%s%d%s%.0f%s%s%s","D:/X/out/DSL/Chi2_Ea",Ea,"Eg",Eg,"Tau",Tau,"slice_vEg_",flag,".dat");
		cout<<b_name<<endl;
		ifstream infile(b_name,ios::in);//The data that need to be fitted
		ii=0;
		while ( getline(infile, line) )
		{
			stringstream(line)>>Energy[i][ii]>>Chi[i][ii];
			cout<<Energy[i][ii]<<'	'<<Chi[i][ii]<<endl;
			ii++;
		}
		peaknum=ii;
	}
	// 	const int Npoint=4;
	// 	double energy_point[Npoint] = {1369,2754,2870,4238};// set location of point for single value
	// 	double err_point[Npoint];  // error on the function at point x0 for single value
	//sprintf(tauflag,"%s","taufree");
	//ofstream outfile("D:/X/out/Si25/outfile/peakcali.dat",ios::out);
	ofstream outfile2("D:/X/out/Si25/outfile/peakcalipara.dat",ios::app);




	for(i=0;i<=0;i++)//which SeGA detector =0<=15
	{
// 		if(i==8||i==11||i==14|i==15) continue;
// 		if(i==5)peaknum=nummax-1;//
// 		if(i==0)peaknum=nummax-4;


		//************ sigma fit by linear func *****************************************

		sprintf(hcali_name,"%s%d%s%d%s%.0f%s%s","Chi2_Ea",Ea,"Eg",Eg,"Tau",Tau,"slice_vEg_",flag);
		canvascali[i]=new TCanvas(hcali_name,hcali_name,900,600);//建立画布
		canvascali[i]->cd();//进入画布
		//canvascali[i]->SetGrid();//显示网格
		graph1[i]=new TGraph(peaknum,Energy[i],Chi[i]);//TGraph *gr1=new TGraph(n,x,y);
		//graph1[i]= new TGraphErrors(peaknum,Energy[i],Chi[i],meanerr[i],Chierr[i]);//画error bars TGraph(n,x,y,ex,ey);
		graph1[i]->SetTitle(hcali_name);
		//graph1[i]->GetXaxis()->SetTitle("Lifetime (fs)");//轴名
		graph1[i]->GetXaxis()->SetTitle("E_{#gamma} (keV)");//轴名
		graph1[i]->GetYaxis()->SetTitle("Chi2");//轴名
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
		if(Eg==1248){ minrange=1234; maxrange=1252; }//DSLmodify
		if(Eg==2234){ minrange=2222; maxrange=2244; }
		if(Eg==3076){ minrange=3067; maxrange=3085; }
		if(Eg==4971){ minrange=4961; maxrange=4979; }
		if(Eg==5156){ minrange=5144; maxrange=5172; }
		if(Eg==6129){ minrange=6125; maxrange=6132; }
		if(Eg==2814){ minrange=2812; maxrange=2816; }
		graph1[i]->GetXaxis()->SetRangeUser(minrange,maxrange);

		//TF1 *pol6=new TF1("pol6","[0]+[1]*x+[2]*x^2+[3]*x^3+[4]*x^4+[5]*x^5+[6]*x^6",minrange,maxrange);//多项式拟合，调节合适的拟合range
		TF1 *pol6=new TF1("pol6","pol6",minrange,maxrange);//多项式拟合，调节合适的拟合range
		pol6->SetNpx((maxrange-minrange)*10);
		pol6->SetParNames("p0","p1","p2","p3","p4","p5","p6");//y=p3*x^3+p2*x^2+p1*x+p0
		
		//Give initial guess. ROOT is not smart enough to do good pol>3 fit without initial guess.
		//DSLmodify
// 		pol6->SetParameter(6,0);
// 		pol6->SetParameter(5,-1.94487760396441E-01);
// 		pol6->SetParameter(4,1.20676089831191E+03);
// 		pol6->SetParameter(3,-2.99500517881828E+06);
// 		pol6->SetParameter(2,3.71647522991241E+09);
// 		pol6->SetParameter(1,-2.30580264626732E+12);
// 		pol6->SetParameter(0,5.72216193814757E+14);
		//1248 projection y = -1.94487760396441E-01x5 + 1.20676089831191E+03x4 - 2.99500517881828E+06x3 + 3.71647522991241E+09x2 - 2.30580264626732E+12x + 5.72216193814757E+14

// 		pol6->SetParameter(6,0);
// 		pol6->SetParameter(5,0);
// 		pol6->SetParameter(4, -7.16833333334762E-01);
// 		pol6->SetParameter(3,3.57431033334042E+03);
// 		pol6->SetParameter(2,-6.68338463467988E+06);
// 		pol6->SetParameter(1,5.55414233197911E+09);
// 		pol6->SetParameter(0,-1.73088291260178E+12);
		//-7.16833333334762E-01x4 + 3.57431033334042E+03x3 - 6.68338463467988E+06x2 + 5.55414233197911E+09x - 1.73088291260178E+12
// 		pol6->SetParLimits(6,2.3086356202952E-03,2.30863562029526E-03);
// 		pol6->SetParLimits(5,-1.7279455679749E+01,-1.72794556797498E+01);
// 		pol6->SetParLimits(4,5.3887808519897E+04,5.38878085198977E+04);
// 		pol6->SetParLimits(3,-8.9628630312873E+07,-8.96286303128730E+07);
// 		pol6->SetParLimits(2,8.3853766054740E+10,8.38537660547409E+10);
// 		pol6->SetParLimits(1,-4.1840239883899E+13,-4.18402398838992E+13);
// 		pol6->SetParLimits(0,8.6986443332496E+15,8.69864433324962E+15);
		//1248 slice Tau=1400 y = 3.11599394933637E-04x6 - 2.34828541084312E+00x5 + 7.37319826797693E+03x4 - 1.23458873643384E+07x3 + 1.16271826293853E+10x2 - 5.83968295468444E+12x + 1.22196074818330E+15
		//1248 slice Tau=1400 y = 2.30863562029526E-03x6 - 1.72794556797498E+01x5 + 5.38878085198977E+04x4 - 8.96286303128730E+07x3 + 8.38537660547409E+10x2 - 4.18402398838992E+13x + 8.69864433324962E+15

// 		pol6->SetParameter(6,0);
// 		pol6->SetParameter(5,-4.30583333325103E-02);
// 		pol6->SetParameter(4,2.67973398432374E+02);
// 		pol6->SetParameter(3,-6.67090190726813E+05);
// 		pol6->SetParameter(2,8.30320658309999E+08);
// 		pol6->SetParameter(1,-5.16744131369674E+11);
// 		pol6->SetParameter(0,1.28636339740528E+14);
		//1248 slice Tau=1400 y = -4.30583333325103E-02x5 + 2.67973398432374E+02x4 - 6.67090190726813E+05x3 + 8.30320658309999E+08x2 - 5.16744131369674E+11x + 1.28636339740528E+14

// 		pol6->SetParameter(6,0);
// 		pol6->SetParameter(5,-9.10468750009841E-03);
// 		pol6->SetParameter(4,5.65382161464470E+01);
// 		pol6->SetParameter(3,-1.40433165459864E+05);
// 		pol6->SetParameter(2,1.74404529626795E+08);
// 		pol6->SetParameter(1,-1.08294469701295E+11);
// 		pol6->SetParameter(0,2.68970828188248E+13);
		//1248 slice Tau=2600 y = -9.10468750009841E-03x5 + 5.65382161464470E+01x4 - 1.40433165459864E+05x3 + 1.74404529626795E+08x2 - 1.08294469701295E+11x + 2.68970828188248E+13

// 		pol6->SetParameter(6,4.63064233713339E-04);
// 		pol6->SetParameter(5,-6.19927835329918E+00);
// 		pol6->SetParameter(4,3.45803467261303E+04);
// 		pol6->SetParameter(3,-1.02876416307665E+08);
// 		pol6->SetParameter(2,1.72156881480732E+11);
// 		pol6->SetParameter(1,-1.53649525775788E+14);
//  	pol6->SetParameter(0,5.71381755372138E+16);

 		//2234 slice Tau=270 y = 4.63064233713339E-04x6 - 6.19927835329918E+00x5 + 3.45803467261303E+04x4 - 1.02876416307665E+08x3 + 1.72156881480732E+11x2 - 1.53649525775788E+14x + 5.71381755372138E+16
// 		pol6->SetParameter(6,-0.00000003962278);
// 		pol6->SetParameter(5,0.00122413264673);
// 		pol6->SetParameter(4,-15.75798079119190);
// 		pol6->SetParameter(3,108186.02957838600000);
// 		pol6->SetParameter(2,-417795857.200374);
// 		pol6->SetParameter(1,860509384185.806);
// 		pol6->SetParameter(0,-738474766394845);

		//y = -0.00000003962278 x6 + 0.00122413264673 x5 - 15.75798079119190 x4 + 108186.02957838600000 x3 - 417795857.20037400000000 x2 + 860509384185.80600000000000 x - 738474766394845.00000000000000 

// 		pol6->SetParLimits(6,9.68031e-014,9.68031e-014);
// 		pol6->SetParLimits(5,-1.52618e-011,-1.52618e-011);
// 		pol6->SetParLimits(4,-2.12054e-006,-2.12054e-006);
// 		pol6->SetParLimits(3,-0.0174038,-0.0174038);
// 		pol6->SetParLimits(2,-79.6153,-79.6153);
// 		pol6->SetParLimits(1,-22982.9,-22982.9);
// 		pol6->SetParLimits(0,5.13151e+009,5.13151e+009);

// 		pol6->SetParameter(6,-3.96227807231142E-08);
// 		pol6->SetParameter(5,1.22413264673238E-03);
// 		pol6->SetParameter(4,-1.57579807911919E+01);
// 		pol6->SetParameter(3,1.08186029578386E+05);
// 		pol6->SetParameter(2,-4.17795857200374E+08);
//		pol6->SetParameter(1,8.60509384185806E+11);
//		pol6->SetParameter(0,-7.38474766394845E+14);
//		y = -3.96227807231142E-08x6 + 1.22413264673238E-03x5 - 1.57579807911919E+01x4 + 1.08186029578386E+05x3 - 4.17795857200374E+08x2 + 8.60509384185806E+11x - 7.38474766394845E+14

// 		pol6->SetParameter(6,4.88903356981041E-05);
// 		pol6->SetParameter(5,-1.45809122344581E+00);
// 		pol6->SetParameter(4,1.81190322994996E+04);
// 		pol6->SetParameter(3,-1.20083655229863E+08);
// 		pol6->SetParameter(2,4.47667116178120E+11);
// 		pol6->SetParameter(1,-8.90071890356739E+14);
// 		pol6->SetParameter(0,7.37366841861574E+17);
//		y = 4.88903356981041E-05x6 - 1.45809122344581E+00x5 + 1.81190322994996E+04x4 - 1.20083655229863E+08x3 + 4.47667116178120E+11x2 - 8.90071890356739E+14x + 7.37366841861574E+17


		//graph1[i]->Fit("pol1");//pol1 can be used directly without TF1 constructor in CINT
		graph1[i]->Fit("pol6","ME","",minrange,maxrange);		graph1[i]->Fit("pol6","ME","",minrange,maxrange);		graph1[i]->Fit("pol6","ME","",minrange,maxrange);		graph1[i]->Fit("pol6","ME","",minrange,maxrange);
		TFitResultPtr r_chi2 = graph1[i]->Fit("pol6","MES","",minrange,maxrange);//"S" means the result of the fit is returned in the TFitResultPtr
		//pol6->Draw();
		//“E” Perform better errors estimation using the Minos technique.
		//“M” Improve fit results, by using the IMPROVE algorithm of TMinuit.
		//If the option "S" is instead used, TFitResultPtr contains the TFitResult and behaves as a smart pointer to it.
		//The fit parameters, error and chi2 (but not covariance matrix) can be retrieved also from the fitted function.
//		TH1D *hint_chi2 = new TH1D("hint_chi2", "Fitted exponential with conf.band", 9000, 0, 9000);//Create a histogram to hold the confidence intervals
//		(TVirtualFitter::GetFitter())->GetConfidenceIntervals(hint_chi2, 0.683);//By default the intervals are corrected using the chi2/ndf value of the fit if a chi2 fit is performed
		//hint_chi2 will contain the CL result that you can draw on top of your fitted graph.
		//where hint_chi2 will hold the errors and could superimpose it on the same canvas where you plot central values.
		//The method TVirtualFitter::GetFitter())->GetConfidenceIntervals(hint_chi2, 0.95) computes the confidence level of your model (fitting) function after having it fitted to an histogram.
//		hint_chi2->SetStats(kFALSE);
//		hint_chi2->SetFillColor(kRed-9);
		TMatrixDSym cov = r_chi2->GetCovarianceMatrix();//useless for now
//		r_chi2->GetConfidenceIntervals(Npoint, 1, 1, energy_point, err_point, 0.683, false);//(Number of x points, 1, 1, x, err, confidence level, false); norm is a flag to control if the intervals need to be normalized to the chi2/ndf value.
// 		for(ii=0;ii<Npoint;ii++)
// 		{
// 			outfile<<"SeGA_"<<i<<"	Eg=	"<<energy_point[ii]<<"	Chi=	"<<pol1->Eval(energy_point[ii])<<"	err_sig=	"<<err_point[ii]<<endl;
// 		}
		graph1[i]->Draw("AP");//"A": Axis are drawn around the graph, "P": The current marker is plotted at each point
//		hint_chi2->Draw("e3 same");
//		graph1[i]->Draw("P same");//draw the points again above error band
		
		p0[i]=pol6->GetParameter(0);
		p0err[i]=pol6->GetParError(0);
		p1[i]=pol6->GetParameter(1);
		p1err[i]=pol6->GetParError(1);
		p2[i]=pol6->GetParameter(2);
		p2err[i]=pol6->GetParError(2);
		p3[i]=pol6->GetParameter(3);
		p3err[i]=pol6->GetParError(3);
		p4[i]=pol6->GetParameter(4);
		p4err[i]=pol6->GetParError(4);
		p5[i]=pol6->GetParameter(5);
		p5err[i]=pol6->GetParError(5);
		p6[i]=pol6->GetParameter(6);
		p6err[i]=pol6->GetParError(6);
		parchi[i]=pol6->GetChisquare();
		parNDF[i]=pol6->GetNDF();
		p_value[i]=pol6->GetProb();//This probability is not the “probability that your fit is good.” If you did many fake experiments (draw many random samples of data points from the assumed distribution (your fit function)), this is the percentage of experiments that would give χ2 values ≥ to the one you got in this experiment.

 		double x,y,z,dx,xlow,xhigh;
// 		x=-p1[i]/(2*p2[i]);
// 		y=p2[i]*x[ii]*x[ii]+p1[i]*x[ii]+p0[i];
// 		z=p0[i]-y-1;
// 		xlow=(-p1[i]-sqrt(p1[i]*p1[i]-4*p2[i]*z))/(2*p2[i]);
// 		xhigh=(-p1[i]+sqrt(p1[i]*p1[i]-4*p2[i]*z))/(2*p2[i]);
// 		dx=x-xlow;

		y=pol6->GetMinimum();//Returns the minimum value of the function on the (xmin, xmax) interval.
		x=pol6->GetMinimumX();//Returns the X value corresponding to the minimum value of the function on the (xmin, xmax) interval.
		xlow=pol6->GetX(y+1,pol6->GetXmin(),x);//Returns the X value corresponding to the function value fy for (xmin<x<xmax)
		xhigh=pol6->GetX(y+1,x,pol6->GetXmax());//Returns the X value corresponding to the function value fy for (xmin<x<xmax)
		//GetXmin is 0; GetXmax is 6000ish.

		TPaveText *textpol6 = new TPaveText(0.30,0.70,0.70,0.90,"brNDC");//left, down, right, up
		textpol6->SetBorderSize(1);
		textpol6->SetFillColor(0);
		textpol6->SetTextAlign(12);//align = 10*HorizontalAlign + VerticalAlign, 12 means水平左对齐、垂直居中对齐
		textpol6->SetTextFont(132);//font = 10 * fontID + precision, 12+2,12 means Symbol; 13+2, 13 means News Times Roman
// 		sprintf(paraprint,"y=p6*x^6+p5*x^5+p4*x^4+p3*x^3+p2*x^2+p1*x+p0");
// 		textpol6->AddText(paraprint);
// 		sprintf(paraprint,"p6=%e+/-%e",p6[i],p6err[i]);
// 		textpol6->AddText(paraprint);
// 		sprintf(paraprint,"p5=%e+/-%e",p5[i],p5err[i]);
// 		textpol6->AddText(paraprint);
// 		sprintf(paraprint,"p4=%e+/-%e",p4[i],p4err[i]);
// 		textpol6->AddText(paraprint);
// 		sprintf(paraprint,"p3=%e+/-%e",p3[i],p3err[i]);
// 		textpol6->AddText(paraprint);
// 		sprintf(paraprint,"p2=%e+/-%e",p2[i],p2err[i]);
// 		textpol6->AddText(paraprint);
// 		sprintf(paraprint,"p1=%e+/-%e",p1[i],p1err[i]);
// 		textpol6->AddText(paraprint);
// 		sprintf(paraprint,"p0=%e+/-%e",p0[i],p0err[i]);
// 		textpol6->AddText(paraprint);
		sprintf(paraprint,"Chisquare=%.3f",parchi[i]);
		textpol6->AddText(paraprint);
		sprintf(paraprint,"NDF=%.0f",parNDF[i]);
		textpol6->AddText(paraprint);
		sprintf(paraprint,"p-value=%f",p_value[i]);
		textpol6->AddText(paraprint);
		sprintf(paraprint,"ymin=%.5f",y);
		textpol6->AddText(paraprint);
		sprintf(paraprint,"Eg=%.5f-%.5f+%.5f keV",x,x-xlow,xhigh-x);
		textpol6->AddText(paraprint);
		if(parchi[i]/parNDF[i]>1)
		{
			sprintf(paraprint,"Eg=%.5f-%.5f+%.5f keV",x,(x-xlow)*sqrt(parchi[i]/parNDF[i]),(xhigh-x)*sqrt(parchi[i]/parNDF[i]));
			textpol6->AddText(paraprint);
		}
		textpol6->Draw();
		//sprintf(hcali_name,"%s.png",hcali_name);
		sprintf(hcali_name,"%s%s%s","D:/X/out/Chi2/",hcali_name,".png");
		canvascali[i]->SaveAs(hcali_name);//存图
		outfile2<<"y=p2*x^2+p1*x+p0	chip2=	"<<setprecision(11)<<p2[i]<<"	+/-	"<<p2err[i]<<"	chip1=	"<<setprecision(8)<<p1[i]<<"	+/-	"<<p1err[i]<<"	chip0=	"<<p0[i]<<"	+/-	"<<p0err[i]<<"	Chi2=	"<<parchi[i]<<"	ndf=	"<<parNDF[i]<<"	p-value=	"<<p_value[i]<<"	Chi2min=	"<<y<<setprecision(9)<<"	Eg=	"<<x<<"	+/-	"<<dx<<endl;

/*




		//************ tau fit by sqrt + const func *****************************************

		sprintf(hcali_name,"%s%d","Tau of SeGA",i);
		canvascali[i]=new TCanvas(hcali_name,hcali_name,900,600);//建立画布
		canvascali[i]->cd();//进入画布
		//canvascali[i]->SetGrid();//显示网格
		//		graph1[i]=new TGraph(peaknum,Energy,tau);//TGraph *gr1=new TGraph(n,x,y);
		graph3[i]= new TGraphErrors(peaknum,Energy[i],tau[i],meanerr[i],tauerr[i]);//画error bars TGraph(n,x,y,ex,ey);
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
		r_tau_sqrt_const->GetConfidenceIntervals(Npoint, 1, 1, energy_point, err_point, 0.683, false);//(Number of x points, 1, 1, x, err, confidence level, norm); norm is a flag to control if the intervals need to be normalized to the chi2/ndf value.
		for(ii=0;ii<Npoint;ii++)
		{
			outfile<<"SeGA_"<<i<<"	Eg=	"<<energy_point[ii]<<"	tau=	"<<func_sqrt_const->Eval(energy_point[ii])<<"	err_tau=	"<<err_point[ii]<<endl;
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
			//x[ii]=energy_point[ii];
			//y[ii]=func_sqrt_const->Eval(energy_point[ii]);//is equal to y[ii]=p1[i]*sqrt(x[ii])+p0[i];
			//dy[ii]=sqrt(sqrt(x[ii])*sqrt(x[ii])*p1err[i]*p1err[i]+p0err[i]*p0err[i]);
			//dy[ii]=err_point[ii];
			//cout<<"	"<<energy_point[ii]<<"	"<<err_point[ii]<<endl;
			errylow[ii]=func_sqrt_const->Eval(energy_point[ii])-err_point[ii];
			erryhigh[ii]=func_sqrt_const->Eval(energy_point[ii])+err_point[ii];
		}
// 		graph3bandlow[i]= new TGraph(Npoint,energy_point,errylow);//画error bars TGraph(n,x,y,ex,ey);
// 		graph3bandlow[i]->SetLineColor(2);
// 		graph3bandlow[i]->SetLineWidth(1);
// 		graph3bandlow[i]->Draw("C");
// 		graph3bandhigh[i]= new TGraph(Npoint,energy_point,erryhigh);//画error bars TGraph(n,x,y,ex,ey);
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
		//		graph1[i]=new TGraph(peaknum,Energy,tau);//TGraph *gr1=new TGraph(n,x,y);
		graph2[i]= new TGraphErrors(peaknum,Energy[i],tau[i],meanerr[i],tauerr[i]);//画error bars TGraph(n,x,y,ex,ey);
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
		TF1 *pol2=new TF1("pol2","pol2",0,60000);//多项式拟合，调节合适的拟合range
		pol2->SetParNames("p0","p1","p2");//y=p2*x^2+p1*x+p0
//		graph2[i]->Fit("pol2");//pol1 can be used directly without TF1 constructor in CINT

		TFitResultPtr r_tau_pol2 = graph2[i]->Fit("pol2","S");//"S" means the result of the fit is returned in the TFitResultPtr
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
		r_tau_pol2->GetConfidenceIntervals(Npoint, 1, 1, energy_point, err_point, 0.683, false);//(Number of x points, 1, 1, x, err, confidence level, false); norm is a flag to control if the intervals need to be normalized to the chi2/ndf value.
		outfile<<"Eg=	"<<energy_point[3]<<"	tau=	"<<pol2->Eval(energy_point[3])<<"	err_tau=	"<<err_point[3]<<endl;
		graph2[i]->Draw("AP");//"A": Axis are drawn around the graph, "P": The current marker is plotted at each point
		hint_tau_pol2->Draw("e3 same");
		graph2[i]->Draw("P same");//draw the points again above error band

		p0[i]=pol2->GetParameter(0);
		p1[i]=pol2->GetParameter(1);
		p2[i]=pol2->GetParameter(2);
		p0err[i]=pol2->GetParError(0);
		p1err[i]=pol2->GetParError(1);
		p2err[i]=pol2->GetParError(2);
		parchi[i]=pol2->GetChisquare();
		parNDF[i]=pol2->GetNDF();
		p_value[i]=pol2->GetProb();//This probability is not the “probability that your fit is good.” If you did many fake experiments (draw many random samples of data points from the assumed distribution (your fit function)), this is the percentage of experiments that would give χ2 values ≥ to the one you got in this experiment.
		double errylow[Npoint],erryhigh[Npoint];
		for(ii=0;ii<Npoint;ii++)
		{
			//x[ii]=energy_point[ii];
			//y[ii]=pol2->Eval(energy_point[ii]);//is equal to y[ii]=p2[i]*x[ii]*x[ii]+p1[i]*x[ii]+p0[i];
			//dy[ii]=sqrt(x[ii]*x[ii]*x[ii]*x[ii]*p2err[i]*p2err[i]+x[ii]*x[ii]*p1err[i]*p1err[i]+p0err[i]*p0err[i]+(2*p2[i]+p1[i])*(2*p2[i]+p1[i])*dx[ii]);
			//dy[ii]=err_point[ii];
			//cout<<"	"<<energy_point[ii]<<"	"<<err_point[ii]<<endl;
			errylow[ii]=pol2->Eval(energy_point[ii])-err_point[ii];
			erryhigh[ii]=pol2->Eval(energy_point[ii])+err_point[ii];
		}
		graph2bandlow[i]= new TGraph(Npoint,energy_point,errylow);//画error bars TGraph(n,x,y,ex,ey);
		graph2bandlow[i]->SetLineColor(2);
		graph2bandlow[i]->SetLineWidth(1);
		graph2bandlow[i]->Draw("C");
		graph2bandhigh[i]= new TGraph(Npoint,energy_point,erryhigh);//画error bars TGraph(n,x,y,ex,ey);
		graph2bandhigh[i]->SetLineColor(2);
		graph2bandhigh[i]->SetLineWidth(1);
		graph2bandhigh[i]->Draw("C");

		TPaveText *textpol6 = new TPaveText(0.10,0.64,0.34,0.89,"brNDC");//left, down, right, up
		textpol6->SetBorderSize(1);
		textpol6->SetFillColor(0);
		textpol6->SetTextAlign(12);//align = 10*HorizontalAlign + VerticalAlign, 12 means水平左对齐、垂直居中对齐
		textpol6->SetTextFont(132);//font = 10 * fontID + precision, 12+2,12 means Symbol; 13+2, 13 means News Times Roman
		sprintf(paraprint,"y=p2*x^2+p1*x+p0");
		textpol6->AddText(paraprint);
		sprintf(paraprint,"p2=%.9f+/-%.9f",p2[i],p2err[i]);
		textpol6->AddText(paraprint);
		sprintf(paraprint,"p1=%.5f+/-%.5f",p1[i],p1err[i]);
		textpol6->AddText(paraprint);
		sprintf(paraprint,"p0=%.5f+/-%.5f",p0[i],p0err[i]);
		textpol6->AddText(paraprint);
		sprintf(paraprint,"Chi2=%.2f",parchi[i]);
		textpol6->AddText(paraprint);
		sprintf(paraprint,"NDF=%.0f",parNDF[i]);
		textpol6->AddText(paraprint);
		sprintf(paraprint,"p-val=%e",p_value[i]);
		textpol6->AddText(paraprint);
		textpol6->Draw();
		//sprintf(hcali_name,"%s.png",hcali_name);
		sprintf(hcali_name,"%s%s%s%s%s","D:/X/out/responseforDoppler/",hcali_name,"_",tauflag,"_band.png");
		canvascali[i]->SaveAs(hcali_name);//存图
		outfile2<<"	SeGA_"<<i<<"	y=p2*x^2+p1*x+p0	taup2=	"<<setprecision(11)<<p2[i]<<"	+/-	"<<p2err[i]<<"	taup1=	"<<setprecision(8)<<p1[i]<<"	+/-	"<<p1err[i]<<"	taup0=	"<<p0[i]<<"	+/-	"<<p0err[i]<<"	Chi2=	"<<parchi[i]<<"	ndf=	"<<parNDF[i]<<"	p-value=	"<<p_value[i];





		//************ tau fit by sqrt func *****************************************

		sprintf(hcali_name,"%s%d","AlSi_tau_sqrt_SeGA_",i);
		canvascali[i]=new TCanvas(hcali_name,hcali_name,900,600);//建立画布
		canvascali[i]->cd();//进入画布
		//canvascali[i]->SetGrid();//显示网格
		//		graph1[i]=new TGraph(peaknum,Energy,tau);//TGraph *gr1=new TGraph(n,x,y);
		graph4[i]= new TGraphErrors(peaknum,Energy[i],tau[i],meanerr[i],tauerr[i]);//画error bars TGraph(n,x,y,ex,ey);
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
		r_tau_sqrt->GetConfidenceIntervals(Npoint, 1, 1, energy_point, err_point, 0.683, false);//(Number of x points, 1, 1, x, err, confidence level, norm); norm is a flag to control if the intervals need to be normalized to the chi2/ndf value.
		outfile<<"Eg=	"<<energy_point[3]<<"	tau=	"<<func_sqrt->Eval(energy_point[3])<<"	err_tau=	"<<err_point[3]<<endl;
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
			//x[ii]=energy_point[ii];
			//y[ii]=func_sqrt->Eval(energy_point[ii]);//is equal to y[ii]=p1[i]*sqrt(x[ii])+p0[i];
			//dy[ii]=sqrt(sqrt(x[ii])*sqrt(x[ii])*p1err[i]*p1err[i]+p0err[i]*p0err[i]);
			//dy[ii]=err_point[ii];
			//cout<<"	"<<energy_point[ii]<<"	"<<err_point[ii]<<endl;
			errylow[ii]=func_sqrt->Eval(energy_point[ii])-err_point[ii];
			erryhigh[ii]=func_sqrt->Eval(energy_point[ii])+err_point[ii];
		}
		graph4bandlow[i]= new TGraph(Npoint,energy_point,errylow);//画error bars TGraph(n,x,y,ex,ey);
		graph4bandlow[i]->SetLineColor(2);
		graph4bandlow[i]->SetLineWidth(1);
		graph4bandlow[i]->Draw("C");
		graph4bandhigh[i]= new TGraph(Npoint,energy_point,erryhigh);//画error bars TGraph(n,x,y,ex,ey);
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