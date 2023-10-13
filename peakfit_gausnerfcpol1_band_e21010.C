#include <iostream>
#include <fstream>
//#include <iomanip.h>
#include <math.h>
#include <time.h>
#include <map>
#include <TROOT.h>
#include <TSystem.h>
#include <TCutG.h>
#include "TChain.h"
#include "TStyle.h"
#include "TLine.h"
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
void peakfit_gausnerfcpol1_band_e21010() // get histogram and EMG fit peaks
{
	const int ID1 = 0;// i=ID1//which detector
	const int ID2 = 1;// i<=ID2// modify which detector
	double binwidth = 1;
	int minrange = 0, maxrange = 0, minbin, maxbin;
	float gaplow = 70., gaphigh = 70.;//fitting range随分辨不同调整
	int i, ii, jj, ibin;
	char paraprint[100], histo_name[200], hfit_name[200];
	TH1F* histo[ID2 + 1];//TH1F peak search+gauss fit,creat histograms
	const int peaknum = 14;//search peak numbers, modify
	TF1* fEMG[peaknum];//creat function
	TF1* p[peaknum], * g[peaknum], * b[peaknum], * p2[peaknum];
	Double_t peakx[peaknum], peakxerr[peaknum];
	Double_t peaky[peaknum], peakyerr[peaknum];
	Double_t sig[peaknum], sigerr[peaknum];
	Double_t tau[peaknum], tauerr[peaknum];
	TH1D* h_confidence_interval[ID2 + 1][peaknum];
	TCanvas* canvaspeak[ID2 + 1][peaknum];
	char pathname[150];
	char filename[150];
	sprintf(pathname, "%s", "F:/e21010/unpacked/");
	sprintf(filename, "%s%s", pathname, "lower_bounds.dat");
	ifstream infilelowerbounds(filename, ios::in);
	sprintf(filename, "%s%s", pathname, "upper_bounds.dat");
	ifstream infileupperbounds(filename, ios::in);
	double par[peaknum][6], par_err[peaknum][6];
	double parChi[peaknum], parNDF[peaknum], p_value[peaknum];
	int lower_bound[ID2 + 1][peaknum], upper_bound[ID2 + 1][peaknum];
	for (i = ID1; i <= ID2; i++)//which detector no need to modify
	{
		for (ii = 0; ii < peaknum; ii++)//which peak
		{
			infilelowerbounds >> lower_bound[i][ii];
			infileupperbounds >> upper_bound[i][ii];
		}
	}
	for (i = ID1; i <= ID2; i++)//which detector no need to modify
	{
		for (ii = 0; ii < peaknum; ii++)//which peak
		{
			cout << lower_bound[i][ii] << "	";
			cout << upper_bound[i][ii] << "	";
		}
		cout << endl;
	}

	sprintf(filename, "%s%s", pathname, "peakpara.dat");
	ofstream outfile(filename, ios::out);
	sprintf(filename, "%s%s", pathname, "sum_122_122_226Ra.root");
	// sum_120_120_gamma_bkg
	// sum_121_121_152Eu
	// sum_122_122_226Ra
	TFile* fin = new TFile(filename);//after this statement, you can use any ROOT command1 for this rootfile
	cout << filename << endl;

	for (i = ID1; i <= ID2; i++)//which detector no need to modify
	{
		sprintf(histo_name, "%s%d", "hG", i);
		histo[i] = (TH1F*)fin->Get(histo_name); //Get spectrum
		histo[i]->Rebin(1);
		histo[i]->SetBinErrorOption(TH1::kPoisson);
	}

	for (i = ID1; i <= ID2; i++) // no need to modify
	{
		for (ii = 13; ii <= 13; ii++)// modify which peak in one detector =0<=13
		{
			sprintf(hfit_name, "%s%d%s%d", "hG", i, "_peak", ii);
			canvaspeak[i][ii] = new TCanvas(hfit_name, hfit_name, 1200, 500);//建立画布
			canvaspeak[i][ii]->cd();//进入画布
			canvaspeak[i][ii]->SetTopMargin(0.02);
			canvaspeak[i][ii]->SetRightMargin(0.03);
			canvaspeak[i][ii]->SetLeftMargin(0.11);
			canvaspeak[i][ii]->SetBottomMargin(0.12);
			histo[i]->SetTitle(hfit_name);//图名
			histo[i]->GetXaxis()->SetTitle("ADC channel");//轴名
			histo[i]->GetYaxis()->SetTitle("Counts per channel");//轴名
			histo[i]->GetXaxis()->CenterTitle();//居中
			histo[i]->GetYaxis()->CenterTitle();//居中
			histo[i]->GetXaxis()->SetLabelFont(132);//坐标字体
			histo[i]->GetYaxis()->SetLabelFont(132);//坐标字体
			histo[i]->GetXaxis()->SetLabelSize(0.05);
			histo[i]->GetYaxis()->SetLabelSize(0.05);
			histo[i]->GetXaxis()->SetTitleFont(132);//轴名字体
			histo[i]->GetYaxis()->SetTitleFont(132);//轴名字体
			histo[i]->GetXaxis()->SetTitleOffset(0.9);//轴名偏移
			histo[i]->GetYaxis()->SetTitleOffset(0.9);//轴名偏移
			histo[i]->GetXaxis()->SetTitleSize(0.06);
			histo[i]->GetYaxis()->SetTitleSize(0.06);
			histo[i]->SetLineWidth(2);
			histo[i]->Draw("e");

			peaky[ii] = 0; peakx[ii] = 0; sig[ii] = 0; tau[ii] = 0;
			histo[i]->GetXaxis()->SetRangeUser(lower_bound[i][ii], upper_bound[i][ii]);
			peaky[ii] = histo[i]->GetMaximum();
			peakx[ii] = histo[i]->GetBinCenter(histo[i]->GetMaximumBin());
			if (ii == 0) { gaplow = 30; gaphigh = 30; }//40K 1461
			if (ii == 1) { gaplow = 50; gaphigh = 60; }//208Tl 2615
			if (ii == 2) { gaplow = 10; gaphigh = 10; }//152Eu 121
			if (ii == 3) { gaplow = 13; gaphigh = 13; }//152Eu	244
			if (ii == 4) { gaplow = 15; gaphigh = 14; }//152Eu	344
			if (ii == 5) { gaplow = 20; gaphigh = 19; }//152Eu	779
			if (ii == 6) { gaplow = 30; gaphigh = 28; }//152Eu	964
			if (ii == 7) { gaplow = 32; gaphigh = 30; }//152Eu	1408
			if (ii == 8) { gaplow = 15; gaphigh = 14; }//226Ra	352
			if (ii == 9) { gaplow = 18; gaphigh = 17; }//226Ra	609
			if (ii == 10) { gaplow = 32; gaphigh = 30; }//226Ra 1120
			if (ii == 11) { gaplow = 28; gaphigh = 32; }//226Ra 1764
			if (ii == 12) { gaplow = 48; gaphigh = 46; }//226Ra 2204
			if (ii == 13) { gaplow = 52; gaphigh = 48; }//226Ra 2448

			histo[i]->GetXaxis()->SetRangeUser(peakx[ii] - gaplow, peakx[ii] + gaphigh);//zoom the axis
			//cout<<"************"<<peakx[ii]<<"	"<<peaky[ii]<<endl;
			minrange = peakx[ii] - gaplow;
			maxrange = peakx[ii] + gaphigh;
			minbin = histo[i]->FindBin(minrange);
			maxbin = histo[i]->FindBin(maxrange);
			ibin = minbin;
			while (histo[i]->GetBinContent(ibin) < (peaky[ii] / 2))
			{
				ibin++;
				if (ibin >= maxbin)break;
			}
			double sigmaguess = 2 * (peakx[ii] - histo[i]->GetBinCenter(ibin)) / 2.355;
			float highcounts = 0, lowcounts = 0;
			for (jj = 0; jj < 10; jj++)
			{
				highcounts += histo[i]->GetBinContent(maxbin - jj);
				lowcounts += histo[i]->GetBinContent(minbin + jj);
			}
			highcounts = highcounts / 10; lowcounts = lowcounts / 10;
			double aguess = (highcounts - lowcounts) / (maxrange - minrange);
			double bguess = lowcounts - minrange * aguess;
			cout << aguess << "	" << bguess << "	" << sigmaguess << "	" << peakx[ii] << "	" << peaky[ii] << endl;
			fEMG[ii] = new TF1("fEMG", "[0]*x+[1]+[2]/2/[3]*exp(0.5*([4]*[4]/([3]*[3]))+(x-[5])/[3])*ROOT::Math::erfc(1/sqrt(2)*([4]/[3]+(x-[5])/[4]))", peakx[ii] - gaplow, peakx[ii] + gaphigh);// Sun PRC2021 low-energy tail
			g[ii] = new TF1("g", "gausn", peakx[ii] - gaplow, peakx[ii] + gaphigh);// The [2]-N parameter in total is equivalent to the Constant in gausn
			p[ii]=new TF1("p","[0]*x+[1]-[0]*x-[1]+[2]/2/[3]*exp(0.5*([4]*[4]/([3]*[3]))+(x-[5])/[3])*ROOT::Math::erfc(1/sqrt(2)*([4]/[3]+(x-[5])/[4]))",peakx[ii]-gaplow,peakx[ii]+gaphigh);//pure peak
			//p2[ii] = new TF1("p2", "[0]*exp(-0.5*((x-[1])/[2])^2) / (sqrt(2*3.141592654)*[2])", peakx[ii] - gaplow, peakx[ii] + gaphigh);//pure peak2
			b[ii] = new TF1("b", "[0]*x+[1]", peakx[ii] - gaplow, peakx[ii] + gaphigh);//pure bkg

			// 			fEMG[ii]=new TF1("total","[0]*x+[1]+sqrt(3.141592654/2)*[2]/[3]*[4]*exp(0.5*([4]*[4]/([3]*[3]))+(x-[5])/[3])*ROOT::Math::erfc(1/sqrt(2)*([4]/[3]+(x-[5])/[4]))",peakx[ii]-gaplow,peakx[ii]+gaphigh);// Glassman PRC2019 low-energy tail
			// 			g[ii]=new TF1("g","gausn",peakx[ii]-gaplow,peakx[ii]+gaphigh);// The [2]-N parameter in total is equivalent to the Constant in gausn
			// 			p[ii]=new TF1("p","[0]*x+[1]-[0]*x-[1]+sqrt(3.141592654/2)*[2]/[3]*[4]*exp(0.5*([4]*[4]/([3]*[3]))+(x-[5])/[3])*ROOT::Math::erfc(1/sqrt(2)*([4]/[3]+(x-[5])/[4]))",peakx[ii]-gaplow,peakx[ii]+gaphigh);//pure peak
			// 			b[ii]=new TF1("b","[0]*x+[1]",peakx[ii]-gaplow,peakx[ii]+gaphigh);//pure bkg
			fEMG[ii]->SetNpx((gaphigh + gaplow) * 10);
			g[ii]->SetNpx((gaphigh + gaplow) * 10);
			p[ii]->SetNpx((gaphigh+gaplow)*10);
			// 			p2[ii]->SetNpx((gaphigh + gaplow) * 10);
			b[ii]->SetNpx((gaphigh + gaplow) * 10);
			//fEMG[ii]->SetParameters(0.1,15,peaky[ii],10,10,peakx[ii]);//initial value [0]-A, [1]-B, [2]-N, [3]-τ, [4]-σ, [5]-μ
			fEMG[ii]->SetParameters(aguess, bguess, peaky[ii], 2, sigmaguess, peakx[ii]);//initial value [0]-A, [1]-B, [2]-N, [3]-τ, [4]-σ, [5]-μ
			// 			fEMG[ii]->SetParLimits(0,-500,500);//Bkg A
			// 			fEMG[ii]->SetParLimits(1,-50000,300000);//Bkg B
			fEMG[ii]->SetParLimits(2, 1000, 500000);//Constant,min,max
			fEMG[ii]->SetParLimits(3, 0.0001, 20);//Tau
			fEMG[ii]->SetParLimits(4, 0.2, 16);//Sigma
			fEMG[ii]->SetParLimits(5, 3028, 3088);//Mean
			fEMG[ii]->SetParNames("BkgA", "BkgB", "Const*bin", "Tau", "Sigma", "Mean");

			TFitResultPtr Fit_result_pointer = histo[i]->Fit("fEMG", "MLES", "", peakx[ii] - gaplow, peakx[ii] + gaphigh);
			//"S" means the result of the fit is returned in the TFitResultPtr
			//“E” Perform better errors estimation using the Minos technique.
			//“M” Improve fit results, by using the IMPROVE algorithm of TMinuit.
			//If the option "S" is instead used, TFitResultPtr contains the TFitResult and behaves as a smart pointer to it.
			//The fit parameters, error and chi2 (but not covariance matrix) can be retrieved also from the fitted function.

			fEMG[ii]->GetParameters(par[ii]);//二维数组的par[ii]是地址,pointer to the TF1, GetParameters的数组得是double类型Obtaining the value of parameters and saving them to par[]; 
			par_err[ii][0] = fEMG[ii]->GetParError(0);//Obtaining the error of the 1st parameter
			par_err[ii][1] = fEMG[ii]->GetParError(1);//Obtaining the error of the 2nd parameter
			par_err[ii][2] = fEMG[ii]->GetParError(2);//Obtaining the error of the 3rd parameter
			par_err[ii][3] = fEMG[ii]->GetParError(3);//Obtaining the error of the 4th parameter
			par_err[ii][4] = fEMG[ii]->GetParError(4);//Obtaining the error of the 5th parameter
			par_err[ii][5] = fEMG[ii]->GetParError(5);//Obtaining the error of the 6th parameter
			parChi[ii] = fEMG[ii]->GetChisquare();
			parNDF[ii] = fEMG[ii]->GetNDF();
			p_value[ii] = fEMG[ii]->GetProb();//This probability is not the “probability that your fit is good.” If you did many fake experiments (draw many random samples of data points from the assumed distribution (your fit function)), this is the percentage of experiments that would give χ2 values ≥ to the one you got in this experiment.
			g[ii]->SetParameters(par[ii][2], par[ii][5], par[ii][4]);//set parameters for drawing gausn
			g[ii]->SetLineColor(4);
			p[ii]->SetParameters(par[ii][0],par[ii][1],par[ii][2],par[ii][3],par[ii][4],par[ii][5]);//set parameters for drawing peak
			p[ii]->SetLineColor(6);
			// 			p2[ii]->SetParameters(par[ii][6], par[ii][7], par[ii][8]);//set parameters for drawing peak
			// 			p2[ii]->SetLineColor(4);
			b[ii]->SetParameters(par[ii][0], par[ii][1]);//set parameters for drawing bkg
			b[ii]->SetLineColor(8);
			fEMG[ii]->SetLineWidth(3);
			// The confidence band is not always properly displayed.
			
			// Uncertainty Band
			sprintf(filename, "%s%s%d%s%d", "h_confidence_interval", "_", i, "_", ii);
			h_confidence_interval[i][ii] = new TH1D(filename, filename, 16000, 10, 4010);//Create a histogram to hold the confidence intervals
			TVirtualFitter* fitter = TVirtualFitter::GetFitter();//The method TVirtualFitter::GetFitter())->Get the parameters of your fitting function after having it fitted to an histogram.
			fitter->GetConfidenceIntervals(h_confidence_interval[i][ii], 0.95);//By default the intervals are inflated using the chi2/ndf value of the fit if a chi2 fit is performed
			//confidence interval for the colored band: 1σ confidence interval: P=0.683, 1σ confidence interval: P=0.95, 3σ confidence interval: P=0.997
			//h_confidence_interval will contain the CL result that you can draw on top of your fitted graph.
			//where h_confidence_interval will hold the errors and could superimpose it on the same canvas where you plot central values.
			h_confidence_interval[i][ii]->SetStats(kFALSE);
			h_confidence_interval[i][ii]->SetFillColor(kRed - 10);
			h_confidence_interval[i][ii]->Draw("e3 same"); // plot the uncertainty band
			fEMG[ii]->Draw("same");
//			g[ii]->Draw("same");
// 			p[ii]->Draw("same");
// 			p2[ii]->Draw("same");
			b[ii]->Draw("same");
			histo[i]->Draw("e same");


			TMatrixD cov = Fit_result_pointer->GetCovarianceMatrix();//error matrix
			TMatrixD cor = Fit_result_pointer->GetCorrelationMatrix();//parameter correlation coefficients
			cov.Print();
			cor.Print();
			double Utau_Utau = fitter->GetCovarianceMatrixElement(3, 3);//(Utau)^2
			cout << sqrt(Utau_Utau) << endl;
			double rou_Utau_Umean = fitter->GetCovarianceMatrixElement(3, 5);//(ρ*Utau*Umean), (3,5) (5,3) doesn't matter.
			cout << rou_Utau_Umean << endl;
			double Umean_Umean = fitter->GetCovarianceMatrixElement(5, 5);//(Umean)^2
			cout << sqrt(Umean_Umean) << endl;
			cout << sqrt(Utau_Utau + Umean_Umean + 2 * rou_Utau_Umean) << endl;

			sig[ii] = par[ii][4]; sigerr[ii] = par_err[ii][4];
			tau[ii] = par[ii][3]; tauerr[ii] = par_err[ii][3];
			peaky[ii] = fEMG[ii]->GetMaximum();
			peakx[ii] = fEMG[ii]->GetMaximumX();
			peakxerr[ii] = sqrt(Utau_Utau + Umean_Umean + 2 * rou_Utau_Umean);
			outfile << "G" << i << "	Constant/binsize" << ii << "=	" << par[ii][2] << "	+/-	" << par_err[ii][2] << "	Mean" << ii << "=	" << par[ii][5] << "	+/-	" << par_err[ii][5] << "	Maximum" << ii << "=	" << peakx[ii] << "	+/-	" << peakxerr[ii] << "	Sigma" << ii << "=	" << par[ii][4] << "	+/-	" << par_err[ii][4] << "	Tau" << ii << "=	" << par[ii][3] << "	+/-	" << par_err[ii][3] << "	A" << ii << "=	" << par[ii][0] << "	+/-	" << par_err[ii][0] << "	B" << ii << "=	" << par[ii][1] << "	+/-	" << par_err[ii][1] << "	Chi2" << ii << "=	" << parChi[ii] << "	NDF" << ii << "=	" << parNDF[ii] << "	Area" << ii << "=	" << par[ii][2] / binwidth << endl;//输出文本查看

			TPaveText* textgaus = new TPaveText(0.7, 0.4, 0.99, 0.97, "brNDC");//加标注left, down, right, up
			textgaus->SetBorderSize(1);//边框宽度
			textgaus->SetFillColor(0);//填充颜色
			textgaus->SetTextAlign(12);//align = 10*HorizontalAlign + VerticalAlign, 12 means水平左对齐、垂直居中对齐
			textgaus->SetTextFont(132);//font = 10 * fontID + precision, 12+2,12 means Symbol; 13+2, 13 means News Times Roman
			//text->SetTextColor(2);//文本颜色
			sprintf(paraprint, "Constant/binsize%d=%.2f%s%.3f", ii, par[ii][2] / binwidth, "+/-", par_err[ii][2] / binwidth);//par数组还保持着刚才的参数
			textgaus->AddText(paraprint);
			sprintf(paraprint, "Mean%d=%.4f%s%.4f", ii, par[ii][5], "+/-", par_err[ii][5]);
			textgaus->AddText(paraprint);
			sprintf(paraprint, "Sigma%d=%.3f%s%.3f", ii, par[ii][4], "+/-", par_err[ii][4]);
			textgaus->AddText(paraprint);
			sprintf(paraprint, "Res%d=%.2f%%", ii, par[ii][4] / par[ii][5] * 2.355 * 100);
			textgaus->AddText(paraprint);
			sprintf(paraprint, "Tau%d=%.2f%s%.3f", ii, par[ii][3], "+/-", par_err[ii][3]);
			textgaus->AddText(paraprint);
			sprintf(paraprint, "A%d(%.4f%s%.4f)*x+B%d(%.3f%s%.3f)", ii, par[ii][0], "+/-", par_err[ii][0], ii, par[ii][1], "+/-", par_err[ii][1]);
			textgaus->AddText(paraprint);
			sprintf(paraprint, "Chisquare%d=%.2f", ii, parChi[ii]);
			textgaus->AddText(paraprint);
			sprintf(paraprint, "NDF%d=%.2f", ii, parNDF[ii]);
			textgaus->AddText(paraprint);
			sprintf(paraprint, "p-val=%e", p_value[ii]);
			textgaus->AddText(paraprint);
			sprintf(paraprint, "Maximum%d=%.4f%s%.4f", ii, peakx[ii], "+/-", peakxerr[ii]);//par数组还保持着刚才的参数
			textgaus->AddText(paraprint);
			textgaus->Draw();
			sprintf(filename, "%s%s%s", pathname, hfit_name, ".png");
			canvaspeak[i][ii]->SaveAs(filename);
		}//for(ii=0;ii<peaknum;ii++)
		outfile << "\n\n" << endl;
	}//for (i=0;i<ID;i++)
}//peakcali main
