#include <iostream>
#include <fstream>
#include <iomanip.h>
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
void Fill2D()// read data from a txt file and fill a TH2D and find the minimum and then output the minimum+2.30 using a grid search
{
	int Eg=3076;//DSLmodify
//1248, 2234, 3076, 4971, 5156, 5141, 4156, 3435, 2186, 2838, 4270, 3541, 5294
//  49,     47,      46,     42,     42,     39,      39,      45,      45,     44,     41,     40,      39
	if (Eg==1248) Ea=49;
	if (Eg==2234) Ea=47;
	if (Eg==3076) Ea=46;
	if (Eg==4971) Ea=42;
	if (Eg==5156) Ea=42;
	if (Eg==5141) Ea=39;
	if (Eg==4156) Ea=39;
	if (Eg==3435) Ea=45;
	if (Eg==2186) Ea=45;
	if (Eg==2838) Ea=44;
	if (Eg==4270) Ea=41;
	if (Eg==3541) Ea=40;
	if (Eg==5294) Ea=39;
	bool flag_smooth=true;//DSLmodify
	char flag[20];
	sprintf(flag,"base");//DSLmodify
//	sprintf(flag,"SP0.9");
// 	sprintf(flag,"SP1.1");
	int Egfirst,Egbinwidth;
	int Taufirst,Taubinwidth;
	double Nbinsx_Tau, Xlow, Xup, Nbinsy_Eg, Ylow, Yup;
	double input_value[100][100];
	double Sum_input_X[100]={0};
	double Sum_input_Y[100]={0};
	int i,j,k=0;
	double c;
	char h_name[300],b_name[300],rootname[300],datname[300];
	sprintf(h_name,"%s%d%s%d%s%s","Chi2_Ea",Ea,"Eg",Eg,"_2D_",flag);
	sprintf(b_name,"%s%s","D:/X/out/DSL/",h_name);
	if (flag_smooth==false) sprintf(rootname,"%s%s",b_name,"histo_nosmooth.root");
	if (flag_smooth==true) sprintf(rootname,"%s%s",b_name,"histo_smooth.root");
	sprintf(datname,"%s%s",b_name,".dat");
	TFile *fout = new TFile(rootname,"RECREATE");
	if (Eg==1248){ Nbinsx_Tau=58, Xlow=0, Xup=5700, Nbinsy_Eg=34, Ylow=1230.5, Yup=1264.5; Taufirst=0; Egfirst=1231; Taubinwidth=100; Egbinwidth=1; }//Eg=1231->1264, 33 keV, 34 bins //low and up should be half-bin broader in order to make bin display correctly.
	if (Eg==2234){ Nbinsx_Tau=43, Xlow=75, Xup=505, Nbinsy_Eg=34, Ylow=2216.5, Yup=2250.5; Taufirst=80; Egfirst=2217; Taubinwidth=10; Egbinwidth=1; }//Eg=2217->2250, 33 keV, 34 bins
	if (Eg==3076){ Nbinsx_Tau=51, Xlow=-0.5, Xup=50.5, Nbinsy_Eg=34, Ylow=3058.5, Yup=3092.5; Taufirst=0; Egfirst=3059; Taubinwidth=1; Egbinwidth=1; }//Eg=3059->3092, 33 keV, 34 bins
	if (Eg==4971){ Nbinsx_Tau=51, Xlow=-0.5, Xup=50.5, Nbinsy_Eg=34, Ylow=4953.5, Yup=4987.5; Taufirst=0; Egfirst=4954; Taubinwidth=1; Egbinwidth=1; }//Eg=4954->4987, 33 keV, 34 bins
	if (Eg==5156){ Nbinsx_Tau=51, Xlow=-0.5, Xup=50.5, Nbinsy_Eg=34, Ylow=5138.5, Yup=5172.5; Taufirst=0; Egfirst=5139; Taubinwidth=1; Egbinwidth=1; }//Eg=5139->5172, 33 keV, 34 bins
	//if (Eg==4045){ Nbinsx_Tau=51, Xlow=-0.5, Xup=50.5, Nbinsy_Eg=34, Ylow=4027.5, Yup=4061.5; Taufirst=0; Egfirst=4028; Taubinwidth=1; Egbinwidth=1; }//Eg=4028->4061, 33 keV, 34 bins
	if (Eg==5141){ Nbinsx_Tau=51, Xlow=-0.5, Xup=50.5, Nbinsy_Eg=34, Ylow=5123.5, Yup=5157.5; Taufirst=0; Egfirst=5124; Taubinwidth=1; Egbinwidth=1; }//Eg=5124->5157, 33 keV, 34 bins
	if (Eg==4156){ Nbinsx_Tau=51, Xlow=-0.5, Xup=50.5, Nbinsy_Eg=34, Ylow=4138.5, Yup=4172.5; Taufirst=0; Egfirst=4139; Taubinwidth=1; Egbinwidth=1; }//Eg=4139->4172, 33 keV, 34 bins
	if (Eg==3435){ Nbinsx_Tau=51, Xlow=-0.5, Xup=50.5, Nbinsy_Eg=34, Ylow=3417.5, Yup=3451.5; Taufirst=0; Egfirst=3418; Taubinwidth=1; Egbinwidth=1; }//Eg=3418->3451, 33 keV, 34 bins
	if (Eg==2186){ Nbinsx_Tau=91, Xlow=-0.5, Xup=90.5, Nbinsy_Eg=34, Ylow=2168.5, Yup=2202.5; Taufirst=0; Egfirst=2169; Taubinwidth=1; Egbinwidth=1; }//Eg=2169->2202, 33 keV, 34 bins
	if (Eg==2838){ Nbinsx_Tau=51, Xlow=-0.5, Xup=50.5, Nbinsy_Eg=34, Ylow=2820.5, Yup=2854.5; Taufirst=0; Egfirst=2821; Taubinwidth=1; Egbinwidth=1; }//Eg=2821->2854, 33 keV, 34 bins
	if (Eg==4270){ Nbinsx_Tau=51, Xlow=-0.5, Xup=50.5, Nbinsy_Eg=34, Ylow=4252.5, Yup=4286.5; Taufirst=0; Egfirst=4253; Taubinwidth=1; Egbinwidth=1; }//Eg=4253->4286, 33 keV, 34 bins
	if (Eg==3541){ Nbinsx_Tau=51, Xlow=-0.5, Xup=50.5, Nbinsy_Eg=34, Ylow=3523.5, Yup=3557.5; Taufirst=0; Egfirst=3524; Taubinwidth=1; Egbinwidth=1; }//Eg=3524->3557, 33 keV, 34 bins
	if (Eg==5294){ Nbinsx_Tau=51, Xlow=-0.5, Xup=50.5, Nbinsy_Eg=34, Ylow=5276.5, Yup=5310.5; Taufirst=0; Egfirst=5277; Taubinwidth=1; Egbinwidth=1; }//Eg=5277->5310, 33 keV, 34 bins
	//TGraph2D *gChisquare2D = new TGraph2D();//create a graph without parameters.
	//TGraph2D(const char* name, const char* title, Int_t n, Double_t* x, Double_t* y, Double_t* z);//if 20 x, 30 y, there should be 600 z, n=600.
	TH2D *hChisquare2D = new TH2D(h_name,h_name, Nbinsx_Tau, Xlow, Xup, Nbinsy_Eg, Ylow, Yup);//create and name a histogram
	sprintf(h_name,"%s%s",h_name,"min");
	TH2D *hChisquare2Dmin = new TH2D(h_name,h_name, Nbinsx_Tau*100, Xlow, Xup, Nbinsy_Eg*100, Ylow, Yup);//create and name a histogram
	//TGraph *gr1=new TGraph(n,x,y);
// 	TH1F *hChisquare1DX = new TH1F("Tau","Tau", Nbinsx_Tau, Xlow, Xup);//create and name a histogram
// 	TH1F *hChisquare1DY = new TH1F("Eg","Eg", Nbinsy_Eg, Ylow, Yup);//create and name a histogram
	//TH2F (const char *name, const char *title, Int_t nbinsx, Double_t xlow, Double_t xup, Int_t nbinsy, Double_t ylow, Double_t yup)
	fstream InputByChar_infile;
	InputByChar_infile.open(datname,ios::in);
	cout<<datname<<endl;
	while(!InputByChar_infile.eof())//read in values from txt
	{
		for(i=0;i<Nbinsx_Tau;i++)//i is row number, j is column number
		{
			for(j=0;j<Nbinsy_Eg;j++)
			{
				InputByChar_infile>>input_value[i][j];
				//cout<<input_value[i][j]<<"	";
			}
			//cout<<endl;
		}
	}
// 	for(i=0;i<Nbinsx_Tau;i++)//i is row number, j is column number
// 	{
// 		for(j=0;j<Nbinsy_Eg;j++)
// 		{
// 			gChisquare2D->SetPoint(k,Taufirst+i*Taubinwidth,Egfirst+j*Egbinwidth,input_value[i][j]);
// 			cout<<k<<"	"<<Taufirst+i*Taubinwidth<<"	"<<Egfirst+j*Egbinwidth<<"	"<<input_value[i][j]<<endl;
// 			//g->SetPoint(n,x,y,z);//n starts from 0, if 20 x, 30 y, there should be 600 z, last k=599.
// 			k++;
// 		}
// 	}
// 	for(i=0;i<Nbinsx_Tau;i++)//i is row, j is column
// 	{
// 		for(j=0;j<Nbinsy_Eg;j++)
// 		{
// 			Sum_input_X[i]+=input_value[i][j];
// 		}
// 	}
// 
// 	for(j=0;j<Nbinsy_Eg;j++)//i is row, j is column
// 	{
// 		for(i=0;i<Nbinsx_Tau;i++)
// 		{
// 			Sum_input_Y[j]+=input_value[i][j];
// 		}
// 	}
// 
// 	for(i=0;i<Nbinsx_Tau;i++)//i is row, j is column
// 	{
// 		for (k=0;k<Sum_input_X[i];k++)
// 		{
// 			hChisquare1DX->Fill(Taufirst+i*Taubinwidth);//Fill histograms z by event
// 		}
// 	}
// 
// 	for(j=0;j<Nbinsy_Eg;j++)//i is row, j is column
// 	{
// 		for (k=0;k<Sum_input_Y[j];k++)
// 		{
// 			hChisquare1DY->Fill(Egfirst+j*Egbinwidth);//Fill histograms z by event
// 		}
// 	}

	for(i=0;i<Nbinsx_Tau;i++)
	{
		for(j=0;j<Nbinsy_Eg;j++)
		{
			//cout<<input_value[i][j]<<"	";
			for (k=0;k<input_value[i][j];k++)
			{
				hChisquare2D->Fill(Taufirst+i*Taubinwidth,Egfirst+j*Egbinwidth);//Fill histograms z by event
			}
		}
		//cout<<endl;
	}

	TCanvas* c1   = new TCanvas("c1","c1",1000,700);
	c1->cd();
	hChisquare2D->GetXaxis()->SetTitle("#tau (fs)");
	hChisquare2D->GetYaxis()->SetTitle("E_{#gamma} (keV)");
	hChisquare2D->GetZaxis()->SetTitle("#chi^{2}");
	hChisquare2D->GetXaxis()->SetTitleOffset(1);
	hChisquare2D->GetYaxis()->SetTitleOffset(1.5);
	hChisquare2D->GetZaxis()->SetTitleOffset(1.5);
	hChisquare2D->GetXaxis()->CenterTitle();
	hChisquare2D->GetYaxis()->CenterTitle();
	hChisquare2D->GetZaxis()->CenterTitle();
	hChisquare2D->SetContour(99);
	hChisquare2D->Draw("cont4");//surf2 for 3D, cont4 for 2D
	hChisquare2D->SetStats(0);
	
//	gStyle->SetPalette(1);
// 	gChisquare2D->SetTitle(h_name);
// 	gChisquare2D->SetName(h_name);
// 	gChisquare2D->SetNpx(Nbinsx_Tau);
// 	gChisquare2D->SetNpy(Nbinsy_Eg);
// 	gStyle->SetNumberContours(99);
// 	gChisquare2D->Draw("surf2");
 	c1->Update();//The axis settings (title, ranges etc ...) can be changed accessing the axis via the GetXaxis GetYaxis and GetZaxis methods. They access the histogram axis created at drawing time only. Therefore they should called after the TGraph2D is drawn
// 	gChisquare2D->GetXaxis()->SetTitle("#tau (fs)");
// 	gChisquare2D->GetYaxis()->SetTitle("E_{#gamma} (keV)");
// 	gChisquare2D->GetZaxis()->SetTitle("#chi^{2}");
// 	gChisquare2D->GetXaxis()->SetTitleOffset(2.0);
// 	gChisquare2D->GetYaxis()->SetTitleOffset(2.0);
// 	gChisquare2D->GetZaxis()->SetTitleOffset(1.5);
// 	gChisquare2D->GetXaxis()->CenterTitle();
// 	gChisquare2D->GetYaxis()->CenterTitle();
// 	gChisquare2D->GetZaxis()->CenterTitle();
	
	//Projects a 2D histogram into 1D histograms
	if (flag_smooth==true) hChisquare2D->Smooth();
	TH1D* hChisquare1DX=hChisquare2D->ProjectionX("ProjectionX");//ProjectionX->Draw();
	TH1D* hChisquare1DY=hChisquare2D->ProjectionY("ProjectionY");//ProjectionY->Draw();

	double Chi2min=hChisquare2D->GetMinimum();//Get z minimum
	cout<<"Chi2min=	"<<Chi2min<<endl;
	double Chi2grid=0;
// 	double Y_Eg_low[1000]={0};
// 	double Y_Eg_high[1000]={0};
// 	double X_Tau_low[1000]={0};
// 	double X_Tau_high[1000]={0};
	ofstream outfile("D:/X/out/Si25/outfile/peakcali.dat",ios::out);
	for(double ii=Taufirst;ii<Taufirst+(Nbinsx_Tau-1)*Taubinwidth;ii=ii+Taubinwidth*0.1)
	{
		for(double jj=Egfirst;jj<Egfirst+(Nbinsy_Eg-1)*Egbinwidth;jj=jj+Egbinwidth*0.1)
		{
			//cout<<ii<<"	"<<jj<<endl;
			Chi2grid=hChisquare2D->Interpolate(ii,jj);//Get z at (x,y);
			//Given a point P(x,y), Interpolate approximates the value via bilinear interpolation based on the four nearest bin centers.
			if (Chi2grid>1&&Chi2grid-Chi2min<2.30)//find the 1¦Ò uncertainty contour
			{
				outfile<<Chi2grid<<"	"<<ii<<"	"<<jj<<endl;
				hChisquare2Dmin->Fill(ii,jj);//Fill 2D histogram with events selected by 1¦Ò uncertainty contour
			}
		}
	}
//	gChisquare2D->SetStats(0);
//	gChisquare2D->Write();
	fout->Write();
	//fout->Close();
}//Fill2D main