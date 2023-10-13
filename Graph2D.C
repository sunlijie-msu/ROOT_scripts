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
#include "TGraph2D.h"
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
#include "TF2.h"
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
void Graph2D()// read data from a txt file and plot a Graph2D
{
	int Eg=1248;//DSLmodify  1248, 2234, 3076, 4971, 5156
	int Ea=49;//DSLmodify         49,     47,      46,     42,     42
	int Egfirst,Egbinwidth;
	int Taufirst,Taubinwidth;
	double Nbinsx_Tau, Xlow, Xup, Nbinsy_Eg, Ylow, Yup;
	double input_value[80][80];
	double Sum_input_X[80]={0};
	double Sum_input_Y[80]={0};
	int i,j,k=0;
	double c;
	char h_name[300],g_name[300],b_name[300],rootname[300],datname[300];
	//sprintf(h_name,"%s%d%s%d%s","Prob_Ea",Ea,"Eg",Eg,"_2D");
	sprintf(h_name,"%s%d%s%d%s","Prob_Ea",Ea,"Eg",Eg,"_2D");//DSLmodify
	sprintf(g_name,"%s%s",h_name,"interpolate");//DSLmodify
	sprintf(b_name,"%s%s","D:/X/out/DSL/",h_name);
	sprintf(rootname,"%s%s",b_name,".root");
	sprintf(datname,"%s%s",b_name,".dat");
	TFile *fout = new TFile(rootname,"RECREATE");
	if (Eg==1248){ Nbinsx_Tau=27, Xlow=100, Xup=5500, Nbinsy_Eg=17, Ylow=1231, Yup=1265; Taufirst=200; Egfirst=1232; Taubinwidth=200; Egbinwidth=2; }//Eg=1232->1264, 32 keV, 17 bins //low and up should be half-bin broader in order to make bin display correctly.
	if (Eg==2234){ Nbinsx_Tau=43, Xlow=75, Xup=505, Nbinsy_Eg=17, Ylow=2217, Yup=2251; Taufirst=80; Egfirst=2218; Taubinwidth=10; Egbinwidth=2; }//Eg=2218->2250, 32 keV, 17 bins
	if (Eg==3076){ Nbinsx_Tau=25, Xlow=-0.5, Xup=24.5, Nbinsy_Eg=17, Ylow=3059, Yup=3093; Taufirst=0; Egfirst=3060; Taubinwidth=1; Egbinwidth=2; }//Eg=3060->3092, 32 keV, 17 bins
	if (Eg==4971){ Nbinsx_Tau=25, Xlow=-0.5, Xup=24.5, Nbinsy_Eg=17, Ylow=4954, Yup=4988; Taufirst=0; Egfirst=4955; Taubinwidth=1; Egbinwidth=2; }//Eg=4955->4987, 32 keV, 17 bins
	if (Eg==5156){ Nbinsx_Tau=25, Xlow=-0.5, Xup=24.5, Nbinsy_Eg=17, Ylow=5139, Yup=5173; Taufirst=0; Egfirst=5140; Taubinwidth=1; Egbinwidth=2; }//Eg=5140->5172, 32 keV, 17 bins
	TGraph2D *gChisquare2D = new TGraph2D();//create a graph without parameters.
	//TGraph2D(const char* name, const char* title, Int_t n, Double_t* x, Double_t* y, Double_t* z);//if 20 x, 30 y, there should be 600 z, n=600.
	//TH2D *hChisquare2D = new TH2D(h_name,h_name, Nbinsx_Tau, Xlow, Xup, Nbinsy_Eg, Ylow, Yup);//create and name a histogram
	TGraph2D *gChisquare2Dinterpolate = new TGraph2D();//create and name a histogram
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
	for(i=0;i<Nbinsx_Tau;i++)//i is row number, j is column number
	{
		for(j=0;j<Nbinsy_Eg;j++)
		{
			gChisquare2D->SetPoint(k,Taufirst+i*Taubinwidth,Egfirst+j*Egbinwidth,input_value[i][j]);
			//void TGraph2D::SetPoint ( Int_t  n,  Double_t  x,  Double_t  y,  Double_t  z ) 
			cout<<k<<"	"<<Taufirst+i*Taubinwidth<<"	"<<Egfirst+j*Egbinwidth<<"	"<<input_value[i][j]<<endl;
			//g->SetPoint(n,x,y,z);//n starts from 0, if 20 x, 30 y, there should be 600 z, last k=599.
			k++;
		}
	}
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

// 	for(i=0;i<Nbinsx_Tau;i++)
// 	{
// 		for(j=0;j<Nbinsy_Eg;j++)
// 		{
// 			cout<<input_value[i][j]<<"	";
// 			for (k=0;k<input_value[i][j];k++)
// 			{
// 				hChisquare2D->Fill(Taufirst+i*Taubinwidth,Egfirst+j*Egbinwidth);//Fill histograms z by event
// 			}
// 		}
// 		cout<<endl;
// 	}

	TCanvas* c1   = new TCanvas("c1","c1",1000,700);
	c1->cd();
// 	hChisquare2D->GetXaxis()->SetTitle("#tau (fs)");
// 	hChisquare2D->GetYaxis()->SetTitle("E_{#gamma} (keV)");
// 	hChisquare2D->GetZaxis()->SetTitle("#chi^{2}");
// 	hChisquare2D->GetXaxis()->SetTitleOffset(2.0);
// 	hChisquare2D->GetYaxis()->SetTitleOffset(2.0);
// 	hChisquare2D->GetZaxis()->SetTitleOffset(1.5);
// 	hChisquare2D->GetXaxis()->CenterTitle();
// 	hChisquare2D->GetYaxis()->CenterTitle();
// 	hChisquare2D->GetZaxis()->CenterTitle();
// 	hChisquare2D->SetContour(99);
// 	hChisquare2D->Draw("surf2");
// 	hChisquare2D->SetStats(0);
	
//	gStyle->SetPalette(1);
	gChisquare2D->SetTitle(h_name);
	gChisquare2D->SetName(h_name);
	gChisquare2D->SetNpx(Nbinsx_Tau);
	gChisquare2D->SetNpy(Nbinsy_Eg);
	gStyle->SetNumberContours(99);
	gChisquare2D->Draw("surf2");
//	gChisquare2D->Draw("cont4");
	c1->Update();//The axis settings (title, ranges etc ...) can be changed accessing the axis via the GetXaxis GetYaxis and GetZaxis methods. They access the histogram axis created at drawing time only. Therefore they should called after the TGraph2D is drawn
	gChisquare2D->GetXaxis()->SetTitle("#tau (fs)");
	gChisquare2D->GetYaxis()->SetTitle("E_{#gamma} (keV)");
	gChisquare2D->GetZaxis()->SetTitle("#chi^{2}");
	gChisquare2D->GetXaxis()->SetTitleOffset(2.0);
	gChisquare2D->GetYaxis()->SetTitleOffset(2.0);
	gChisquare2D->GetZaxis()->SetTitleOffset(1.5);
	gChisquare2D->GetXaxis()->CenterTitle();
	gChisquare2D->GetYaxis()->CenterTitle();
	gChisquare2D->GetZaxis()->CenterTitle();
	//create 1D histogram from 2D graph. Projects a 2D graph into 1D histograms
	TH1D *hChisquare1DX=gChisquare2D->Project("x");//Prob_Ea49Eg1248_2D_x->Draw();
	TH1D *hChisquare1DY=gChisquare2D->Project("y");//Prob_Ea49Eg1248_2D_y->Draw();
	
	double Chi2min=gChisquare2D->GetZmax();
	cout<<"Probmax=	"<<Chi2min<<endl;
	double Chi2grid=0;
	k=0;
	for(double ii=Taufirst;ii<Taufirst+(Nbinsx_Tau)*Taubinwidth;ii=ii+Taubinwidth*0.1)
	{
		for(double jj=Egfirst;jj<Egfirst+(Nbinsy_Eg)*Egbinwidth;jj=jj+Egbinwidth*0.1)
		{
			//cout<<ii<<"	"<<jj<<endl;
			Chi2grid=gChisquare2D->Interpolate(ii,jj);//Get z value, Double_t Interpolate(Double_t x, Double_t y)
			//Given a point P(x,y), Interpolate approximates the value via bilinear interpolation based on the four nearest bin centers.
			//if (Chi2grid>1&&Chi2grid-Chi2min<2.30)
			
				//outfile<<Chi2grid<<"	"<<ii<<"	"<<jj<<endl;
				gChisquare2Dinterpolate->SetPoint(k++,ii,jj,Chi2grid);
		}
	}
	gChisquare2Dinterpolate->SetTitle(g_name);
	gChisquare2Dinterpolate->SetName(g_name);
	gChisquare2Dinterpolate->SetNpx(Nbinsx_Tau*10);
	gChisquare2Dinterpolate->SetNpy(Nbinsy_Eg*10);
	//create 1D histogram from 2D histogram
	TH1D *hChisquare1DinterpolateX=gChisquare2Dinterpolate->Project("x");//Prob_Ea49Eg1248_2Dinterpolate_x->Draw();
	TH1D *hChisquare1DinterpolateY=gChisquare2Dinterpolate->Project("y");//Prob_Ea49Eg1248_2Dinterpolate_y->Draw();

//	gChisquare2D->SetStats(0);
	gChisquare2D->Write();
	gChisquare2Dinterpolate->Write();
	fout->Write();
	//fout->Close();
}//Fill2D main