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
void Fill1D_DSL()// read data from a txt file and fill a TH1F and a graph and save as a root file
 // million MCMC sample points → likelihood and posterior (normalized) distributions
{
	int Eg,Egbinwidth;
	int Taufirst,Taubinwidth;
	int Nbinsx_Tau, Nbinsy_E;
	double Xlow, Xup, Ylow, Yup;
	double input_value[80][80];
	double Sum_input_X[80]={0};
	double Sum_input_Y[80]={0};
	int i,j,k, ;
	double c, d, istep;
	double hight90, lowt, hight, centralt;
	TCanvas *hcanvas[30];
	TCanvas *gcanvas[30];
	TH1D *h1DX [30];
	TGraph *g1DX[30];
	double scale_Y[6000]={0};
	double scale_X[6000]={0};
	char paraprint[300], h_name[300],g_name[300],b_name[300],rootname[300],datname[300], jointname[300];

	int un_flag = 0; //modify choose 0 or 1
	int ipeak = 11;
	if (ipeak==0) { Eg=1248; Xlow=0; Xup=30; Nbinsx_Tau=(Xup-Xlow)*100; Nbinsy_E=50; Ylow=3075.44; Yup=3077.04;  }
	if (ipeak==1) { Eg=2234; Xlow=0; Xup=30; Nbinsx_Tau=(Xup-Xlow)*100; Nbinsy_E=50; Ylow=3075.44; Yup=3077.04;  }
	if (ipeak==2) { Eg=3076; Xlow=0; Xup=30; Nbinsx_Tau=(Xup-Xlow)*100; Nbinsy_E=50; Ylow=3075.44; Yup=3077.04;  }
	if (ipeak==3) { Eg=4970; Xlow=0; Xup=20; Nbinsx_Tau=(Xup-Xlow)*100; Nbinsy_E=50; Ylow=3075.44; Yup=3077.04;  }
	if (ipeak==4) { Eg=5156; Xlow=0; Xup=20; Nbinsx_Tau=(Xup-Xlow)*100; Nbinsy_E=50; Ylow=3075.44; Yup=3077.04;  }
	if (ipeak==5) { Eg=5141; Xlow=0; Xup=60; Nbinsx_Tau=(Xup-Xlow)*100; Nbinsy_E=50; Ylow=3075.44; Yup=3077.04;  }
	if (ipeak==6) { Eg=4156; Xlow=0; Xup=60; Nbinsx_Tau=(Xup-Xlow)*100; Nbinsy_E=50; Ylow=3075.44; Yup=3077.04;  }
	if (ipeak==7) { Eg=3435; Xlow=0; Xup=30; Nbinsx_Tau=(Xup-Xlow)*100; Nbinsy_E=50; Ylow=3075.44; Yup=3077.04;  }
	if (ipeak==8) { Eg=2186; Xlow=0; Xup=30; Nbinsx_Tau=(Xup-Xlow)*100; Nbinsy_E=50; Ylow=3075.44; Yup=3077.04;  }
	if (ipeak==9) { Eg=2838; Xlow=0; Xup=60; Nbinsx_Tau=(Xup-Xlow)*100; Nbinsy_E=50; Ylow=3075.44; Yup=3077.04;  }
	if (ipeak==10) { Eg=4270; Xlow=0; Xup=60; Nbinsx_Tau=(Xup-Xlow)*100; Nbinsy_E=50; Ylow=3075.44; Yup=3077.04;  }
	if (ipeak==11) { Eg=3541; Xlow=0; Xup=60; Nbinsx_Tau=(Xup-Xlow)*100; Nbinsy_E=50; Ylow=3075.44; Yup=3077.04;  }
	if (ipeak==12) { Eg=5293; Xlow=0; Xup=60; Nbinsx_Tau=(Xup-Xlow)*100; Nbinsy_E=50; Ylow=3075.44; Yup=3077.04;  }

	if (ipeak==13) { Eg=41565141; Xlow=0; Xup=60; Nbinsx_Tau=(Xup-Xlow)*100; Nbinsy_E=50; Ylow=3075.44; Yup=3077.04;  }
	if (ipeak==14) { Eg=21863435; Xlow=0; Xup=30; Nbinsx_Tau=(Xup-Xlow)*100; Nbinsy_E=50; Ylow=3075.44; Yup=3077.04;  }

	sprintf(h_name,"%s%d","Likelihood",Eg);
	sprintf(g_name,"%s%d","Posterior",Eg);
	sprintf(b_name,"%s%s","D:/X/out/DSL_Final_Tau/",h_name);

	if (un_flag==0)
	{
		sprintf(rootname,"%s%s",b_name,"n.root");
		sprintf(datname,"%s%s",b_name,"n.txt");// Use ANSI txt. Using UTF or Unicode is fucking slow.
		sprintf(jointname,"%s%s",b_name,"forjointn.txt");
	}
	
	if (un_flag==1)
	{
		sprintf(rootname,"%s%s",b_name,".root");
		sprintf(datname,"%s%s",b_name,".txt");// Use ANSI txt. Using UTF or Unicode is fucking slow.
		sprintf(jointname,"%s%s",b_name,"forjoint.txt");
	}
	ofstream outfile(jointname,ios::out); // output p.d.f.
	TFile *fout = new TFile(rootname,"RECREATE");

	hcanvas[ipeak]=new TCanvas(h_name,h_name,900,600);
	hcanvas[ipeak]->cd();
	h1DX[ipeak] = new TH1D(h_name,h_name, Nbinsx_Tau, Xlow, Xup);//create and name a histogram
//	TH2F *h2DX = new TH2F("TauEg","TauEg", Nbinsx_Tau, Xlow, Xup, Nbinsy_E, Ylow, Yup);//create and name a histogram
	//TH1F *hChisquare1DY = new TH1F("Eg","Eg", Nbinsy_E, Ylow, Yup);//create and name a histogram
	//TH2F (const char *name, const char *title, Int_t nbinsx, Double_t xlow, Double_t xup, Int_t nbinsy, Double_t ylow, Double_t yup)
	fstream InputByChar_infile;
	InputByChar_infile.open(datname,ios::in);
	int ievent=0;
	while(!InputByChar_infile.eof())
	{
		//InputByChar_infile>>c>>d;
		InputByChar_infile>>c; //read in million MCMC samples
		ievent++;
		h1DX[ipeak]->Fill(c);
		//h1DX->SetBinContent(i+1,x);//SetBinContent(i starts from 1, bincontent);
	}

	double totalarea=0, ULarea=0;
	for (int ibin = 1; ibin <= h1DX[ipeak]->GetNbinsX(); ibin++)
	{
		totalarea += h1DX[ipeak]->GetBinContent(ibin); //Get total area under pdf curve
	}

	for (int ibin=1; ibin<=h1DX[ipeak]->GetNbinsX(); ibin++)
	{
		ULarea+=h1DX[ipeak]->GetBinContent(ibin);
		if (ULarea>=totalarea*0.683) break; // Get 68% CL
	}
	cout<<"68%UL=	"<<h1DX[ipeak]->GetBinCenter(ibin)<<endl;
	
	ULarea=0;
	for (int ibin=1; ibin<=h1DX[ipeak]->GetNbinsX(); ibin++)
	{
		ULarea+=h1DX[ipeak]->GetBinContent(ibin);
		if (ULarea>=totalarea*0.90) break; // Get 90% CL
	}
	hight90=h1DX[ipeak]->GetBinCenter(ibin);
	cout<<"90%UL=	"<<hight90<<endl;

	gPad->SetRightMargin(0.03);
	h1DX[ipeak]->Draw();
	hcanvas[ipeak]->Update();


	for (int ibin = 1; ibin <= h1DX[ipeak]->GetNbinsX(); ibin++)
	{
		scale_Y[ibin-1] = h1DX[ipeak]->GetBinContent(ibin) / totalarea; //normalization
		scale_X[ibin-1] = h1DX[ipeak]->GetBinCenter(ibin);
		if(ibin<h1DX[ipeak]->GetNbinsX()) outfile<<scale_X[ibin-1]<<"	"<<setprecision(14)<<scale_Y[ibin-1]<<endl;
		if(ibin==h1DX[ipeak]->GetNbinsX()) outfile<<scale_X[ibin-1]<<"	"<<setprecision(14)<<scale_Y[ibin-1];
	}
	g1DX[ipeak] = new TGraph(Nbinsx_Tau, scale_X, scale_Y);//TGraph *gr1=new TGraph(n,x,y);
	g1DX[ipeak]->SetTitle(g_name);
	g1DX[ipeak]->SetName(g_name);
	//g1DX[ipeak]->GetXaxis()->SetTitle("Lifetime (fs)");//轴名
	g1DX[ipeak]->GetXaxis()->SetTitle("Parameter #tau (fs)");//轴名
	g1DX[ipeak]->GetYaxis()->SetTitle("Likelihood");//轴名
	g1DX[ipeak]->GetXaxis()->CenterTitle();//居中
	g1DX[ipeak]->GetYaxis()->CenterTitle();//居中
	g1DX[ipeak]->GetXaxis()->SetLabelFont(132);//坐标字体
	g1DX[ipeak]->GetYaxis()->SetLabelFont(132);//坐标字体
	g1DX[ipeak]->GetXaxis()->SetTitleFont(132);//轴名字体
	g1DX[ipeak]->GetYaxis()->SetTitleFont(132);//轴名字体
	//g1DX[ipeak]->GetYaxis()->SetLabelSize(0.05);//坐标字号
	//g1DX[ipeak]->GetYaxis()->SetTitleSize(0.05);//轴名字号
	g1DX[ipeak]->GetXaxis()->SetTitleOffset(1.2);//轴名偏移
	g1DX[ipeak]->GetYaxis()->SetTitleOffset(1.3);//轴名偏移
	//g1DX[ipeak]->SetMarkerStyle(21);
	//g1DX[ipeak]->SetMarkerColor(1);

	gcanvas[ipeak]=new TCanvas(g_name,g_name,900,600);
	gcanvas[ipeak]->cd();
	gPad->SetRightMargin(0.03);

	TSpline3 *s3 = new TSpline3("s3",scale_X,scale_Y,Nbinsx_Tau);
	s3->SetLineColor(kRed);
	g1DX[ipeak]->SetLineWidth(2);
	g1DX[ipeak]->Draw("AL");
	s3->SetLineWidth(2);
	s3->Draw("l same");

	TPaveText *textpol6 = new TPaveText(0.50,0.75,0.97,0.90,"brNDC");//left, down, right, up
	textpol6->SetBorderSize(1);
	textpol6->SetFillColor(0);
	textpol6->SetTextAlign(12);//align = 10*HorizontalAlign + VerticalAlign, 12 means水平左对齐、垂直居中对齐
	textpol6->SetTextFont(132);//font = 10 * fontID + precision, 12+2,12 means Symbol; 13+2, 13 means News Times Roman
	sprintf(paraprint,"90%% UL=%.1f%s",hight90," fs"); // modify
	textpol6->AddText(paraprint);
	textpol6->Draw();

	gcanvas[ipeak]->Update();
	//h2DX->Draw("colz");
	g1DX[ipeak]->Write();
	fout->Write();
	//fout->Close();
}//Fill2D main