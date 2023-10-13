#include "TH1F.h"
#include "TF1.h"
#include <cmath>
#include <stdlib.h>
#include "TMinuit.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TStyle.h"
#include "TLine.h"
#define Nbpk 10 //Number of bins per keV
#define minrange 1275 //lower bound of fit
#define maxrange 1320 //upper bound of fit
using namespace std;

TFile *inFile,*inFile2,*outFile;
TH1F *hFit;

TH1F *p1,*p2,*p3,*p4,*p5,*p6,*hSega,*hSega1,*background;//>
TCanvas *c1;
int Nbins=(maxrange-minrange)*Nbpk;//make sure this is an integer

void readHists()
	//generate example histograms
{ 


	inFile = new TFile("/mnt/analysis/e14066/BEG/DopplerB/simulation.root"); //Rootfile with Simulation histograms
	inFile2 = new TFile("/mnt/analysis/e14066/BEG/DopplerShift/ExperimentData.root"); //Rootfile with Experimental Data Histogram

	//Get Histograms from Simulation
	p4= (TH1F*)inFile->Get("1034keVproton");
	p3= (TH1F*)inFile->Get("1878keVproton");
	p2= (TH1F*)inFile->Get("2546keVproton");
	p1= (TH1F*)inFile->Get("2769keVproton");//IAS
	p5= (TH1F*)inFile->Get("3194keVproton");
	p6= (TH1F*)inFile->Get("3714keVproton");





	hSega1= (TH1F*)inFile2->Get("hName"); //Get Gamma ray spectrum

	hSega=new TH1F("hSega","hSega",Nbins,minrange,maxrange);
	background=new TH1F("background","Background",Nbins,minrange,maxrange);
	for(int i=0;i<Nbins;i++){

		int power=1;//change this to whatever you want. 3 for high statistics
		TF1 *f1=new TF1("f1","[0]+x*[1]",minrange,1288);//range of left half of hist you want to fit the background to
		TF1 *f2=new TF1("f2","[0]+x*[1]",1309,maxrange);//range of right half of hist you want to fit the background to
		hSega1->Fit("f1","0qR");
		hSega1->Fit("f2","0qR");

		float a = (f1->GetParameter(0)+(minrange+(float)i/(float)Nbpk)*f1->GetParameter(1));//Lower BG
		float c = (f2->GetParameter(0)+(minrange+(float)i/(float)Nbpk)*f2->GetParameter(1));//Upper BG
		float f=((float)Nbins-(float)i)/(float)Nbins;
		float g =(float)i/(float)Nbins;
		float h =pow(f,power)+pow(g,power);

		float x = (a*pow(f,power)+c*pow(g,power))/h;
		int b = x;
		hSega->AddBinContent(i+1,hSega1->GetBinContent(minrange*Nbpk+i));
		background->AddBinContent(i+1,b);
	}

	p1->SetStats(0);
	p2->SetStats(0);
	p3->SetStats(0);
	p4->SetStats(0);

	background->SetStats(0);
	hSega->SetStats(0);

	p1->SetLineColor(4); 
	p2->SetLineColor(3); 
	p3->SetLineColor(6); 
	p4->SetLineColor(12); 

	background->SetLineColor(8);
	hSega->SetLineColor(1);


	p1->SetMarkerColor(4); 
	p2->SetMarkerColor(3); 
	p3->SetMarkerColor(6); 
	p4->SetMarkerColor(12); 

	background->SetMarkerColor(8);
	hSega->SetMarkerColor(1);

	int minim=0; //Only change this if you need to shift the histogram some amount of bins. Wouldn't recomment for bins larger than 0.1 keV
	if(minim>=0){
		for(int i=0;i<(p1->GetNbinsX()-minim);i++){
			p1->SetBinContent(p1->GetNbinsX()-i,p1->GetBinContent(p1->GetNbinsX()-i-minim));
			p2->SetBinContent(p1->GetNbinsX()-i,p2->GetBinContent(p1->GetNbinsX()-i-minim));
			p3->SetBinContent(p1->GetNbinsX()-i,p3->GetBinContent(p1->GetNbinsX()-i-minim));
			p4->SetBinContent(p1->GetNbinsX()-i,p4->GetBinContent(p1->GetNbinsX()-i-minim));
			p5->SetBinContent(p1->GetNbinsX()-i,p5->GetBinContent(p1->GetNbinsX()-i-minim));
			p6->SetBinContent(p1->GetNbinsX()-i,p6->GetBinContent(p1->GetNbinsX()-i-minim));
		}
	}
	else{
		for(int i=0;i<(p1->GetNbinsX()+minim);i++){
			p1->SetBinContent(i,p1->GetBinContent(i-minim));
			p2->SetBinContent(i,p2->GetBinContent(i-minim));
			p3->SetBinContent(i,p3->GetBinContent(i-minim));
			p4->SetBinContent(i,p4->GetBinContent(i-minim));
		}

	}
	hSega->SetLineWidth(2);



}

double CompareHists(TH1F* his1, TH1F* his2, int binMin = 1, int binMax = 1000000)
	//compare two histograms. returns chi-square
	//Two histograms must be with the same length
{
	double chiSquare = 0;                                   //chi square 
	int N = his1->GetNbinsX();                         //get number of bins in his1
	//   return 0;

	//if not identical to his2, exit with error.
	if (N!=his2->GetNbinsX())
	{
		printf("ERROR <CompareHists>: Different number of bins in histograms!%d\t%d\n",N,his2->GetNbinsX());
		exit(1);
	}

	if (N>binMax)
		N = binMax;

	double y[2],dy[2];                                         //bin information
	double delta,err;                                            //bin difference and relevant uncertainty
	for (int i=binMin;i<=N;i++)                         //loop over all bins   
	{
		//if(i<251||i>301){
		//get bin information
		y[0] = his1->GetBinContent(i); 
		dy[0] = his1->GetBinError(i);
		y[1] = his2->GetBinContent(i);
		// dy[1] = his2->GetBinError(i);
		dy[1] =0; //Set this to zero because simulation uncertainty double counts uncertainty
		delta = y[1]-y[0];
		err = sqrt(dy[0]*dy[0]+dy[1]*dy[1]);   //common uncertainty
		chiSquare += delta*delta/err/err;        
		//}      //add to chi square    
	} //for 
	return chiSquare;
}



void fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
	if (hFit) delete hFit;

	hFit = (TH1F*)p1->Clone("hFit");//2770 keV proton - IAS 
	hFit->Sumw2();
	hFit->Scale(par[0]);  
	//Set your different feeding intensities here
	hFit->Add(p2,0.1960*par[0]);//PIE-2547 from 6.266
	hFit->Add(p4,1.3730*par[0]);//PIE-1056 from 4.8 state
	hFit->Add(p5,0.0392*par[0]);//PIE-3194 from 6.92
	hFit->Add(p6,0.0196*par[0]);//PIE-3714 from 7.44

	hFit->Add(background,par[1]);//Make sure binning is 100 eV



	f = CompareHists(hFit,hSega,1,Nbins);//Nbpk

	cout<<"ChiSquare is " <<f<<endl;


}

void fit1()
{
	readHists();
	//   elasticH->Draw();
	//   int i;
	//   cin >> i;
	int nParams=2;
	TMinuit *gMin = new TMinuit(nParams);  //initialize TMinuit with a maximum of 4 params
	gMin->SetFCN(fcn);
	Double_t arglist[10];
	Int_t ierflg = 0;
	arglist[0] = 1;
	static Double_t vstart[2] = {0.07,1.0};//initial guess
	static Double_t step[2] = {0.0007,0.01};                //fitting step

	gMin->mnparm(0, "a1", vstart[0], step[0], 0,25,ierflg);  //parID, parName, initialGuess, step, lowLimit, highLimit, irrelevant
	gMin->mnparm(1, "a2", vstart[1], step[1], 0.9,1.5,ierflg);
	//gMin->mnparm(2, "a3", vstart[2], step[2], 0,25,ierflg);  

	arglist[0] = 1000;                             //number if fitting iterations
	arglist[1] = 1.;
	gMin->mnexcm("MIGRAD", arglist ,1,ierflg);  

	Double_t amin,edm,errdef;
	Int_t nvpar,nparx,icstat;
	gMin->mnstat(amin,edm,errdef,nvpar,nparx,icstat);

	double pars[nParams],errs[nParams];
	for (int i=0;i<nParams;i++)
	{
		gMin->GetParameter(i,pars[i],errs[i]);
		printf("%f %f\n",pars[i],errs[i]);    
	}


	//visualization

	hFit->SetLineColor(2);
	hFit->SetMarkerColor(2);
	hFit->SetLineWidth(2);




	p1->Scale(pars[0]);
	p2->Scale(0.1960*pars[0]);
	p4->Scale(1.3730*pars[0]);
	p5->Scale(0.0392*pars[0]);
	p6->Scale(0.0196*pars[0]);
	p1->Add(background,pars[1]);
	p2->Add(background,pars[1]);
	p4->Add(background,pars[1]);
	p5->Add(background,pars[1]);
	p6->Add(background,pars[1]);
	background->Scale(pars[1]);



	TCanvas* c1   = new TCanvas("c1","multipads",1100,900);  
	TPad *pad1 = new TPad("pad1", "The pad 80% of the height",0.0,0.35,1.0,1.0);
	TPad *pad2 = new TPad("pad2", "The pad 20% of the height",0.0,0.1,1.0,0.4);

	pad1->Draw();
	pad2->Draw();
	pad1->cd();
	gStyle->SetOptTitle(0);

	hSega->Draw();
	hSega->GetXaxis()->SetTitle("Energy (keV)");
	hSega->GetYaxis()->SetTitle("counts/0.1 keV");
	//hSega->GetXaxis()->SetTitle("Energy (keV)");
	hSega->GetYaxis()->SetTitle("counts/0.1 keV");
	hSega->GetXaxis()->SetLabelSize(0);
	//hSega->GetXaxis()->CenterTitle();
	hSega->GetYaxis()->CenterTitle();
	//hSega->GetXaxis()->SetTitleOffset(0.8);
	hSega->GetYaxis()->SetTitleOffset(0.9);
	//hSega->GetXaxis()->SetTitleSize(0.06);
	hSega->GetYaxis()->SetTitleSize(0.1);
	hSega->GetYaxis()->SetLabelSize(0.07);
	hSega->GetYaxis()->SetNdivisions(505);
	hSega->GetYaxis()->SetRangeUser(0,229000);
	hFit->Draw("samehist");
	p1->Draw("sameh");
	p2->Draw("sameh");
	p4->Draw("sameh");
	p5->Draw("sameh");
	p6->Draw("sameh");
	background->Draw("sameh");






	pad2->cd();
	int rebin=1;//1,2,3,5,6,9,10,15,18,25,30,45,50,75,90,150,225,450
	int Nbins=45*Nbpk/rebin;
	int minbin=1275;
	int maxbin=minbin+45;
	TH1F *Resid=new TH1F("Resid","Residual",Nbins,minbin,maxbin);
	double Resids[Nbins]={};
	double ResidsUncert[Nbins]={};
	double Xvals[Nbins]={};
	double XvalsUncert[Nbins]={};
	int count=0;
	for(int i=0;i<Nbins;i++){
		Xvals[i]=float(minbin)+1.0/Nbpk*(i*rebin)+0.5*rebin/Nbpk;
		// Xvals[i]=Xvals[i]*(1.+(float(maxbin-minbin)/2.-1.0/Nbpk*(i*rebin))/(float(maxbin-minbin)/2.)/(1.95*float(minbin)));
		//  cout<<Xvals[i]<<endl;
		ResidsUncert[i]=0;
		XvalsUncert[i]=0.0;
		for(int j=0;j<rebin;j++){
			Resid->AddBinContent(i+1,hFit->GetBinContent(i*rebin+j)-hSega->GetBinContent(i*rebin+j));
			Resids[i] += hFit->GetBinContent(i*rebin+j)-hSega->GetBinContent(i*rebin+j);
			ResidsUncert[i]+=hSega->GetBinContent(i*rebin+j);
		}
		ResidsUncert[i]=pow(ResidsUncert[i],0.5);
		if(abs(Resids[i])<ResidsUncert[i]) count++;
		//cout<< abs(Resids[i])<<"         "<<ResidsUncert[i]<<endl;
	}
	cout<<endl<<count<<" of "<<Nbins<<" error bars go through the fit."<<endl;
	//TGraphErrors *g;
	//g=new TGraphErrors(Nbins,Xvals,Resids,XvalsUncert,ResidsUncert); 
	TGraph *g;
	g=new TGraph(Nbins,Xvals,Resids); 
	g->GetXaxis()->SetTitle("Energy (keV)");
	Char_t axisName[256];
	sprintf(axisName,"#splitline{ Residual}{cts/%0.1f keV}",Nbins);
	g->GetYaxis()->SetTitle(axisName);
	g->GetXaxis()->CenterTitle();
	g->GetXaxis()->SetRangeUser(minbin,maxbin);
	g->GetYaxis()->CenterTitle();
	g->GetXaxis()->SetTitleOffset(1.0);
	g->GetYaxis()->SetTitleOffset(0.65);
	g->GetXaxis()->SetTitleSize(0.18);
	g->GetXaxis()->SetLabelOffset(.015);
	g->GetYaxis()->SetTitleSize(0.135);
	g->GetYaxis()->SetNdivisions(506);
	g->GetXaxis()->SetNdivisions(505);
	g->GetYaxis()->SetLabelSize(0.14);
	g->GetXaxis()->SetLabelSize(0.14);
	g->Draw("AP*same");
	TLine *T1=new TLine(minbin,0,maxbin,0);
	T1->Draw("R");

	outFile = new TFile("out.root","RECREATE");



	hSega->Write("hSega");
	hFit->Write("hFit");
	g->Write("Residual Plot");
	p1->Write("p_2770");
	p2->Write("p_2547");
	p3->Write("p_1878");
	p4->Write("p_1056");
	p5->Write("p_3194");
	p6->Write("p_3714");
	background->Write("background");
	delete outFile;
}

int main()
{
	fit1();
	return 0;
}
