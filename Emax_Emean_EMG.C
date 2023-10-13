#include <iostream>
#include <fstream>
#include <math.h>
#include <time.h>
#include <map>
#include <TROOT.h>
#include <TSystem.h>
#include <TStyle.h>
#include <TCutG.h>
#include <TChain.h>
#include <TMinuit.h>
#include "TChain.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TEllipse.h"
#include "TString.h"
#include "TGraph.h"
#include "TSpectrum.h"
#include "TVirtualFitter.h"
#include "TRandom.h"
#include <TRandom3.h>
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TF1.h"
#include "TMath.h"
#include "TPaveText.h"
#include "TPad.h"
#include "string.h"
#include "stdlib.h"
#include "stdio.h"
#include "TLegend.h"
#include "TPaletteAxis.h"

#define PI 3.141592653589793

using namespace std;

void Emax_Emean_EMG()//Find out the energy difference between Emax and Emean(米) for the EMG function used by DSL simulation
{
	double E0_gamma[6]={1248, 2234, 3076, 4971, 5156, 4045};
	double range=50;
	TCanvas *canvas1[6];
	TCanvas *canvas2[6];
	for (int ii=5;ii<6;ii++)
	{
		TRandom3* ran = new TRandom3(61520154+49092137+200000);
		TH1F *hGammaspec1=new TH1F("hf1","hf1",2*range*1000,E0_gamma[ii]-range,E0_gamma[ii]+range);
		TH1F *hGammaspec2=new TH1F("hf2","hf2",2*range*1000,E0_gamma[ii]-range,E0_gamma[ii]+range);
		/*Probability Density Function*/
		TF1 *f1=new TF1("f1","[0]*[3]/[1]*1.253314137*TMath::Exp(0.5*([3]*[3]/[1]/[1])+(x-[2])/[1])*TMath::Erfc(1/1.41421356*([3]/[1]+(x-[2])/[3]))",0,8000);//Glassman_PRC2019
		f1->SetNpx(1000);//Set the number of points used to draw the function. [0]-N, [1]-而, [2]-米, [3]-考,
		/*normalized cumulative Probability Density Function*/
		TF1 *f2=new TF1("f2","[0]/2*(TMath::Exp((0.5*[3]*[3]-[1]*[2]+[1]*x)/([1]*[1]))*TMath::Erfc(1/1.41421356*([3]/[1]+(x-[2])/[3]))-TMath::Erfc((x-[2])/(1.41421356*[3])))+1",0,8000);//Glassman_PRC2019
		f2->SetNpx(1000);//Set the number of points used to draw the function.

		double Ef_gamma_det1,Ef_gamma_det2;
		double sig,tau,errsig,errtau;
		double sigp1=	0.00011508152;
		double sigp0=	1.142295;
		double taup1=	0.000823321;
		double taup0=	-0.354105;
		sig = sigp1*E0_gamma[ii] + sigp0;
		tau = taup1*E0_gamma[ii] + taup0;

		for (int i=0;i<1;i++)
		{
			f1->SetParameter(0,1.);//N
			f1->SetParameter(1,tau);//而
			f1->SetParameter(2,E0_gamma[ii]);//米
			f1->SetParameter(3,sig);//考
			Ef_gamma_det1 = f1->GetMaximumX();
			//Ef_gamma_det1 = f1->GetRandom(E0_gamma[ii]-range,E0_gamma[ii]+range);//slow
			cout<<setprecision(11)<<Ef_gamma_det1<<endl;

			f2->SetParameter(0,1.);//N
			f2->SetParameter(1,tau);//而
			f2->SetParameter(2,E0_gamma[ii]);//米
			f2->SetParameter(3,sig);//考
			Ef_gamma_det2 = f2->GetX(ran->Rndm(), E0_gamma[ii]-range,E0_gamma[ii]+range);//fast

			hGammaspec1->Fill(Ef_gamma_det1);
			hGammaspec2->Fill(Ef_gamma_det2);
		}

		gStyle->SetStatFormat("9.8g");
		canvas1[ii]=new TCanvas("Individual_gamma1","Individual_gamma1",800,600);
		canvas1[ii]->cd();
		hGammaspec1->Draw();
		canvas2[ii]=new TCanvas("Individual_gamma2","Individual_gamma2",800,600);
		canvas2[ii]->cd();
		hGammaspec2->Draw();
	}
}