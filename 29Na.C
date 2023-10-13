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
#include "TF1.h"
#include "TMath.h"
#include "TPaveText.h"
#include "TPad.h"
#include "string.h"
#include "stdlib.h"
#include "stdio.h"
#include "TLegend.h"
#include "TPaletteAxis.h"

#define PI 3.14159265

using namespace std;

// Step 1: Determine which variable we want to solve for i.e. [S(E) stopping power OR proton feedings]
/* Step 2:

IF determining feedings: How many states in 29Mg (N) can feed 28Mg states?
1.Create random number to determine how long the recoiled 28Mg will travel.
2.Stopping power is known -- Determine final COM energy of 28Mg from travel time in material(EJ-200 scintillator)
3.Emit photon in random direction and calculate resulting photon energy
4.Determine ChiSquare which minimized the feedings from N states.

IF determining S(E)
1.Create random number to determine how long the recoiled 19Ne will travel.
2.Final energy of 
3.Emit photon in random direction and calculate resulting photon energy

*/

void DopplerBroad_29Na()//Doppler broadening simulation for 29Na
{
	time_t start,tim;
	struct tm *at;
	char now[80];
	float speed;
	/******************* Change These Numbers **********************/
	float MassDaughter = 28;             // Mass of daughter in nucleons. Multiplied by 938.5 to get total mass
	float T_Lifetime = 300;                   // Excited State LIFETIME (fs)
	float E0_gamma=3085;                  // Energy of Gamma-ray (keV)
	const int nBranches=4;                  // Number of proton Branches feeding daughter
	float E_CoM[nBranches] = {1,2,3,4.0};         //Center of Mass energy (MeV)
	int Binsperkev = 10;                        //(10 bins per keV is 0.1 keV bins) etc...
	unsigned long CountsPerBranch = 1e3;         // Number of counts simulated per proton branch. Recommend 1e6 or more
	/***************************************************************/
	int binmin=E0_gamma-25,binmax=E0_gamma+25;
	int binnum = (binmax-binmin)*Binsperkev;
	char simurootname[80];
	char hname[300];
	int ic,ii;
	long i;
	sprintf(simurootname,"%s","D:/X/out/Sim_28Mg_Output.root");
	TFile *fout=new TFile(simurootname,"RECREATE");
	TCanvas *canvas[30];
	TH1F *hDecaytime_distribution=new TH1F("Decaytime","Decaytime",200,0,20*T_Lifetime);
	TH1F *hv_c=new TH1F("beta","beta",1000,0,0.01);
	TH1F *hTotalGamma=new TH1F("TotalGamma_AllneutronFeedings","TotalGamma_AllneutronFeedings",binnum,binmin,binmax);
	TH1F *hAngulardist=new TH1F("Angular_Dist","Angular_Dist",100,-1,1);
	TH1F *hProton[nBranches];
	for(ii=0;ii<nBranches;ii++)
	{
		sprintf(hname,"Ecm%dkeVneutron",int(E_CoM[ii]*1000));
		hProton[ii]=new TH1F(hname,hname,binnum,binmin,binmax);
	}

	long TotalCounts=CountsPerBranch*nBranches;
	
	/*Probability Density Function*/
	//TF1 *f1=new TF1("f1","[0]*[3]/[1]*1.2533*TMath::Exp(0.5*([3]*[3]/[1]/[1])+(x-[2])/[1])*TMath::Erfc(1/1.414*([3]/[1]+(x-[2])/[3]))",binmin,binmax);
	/*Probability Density Function Integral*/
	TF1 *f2=new TF1("f2","[0]/2*(TMath::Exp((0.5*[3]*[3]-[1]*[2]+[1]*x)/([1]*[1]))*TMath::Erfc(1/1.414*([3]/[1]+(x-[2])/[3]))-TMath::Erfc((x-[2])/(1.414*[3])))+1",binmin,binmax);
	f2->SetNpx(100);//函数画图时的取样点

	float sigSQ0[13] = {0.0145,            0.0161,  0.0179,  0.0153,  0.0179,  0.0164,  0.0139,  0.0146,  0.0144,  0.0147,  0.0164,                0.0135,  0.0137};
	float sigSQ1[13] = {0.5484,            0.5441,  0.57,      0.5789,  0.5481,  0.4199,  0.5961,  0.4957,  0.7607,  0.5544,  0.4591,                0.7276,  0.8652};
	//σ

	float integ[13] = {
		0.0822727070,    0.162619677,      0.23000568,        0.291370840,      0.371314436,      0.437634883,      0.530582334,      
		0.619391816,      0.696641487,      0.768618328,      0.842208471,      0.918817925,      1};//efficiency

	float feed[nBranches];
	for(int m=0;m<nBranches;m++){feed[m]=(m+1)*1.0/float(nBranches);}//if nBranches = 4, feed[0,1,2,3]=0.25,0.5,0.75,1;

	float SP_Energy[72] = {0.01,        0.011,    0.012,    0.013,    0.014,    0.015,    0.016,    0.017,    0.018,    0.02,      0.0225,  0.025,                0.0275,  0.03,      0.0325,  0.035,    0.0375,  0.04,      0.045,    0.05,      0.055,    0.06,      0.065,    0.07,      0.08,      0.09,                0.1,         0.11,      0.12,      0.13,      0.14,      0.15,      0.16,      0.17,      0.18,      0.2,         0.225,    0.25,      0.275,    0.3,                0.325,    0.35,      0.375,    0.4,         0.45,      0.5,         0.55,      0.6,         0.65,      0.7,         0.8,         0.9,         1,            1.1,                1.2,         1.3,         1.4,         1.5,         1.6,         1.7,         1.8,         2,            2.25,      2.5,         2.75,      3,            3.25,      3.5,                3.75,      4,            4.5,         5};//MeV
	//19Ne Stopping power
	//float SP_dEdx[72]={273.32,         274.24,  274.92,  275.46,  275.88,  276.2,    276.52,  276.65,  276.9,  277.2,    277.6,    278,                278.5,    279,        279.7,    280.4,    281.2,    286.3,    295.4,    301.7,    306.1,    309.3,    311.6,    313.5,    316.2,                318.35,  320.5,    322.79,  325.42,  328.39,  331.73,  335.28,  339.11,  343.17,  347.44,  356.39,  368.09,  380.15,                392.44,  405.05,  417.79,  430.61,  443.66,  456.61,  482.57,  508.02,  532.76,  556.79,  579.97,  602.33,  644.81,                684.35,  721.63,  756.75,  790.18,  822.16,  852.86,  882.37,  910.86,  938.43,  965.16,  1016.46,               1075.44,                1130.593,             1181.886,             1228.284,             1270.764,             1310.31,               1346.91,               1379.555,                1436.951,             1484.457};//  keV/micron
	//28Mg Stopping power
	float SP_dEdx[72] = {387.47,	388.6,	389.42,	389.74,	389.88,	389.65,	389.35,	388.89,	388.37,	386.98,	384.85,	382.69,	380.33,	378.1,	375.82,	373.71,	371.58,	369.64,	366,	362.7,	359.8,	360.7,	361,	360.6,	358.2,	355.2,	352.2,	349.6,	347.5,	345.9,	344.9,	344.5,	344.5,	345,	345.9,	348.9,	354.3,	361,	368.8,	377.1,	385.99,	395.17,	404.6,	414.22,	433.99,	454.37,	475.12,	496.29,	517.55,	538.81,	581.4,	623.46,	664.84,	705.46,	745.31,	784.5,	822.87,	860.57,	897.66,	934.11,	969.9,	1040.22,	1123.85,	1204.9,	1282.25,	1356.85,	1427.63,	1494.57,	1558.63,	1619.79,	1731.37,	1830.2};//  keV/micron
	ic=0;
	canvas[ic]=new TCanvas("Stopping_power","Stopping_power",640,480);//建立画布
	canvas[ic]->cd();//进入画布
	TGraph *g=new TGraph(72,SP_Energy,SP_dEdx);//TGraph(n,x,y);
	g->GetXaxis()->SetTitle("SP_Energy/MeV");//轴名
	g->GetYaxis()->SetTitle("SP_dEdx (keV/um)");//轴名
	g->GetXaxis()->CenterTitle();//居中
	g->GetYaxis()->CenterTitle();//居中
	g->SetMarkerStyle(21);
	g->SetMarkerColor(1);
	g->Draw("AP");
	canvas[ic]->SaveAs("D:/X/out/Sim_28Mg_Output.png");//存图
	ic++;
	float E0=0;
	float Ef_gamma;
	int TimeSteps=(int)T_Lifetime/300+1;//lifetime<300 fs,step=1; lifetime=300~600 fs, step=2; lifetime=600~900 fs, step=3, etc.
	//Divide the entire lifetime into several steps, short lifetime will have fewer steps, and long lifetime will need more steps, of course.
	float c = 0.2998;//0.3 is speed of light in microns/femtosecond
	int itimes=0;
	gRandom->SetSeed(0);
	//TRandom3* ran = new TRandom3(1e6*TStopwatch::GetCPUTime());
	TRandom3* ran = new TRandom3(61520154+49092137+200000);
	for(i=0;i<TotalCounts;i++)
	{
// 		if(TotalCounts>100)
// 		{
// 			if(i%(TotalCounts/100)==0)
// 			{
// 				cout<<j<<"%"<<" done"<<endl;j++;
// 			}
// 		}
		if(i<(feed[0]*TotalCounts))
		{ 
			E0 = E_CoM[0];//E0= E_CoM[0,1,2,3] according to i
		}
		else
		{
			for(int n=0;n<nBranches-1;n++)
			{
				if(i>(feed[n]*TotalCounts)&&i<(feed[n+1]*TotalCounts))
					E0 = E_CoM[n+1];
			}
		}
		
		//sample decay time, calculate the distance using the velocity, calculate the energy loss using the stopping power and the distance
		//get the residual recoil energy, calculate the Dopper shift for gamma energy.
		float rand2 = ran->Rndm();//same as Uniform(0,1)
		float decaytime = T_Lifetime*TMath::Log(1.0/rand2);//ln(1/0~1)
		hDecaytime_distribution->Fill(decaytime);
		float range=0;
		float ERecoil=E0/(MassDaughter+1);
		float beta = pow(2.*ERecoil/(MassDaughter*938.5),0.5);//v0/c
		for(ii=0;ii<TimeSteps;ii++)//before the decay occurs, how far is the recoil move? Using the v0 obtain the distance.
		{
			float dx = beta * c * decaytime/float(TimeSteps);
			range+=dx;
			float dE = 0.001 * dx * g->Eval(ERecoil);//Get the stopping power for this Erecoil, combined with the distance, obtain the energy loss (MeV).
			if(dE<ERecoil)
				ERecoil-=dE;
			else
				ERecoil=0;
			beta = pow(2.*ERecoil/(MassDaughter*938.5),0.5);//v/c
		}
		hv_c->Fill(beta);//Speed of recoiling Atom when decay occurs

		/*Calculate Final energy of gamma*/
		float rand3 = ran->Uniform(-1,1);//random angle between photon and recoil
		float top=1-rand3*beta;
		float bot=1+rand3*beta;
		hAngulardist->Fill(rand3);
		Ef_gamma=E0_gamma*pow(top/bot,0.5);
		float randx = ran->Rndm();//same as Uniform(0,1)
		int det=-5;
		if(randx<integ[0]) {det=0;}
		else
			for(int k=0;k<12;k++)
			{
				if(randx>=integ[k]&&randx<integ[k+1])
				{
					det=k+1;
						break;//Randomly choose a detector according to the efficiency.
				}
			}
		if (det==-5) cout<<randx<<endl;
		float sSq = sigSQ0[det]*sqrt(Ef_gamma) + sigSQ1[det];
		f2->SetParameter(0,1.);//N
		f2->SetParameter(1,0.7);//λ
		f2->SetParameter(2,Ef_gamma);//μ
		f2->SetParameter(3,sSq);//σ
		Double_t Energy_det = f2->GetX(ran->Rndm(), binmin,binmax);

		if(i<(feed[0]*TotalCounts))
		{ 
			hProton[0]->Fill(Energy_det);
		}
		else
		{
			for(int n=0;n<nBranches-1;n++)
			{
				if(i>(feed[n]*TotalCounts)&&i<(feed[n+1]*TotalCounts))
					hProton[n+1]->Fill(Energy_det);
			}
		}
		hTotalGamma->Fill(Energy_det);

		if(i==0)
		{
			time(&start);
			at=localtime(&start);
			strftime(now,79,"%Y-%m-%d %H:%M:%S",at);
			cout<<now<<" started."<<endl;
		}
		if(i%10000==0&&i!=0)
		{
			time(&tim);
			at=localtime(&tim);
			strftime(now,79,"%Y-%m-%d %H:%M:%S",at);
			cout<<now;
			speed=(float)i/(float)difftime(tim,start);
			printf(" %.1f%s%.0f%s%.0f%s\n",100*(float)i/(float)TotalCounts,"%. Speed is ",speed," e/s. Still need ",float(TotalCounts-i)/speed," s.");
		}
 	}//for(i=0;i<nentries;i++)
	
	canvas[ic]=new TCanvas("gamma","gamma",640,480);//建立画布
	//canvas[ic++]->cd();//进入画布
	hProton[nBranches-1]->Draw();
	for (ii=0;ii<(nBranches-1);ii++)
	{
		hProton[ii]->Draw("same");
	}

	canvas[ic]=new TCanvas("Decaytime_distribution","Decaytime_distribution",640,480);//建立画布
	//canvas[ic]->cd();//进入画布
	hDecaytime_distribution->Draw();
	ic++;
	canvas[ic]=new TCanvas("beta","beta",640,480);//建立画布
	//canvas[ic]->cd();//进入画布
	//canvas[ic++]->SetLogy();
	hv_c->Draw();

	//canvas[ic]=new TCanvas("TotalGamma","TotalGamma",640,480);//建立画布
	//canvas[ic]->cd();//进入画布
	hTotalGamma->Draw();
	ic++;
	//canvas[ic]=new TCanvas("Angulardist","Angulardist",640,480);//建立画布
	//canvas[ic]->cd();//进入画布
	hAngulardist->GetYaxis()->SetRangeUser(0,TotalCounts/50);
	hAngulardist->Draw();
	
	for(ii=0;ii<nBranches;ii++)
	{
		hProton[ii]->Write();
	}
	hDecaytime_distribution->Write();
	hv_c->Write();
	hTotalGamma->Write();
	hAngulardist->Write();
	fout->Close();
}
