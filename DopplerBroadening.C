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

#define PI 3.141592653589793

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

void DopplerBroadening()//Doppler broadening simulation for 29Na
{
	time_t start,tim;
	struct tm *at;
	char now[80];
	float speed;
	/******************* Change These Numbers **********************/
	float MassDaughter = 28;             // Mass of daughter in nucleons. Multiplied by 938.5 to get total mass
	float T_Lifetime = 100;                   // This is Excited State Mean Lifetime (fs), Mean Lifetime τ=T/ln(2); half-life T=τ*ln(2)
	float E0_gamma=2000;                  // Energy of Gamma-ray (keV)
	const int nBranches=1;                  // Number of proton Branches feeding daughter
	float E_CoM[nBranches] = {3};         //Center of Mass energy (MeV) in descending order
	int Binsperkev = 10;                        //(10 bins per keV is 0.1 keV bins) etc...
	long CountsPerBranch = 10000;         // Number of counts simulated per proton branch. Recommend 1e6 or more
	bool Flag_SP = true;                         //Stopping power
	bool Flag_Det = true;                      //Detector resolution
	bool Flag_Ang = true;                       //Tests random cos vs random between -1 and 1 for angular dist (false tests cos and true tests random X(-1,1)
	/***************************************************************/
	int binmin=E0_gamma-25,binmax=E0_gamma+25;
	int binnum = (binmax-binmin)*Binsperkev;
	float Ecm_neutron=0;
	float Ef_gamma;
	float decaytime;
	int NumberofSteps=30;	//Divide the entire lifetime into several steps
	float c = 0.299792458;//0.3 is speed of light in microns/femtosecond
	float range=0;
	float ERecoil=0;//MeV
	float beta = pow(2.*ERecoil/(MassDaughter*938.5),0.5);//v0/c
	float dx;
	float dE;
	int det=-5;
	Double_t Ef_gamma_det1,Ef_gamma_det2;
	char simurootname[100];
	char hname[300];
	int ic,ii;
	long i;
	ofstream outfile("D:/X/out/Dopplershow.dat",ios::out);
	TCanvas *canvas[300];
	TH1F *hGammaspec1[nBranches];
	TH1F *hGammaspec2[nBranches];
	for(ii=0;ii<nBranches;ii++)
	{
		sprintf(hname,"Ecm%dkeVneutronf1",int(E_CoM[ii]*1000));
		hGammaspec1[ii]=new TH1F(hname,hname,binnum,binmin,binmax);
		sprintf(hname,"Ecm%dkeVneutronf2",int(E_CoM[ii]*1000));
		hGammaspec2[ii]=new TH1F(hname,hname,binnum,binmin,binmax);
	}
	TH1F *hDecaytime_distribution=new TH1F("Decaytime","Decaytime",200,0,20*T_Lifetime);
	TH1F *hv_c=new TH1F("beta","beta",1000,0,0.01);
	TH1F *hTotalGamma=new TH1F("TotalGamma_AllneutronFeedings","TotalGamma_AllneutronFeedings",binnum,binmin,binmax);
	TH1F *hAngulardist=new TH1F("Angular_Dist","Angular_Dist",100,-1,1);
	sprintf(simurootname,"%s","D:/X/out/Sim_script_Output_sqrt.root");
	TFile *fout=new TFile(simurootname,"RECREATE");//输出文件。It's better to define histograms and then define fout, in case of draw bugs.
	long TotalCounts=CountsPerBranch*nBranches;
	
	/*Probability Density Function*/
	TF1 *f1=new TF1("f1","[0]*[3]/[1]*1.2533*TMath::Exp(0.5*([3]*[3]/[1]/[1])+(x-[2])/[1])*TMath::Erfc(1/1.414*([3]/[1]+(x-[2])/[3]))",binmin,binmax);
	f1->SetNpx(500);//Set the number of points used to draw the function.
	/*Probability Density Function Integral*/
	TF1 *f2=new TF1("f2","[0]/2*(TMath::Exp((0.5*[3]*[3]-[1]*[2]+[1]*x)/([1]*[1]))*TMath::Erfc(1/1.414*([3]/[1]+(x-[2])/[3]))-TMath::Erfc((x-[2])/(1.414*[3])))+1",binmin,binmax);
	f2->SetNpx(500);//Set the number of points used to draw the function.

	float sigSQ0[13] = {0.0145,            0.0161,  0.0179,  0.0153,  0.0179,  0.0164,  0.0139,  0.0146,  0.0144,  0.0147,  0.0164,                0.0135,  0.0137};
	float sigSQ1[13] = {0.5484,            0.5441,  0.57,      0.5789,  0.5481,  0.4199,  0.5961,  0.4957,  0.7607,  0.5544,  0.4591,                0.7276,  0.8652};
	float sSq;
	//σ

	float integ[13] = {
		0.0822727070,    0.162619677,      0.23000568,        0.291370840,      0.371314436,      0.437634883,      0.530582334,      
		0.619391816,      0.696641487,      0.768618328,      0.842208471,      0.918817925,      1};//efficiency

	float feed[nBranches];
	for(ii=0;ii<nBranches;ii++)		feed[ii]=(ii+1)*1/(float)nBranches;//if nBranches = 4, feed[0,1,2,3]=0.25,0.5,0.75,1;
	float SP_Energy[72] = {0.01,        0.011,    0.012,    0.013,    0.014,    0.015,    0.016,    0.017,    0.018,    0.02,      0.0225,  0.025,                0.0275,  0.03,      0.0325,  0.035,    0.0375,  0.04,      0.045,    0.05,      0.055,    0.06,      0.065,    0.07,      0.08,      0.09,                0.1,         0.11,      0.12,      0.13,      0.14,      0.15,      0.16,      0.17,      0.18,      0.2,         0.225,    0.25,      0.275,    0.3,                0.325,    0.35,      0.375,    0.4,         0.45,      0.5,         0.55,      0.6,         0.65,      0.7,         0.8,         0.9,         1,            1.1,                1.2,         1.3,         1.4,         1.5,         1.6,         1.7,         1.8,         2,            2.25,      2.5,         2.75,      3,            3.25,      3.5,                3.75,      4,            4.5,         5};//MeV
	//19Ne Stopping power
	float SP_dEdx[72]={273.32,         274.24,  274.92,  275.46,  275.88,  276.2,    276.52,  276.65,  276.9,  277.2,    277.6,    278,                278.5,    279,        279.7,    280.4,    281.2,    286.3,    295.4,    301.7,    306.1,    309.3,    311.6,    313.5,    316.2,                318.35,  320.5,    322.79,  325.42,  328.39,  331.73,  335.28,  339.11,  343.17,  347.44,  356.39,  368.09,  380.15,                392.44,  405.05,  417.79,  430.61,  443.66,  456.61,  482.57,  508.02,  532.76,  556.79,  579.97,  602.33,  644.81,                684.35,  721.63,  756.75,  790.18,  822.16,  852.86,  882.37,  910.86,  938.43,  965.16,  1016.46,               1075.44,                1130.593,             1181.886,             1228.284,             1270.764,             1310.31,               1346.91,               1379.555,                1436.951,             1484.457};//  keV/micron
	//28Mg Stopping power
	//float SP_dEdx[72] = {387.47,	388.6,	389.42,	389.74,	389.88,	389.65,	389.35,	388.89,	388.37,	386.98,	384.85,	382.69,	380.33,	378.1,	375.82,	373.71,	371.58,	369.64,	366,	362.7,	359.8,	360.7,	361,	360.6,	358.2,	355.2,	352.2,	349.6,	347.5,	345.9,	344.9,	344.5,	344.5,	345,	345.9,	348.9,	354.3,	361,	368.8,	377.1,	385.99,	395.17,	404.6,	414.22,	433.99,	454.37,	475.12,	496.29,	517.55,	538.81,	581.4,	623.46,	664.84,	705.46,	745.31,	784.5,	822.87,	860.57,	897.66,	934.11,	969.9,	1040.22,	1123.85,	1204.9,	1282.25,	1356.85,	1427.63,	1494.57,	1558.63,	1619.79,	1731.37,	1830.2};//  keV/micron
	ic=0;
	canvas[ic]=new TCanvas("dEdx","dEdx",800,600);//建立画布
	canvas[ic]->cd();//进入画布
	TGraph *g=new TGraph(72,SP_Energy,SP_dEdx);//TGraph(n,x,y);
	g->SetNameTitle("Stopping_power","Stopping_power");//Set graph name and title.
	g->GetXaxis()->SetTitle("SP_Energy/MeV");//轴名
	g->GetYaxis()->SetTitle("SP_dEdx (keV/um)");//轴名
	g->GetXaxis()->CenterTitle();//居中
	g->GetYaxis()->CenterTitle();//居中
	g->SetMarkerStyle(7);
	g->SetMarkerColor(1);
	g->Draw("AP");
	canvas[ic]->SaveAs("D:/X/out/Sim_28Mg_Output.png");//存图
	ic++;

	gRandom->SetSeed(0);
	//TRandom3* ran = new TRandom3(1e6*TStopwatch::GetCPUTime());
	TRandom3* ran = new TRandom3(61520154+49092137+200000);
	float rand2,rand3,rand4;
	
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
			Ecm_neutron = E_CoM[0];//Ecm_neutron= E_CoM[0,1,2,3] according to i
		}
		else
		{
			for(ii=0;ii<nBranches-1;ii++)
			{
				if(i>=(feed[ii]*TotalCounts)&&i<(feed[ii+1]*TotalCounts))
					Ecm_neutron = E_CoM[ii+1];
			}
		}
		
		ERecoil=Ecm_neutron/(MassDaughter+1);//MeV E0_recoil
		beta = pow(2.*ERecoil/(MassDaughter*938.5),0.5);//v0/c
		range=0;//range μm
		dx=0; dE=0;
		//(1) sample decay time, divided into 300 steps
		//(2) In every step: calculate the distance (dx) using the initial velocity, then calculate the energy loss (dE) using the stopping power from graph (E_recoil → dE/dx) and the distance (dx). Calculate the E_recoil after this time step.
		//(3) Until the decay occurs or the Recoil stops. After steps, calculate the final velocity using the final E_recoil.
		//(4) Sample an angle, then calculate the Dopper shift for γ energy using final E_recoil and angle.
		//(5) Choose a detector, sample its response function to obtain the measured γ energy.
		rand2 = ran->Rndm();//same as Uniform(0,1)
		decaytime = T_Lifetime*TMath::Log(1.0/rand2);//ln(1/0~1)
		//float decaytime = T_Lifetime/0.69314718*TMath::Log(1.0/rand2);//ln(1/0~1) if using half-life at the beginning
		//outfile<<decaytime<<"	"<<rand2<<endl;
		hDecaytime_distribution->Fill(decaytime);
		
		if (Flag_SP)
		{
			for(ii=0;ii<NumberofSteps;ii++)//before the decay occurs, how far is the recoil move? Using the v0 obtain the first distance dx.
			{//then using dx obtain the dE, calculating the E_recoil after this time step.
				dx = beta*c*decaytime/(float)NumberofSteps;//dx μm
				range+=dx;//range μm
				dE = 0.001*dx*g->Eval(ERecoil);//a linear interpolation between the two points close to x is computed. If x is outside the graph range, a linear extrapolation is computed.
				//dE MeV Get the stopping power for this Erecoil, combined with the distance, obtain the energy loss in this distance (MeV).
				if(dE<ERecoil)
					ERecoil-=dE;//ERecoil MeV
				else
					ERecoil=0;
				beta = pow(2.*ERecoil/(MassDaughter*938.5),0.5);//v/c
				//outfile<<i<<"	Timesteps=	"<<NumberofSteps<<"	decaytime=	"<<decaytime<<"	dx=	"<<dx<<"	range=	"<<range<<"		dE=	"<<dE<<"	ERecoil=	"<<ERecoil<<endl;
			}//until the decay occurs or the recoil stops, calculate the final velocity using the final E_recoil.
		}
		
		hv_c->Fill(beta);//velocity of recoiling Atom when decay occurs

		/*Calculate Final energy of gamma*/
		if(Flag_Ang)
		{
			rand3 = ran->Uniform(-1,1);//random angle between photon and recoil
		}
		else
		{
			rand3 = 0;
		}
		float top=1-rand3*beta;
		float bot=1+rand3*beta;
		hAngulardist->Fill(rand3);
		Ef_gamma=E0_gamma*sqrt(top/bot);
		rand4 = ran->Rndm();//same as Uniform(0,1)
		det=-5;
		if(rand4<integ[0]) {det=0;}
		else
			for(ii=0;ii<12;ii++)
			{
				if(rand4>=integ[ii]&&rand4<integ[ii+1])
				{
					det=ii+1;
						break;//Randomly choose a detector according to the efficiency.
				}
			}

		if(Flag_Det)
		{
			sSq = sigSQ0[det]*sqrt(Ef_gamma) + sigSQ1[det];
			//cout<<det<<"	"<<sSq<<endl;
			f1->SetParameter(0,1.);//N
			f1->SetParameter(1,0.7);//λ
			f1->SetParameter(2,Ef_gamma);//μ
			f1->SetParameter(3,sSq);//σ
			//Ef_gamma_det1 = f1->GetRandom();
			//outfile<<Ef_gamma<<"	"<<Ef_gamma_det1<<endl;
			Ef_gamma_det1 = f1->GetRandom(binmin,binmax);//Sample a detector, sample its response function to obtain the measured γ energy.
			//outfile<<Ef_gamma<<"	"<<Ef_gamma_det1<<endl;
			f2->SetParameter(0,1.);//N
			f2->SetParameter(1,0.7);//λ
			f2->SetParameter(2,Ef_gamma);//μ
			f2->SetParameter(3,sSq);//σ
			Ef_gamma_det2 = f2->GetX(ran->Rndm(), binmin,binmax);
			//cout<<Ef_gamma<<"	"<<Ef_gamma_det2<<endl;
			//sprintf(cname,"c1_%d",i);
			//canvas[ic]=new TCanvas(cname,cname,600,400);//建立画布
			//canvas[ic]->cd();//进入画布
			//f1->Draw();
			//canvas[ic++]->Update();//script里必须加，命令行里不必须，最好一个Canvas里cd、Draw然后都Update一下
			//sprintf(cname,"c2_%d",i);
			//canvas[ic]=new TCanvas(cname,cname,600,400);//建立画布
			//canvas[ic]->cd();//进入画布
			//f2->Draw();
			//canvas[ic++]->Update();//script里必须加，命令行里不必须，最好一个Canvas里cd、Draw然后都Update一下
		}
		else
		{
			Ef_gamma_det1 = Ef_gamma;
			Ef_gamma_det2 = Ef_gamma;
		}
		
		if(i<(feed[0]*TotalCounts))
		{
			hGammaspec1[0]->Fill(Ef_gamma_det1);
			hGammaspec2[0]->Fill(Ef_gamma_det2);
		}
		else
		{
			for(ii=0;ii<nBranches-1;ii++)
			{
				if(i>=(feed[ii]*TotalCounts)&&i<(feed[ii+1]*TotalCounts))
				{
					hGammaspec1[ii+1]->Fill(Ef_gamma_det1);
					hGammaspec2[ii+1]->Fill(Ef_gamma_det2);
				}
			}
		}
		hTotalGamma->Fill(Ef_gamma_det2);

		if(i==0)
		{
			time(&start);
			at=localtime(&start);
			strftime(now,79,"%Y-%m-%d %H:%M:%S",at);
			cout<<now<<" started."<<endl;
		}
		if(i%(TotalCounts/100)==0&&i!=0)
		{
			time(&tim);
			at=localtime(&tim);
			strftime(now,79,"%Y-%m-%d %H:%M:%S",at);
			cout<<now;
			speed=(float)i/(float)difftime(tim,start);
			printf(" %.1f%s%.0f%s%.0f%s\n",100*(float)i/(float)TotalCounts,"%. Speed is ",speed," e/s. Still need ",float(TotalCounts-i)/speed," s.");
		}
 	}//for(i=0;i<nentries;i++)
	
	canvas[ic]=new TCanvas("Individual_gamma1","Individual_gamma1",600,400);//建立画布
	canvas[ic++]->cd();//进入画布
	hGammaspec1[nBranches-1]->Draw();
	for (ii=0;ii<(nBranches-1);ii++)
	{
		hGammaspec1[ii]->SetLineColor((Color_t)(ii+1));
		hGammaspec1[ii]->Draw("same");
	}
	canvas[ic]=new TCanvas("Individual_gamma2","Individual_gamma2",600,400);//建立画布
	canvas[ic++]->cd();//进入画布
	hGammaspec2[nBranches-1]->Draw();
	for (ii=0;ii<(nBranches-1);ii++)
	{
		hGammaspec2[ii]->SetLineColor((Color_t)(ii+1));
		hGammaspec2[ii]->Draw("same");
	}
	
	canvas[ic]=new TCanvas("Decaytime_distribution","Decaytime_distribution",600,400);//建立画布
	canvas[ic++]->cd();//进入画布
	//hDecaytime_distribution->SetLineColor((Color_t)ic);
 	hDecaytime_distribution->Draw();
	//canvas[ic]->Update();
	canvas[ic]=new TCanvas("beta","beta",600,400);//建立画布
	canvas[ic]->cd();//进入画布
	canvas[ic++]->SetLogy();
 	hv_c->Draw();
	//canvas[ic]->Update();
	canvas[ic]=new TCanvas("Total_Gamma","Total_Gamma",600,400);//建立画布
	canvas[ic++]->cd();//进入画布
 	hTotalGamma->Draw();
	canvas[ic]=new TCanvas("Angular_dist","Angular_dist",600,400);//建立画布
	canvas[ic++]->cd();//进入画布
	hAngulardist->GetYaxis()->SetRangeUser(0,TotalCounts/50);
	hAngulardist->Draw();

	for (ii=0;ii<nBranches;ii++)
	{
		hGammaspec1[ii]->Write();
		hGammaspec2[ii]->Write();
	}
	hDecaytime_distribution->Write();
	hv_c->Write();
	hTotalGamma->Write();
	hAngulardist->Write();
	
	//fout->Write();
	fout->Close();
}
