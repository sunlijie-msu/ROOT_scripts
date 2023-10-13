#include <TH2.h>
#include <TCanvas.h>
#include <iostream>
#include <TFile.h>
#include <stdio.h>
#include <TROOT.h>
#include <TRandom3.h>
#include <TMath.h>
#include <TGraph.h>
#include <TF1.h>
#include <TMinuit.h>

#define PI 3.14159265

using namespace std;

// Step 1: Determine which variable we want to solve for i.e. [S(E) stopping power OR proton feedings]
/* Step 2:

IF determining feedings: How many states in 20Na (N) can feed 19Ne states?
1.Create random number to determine how long the recoiled 19Ne will travel.
2.Stopping power is known -- Determine final COM energy of 19Ne from travel time in material(EJ-200 scintillator)
3.Emit photon in random direction and calculate resulting photon energy
4.Determine ChiSquare which minimized the feedings from N states.

IF determining S(E)
1.Create random number to determine how long the recoiled 19Ne will travel.
2.Final energy of 
3.Emit photon in random direction and calculate resulting photon energy

*/

void Brentcos(){

	/******************* Change These Numbers **********************/
	float MassDaughter = 28;             // Mass of daughter in nucleons. Multiplied by 938.5 to get total mass
	float T_Lifetime = 100;                   // Excited State LIFETIME (fs)
	float E0_gamma=2000;                  // Energy of Gamma-ray (keV)
	const int nBranches=1;                  // Number of proton Branches feeding daughter
	float E_CoM[nBranches] = {3};     //Center of Mass energy (MeV)
	int BinsperKev = 10;                        //(10 bins per keV is 0.1 keV bins) etc...
	int CountsPerBranch = 1e4;         // Number of counts simulated per proton branch. Recommend 1e6 or more
	bool SP_on = true;                         //Stopping power
	bool Det_res = true;                      //Detector resolution
	bool CosVsX = false;                       //Tests random cos vs random between -1 and 1 for angular dist (false tests cos and true tests random X(-1,1)
	/***************************************************************/

	int nbins = 50*BinsperKev;
	int minbin=E0_gamma-25;
	int maxbin=E0_gamma+25;

	TH1F *h5=new TH1F("time","time",200,0,20*T_Lifetime);
	TH1F *h6=new TH1F("beta","beta",1000,0,0.01);
	TH1F *h7=new TH1F("SummedAllProtonFeedings","SummedAllProtonfeedings",nbins,minbin,maxbin);
	TH1F *h8=new TH1F("AngularDist","Angular Dist",100,-1,1);


	TH1F *hProton[nBranches];
	for(int i=0;i<nBranches;i++){
		Char_t hname[256];
		sprintf(hname,"%dkeVproton",int(E_CoM[i]*1000));
		hProton[i]=new TH1F(hname,hname,nbins,minbin,maxbin);
	}

	int nCounts=CountsPerBranch*nBranches;

	/*Probability Density Function*/
	//TF1 *f1=new TF1("f1","[0]*[3]/[1]*1.2533*TMath::Exp(0.5*([3]*[3]/[1]/[1])+(x-[2])/[1])*TMath::Erfc(1/1.414*([3]/[1]+(x-[2])/[3]))",minbin,maxbin);
	/*Probability Density Function Integral*/
	TF1 *f2=new TF1("f2","[0]/2*(TMath::Exp((0.5*[3]*[3]-[1]*[2]+[1]*x)/([1]*[1]))*TMath::Erfc(1/1.414*([3]/[1]+(x-[2])/[3]))-TMath::Erfc((x-[2])/(1.414*[3])))+1",minbin,maxbin);
	f2->SetNpx(100);


	float sigSQ0[13] = {0.0145,            0.0161,  0.0179,  0.0153,  0.0179,  0.0164,  0.0139,  0.0146,  0.0144,  0.0147,  0.0164,                0.0135,  0.0137};
	float sigSQ1[13] = {0.5484,            0.5441,  0.57,      0.5789,  0.5481,  0.4199,  0.5961,  0.4957,  0.7607,  0.5544,  0.4591,                0.7276,  0.8652};


	float integ[13] = {
		0.0822727070,    0.162619677,      0.23000568,        0.291370840,      0.371314436,      0.437634883,      0.530582334,      
		0.619391816,      0.696641487,      0.768618328,      0.842208471,      0.918817925,      1};


		TCanvas *c2=new TCanvas("c2");
		float Efeed1536[nBranches] = {};
		for(int m=0;m<nBranches;m++){Efeed1536[m]=(m+1)*1.0/float(nBranches);}

		float SP_Energy[72] = {0.01,        0.011,    0.012,    0.013,    0.014,    0.015,    0.016,    0.017,    0.018,    0.02,      0.0225,  0.025,                0.0275,  0.03,      0.0325,  0.035,    0.0375,  0.04,      0.045,    0.05,      0.055,    0.06,      0.065,    0.07,      0.08,      0.09,                0.1,         0.11,      0.12,      0.13,      0.14,      0.15,      0.16,      0.17,      0.18,      0.2,         0.225,    0.25,      0.275,    0.3,                0.325,    0.35,      0.375,    0.4,         0.45,      0.5,         0.55,      0.6,         0.65,      0.7,         0.8,         0.9,         1,            1.1,                1.2,         1.3,         1.4,         1.5,         1.6,         1.7,         1.8,         2,            2.25,      2.5,         2.75,      3,            3.25,      3.5,                3.75,      4,            4.5,         5};//MeV
		float SP_dEdx[72]={273.32,         274.24,  274.92,  275.46,  275.88,  276.2,    276.52,  276.65,  276.91,  277.2,    277.6,    278,                278.5,    279,        279.7,    280.4,    281.2,    286.3,    295.4,    301.7,    306.1,    309.3,    311.6,    313.5,    316.2,                318.35,  320.5,    322.79,  325.42,  328.39,  331.73,  335.28,  339.11,  343.17,  347.44,  356.39,  368.09,  380.15,                392.44,  405.05,  417.79,  430.61,  443.66,  456.61,  482.57,  508.02,  532.76,  556.79,  579.97,  602.33,  644.81,                684.35,  721.63,  756.75,  790.18,  822.16,  852.86,  882.37,  910.86,  938.43,  965.16,  1016.46,               1075.44,                1130.593,             1181.886,             1228.284,             1270.764,             1310.31,               1346.91,               1379.555,                1436.951,             1484.457};//  keV/micron
		float SPmult=1.0;
		TGraph *g=new TGraph(72,SP_Energy,SP_dEdx);
		g->Draw();
		float E0=0;
		float Ef_gamma;
		float SP_mult=SPmult*1.0e-3;
		int TimeSteps=30;
		float c = 0.2998;//0.3 is speed of light in microns/femtosecond
		int j=0;
		gRandom->SetSeed(0);
		//TRandom3* ran = new TRandom3(1e6*TStopwatch::GetCPUTime());
		TRandom3* ran = new TRandom3(61520154+49092137+200000);
		for(int i=0;i<nCounts;i++){
			if(nCounts>100){
				if(i%(nCounts/100)==0){cout<<j<<"%"<<" done"<<endl;j++;}}
			if(i<(Efeed1536[0]*nCounts)){ 
				E0 = E_CoM[0];
			}
			else{
				for(int n=0;n<nBranches-1;n++){
					if(i>(Efeed1536[n]*nCounts)&&i<(Efeed1536[n+1]*nCounts)) E0 = E_CoM[n+1];}
			}

			float rand2 = ran->Rndm();
			float time = T_Lifetime*TMath::Log(1.0/rand2);
			h5->Fill(time);
			float range=0;
			float ERecoil=E0/(MassDaughter+1);
			float beta = pow(2.*ERecoil/(MassDaughter*938.5),0.5);
			if(SP_on){
				for(int jj=0;jj<TimeSteps;jj++){
					float dx = beta * c * time/float(TimeSteps);
					range+=dx;
					float dE = SP_mult * dx * g->Eval(ERecoil);
					if(dE<ERecoil){ERecoil-=dE;}
					else {ERecoil=0;}
					beta = pow(2.*ERecoil/(MassDaughter*938.5),0.5);
				}
			}

			h6->Fill(beta);//Speed of recoiling Atom

			/*Calculate Final energy of photon*/
			float rand3;
			if(!CosVsX){
				rand3 = ran->Uniform(-1,1);//random angle between photon and recoil
			}
			else{
				rand3 = TMath::Cos(ran->Uniform(-PI,PI));
			}
			float top=1-rand3*beta;
			float bot=1+rand3*beta;
			h8->Fill(rand3);
			Ef_gamma=E0_gamma*pow(top/bot,0.5);
			float randx = ran->Rndm();
			int det=-5;
			if(randx<integ[0]) {det=0;}
			else{
				for(int k=0;k<12;k++){
					if(randx>integ[k]&&randx<integ[k+1]){
						det=k+1;
						break;
					}
				}
			}
			Double_t xx;
			if(Det_res){
				float sSq = sigSQ0[det]*sqrt(Ef_gamma) + sigSQ1[det];
				f2->SetParameter(0,1.);
				f2->SetParameter(1,0.7);
				f2->SetParameter(2,Ef_gamma);
				f2->SetParameter(3,sSq);
				xx = f2->GetX(ran->Rndm(), minbin,maxbin);
			}
			else{
				xx = Ef_gamma;
			}

			if(i<(Efeed1536[0]*nCounts)){ 
				hProton[0]->Fill(xx);
			}
			else{
				for(int n=0;n<nBranches-1;n++){
					if(i>(Efeed1536[n]*nCounts)&&i<(Efeed1536[n+1]*nCounts)) hProton[n+1]->Fill(xx);}
			}
			h7->Fill(xx);        

		}

		TCanvas *c1=new TCanvas("c1");
		hProton[nBranches-1]->Draw();
		for (int i=0;i<(nBranches-1);i++){
			hProton[i]->Draw("same");
		}
		TCanvas *c5=new TCanvas("c5");
		h5->Draw();
		TCanvas *c6=new TCanvas("c6");
		c6->SetLogy();
		h6->Draw();
		TCanvas *c7=new TCanvas("c7");
		h7->Draw();
		TCanvas *c8=new TCanvas("c8");
		h8->GetYaxis()->SetRangeUser(0,nCounts/50);
		h8->Draw();




		TFile *out=new TFile("Sim_29Mg_Output.root","RECREATE");

		for(int i=0;i<nBranches;i++){
			hProton[i]->Write();
		}
		h5->Write();
		h6->Write();
		h7->Write();
		h8->Write();


		out->Close();














}
