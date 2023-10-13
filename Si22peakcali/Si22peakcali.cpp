#include <iostream>
#include <fstream>
#include <math.h>
#include <map>
#include <TROOT.h>
#include "TChain.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TEllipse.h"
#include "TString.h"
#include "TGraph.h"
#include "TSpectrum.h"
#include "TVirtualFitter.h"
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
#include <time.h>
using namespace std;
void main()
{
	time_t tim;
	struct tm *at;
	char now[80];
	int irawroot;
	char rawrootname[80];
	double slope[300];
	double intercept[300];
	int ID=2;
	int binmin=500,binmax=3900;
	double sigma=12,thersh=0.4;
	float gaplow=60.,gaphigh=60.;//��ֱ治ͬ����
	unsigned long i;
	int j,ii;
	char paraprint[30],h_name[80],histo_name[80],b_name[80],hfit_name[80],hcali_name[80];
	TCanvas *canvaspeak[300];
	TCanvas *canvascali[300];
	TH1D *histo[16];//TH1F peak search+gauss fit,creat histograms, 
	//TH1D *h300AL[16] cannot be TH1D *h300AL[ii]= new TH1D(h_name,h_name,4000,0,4000);in loop
	TF1 *g[10];//creat function
	TGraph *graph[300];//TGraph�̶�
	float *peakch;
	double par[10][3];//�������ͬ�Ķ�
	//float energy[3]={1320,1860,2030};
	//float energy[3]={5157,5486,5805};
	float energy[3]={797,1670,2700};
	int peaknum;
	int runstart,runstop;
	cout<<"input runstart:";
	cin>>runstart;
	cout<<"input runstop:";
	cin>>runstop;
	for(irawroot=runstart; irawroot<=runstop; irawroot++)
	{
		//if()continue;
		sprintf(rawrootname,"%s%d%s","V:/RIBLL2015/Si22_0347_0472Mg",irawroot,".root");
		//sprintf(calirootname,"%s%04d%s","X:/RIBLL2015/cali",irawroot,".root");
		TFile *fin = new TFile(rawrootname);
		TTree *T999 = (TTree*)fin->Get("T999");
		unsigned long nentries=T999->GetEntries();//���¼���
		cout<<rawrootname;
		cout<<"  Entries="<<nentries<<endl;
		Float_t D300ALch[16];//define the variables to hold the read values
		Float_t D300BLch[16];//define the variables to hold the read values
		Float_t D60ALch[16];//define the variables to hold the read values
		Float_t D60BLch[16];//define the variables to hold the read values
		T999->SetBranchAddress("D300ALch", D300ALch);//other channels can be fitted by rename D300ALch
		T999->SetBranchAddress("D300BLch", D300BLch);//other channels can be fitted by rename D300ALch
		T999->SetBranchAddress("D60ALch", D60ALch);//other channels can be fitted by rename D300ALch
		T999->SetBranchAddress("D60BLch", D60BLch);//other channels can be fitted by rename D300ALch
		cout<<"input h_name for cali:";
		cin>>h_name;
		for(ii=0;ii<16;ii++)
		{
			sprintf(histo_name,"%s%d",h_name,ii);
			histo[ii] = new TH1D(histo_name,histo_name,3900,100,4000);//name a histogram
		}//name sequence decide the TH1F sequence in rootfile, name histograms if use Fill(branch)
		for(i=0;i<nentries;i++)
		{
			memset(D300ALch,0,sizeof(D300ALch));
			memset(D300BLch,0,sizeof(D300BLch));
			memset(D60ALch,0,sizeof(D60ALch));
			memset(D60ALch,0,sizeof(D60BLch));
			for(ii=0;ii<16;ii++)
			{
				T999->GetEntry(i);//Read the entry i in the tree (the i th event in the tree)
				sprintf(b_name,"%s%s%d%s",h_name,"ch[",ii,"]");
				histo[ii]->Fill(b_name);//no need to Fill Branch, TH1F is enough for cali
			}
		}
		//��tree Fill TH1F���ߴ�root�ļ�Get TH1F�γ��µ�TH1F����Ѱ�����
		//for(ii=0;ii<16;ii++)
		//{
		//sprintf(h_name,"%s%d","h300ALch",ii);
		//h300AL[ii] = (TH1D*)fin->Get(h_name);
		//}//read histograms from a file, no need to name histograms, also cannot rebin

		// 	h300AL[0]=(TH1D*)fin->Get("D300ALch0");
		// 	h300AL[1]=(TH1D*)fin->Get("D300ALch1");
		// 	h300AL[2]=(TH1D*)fin->Get("D300ALch2");
		// 	h300AL[3]=(TH1D*)fin->Get("D300ALch3");
		// 	h300AL[4]=(TH1D*)fin->Get("D300ALch4");
		// 	h300AL[5]=(TH1D*)fin->Get("D300ALch5");
		// 	h300AL[6]=(TH1D*)fin->Get("D300ALch6");
		// 	h300AL[7]=(TH1D*)fin->Get("D300ALch7");
		// 	h300AL[8]=(TH1D*)fin->Get("D300ALch8");
		// 	h300AL[9]=(TH1D*)fin->Get("D300ALch9");
		// 	h300AL[10]=(TH1D*)fin->Get("D300ALch10");
		// 	h300AL[11]=(TH1D*)fin->Get("D300ALch11");
		// 	h300AL[12]=(TH1D*)fin->Get("D300ALch12");
		// 	h300AL[13]=(TH1D*)fin->Get("D300ALch13");
		// 	h300AL[14]=(TH1D*)fin->Get("D300ALch14");
		// 	h300AL[15]=(TH1D*)fin->Get("D300ALch15");

		ofstream outfile("V:/Si23peakcali/Sipeakcali.dat",ios::out);
		for(i=0;i<ID;i++)
		{
			sprintf(hfit_name,"peak_%s%d",h_name,i);
			histo[i]->GetXaxis()->SetRangeUser(binmin,binmax);//zoom the axis
			canvaspeak[i]=new TCanvas(hfit_name,hfit_name,600,480);//��������
			canvaspeak[i]->cd();//���뻭��
			histo[i]->SetTitle(hfit_name);//ͼ��
			histo[i]->GetXaxis()->SetTitle("Channel");//����
			histo[i]->GetYaxis()->SetTitle("Counts");//����
			histo[i]->GetXaxis()->CenterTitle();//����
			histo[i]->GetYaxis()->CenterTitle();//����
			histo[i]->GetXaxis()->SetLabelFont(132);//��������
			histo[i]->GetYaxis()->SetLabelFont(132);//��������
			histo[i]->GetXaxis()->SetTitleFont(132);//��������
			histo[i]->GetYaxis()->SetTitleFont(132);//��������
			histo[i]->GetXaxis()->SetTitleOffset(1.2);//����ƫ��
			histo[i]->GetYaxis()->SetTitleOffset(1.3);//����ƫ��
			histo[i]->Draw();
			TSpectrum *s = new TSpectrum();
			//TSpectrum contains advanced spectra processing functions for 1- and 2-dimensional background estimation, smoothing, deconvolution, peak search and fitting, and orthogonal transformations.
			peaknum = s->Search(histo[i],sigma,"",thersh);//��Ѱ�����
			// Int_t Search(const TH1* hist, Double_t sigma = 2, Option_t* option = "", Double_t threshold = 0.05)
			//hin: pointer to the histogram of source spectrum
			//sigma: sigma of searched peaks, for details we refer to manual
			//threshold: (default=0.05) peaks with amplitude less than threshold*highest_peak are discarded. 0 By default, the background is removed before deconvolution. Specify the option "nobackground" to not remove the background.
			//printf("Found %d candidate peaks to fit in DSSD%d.\n",peaknum,i);
			peakch = s->GetPositionX();//��Ѱ�����
			for (j=0;j<peaknum;j++)
			{
				outfile<<hfit_name<<" canvaspeak"<<j<<" "<<peakch[j]<<endl;//����ı��鿴
			}
			for(ii=0;ii<peaknum;ii++)
			{
				g[ii]=new TF1("g","gaus",peakch[ii]-gaplow,peakch[ii]+gaphigh);//��˹��ϣ����ں��ʵ����range
				histo[i]->Fit("g","R+");//restrict the fit to the range specified in the TF1 constructor
				g[ii]->GetParameters(par[ii]);//g[ii],pointer to the TF1, GetParameters���������double����
			}
			//TF1 *g2=new TF1("g2","gaus",peakch[1]-gaplow/2,peakch[1]+gaphigh/2);//��ͬ��С�ķ��в�ͬ�����range�����ʺ�д�Զ����ѭ��
			//gStyle->SetFitColor(3);
			//h300AL[i]->Fit("g2","R+");//
			//TF1 *gt=new TF1("gt","gaus(0)+gaus(3)",peakch[0]-gaplow,peakch[1]+gaphigh);
			//g2->GetParameters(&par[3]);//Get��par���飬ͬʱҲ��gt��Set���飬���Զ�
			//gt->SetParameters(par);//�Զ���������˹�Ĳ�����gt�ĳ�ʼ����In the more complicated case of the sum of 3 Gaussian functions, the initial values of parameters must be set. In this particular case, the initial values are taken from the result of the individual fits.
			//h300AL[i]->SetLineColor(1);
			//h300AL[i]->Fit("gt","R+");//+ means adding this new fitted function to the list of fitted functions (by default, the previous function is deleted and only the last one is kept)
			//double Sigma1=gt->GetParameter(2);//Obtaining the value of the 3rd parameter (Sigma)
			//double Mean1=gt->GetParameter(1);//Obtaining the value of the 2nd parameter (Mean)
			//double Sigma2=gt->GetParameter(5);//Obtaining the value of the 3rd parameter (Sigma)
			//double Mean2=gt->GetParameter(4);//Obtaining the value of the 2nd parameter (Mean)
			TPaveText *textgaus = new TPaveText(0.13,0.50,0.31,0.85,"brNDC");//�ӱ�ע
			textgaus->SetBorderSize(1);//�߿���
			textgaus->SetFillColor(0);//�����ɫ
			textgaus->SetTextAlign(12);//align = 10*HorizontalAlign + VerticalAlign, 12 meansˮƽ����롢��ֱ���ж���
			textgaus->SetTextFont(132);//font = 10 * fontID + precision, 12+2,12 means Symbol; 13+2, 13 means News Times Roman
			//text->SetTextColor(2);//�ı���ɫ
			for(ii=0;ii<peaknum;ii++)
			{
				sprintf(paraprint,"Mean%d=%.1f",ii,par[ii][1]);//par���黹������g1��g2�Ĳ���
				textgaus->AddText(paraprint);
				sprintf(paraprint,"Sigma%d=%.1f",ii,par[ii][2]);
				textgaus->AddText(paraprint);
				sprintf(paraprint,"Res%d=%.2f%%",ii,par[ii][2]/par[ii][1]*2.355*100);
				textgaus->AddText(paraprint);
			}
			//sprintf(paraprint,"Mean1=%4.2f",Mean1);//gt��ϲ������ƺ�Ѱ�岻̫׼����ϵ�Mean����λ��׼
			textgaus->Draw();
			sprintf(hfit_name,"peak_%s%d.png",h_name,i);
			canvaspeak[i]->SaveAs(hfit_name);//��ͼ
			for(ii=0;ii<peaknum;ii++)
			{
				outfile<<"h300AL"<<i<<" Mean"<<ii<<"="<<par[ii][1]<<endl;//����ı��鿴
			}
			sprintf(hcali_name,"cali_%s%d",h_name,i);
			canvascali[i]=new TCanvas(hcali_name,hcali_name,600,480);//��������
			canvascali[i]->cd();//���뻭��
			//canvascali[i]->SetGrid();//��ʾ����
			graph[i]=new TGraph(peaknum,peakch,energy);//TGraph *gr1=new TGraph(n,x,y);
			//float eenergy[2]={200,280};
			//float epeakch[2]={0,0};
			//graph[i]= new TGraphErrors(peaknum,peakch,energy,epeakch,eenergy);//��error bars
			graph[i]->SetTitle(hcali_name);
			graph[i]->GetXaxis()->SetTitle("Channel");//����
			graph[i]->GetYaxis()->SetTitle("Energy");//����
			graph[i]->GetXaxis()->CenterTitle();//����
			graph[i]->GetYaxis()->CenterTitle();//����
			graph[i]->GetXaxis()->SetLabelFont(132);//��������
			graph[i]->GetYaxis()->SetLabelFont(132);//��������
			graph[i]->GetXaxis()->SetTitleFont(132);//��������
			graph[i]->GetYaxis()->SetTitleFont(132);//��������
			//graph[i]->GetYaxis()->SetLabelSize(0.05);//�����ֺ�
			//graph[i]->GetYaxis()->SetTitleSize(0.05);//�����ֺ�
			graph[i]->GetXaxis()->SetTitleOffset(1.2);//����ƫ��
			graph[i]->GetYaxis()->SetTitleOffset(1.3);//����ƫ��
			graph[i]->SetMarkerStyle(21);
			graph[i]->SetMarkerColor(1);
			TF1 *pol1=new TF1("pol1","pol1",200,4000);//����ʽ��ϣ����ں��ʵ����range
			graph[i]->Fit("pol1");//pol1 can be used directly without TF1 constructor in CINT
			intercept[i]=pol1->GetParameter(0);
			slope[i]=pol1->GetParameter(1);
			graph[i]->Draw("AP"); //A-Axis around the graph,AP is suitable
			TPaveText *textpol1 = new TPaveText(0.13,0.70,0.33,0.85,"brNDC");
			textpol1->SetBorderSize(1);
			textpol1->SetFillColor(0);
			textpol1->SetTextAlign(12);//align = 10*HorizontalAlign + VerticalAlign, 12 meansˮƽ����롢��ֱ���ж���
			textpol1->SetTextFont(132);//font = 10 * fontID + precision, 12+2,12 means Symbol; 13+2, 13 means News Times Roman
			sprintf(paraprint,"slope=%.2f",slope[i]);
			textpol1->AddText(paraprint);
			sprintf(paraprint,"intercept=%.2f",intercept[i]);
			textpol1->AddText(paraprint);
			textpol1->Draw();
			sprintf(h_name,"cali_%s%d.png","h300AL",i+1);
			canvascali[i]->SaveAs(hcali_name);//��ͼ
			outfile<<"h300AL"<<i<<" slope="<<slope[i]<<endl;//����ı��鿴
			outfile<<"h300AL"<<i<<" intercept="<<intercept[i]<<endl;//����ı��鿴
		}//for (i=0;i<ID;i++)
	}//for(irawroot=runstart; irawroot<=runstop; irawroot++)
}//peakcali main