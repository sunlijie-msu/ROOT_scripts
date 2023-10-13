#include<TChain.h>
#include "TH1.h"
#include "TMath.h"
#include "TF1.h"
#include "TLegend.h"
#include "TCanvas.h"


#define R 100
#define halfd 200
#define halfd_delay 9e-7
#define halfgd 179.4
#define halfgd_delay 486
#define halfggd 100000
#define P 1.0

#define neutronfrac 0.0
void decay_fit()
{
 chain();
 TH1F *decayhist = new TH1F("decayhist","decay curve of Tb168",100,0.0,R);
TAxis *xax = decayhist->GetXaxis();
TAxis *yay = decayhist->GetYaxis();
xax->SetNdivisions(5);
xax->SetLabelSize(.06);

 decay->Draw("decay_time/1.0e8>>decayhist","deltaxy<1.0 && ADBK_T_cali[][]>-600 && ADBK_T_cali[][]<-200","E");
//   decay->Project("decayhist","decay_time/1.0e8","deltaxy<1.0 && ADBK_T_cali[][]>-600 && ADBK_T_cali[][]<-200 ");

   double t;
   TF1 *f = new TF1("f",total,0,R,3); 
   f->SetLineColor(kMagenta); 
   f->SetNpx(500);
   f->SetLineWidth(4);

  f->SetParameters(1000,10,400);
  f->SetParLimits(0,0,300000);
  f->SetParLimits(1,0,1000);
  f->SetParLimits(2,0,500000);
  
   decayhist->Fit("f","L","");
   decayhist->Draw("E");
   f->Draw("same");

 decayhist->GetYaxis()->SetRangeUser(0.01,1000);

TF1 *f1 = new TF1("f1",parent,0,R,3);
f1->SetLineColor(kRed);

TF1 *f2 = new TF1("f2",daughter,0,R,3);
f2->SetLineColor(kBlue);
//f2->SetNpx(500);


TF1 *f3 = new TF1("f3",granddaughter,0,R,3);
f3->SetLineColor(kYellow);


TF1 *f4 = new TF1("f4",bkg,0,R,3);
f4->SetLineColor(kGreen);
Double_t par[3];

par[0]=f->GetParameter(0);
par[1]=f->GetParameter(1);
par[2]=f->GetParameter(2);
 double wu = f->GetParError(1);
// cout<< "T1/2 = "<<693/par[1]<<" ("<<693/par[1]/par[1]*wu<<") ms"<< endl;
 cout<< "T1/2 = "<<par[1]<<"("<<wu<<")"<<endl;
 cout<< "X0 = "<<par[0]<<endl;


f1->SetParameters(par);
f1->Draw("same");
f2->SetParameters(par);
f2->Draw("same");
f3->SetParameters(par);
f3->Draw("same");
f4->SetParameters(par);
f4->Draw("same");


}


 void chain()
 {
   char str[100];
   TChain *decay = new TChain("decay","f"); 
   for(int i=6;i<27;i++){
   if(i!=8 && i!= 26){
     sprintf(str,"../decay_100s_%04d.root",i);
     decay->Add(str);
   }
  }
 }

	Long64_t parent(Double_t *t, Double_t *par)
	{
	return par[0]*0.693/par[1]*exp(-0.693/par[1]*t[0]);
	}
	
	Long64_t daughter(Double_t *t, Double_t *par)
	{
	
	return  P*par[0]*0.693/par[1]*(0.693/halfd/(0.693/halfd-0.693/par[1]))*(exp(-0.693/par[1]*t[0])-exp(-0.693/halfd*t[0]));
	}
	Long64_t daughterdelay(Double_t *t, Double_t *par)
	{
	//Double_t halfd=1.26, neutronfrac=0.0043;
	return  (1-P)*par[0]*0.693/par[1]*(0.693/halfd_delay/(0.693/halfd_delay-0.693/par[1]))*(exp(-0.693/par[1]*t[0])-exp(-0.693/halfd_delay*t[0]));
	}
	 Long64_t granddaughter(Double_t *t, Double_t *par)
	{
	
	  return P*par[0]*0.693/halfgd*(exp(-0.693/par[1]*t[0])*0.693/par[1]*0.693/halfd/((0.693/halfd-0.693/par[1])*(0.693/halfgd-0.693/par[1]))+exp(-0.693/halfd*t[0])*0.693/par[1]*0.693/halfd/((0.693/par[1]-0.693/halfd)*(0.693/halfgd-0.693/halfd))+exp(-0.693/halfgd*t[0])*0.693/par[1]*0.693/halfd/((0.693/par[1]-0.693/halfgd)*(0.693/halfd-0.693/halfgd)));
	}
	
	
	Long64_t bkg(Double_t *t, Double_t *par)
	{
	return par[2];
	} 
	
	
	Long64_t total(Double_t *t, Double_t *par){
	
	//return parent(x,par)+bkg(x,par)+daugter(x,par)+granddaugter(x,par)+ggranddaugter(x,par); //+granddaugter(t,par)+daugter(t,par)
	return parent(t,par)+bkg(t,par)+daughter(t,par)+granddaughter(t,par)+daughterdelay(t,par); //+granddaugter(t,par)+daugter(t,par)
	}


