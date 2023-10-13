#include<TChain.h>
#include "TH1.h"
#include "TMath.h"
#include "TF1.h"
#include "TLegend.h"
#include "TCanvas.h"

#define R 1350//T正式结果2700,B正式结果1350
#define halfd 447.9// half-life of daughter
#define halfd_delay 9e-7// half-life of daughter isomer
#define halfgd 179.4// half-life of granddaughter
#define halfgd_delay 486// half-life of granddaughter isomer
#define halfggd 100000// half-life of great-granddaughter
#define P 1//置为0即一级衰变
#define neutronfrac 0.0

void Mg20B()//Branching ratio of each peak, Decay spectrum
{
	//chain();
	TFile *fin = new TFile("V:/RIBLL2015/data20/Mg20calabcd_0347_0622cl750+500+t1200+p300+T1_1350T2_1350Q1repBr.root");//B正式结果
	//TFile *fin = new TFile("V:/RIBLL2015/data20/Mg20calabcd_0347_0622cl750+500+p500+g100+T1T2_2700T.root");//T正式结果
	ofstream outfile("C:/Si24/Si22peakcali/Tdfit.dat",ios::out);
	TH1F *T36 = new TH1F("T36","T36",9000,1,9001);//正式结果
	TH1F *T3_A = new TH1F("T3_A","T3_A",9000,1,9001); //与质子谱符合，T正式结果
	TH1F *T6_A = new TH1F("T6_A","T6_A",9000,1,9001); //与质子谱符合，T正式结果
	TH1F *T3_B = new TH1F("T3_B","T3_B",9000,1,9001); //与质子谱符合，T正式结果
	TH1F *T6_B = new TH1F("T6_B","T6_B",9000,1,9001); //与质子谱符合，T正式结果
 	//T999->Draw("T300>>T3_A","D300Ane+D300Bne>3110&&D300Ane+D300Bne<3600"); //与质子谱符合，T正式结果
	//T999->Draw("T60>>T6_A","D60Ane+D60Bne>3110&&D60Ane+D60Bne<3600"); //与质子谱符合，T正式结果
 	//T999->Draw("T300>>T3_B","D300Ane+D300Bne>1430&&D300Ane+D300Bne<1830"); //与质子谱符合，T正式结果
 	//T999->Draw("T60>>T6_B","D60Ane+D60Bne>1430&&D60Ane+D60Bne<1830"); //与质子谱符合，T正式结果

	//T999->Draw("T300>>T3_B","D300Ane+D300Bne>(712*2)&&D300Ane+D300Bne<(859*2)&&(QSD2[0]>80||QSD2[1]>80||QSD2[2]>100||QSD2[3]>80)"); //808
	//T999->Draw("T300>>T3_B","D300Ane+D300Bne>(1016*2)&&D300Ane+D300Bne<(1080*2)&&(QSD2[0]>80||QSD2[1]>80||QSD2[2]>100||QSD2[3]>80)"); //1048Q2/1071
	//T999->Draw("T300>>T3_B","D300Ane+D300Bne>(1366*2)&&D300Ane+D300Bne<(1422*2)&&(QSD2[0]>80||QSD2[1]>80||QSD2[2]>100||QSD2[3]>80)"); //1394Q2/1416
	//T999->Draw("T300>>T3_B","D300Ane+D300Bne>(1565*2)&&D300Ane+D300Bne<(1726*2)&&(QSD2[0]>80||QSD2[1]>80||QSD2[2]>100||QSD2[3]>80)"); //1679
	T999->Draw("T300>>T3_B","D300Ane+D300Bne>(1844*2)&&D300Ane+D300Bne<(1945*2)&&(QSD2[0]>80||QSD2[1]>80||QSD2[2]>100||QSD2[3]>80)"); //1899
	//T999->Draw("T300>>T3_B","D300Ane+D300Bne>(2200*2)&&D300Ane+D300Bne<(2258*2)&&(QSD2[0]>80||QSD2[1]>80||QSD2[2]>100||QSD2[3]>80)"); //2234Q2/2256
	//T999->Draw("T300>>T3_B","D300Ane+D300Bne>(2309*2)&&D300Ane+D300Bne<(2365*2)&&(QSD2[0]>80||QSD2[1]>80||QSD2[2]>100||QSD2[3]>80)"); //2337Q2/2359
	//T999->Draw("T300>>T3_B","D300Ane+D300Bne>(3983*2)&&D300Ane+D300Bne<(4131*2)&&(QSD2[0]>80||QSD2[1]>80||QSD2[2]>100||QSD2[3]>80)"); //4079
	//T999->Draw("T300>>T3_B","D300Ane+D300Bne>(4271*2)&&D300Ane+D300Bne<(4370*2)&&(QSD2[0]>80||QSD2[1]>80||QSD2[2]>100||QSD2[3]>80)"); //4333
	//T999->Draw("T300>>T3_B","D300Ane+D300Bne>(711*2)&&D300Ane+D300Bne<(904*2)"); //808
	//T999->Draw("T300>>T3_B","D300Ane+D300Bne>(1568*2)&&D300Ane+D300Bne<(1789*2)"); //1679
	//T999->Draw("T300>>T3_B","D300Ane+D300Bne>(1829*2)&&D300Ane+D300Bne<(1969*2)"); //1899
	//T999->Draw("T300>>T3_B","D300Ane+D300Bne>(2544*2)&&D300Ane+D300Bne<(2564*2)"); //2576
	//T999->Draw("T300>>T3_B","D300Ane+D300Bne>(3821*2)&&D300Ane+D300Bne<(3901*2)"); //3861
	//T999->Draw("T300>>T3_B","D300Ane+D300Bne>(3995*2)&&D300Ane+D300Bne<(4163*2)"); //4079
	//T999->Draw("T300>>T3_B","D300Ane+D300Bne>(4233*2)&&D300Ane+D300Bne<(4433*2)"); //4333
	
	//T999->Draw("T60>>T6_B","D60Ane+D60Bne>(716*2)&&D60Ane+D60Bne<(900*2)"); //808
	//T999->Draw("T60>>T6_B","D60Ane+D60Bne>(1560*2)&&D60Ane+D60Bne<(1768*2)"); //1664
	//T999->Draw("T60>>T6_B","D60Ane+D60Bne>(1817*2)&&D60Ane+D60Bne<(1977*2)"); //1897
	//T999->Draw("T60>>T6_B","2*D60Bne>(3792*2)&&2*D60Bne<(3872*2)"); //3832
	//T999->Draw("T60>>T6_B","2*D60Bne>(3792*2)&&2*D60Bne<(3872*2)"); //3832
	//T999->Draw("T60>>T6_B","2*D60Bne>(3999*2)&&2*D60Bne<(4131*2)"); //4065
	//T999->Draw("T60>>T6_B","2*D60Bne>(4248*2)&&2*D60Bne<(4400*2)"); //4324
	
	//T999->Draw("TAB60>>T6_B","D60Ane+D60Bne>(437*2)&&D60Ane+D60Bne<(477*2)"); //457
	//T999->Draw("TAB300>>T3_B","D300Ane+D300Bne>(437*2)&&D300Ane+D300Bne<(477*2)"); //457
// 	TSpectrum *s = new TSpectrum();
// 	TH1F *hb = s->Background(h6, 12, "same");
// 	TH1F *h2 = new TH1F(*h6);
// 	h2->Add(h6,hb,1,-1);
// 	h1->Draw();

	//T36->Add(T36,T3_A); //与质子谱符合，T正式结果
	//T36->Add(T36,T6_A); //与质子谱符合，T正式结果
	T36->Add(T36,T3_B); //与质子谱符合，T正式结果，B正式结果
	//T36->Add(T36,T6_B); //与质子谱符合，T正式结果，B正式结果部分
	//T36->Add(T36,Td300);
	//T36->Add(T36,Td60);
	T36->Rebin(10);//little influence for log likelihood method，T正式结果 10 ms bin
	//TAxis *xax = T36->GetXaxis();
	//TAxis *yay = T36->GetYaxis();
	//xax->SetNdivisions(5);
	//xax->SetLabelSize(.06);

	//decay->Draw("decay_time/1.0e8>>T36","deltaxy<1.0 && ADBK_T_cali[][]>-600 && ADBK_T_cali[][]<-200","E");
	//   decay->Project("T36","decay_time/1.0e8","deltaxy<1.0 && ADBK_T_cali[][]>-600 && ADBK_T_cali[][]<-200 ");

	double t;
	TF1 *SiDEC = new TF1("SiDEC",total,0,R,3); 
	SiDEC->SetLineColor(kRed); 
	SiDEC->SetNpx(5000);//函数画图时的取样点
	//SiDEC->SetLineWidth(4);
	SiDEC->SetParNames("N","T","B");//N: the number of decaying times time unit, here, the unit is 10 ms.
	SiDEC->SetParameters(10000,100,2);
	//SiDEC->SetParLimits(0,0,300000);//N
	SiDEC->SetParLimits(1,87,94);//T
	SiDEC->SetParLimits(2,0,5000);//B
	T36->SetLineColor(1);
	T36->Draw();
	//T36->GetYaxis()->SetRangeUser(0.01,1000);
	gPad->SetLogy(1);
	T36->Fit("SiDEC","L","",1,1350);//T正式结果2700，B正式结果1350
	//SiDEC->Draw("same");
	Double_t par[3]; //get the fit results to draw lines
	par[0]=SiDEC->GetParameter(0);//N
	par[1]=SiDEC->GetParameter(1);//T
	par[2]=SiDEC->GetParameter(2);//B
	double fiterror = SiDEC->GetParError(1);
	cout<<"Chisquare/NDF="<<SiDEC->GetChisquare()<<"/"<<SiDEC->GetNDF()<<"="<<SiDEC->GetChisquare()/SiDEC->GetNDF();
	TF1 *f1 = new TF1("f1",parent,0,R,3);
	f1->SetLineColor(4);

	TF1 *f2 = new TF1("f2",daughter,0,R,3);
	f2->SetLineColor(6);
	//f2->SetNpx(500);

	TF1 *f3 = new TF1("f3",granddaughter,0,R,3);
	f3->SetLineColor(5);

	TF1 *f4 = new TF1("f4",bkg,0,R,3);
	f4->SetLineColor(3);

	// cout<< "T1/2 = "<<693/par[1]<<" ("<<693/par[1]/par[1]*fiterror<<") ms"<< endl;
	//cout<< "T1/2 = "<<par[1]<<"("<<fiterror<<")"<<endl;
	//cout<< "X0 = "<<par[0]<<endl;

	f1->SetParameters(par);
	f1->Draw("same");
	f2->SetParameters(par);
	//f2->Draw("same");
	f3->SetParameters(par);
	//f3->Draw("same");
	f4->SetParameters(par);
	f4->Draw("same");
}


// void chain()
// {
// 	char str[100];
// 	TChain *decay = new TChain("decay","SiDEC"); 
// 	for(int i=6;i<27;i++){
// 		if(i!=8 && i!= 26){
// 			sprintf(str,"../decay_100s_%04d.root",i);
// 			decay->Add(str);
// 		}
// 	}
// }

Long64_t parent(Double_t *t, Double_t *par)
{
	return par[0]*0.693/par[1]*exp(-0.693/par[1]*t[0]);//[0]*exp(x/(-[1]/0.693147))+[2], [0]=A, par[0]=N, A=λN, par[1]=T, t[0]=x, [2]=bkg
	//return par[0]*exp(-0.693/par[1]*t[0]);//这样返回par[0]=A, A=λN
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

Long64_t total(Double_t *t, Double_t *par)
{
	//return parent(x,par)+bkg(x,par)+daugter(x,par)+granddaugter(x,par)+ggranddaugter(x,par); //+granddaugter(t,par)+daugter(t,par)
	return parent(t,par)+bkg(t,par);//+daughter(t,par);//+granddaughter(t,par)+daughterdelay(t,par); //+granddaugter(t,par)+daugter(t,par)
}