Ctrl + R can search history commands in ROOT

//在包含C++头文件时一般不用后缀。如果用户自己编写头文件，可以用.h为后缀。
//在VS里用Ctrl+F5执行cpp文件
g++ -o ExampleMacro ExampleMacro.C `root-config --cflags --libs`

pad1->RedrawAxis();//Redraw the axis on this pad.

TGraph* graph = new TGraphErrors(n, x, y, err_x, err_y);//error bars TGraph(n,x,y,ex,ey);
tree->Draw("south_e:lege_e>>hsouth_lege", "south_e>50&&south_e<2000&&lege_e>5&&lege_e<400", "colz"); // y:x
TH2D* hmsd_lege = new TH2D("hmsd_lege", "hmsd_lege", 2000, 0, 100, 6000, 0, 6000); // x,y


#chi^{2} = %.2f / %d = %.5f   E_{#gamma} = 1 keV

TString rootfilename("Si25ET16m.root");

TCutG *cutg = new TCutG("CUTG",6);
cutg->SetVarX("TOF");
cutg->SetVarY("DE1max");
cutg->SetTitle("Graph");
cutg->SetFillColor(1);
cutg->SetPoint(0,185.595,197.372);
cutg->SetPoint(1,184.935,192.492);
cutg->SetPoint(2,186.783,174.85);
cutg->SetPoint(3,188.037,186.862);
cutg->SetPoint(4,186.981,200.751);
cutg->SetPoint(5,185.595,197.372);
cutg->Draw("");
if (cutg->IsInside(TOF,DE1max))
{
	DETOFSivalid=1; SiDETOF++;
}

h304px->SetName("P26");//set Name of a histogram

void Tfit()
{
//★Random
//gRandom is a pointer to the current random number generator
Float_t phi;
phi=gRandom->Uniform(0,2*pi);
phi=gRandom->Gaus(1,0.01);//mu=1, sigma=0.01

TH1F *h1=new TH1F("h1","h1",300, -15,15);
for (int i=0;i<500000;i++){
float myvalue = gRandom->Gaus(0,3);
h1->Fill(myvalue);//Add one x-value at a time
}
h1->Draw();

TF1 *gaus = new TF1("gaus","gaus",-10,10);// predefined 1D Gaussian function
gaus->SetParameters(100,0,3); //amplitude, mean_x, sigma_x, mean_y, sigma_y
gaus->Draw();
double x,y;
TH1F *h1=new TH1F("h1","h1",300, -15,15);
for (int i=0;i<500000;i++){
x=gaus->GetRandom();//Return a random number following this function shape.
h1->Fill(x);//Add one x-value at a time
}
h1->Draw();

// 2D Gaussian beam spot
TF2 *xygaus = new TF2("xygaus","xygaus",-15,15,-15,15);// predefined 2D Gaussian function
xygaus->SetParameters(1,0,3,0,3); //amplitude, mean_x, sigma_x, mean_y, sigma_y
xygaus->SetNpx(200);
xygaus->SetNpy(200);
xygaus->Draw("colz");
double x,y;
TH2F *h2 = new TH2F("h2","h2",300,-15.,15.,300,-15.,15.);
for (int i=0;i<500000;i++){
xygaus->GetRandom2(x,y); // Return two random numbers to x and y following this function shape.
h2->Fill(x,y);//
one x-y-value at a time
}
h2->Draw("colz");
h2->FillRandom("xygaus");

// 2D Uniform beam spot
double x, y, r, theta;
double Rmax = 15;
double Diameter = 2*Rmax;
TH2F *h2 = new TH2F("h2","h2",90,-15.,15.,90,-15.,15.);
for (int i=0;i<500000;i++){
	r = sqrt(gRandom->Uniform(0,1))*Rmax;
	theta = gRandom->Uniform(0,2.*3.14159);
	x= r*cos(theta);
	y = r*sin(theta);
h2->Fill(x,y);//Add one x-y-value at a time
}
h2->Draw("surf2");

//Fill histogram following distribution in function SiDEC.
//The distribution contained in the function SiDEC (TF1) is integrated over the channel contents for the bin range of this histogram. It is normalized to 1.
//Fill h4 channel 1000000 random numbers are generated
TH1F *h4=new TH1F("Fill_Histogram_Randomfrom_Function","Fill_Histogram_Randomfrom_Function",1000, 0,1000);
TF1 *SiDEC=new TF1("SiDEC","[0]*exp(x/(-[1]/0.693147))+[2]",1,1000);//自定义拟合函数
TF1 *SiDEC=new TF1("SiDEC","[0]*exp(x/(-[1]/0.693147))+[2]", 0,1000);//自定义拟合函数
SiDEC->SetNpx(1000);
SiDEC->SetParNames("A","T","B");
SiDEC->SetParameters(500,100,0);//自定义的拟合函数必须赋初值
SiDEC->SetParLimits(2,0,20);
h4->FillRandom("SiDEC",1000000);
h4->Draw();
//Fill histogram5 following distribution in histogram4 //fill a histogram following the distribution in an existing histogram
TH1F *h5=new TH1F("Fill_Histogram_Randomfrom_Histogram","Fill_Histogram_Randomfrom_Histogram",1000, 0,1000);
h5->FillRandom(h4,1000000);
h5->Draw();
//get a random number distributed according the contents of a histogram or a function
TH1F *h6=new TH1F("GetRandomfrom_Histogram","GetRandomfrom_Histogram",1000, 0,1000);
TH1F *h7=new TH1F("GetRandomfrom_Function","GetRandomfrom_Function",1000, 0,1000);
TH1F *h8=new TH1F("GetRandomfrom_Function_Range","GetRandomfrom_Function_Range",1000, 0,1000);
for(Int_t i=0;i<1000000;i++)//If type for loop in root, it just execute once, instead of 1000 times. It's weird.
{
	Double_t n=h5->GetRandom();//Return a random number distributed according the histogram bin contents.
	h6->Fill(n);
	n=SiDEC->GetRandom();//Return a random number following this function shape.
	h7->Fill(n);
	n=SiDEC->GetRandom(200,500);//Return a random number following this function shape in [xmin,xmax].
	h8->Fill(n);//Random numbers distributed according to a user defined function in a limited interval,
}
h6->Draw();
h7->Draw();
h8->Draw();
// or to a user defined histogram, can be generated in a very efficient way using
// TF1::GetRandom() or TH1::GetRandom(). TH1::GetRandom() method which can be used to get a random number distributed according the contents of a histogram.
// randomly fill an histogram using the contents of an existing TF1 function or another TH1 histogram (for all dimensions).

TMath::Log(2.71828)==1//Ln(x)
TMath::Log10(100)==2//Lg(x)

for (Int_t i=1;i<=h2->GetNbinsX();i++)
	h2->SetBinContent(i,value[i]);//Set the content of each bin /SetBinContent(i starts from 1, bincontent);
Double_t val[i] = h2->GetBinContent(i);//Retrieve the content of each bin

hChisquare2D->Draw("cont4");//surf2 for 3D, cont4 for 2D

//★TFit
total=new TF1("total","[0]*x+[1]+sqrt(3.141592654/2)*[2]/[3]*[4]*exp(0.5*([4]*[4]/([3]*[3]))+(x-[5])/[3])*ROOT::Math::erfc(1/sqrt(2)*([4]/[3]+(x-[5])/[4]))",3395,3575);
total->SetParameters(-0.27,1425.2,1482.8,4.21,6.611,3487.2);
total->Draw();
total->SetNpx(2000);
total->DrawDerivative();
total->Derivative(3483.789815422,0,0.001);//Returns the first derivative of the function at point x. x,params=0, eps=0.001 should be fine
total->GetMaximum();//Returns the maximum value of the function.
total->GetMaximumX();//Returns the X value corresponding to the maximum value of the function.
total->GetMinimum();//Returns the minimum value of the function.
total->GetMinimumX();//Returns the X value corresponding to the minimum value of the function.

y=pol6->GetMinimum(xmin,xmax);//Returns the minimum value of the function on the (xmin, xmax) interval.
x=pol6->GetMinimumX(xmin,xmax);//Returns the X value corresponding to the minimum value of the function on the (xmin, xmax) interval.
xlow=pol6->GetX(y+1,pol6->GetXmin(),x);//Returns the X value corresponding to the function value fy for (xmin<x<xmax)
xhigh=pol6->GetX(y+1,x,pol6->GetXmax());//Returns the X value corresponding to the function value fy for (xmin<x<xmax)
//GetXmin is 0; GetXmax is 6000ish. useful for Chisquare minimization

peaky[ii]=histo[i]->GetMaximum();
peakx[ii]=histo[i]->GetBinCenter(histo[i]->GetMaximumBin())
h->GetBinError(bin);
if you have a bin number of i
h1->GetBinCenter(i) returns the real x value.
If you have an x value.
h1->FindBin(x) returns the bin number corresponding to x.

TFile *fin = new TFile("Si25ET16m.root");
TTree *T999 = (TTree*)fin->Get("T999");
gStyle->SetOptFit(1111);//gStyle->SetOptFit(kTRUE)效果相同;
gStyle->SetFitFormat("6.5g");
gStyle->SetOptStat(0);
gStyle->SetOptFit();

TMath::Prob(χ2,NDF);//calculate p-value from Chisquare
TMath::Prob(500,498);//=0.466

TH1F *TSi = new TH1F("TSi","TSi",1199,5,6600);
T999->Draw("T>>TSi","T>0");
TF1 *SiDEC=new TF1("SiDEC","[0]*exp(x/(-[1]/0.693147))+[2]",5,6600);//自定义拟合函数
SiDEC->SetParNames("A","T","B");
SiDEC->SetParameters(300,220,1);//自定义的拟合函数必须赋初值
SiDEC->SetParLimits(1,210,230);
TSi->Fit("SiDEC","","",5,4400);//specify a range in the Fit, recommended
//TSi->Fit("SiDEC","R");//restrict the fit to the range specified in the TF1 constructor.
//void Fit(const char *fname, Option_t *option, Option_t *goption, Axis_t xxmin, Axis_t xxmax)
//option是L,W等，第二个goption跟draw的option一样，比如same

TF1 *SiDEC=new TF1("SiDEC","[0]*exp(x/(-[1]/0.693147))+[2]*x+[3]",5,6600);//自定义拟合函数
SiDEC->SetParNames("A","T","Ba","Bb");
SiDEC->SetParameters(1000,220,-0.1,5);//自定义的拟合函数必须赋初值
SiDEC->SetParLimits(1,210,230);

TH1F *hT = new TH1F("hT","hT",280,1,1401);
T999->Draw("TAB60>>Td","TAB60>0");
TF1 *SiDEC=new TF1("SiDEC","[0]*exp(x/(-[1]/0.693147))+[2]",1,6600);//自定义拟合函数
SiDEC->SetParNames("A","T","B");
SiDEC->SetParameters(300,220,1);//自定义的拟合函数必须赋初值
hT->Fit("SiDEC","","",1,1400);//specify a range in the Fit, recommended
SiDEC->SetParLimits(1,110,130);

SiDEC->SetLineColor(2);//fitting line color
SiDEC->SetLineWidth(2);//fitting line width
TSi->SetLineColor(4);
Double_t par[3];
SiDEC->GetParameters(par);
TF1 *Gausmod=new TF1("Gausmod","[0]*exp(x/(-[1]/0.693147))+[2]*x+[3]",5,6600);//自定义拟合函数
Gausmod->SetParameters(par);//自定义的拟合函数用上一个函数拟合后获得的参数作为初值
cout<<SiDEC->GetParError(1);//Obtaining the error of the 2nd parameter (T)
cout<<SiDEC->GetChisquare();
cout<<SiDEC->GetNDF();
cout<<SiDEC->GetProb();//This probability is not the “probability that	your fit is good.” If you did many fake experiments (draw many random samples of data points from the assumed distribution (your fit function)), this is the percentage of experiments that would give χ2 values ≥ to the one you got in this experiment.
cout<<SiDEC->GetChisquare()/SiDEC->GetNDF();
cout<<TSi->GetBinError(1000);//bin number=1000, the error is set equal to the sqrt(bin content)
//SiDEC->GetParameter(1);//Obtaining the value of the 2nd parameter (T)
//SiDEC->GetParameter("T");//Obtaining the value of the 2nd parameter (T)

TF1 *G=new TF1("G","1/(sqrt(2*3.14)*[0])*exp(-(x-[1])*(x-[1])/(2*[0]*[0]))",500,2500);
//★Gaussian Fit
ofstream outfile("fitg.dat",ios::out);//test multiproton 定义输出文件流对象outfile，以输出方式打开磁盘文件fitg.dat
TF1 *g=new TF1("g","gaus",200,3500);//高斯拟合，range可以大一些
h1->Fit("g","","",1980,2040);//调节合适的拟合range，h1是拟合的一维谱名
Double_t Sigma=g->GetParameter(2);//Obtaining the value of the 3rd parameter (Sigma)
Double_t Mean=g->GetParameter(1);//Obtaining the value of the 2nd parameter (Mean)
Double_t Sigmaerr=g->GetParError(2);//Obtaining the error of the 3rd parameter (Sigmaerr)
Double_t Chi=g->Chi2();//Obtaining the fit CHi2
outfile<<Mean<<" "<<Sigma*2.36/Mean*100<<"%"<<endl;
Double_t par[3];//GetParameters的数组得是double类型
g->GetParameters(par);//par[0]=constant; par[1]=mean; par[2]=sigma;

TF1 *g=new TF1("g","gausn",200,3500); //gausn: constant is the sum of the bin contents multiplied by the bin width in x

ofstream outfile("fitg.dat",ios::out);//test multiproton 定义输出文件流对象outfile，以输出方式打开磁盘文件fitg.dat
TH1F *E30060 = new TH1F("E30060","E30060",401,0,8020);
TF1 *g1=new TF1("g1","gaus",1820,2070);//高斯拟合，调节合适的拟合range
g1->SetParNames("Cons1","Mean1","Sig1");
E30060->Fit("g1","R");//restrict the fit to the range specified in the TF1 constructor
TF1 *g2=new TF1("g2","gaus",2045,2105);
g2->SetParNames("Cons2","Mean2","Sig2");
E30060->Fit("g2","R+");//
TF1 *g3=new TF1("g3","gaus",2105,2170);
g3->SetParNames("Cons3","Mean3","Sig3");
E30060->Fit("g3","R+");//
TF1 *g4=new TF1("g4","gaus",2160,2240);
g4->SetParNames("Cons4","Mean4","Sig4");
E30060->Fit("g4","R+");//
TF1 *g5=new TF1("g5","gaus",2240,2295);
g5->SetParNames("Cons5","Mean5","Sig5");
E30060->Fit("g5","R+");//
TF1 *gt=new TF1("gt","gaus(0)+gaus(3)+gaus(6)+gaus(9)+gaus(12)",1820,2295);
g1->SetParameters(45,1941,50);//填入上面各单高斯拟合找的参数
g2->SetParameters(20,2040,40);
g3->SetParameters(20,2134,40);
g4->SetParameters(12,2201,40);
g5->SetParameters(7,2263,40);
gt->SetParameters(2400,1530,11,1600,1535,14);//手动，各个高斯的参数做gt的初始化，双高斯拟合可不做初始化，In the more complicated case of the sum of 3 Gaussian functions, the initial values of parameters must be set. In this particular case, the initial values are taken from the result of the individual fits.
E30060->Fit("gt","R+");//+ means adding this new fitted function to the list of fitted functions (by default, the previous function is deleted and only the last one is kept)
//下面的参数设置更方便
double par[15];
g1->GetParameters(&par[0]);//GetParameters的数组得是double类型
g2->GetParameters(&par[3]);//Get到par数组，同时也做gt的Set数组，更自动
g3->GetParameters(&par[6]);
g4->GetParameters(&par[9]);
g5->GetParameters(&par[12]);
gt->SetParameters(par);//两个高斯的参数做gt的初始化
gt->SetParNames("Cons1","Mean1","Sig1","Cons2","Mean2","Sig2","Cons3","Mean3","Sig3","Cons4","Mean4","Sig4","Cons5","Mean5","Sig5");
E30060>Fit("gt","R+");//+ means adding this new fitted function to the list of fitted functions (by default, the previous function is deleted and only the last one is kept)
double Mean1=gt->GetParameter(1);//Obtaining the value of the 2nd parameter (Mean)
sprintf(paraprint,"Mean1=%4.2f",par[1]);//par数组还保持着g1、g2的参数
sprintf(paraprint,"Mean1=%4.2f",Mean1);//gt拟合参数，似乎寻峰不太准，拟合的Mean作峰位更准
//Draw
gStyle->SetOptFit(1);//显示Fit parameters
gStyle->SetOptStat("nemr")//仅对TH2F管用，对Tree->Draw()不管用，Tree->Draw()可以不加分号，返回值即为Entries
gStyle->SetOptStat("");//remove statistics.
hFit->Sumw2();
double n26 = hmsd26_e->Integral(47128, 47483);
double n31 = hmsd26_e->Integral(47128, 47483);
cout << n26 / n31;
hmsd26_e->Scale(n26 / n31);
hmsd26_e->Draw("same");


TH2(name, title, nbinsx, xlow, xup, nbinsy, ylow, yup);
TTree::Draw(varexp, selection, option = "", nentries = kMaxEntries, firstentry = 0); // y:x
T999->Draw("DE1:TOF", "TOF<199", "COLZ", 1e4, 283);//只画10000个事件,从第283个事件开始画

hFit->Scale(par[0]);//乘par[0]
h1->Scale(factor)// to normalize a histogram.

h3->Add(h1,h2,1,1)// to add two histograms together with the scaling factors in the command, for example, 1,1 means add, 1,-1 means subtract.
TH1F *h3=new TH1F("h3","h3",500,0,100);

// return a smart pointer to TFitResult using option “S”
TFitResultPtr r = h1->Fit("gaus","S");
double chi2 = r->Chi2(); // chi2 of fit
double fmin = r->MinFcnValue(); // minimum of fcn function
const double * par = r->GetParams(); // get fit parameters
const double * err = r->GetErrors(); // get fit errors
TMatrixDSym covMat = r->GetCovarianceMatrix();
TMatrixDSym corMat = r->GetCorrelationMatrix();
r->Print("V"); // full printout of result

h1->Integral()//返回积分值，loca ! i
h1->Integral(binx1,binx2)//include binx1 and binx2
TH1F *ESi0a = new TH1F("ESi0a","ESi0a",275,400,6000);
T999->Draw("Si0A>>ESi0a","Si0A>0");
TH1F *ESi0b = new TH1F("ESi0b","ESi0b",275,400,6000);
T999->Draw("Si0B>>ESi0b","Si0B>0");
T999->Draw("Si0B>>(275,400,6000)");//效果同上
ESi0a->Draw();
ESi0b->SetLineColor(2);
ESi0b->Draw("sames");//sames是statistics也re-drawn，same是只画图
ESi0b->Draw("C");//Draw a smooth curve through the histogram bins
ESi0b->Draw("E");//Draw the error bars
ESi0b->Draw("E0");//Draw the error bars including bins with 0 contents
ESi0b->Draw("E1");//Draw the error bars with perpendicular lines at the edges
ESi0b->Draw("HIST")//When a histogram has errors, it is visualized by default with error bars. To visualize it without errors use HIST
ESi0b->Draw("HIST C SAME");//以上命令可以叠加

T999->Draw("Si0A","Si0A>10&&Si0A<6000");
T999->Draw("Si0A:Si0B","Si0A>10&&Si0A<6000&&Si0B>10&&Si0B<6000");

char command1[100],command2[100];
sprintf(D_command1,"%s%d%s%d%s","adc[",i,"][",j,"]");
puts(D_command1);
sprintf(D_command2,"%s%d%s%d%s","adc[",i,"][",j,"]>0");
puts(D_command2);
tree->Draw(D_command1,D_command2);

T999->Draw("DE1:TOF","DE1>115&&DE1<255&&TOF>180&&TOF<215");
T999->SetMarkerColor(2);
T999->SetMarkerStyle(6);
T999->Draw("DE1:TOF","DE1>174&&DE1<220&&TOF>193&&TOF<198","same");//140
T999->Draw("DE1:TOF","DE1>175&&DE1<220&&TOF>193&&TOF<199&&CUTG1");//CUT,双击结束，右键setname，saveas可得到各点坐标值。
//gcut，如果选Ellipse，saveas可得到椭圆方程。
T999->Draw("DE1:TOF","DE1>175&&DE1<220&&TOF>193&&TOF<199&&CUTG1","same");//不加分号可得到返回值Entries
T999->Draw("DE1:TOF","same");//incorrect!
T999->Draw("DE1:TOF","","same");//correct!
T888->SetScanField(0);//dump numbers before <CR>, 0 means all numbers
T888->Scan("E","E>500&&E<6000&&T<54000"); >ca37.log//dump contents
T888->Scan("SCA[0]:SCA[1]",""); >si.log
T999->Scan("mul:FE:id", "", "", 1, 283)//只输出第283个event开始1个event的mul,FE,id
h1->SetBinErrorOption(TH1::kPoisson);
h1->Print("all"); >si.log//dump contents. This command doesn't work for new version ROOT
// use this instead:
void print_histogram() {
	TFile* fin = new TFile("F:/e21010/pxct/run0012-00_LEGe_55Fe_217min_cal.root");
	TH1F* histo = (TH1F*)fin->Get("hlege_e"); //Get spectrum
	histo->Print("all");
}
//usage: root -l -q print_histogram.C > h1.log

T888->GetEntries()
T888->Show(283)
for( Int_t i=1;i<h1->GetNbinsX();i++)
{
	Double_t Ch[i]=h1->GetBinContent(i);//retrieve the content of each bin
}

for (Int_t i=1;i<=hSega->GetNbinsX();i++){hSega->SetBinContent(i,1000);}//Set all bin contents to be 1000
for (Int_t i=1;i<=hSega->GetNbinsX();i++){hSega->AddBinContent(i,500);}//Add 500 to all bin content, we'll get 1500

for (Int_t i=1;i<=100;i++) {
	sum += h1->GetBinContent(i);
	hint1->SetBinContent(i,sum);
}
GetBinContent() [1/3]
Double_t TH1::GetBinContent 	(	Int_t 	bin	)
	const
	virtual 
	Return content of bin number bin. 
	Implemented in TH1C,S,F,D
	Convention for numbering bins
	For all histogram types: nbins, xlow, xup
	bin = 0; underflow bin
	bin = 1; first bin with low-edge xlow INCLUDED
	bin = nbins; last bin with upper-edge xup EXCLUDED
	bin = nbins+1; overflow bin

his2->Sumw2(kFALSE);
his2->SetBinErrorOption(TH1::kPoisson);//TH1::kNormal or TH1::kPoisson
y[0] = his2->GetBinContent(i); //hSega
if(y[0]>10)	dy[0] = his2->GetBinError(i);
if(y[0]<=10&&y[0]>0)	dy[0] = (his2->GetBinErrorLow(i)+his2->GetBinErrorUp(i))/2;
if(y[0]==0)	dy[0] = his2->GetBinErrorUp(i);//empty bin ErrorUp=1.8, ErrorLow=0, Error=0;

const float centroid=1000;
const float minrange=950;
const float maxrange=1050;
const float minrangeb=900;
const float maxrangeb=1100;
const int Nbins=1000;
const int Nbinsb=2000;
const int Binsperkev=10;
background=new TH1F("background","background",Nbinsb,minrangeb,maxrangeb);
bkgshort=new TH1F("bkgshort","bkgshort",Nbins,minrange,maxrange);
hSega=new TH1F("hSega","hSega",Nbins,minrange,maxrange);
hFit=new TH1F("hFit","hFit",Nbins,minrange,maxrange);

for (Int_t i=1;i<=hSega->GetNbinsX();i++){hSega->SetBinContent(i,i);}(注意bin的操作是=1到<=N)
hSega->GetNbinsX(); =1000;
hSega->Draw();
bin 0=0(不存在); bin 1=1; bin 2=2; bin i=i; ……, bin 1000=1000; bin 1001=0(不存在);

for (Int_t i=1;i<=background->GetNbinsX();i++){background->SetBinContent(i,i+100);}
background->GetNbinsX(); =2000;
background->Draw();
bin 0=0(不存在); bin 1=101; bin 2=102; bin i=100+i; ……, bin 2000=2100; bin 2001=0(不存在);

for (Int_t i=1;i<=bkgshort->GetNbinsX();i++){bkgshort->SetBinContent(i,background->GetBinContent((minrange-minrangeb)*Binsperkev+i));}（截取一段）
	bkgshort->GetNbinsX(); =1000;
(minrange-minrangeb)*Binsperkev;=500;
bkgshort->Draw();
bin 0=0(不存在); bin 1=601; bin 2=602; bin i=600+i; ……, bin 1000=1600; bin 1001=0(不存在);

hFit->Add(hSega,1);
hFit->Add(bkgshort,1);
hFit->Draw();
bin 0=0(不存在); bin 1=602; bin 2=604; bin i=600+2*i; ……, bin 1000=2600; bin 1001=0(不存在);

T888->Draw("Si0APC[8]","Si0APC[8]>0&&Si0APC[8]<6000&&DE1C<2&&TOFC>615");
T888->Draw("Si3AC[8]","Si3AC[8]>0&&Si3AC[8]<6000&&&&DE1C<2&&TOFC>615");
T888->Draw("Si3BC[8]","Si3BC[8]>0&&Si3BC[8]<6000&&DE1C<2&&TOFC>615");
T888->Draw("Si3AC[8]:Si3BC[8]","Si3AC[8]>0&&Si3AC[8]<6000&&Si3BC[8]>0&&Si3BC[8]<6000&&DE1C<2&&TOFC>615");

TH1F *ESi3a = new TH1F("ESi3a","ESi3a",190,500,6000);
T999->Draw("Si3A>>ESi3a","Si3A>0");
TH1F *ESi3b = new TH1F("ESi3b","ESi3b",275,500,6000);
T999->Draw("Si3B>>ESi3b","Si3B>0");
T999->Draw("Si3A:Si3B","Si3A>0&&Si3B>0&&Si3A<6000&&Si3B<6000");
ESi3a->Draw();
ESi3b->SetLineColor(2);
ESi3b->Draw("sames");
ESi3b->Draw("e1");//with error bars
ESi3b->DrawNormalized(”same”, 200);//ESi3b->DrawNormalized(”option_t”, norm);
h1->Rebin(5); //merges 5 bins in one in h1
TH1F *h2=h1->Rebin(5,"h2");// creat a new histogram h2 for merging h1

TH1F *ESi0A3Atot = new TH1F("ESi0A3Atot","ESi0A3Atot",230,300,7200);//ESi0+3A
T999->Draw("Si0A:Si3A","Si0A>10&&Si0A<6000&&Si3A>10&&Si3A<6000");
T999->Draw("Si0A+Si3A>>ESi0A3Atot","Si0A>10&&Si0A<6000&&Si3A>10&&Si3A<6000");
TH2F *ESi0ASi3A = new TH2F("ESi0ASi3A","ESi0ASi3A",190,300,6000,280,400,6000);//定义输出一二维谱
T999->Draw("Si0A:Si3A>>ESi0ASi3A","Si0A>10&&Si0A<6000&&Si3A>10&&Si3A<6000");

T888->Draw("DE1C:TOFC","DE1C>0&&DE1C<350&&TOFC>120&&TOFC<380");//分母
T888->Draw("DE1C:TOFC","DE1C>115&&DE1C<255&&TOFC>180&&TOFC<215");//5 nuclei
T888->SetMarkerColor(2);
T888->Draw("DE1C:TOFC","DE1C>176&&DE1C<220&&TOFC>193&&TOFC<198","same");//25Si
T888->Draw("Si0BHC[0]+Si0BHC[1]+Si0BHC[2]+Si0BHC[3]+Si0BHC[4]+Si0BHC[5]+Si0BHC[6]+Si0BHC[7]+Si0BHC[8]+Si0BHC[9]+Si0BHC[10]+Si0BHC[11]+Si0BHC[12]+Si0BHC[13]+Si0BHC[14]+Si0BHC[15]","Si0BHC[0]+Si0BHC[1]+Si0BHC[2]+Si0BHC[3]+Si0BHC[4]+Si0BHC[5]+Si0BHC[6]+Si0BHC[7]+Si0BHC[8]+Si0BHC[9]+Si0BHC[10]+Si0BHC[11]+Si0BHC[12]+Si0BHC[13]+Si0BHC[14]+Si0BHC[15]>1000&&DE1C>218&&DE1C<243&&TOFC>183&&TOFC<188"); //27S 150

GetBinContent() [1/3]
Double_t TH1::GetBinContent 
	(
	Int_t 
	bin
	)
	const
	virtual 
	Return content of bin number bin. 
	Implemented in TH1C,S,F,D
	Convention for numbering bins
	For all histogram types: nbins, xlow, xup
	bin = 0; underflow bin
	bin = 1; first bin with low-edge xlow INCLUDED
	bin = nbins; last bin with upper-edge xup EXCLUDED
	bin = nbins+1; overflow bin

T999->Draw("DE1:TOF","DE1>115&&DE1<255&&TOF>180&&TOF<215","COLZ");//5 nuclei
T999->SetMarkerColor(2);
T999->Draw("DE1:TOF","DE1>175&&DE1<220&&TOF>193&&TOF<199","COLZ");//25Si
T999->Draw("DE1:TOF","TOF<199","COLZ",1e4);//只画前10000个事件
T999->Draw("DE1:TOF","TOF<199","COLZ",1e4,283);//只画10000个事件,从第283个事件开始画
T999->Scan("mul:FE:id", "", "", 1, 283)//只输出第283个事件开始1个事件的mul,FE,id
T999->Draw("D300A[5]:D300B[6]","","");
int i=5; int j=i+1;
T999->Draw("D300A[5]:D300B[6]","D300A[5]>1000&&D300A[5]<1013","")
T999->Draw(Form("D300A[%d]:D300B[%d]",i,j),Form("D300A[%d]>1000&&D300A[%d]<1013",i,i),"")
T777->Draw(Form("DE2[%d]>>(600,10,3010)",i),"","")

TCanvas *myc_1 =new TCanvas("myc_1","myc_1",10,10,1200,600);//Canvas,画布
myc_1->Divide(1,2);
myc_1->cd(1);
myHisto_3->Draw();
myc_1->Update();//script里必须加，命令行里不必须，最好一个Canvas里cd、Draw然后都Update一下
TCanvas *canvascali[300];
sprintf(name,"cali_%s%d","h",i+1);
canvascali[i]=new TCanvas(name,name,600,480);//建立画布
canvascali[i]->cd();//进入画布

TGraph *graph[300];
graph[i]=new TGraph(peaknum,peakch,energy);//TGraph *gr1=new TGraph(n,x,y);建立曲线图、散点图
graph[i]->SetNameTitle("Stopping_power","Stopping_power");//Set graph name and title.
graph[i]->Draw("ap");
graph[i]->Eval(peakch);//a linear interpolation between the two points close to x is computed. If x is outside the graph range, a linear extrapolation is computed.
graph[i]->Fit("pol1");//pol1 can be used directly without TF1 constructor in CINT

sprintf(name,"cali_%s%d.png","h",i+1);
canvascali[i]->SetLogy();
canvascali[i]->SaveAs(name);//存png图
c1->Print(pdffilename);//存pdf
pad1->SetLogy(1);//setlogx, setlinx, setliny
pad1->SetLogy(0);//setlogx, setlinx, setliny
TH1D *histo=(TH1D*)fin->Get("h1");//直方图
canvaspeak->SaveAs("F:/e21010/pxct/alphagamma_241Am_PID_run70_71.png");

float Channel[72]={};
float Energy[72]={};
float E;
TGraph *g=new TGraph(72,Channel,Energy);//TGraph(n,x,y);
E=g->Eval(Channel);//Computes the Energy value of the given Channel value.

h1->Integral(200,300)
gPad->SetLogz(1) setlogz
gPad->SetLogz(0) setlinz
gPad->SetFrameLineWidth(2)
pad1->SetTopMargin(0.02);
TPad *pad1 = new TPad("pad1", "The pad 70% of the height",0.0,0.3,1.0,1.0);// Double_t xlow, Double_t ylow, Double_t xup, Double_t yup,
pad1->SetRightMargin(0.01);
pad1->SetLeftMargin(0.08);
pad2->SetRightMargin(0.01);
pad2->SetLeftMargin(0.08);
pad2->SetBottomMargin(0.2);
//Graphics
TH1F*h1=new TH1F(*hpG3001);//copy any histogram
h1->Add(hpG3001,hpG3002);//
h1->Add(h1,hpG3003);
h1->Add(h1,hpG3004);
h1->Draw();
TH1F*h2=new TH1F(*hpG600);//copy any histogram
h2->Add(h2,hpG601);
h2->Add(h2,hpG602);//
h2->Add(h2,hpG603);
h2->Add(h2,hpG604);
h2->Draw();
h1->Add(hG3001,hG3002,1,1);//1,1是两者相加，1,-1就是两者相减
TH1F*hD=new TH1F(*E60Bne);//copy any histogram
hD->Add(hD,E300Ane);//
hD->Draw();
TH1F *hT=(TH1F*)Td60->Clone("hT");//clone any histogram
TH1F*hT=new TH1F(*Td60);//copy any histogram
hT->Add(hT,Td300);//note that the bins should be same
hT->Draw();
// Create a new histogram based on the binning of an existing histogram
const Int_t numBins = hnorth_e->GetNbinsX();
const Double_t xMin = hnorth_e->GetXaxis()->GetXmin();
const Double_t xMax = hnorth_e->GetXaxis()->GetXmax();
TH1D* hnorth_e_tbkg = new TH1D("hnorth_e_tbkg", "hnorth_e_tbkg", numBins, xMin, xMax);

TSpectrum *s123 = new TSpectrum();
TH1F *G304bgb = s123->Background(G304bg, 20, "same");
TH1F *G304n = new TH1F(*G304bg);
G304n->Add(G304bg,G304bgb,1,-1);
G304n->Draw();

TH1F*hlineshape_high=new TH1F(*hlineshape_high);
hlineshape_high->SetFillColorAlpha(kAzure-4,1.0);
//hlineshape_high->SetFillStyle(3001);
hlineshape_high->SetLineColor(0);
hlineshape_high->Draw("samec");
TH1F*hlineshape_low=new TH1F(*hlineshape_low);
hlineshape_low->SetFillColor(10);
//hlineshape_low->SetFillStyle(3001);
hlineshape_low->SetLineColor(0);
hlineshape_low->Draw("samec");

In case multiple color filled histograms are drawn on the same pad, the fill area may hide the axis tick marks. One can force a redraw of the axis over all the histograms by calling:
gPad->RedrawAxis();

TF1 *SiDEC=new TF1("SiDEC","[0]*exp(x/(-[1]/0.693147))+[2]",1,6600);//自定义拟合函数
SiDEC->SetParNames("A","T","B");
SiDEC->SetParameters(300,220,1);//自定义的拟合函数必须赋初值
hT->Fit("SiDEC","","",1,1400);//specify a range in the Fit, recommended

T999->Draw("D60Ane>>h1(150,0,3000)","");
TSpectrum *s = new TSpectrum();
TH1F *hb = s->Background(h1, 12, "same");
TH1F *h2 = new TH1F(*h1);//copy any histogram
h2->Add(h1,hb,1,-1);
h1->Draw();
h2->Draw("same");

E60B->SetLineColor(2)
E60B->Draw("sames")
E300B->SetLineColor(3)
E300B->Draw("sames")
TAxis *axis = DETOF->GetXaxis();
DETOF->GetXaxis()->SetTitle("TOF");
DETOF->GetYaxis()->SetTitle("DE");
DETOF->GetXaxis()->CenterTitle();
DETOF->GetYaxis()->CenterTitle();
DETOF->GetYaxis()->SetTitle("#Delta#it{E} (MeV)");//ΔE,LaTex里的\都变成#
DETOF->SetGrid();//网格虚线

AC->GetYaxis()->SetRangeUser(0,AC->GetMaximum()*1.1);
AC->GetXaxis()->SetLabelSize(0.06);
AC->GetXaxis()->SetNdivisions(1002);
Posterior->Divide(5,5);
//Posterior->SetBorderMode(1);
gStyle->SetOptStat(0);
gStyle->SetOptTitle(0);
gStyle->SetPalette(51);
gStyle->SetLineWidth(2); // set the width of the axis lines, not frame lines
TColor::InvertPalette();
Posterior->Update();


TCanvas* canvaspeak = new TCanvas("LEGe", "LEGe", 1400, 540);//建立画布
canvaspeak->cd();//进入画布
canvaspeak->SetTopMargin(0.025);
canvaspeak->SetRightMargin(0.02);
canvaspeak->SetLeftMargin(0.08);
canvaspeak->SetBottomMargin(0.13);
hTauOmegaGamma->GetXaxis()->SetTitle("#it{#tau} (fs)");//斜体 Italic
hTauOmegaGamma->SetContour(100);
hlege_e->GetXaxis()->SetTitle("Energy (keV)");//轴名
hlege_e->GetYaxis()->SetTitle("Counts per 7 eV");//轴名
hlege_e->GetXaxis()->CenterTitle();//居中
hlege_e->GetYaxis()->CenterTitle();//居中
hlege_e->GetXaxis()->SetLabelFont(132);//坐标字体
hlege_e->GetYaxis()->SetLabelFont(132);//坐标字体
hlege_e->GetXaxis()->SetLabelSize(0.05);
hlege_e->GetYaxis()->SetLabelSize(0.05);
hlege_e->GetXaxis()->SetTitleFont(132);//轴名字体
hlege_e->GetYaxis()->SetTitleFont(132);//轴名字体
hlege_e->GetXaxis()->SetTitleOffset(1.0);//轴名偏移
hlege_e->GetYaxis()->SetTitleOffset(0.66);//轴名偏移
hlege_e->GetXaxis()->SetTitleSize(0.06);
hlege_e->GetYaxis()->SetTitleSize(0.06);
hlege_e->GetXaxis()->SetNdivisions(505);//n = n1 + 100*n2 + 10000*n3
hlege_e->GetYaxis()->SetNdivisions(505);//n = n1 + 100*n2 + 10000*n3
hlege_e->GetYaxis()->SetTickLength(0.02);
hlege_e->GetXaxis()->SetRangeUser(1, 435.7);
hlege_e->SetLineWidth(1);
hlege_e->Draw("hist");

tree->SetLineColor(kOrange + 7);
tree->SetLineColor(kGreen + 1);
tree->SetLineColor(kViolet);


axis->SetTitle("TOF");
axis->CenterTitle();
axis->SetTitleSize(0.05);//轴名字号
axis->SetTitleOffset(0.1);//轴名偏移
axis->SetTitleFont(22);//坐标字体, 22 means Times News Roman bold
axis->SetLabelFont(22);//坐标字体, 22 means Times News Roman bold
axis->SetLabelSize(0.04);//坐标字体
axis->SetRangeUser(193,198);//zoom the axis

TPaveText *text = new TPaveText(0.13,0.70,0.33,0.85,"brNDC");
text->SetBorderSize(1);
text->SetFillColor(0);
text->SetTextAlign(12);//align = 10*HorizontalAlign + VerticalAlign, 12 means水平左对齐、垂直居中对齐
text->SetTextFont(132);//font = 10 * fontID + precision, 12+2,12 means Symbol; 13+2, 13 means News Times Roman

colz图
palette右键->SetLabelFont(22);//色板字体, 22 means Times News Roman bold; 132 means News Times Roman
palette右键->SetLabelOffset(0.004);//色板坐标偏移量
palette右键-->SetLabelColor(1);
palette右键-->SetLabelFont(132);
palette右键-->SetLabelSize(0.04);

cDE1TOF->cd();
hDE1TOF->GetXaxis()->SetTitle("T1-T2");
hDE1TOF->GetYaxis()->SetTitle("DE1");
hDE1TOF->GetXaxis()->CenterTitle();
hDE1TOF->GetYaxis()->CenterTitle();
hDE1TOF->GetYaxis()->SetTitleOffset(1.5);
hDE1TOF->Draw();
TEllipse *Mgellipse = new TEllipse(x1Mg,y1Mg,aMg,bMg,0,360,0);
Mgellipse->SetFillStyle(0);
Mgellipse->SetLineColor(2);
Mgellipse->SetLineWidth(1);
Mgellipse->Draw("same");
TEllipse *Siellipse = new TEllipse(x1Si,y1Si,aSi,bSi,0,360,0);
Siellipse->SetFillStyle(0);
Siellipse->SetLineColor(4);
Siellipse->SetLineWidth(1);
Siellipse->Draw("same");
sprintf(beamlogname,"%s%04d%s","X:/RIBLL2015/Si22beam",irawroot,".png");
cDE1TOF->SaveAs(beamlogname);

T999->GetEntries("D300A>0")  Total number of events in the branch D300A, D300A>0
T999->Draw("D300A","D300A>0") return the total number of events in the branch D300A, D300A>0
T999->Draw("D300A","") return the total number of events in the branch D300A
T999->Draw("D300A") return nothing
outTree->Branch("SeGAenergy",&SeGAenergy,Form("SeGAenergy[%d]/D",SEGACHANNELS)); //energy branch

T777->MakeClass("RIBLL")//get RIBLL.h and RIBLL.C，类的定义以及Branch地址设定、分析框架都已经自动完成
}//void Tfit()

cout<<"RMS= "<<setiosflags(ios::fixed)<<setprecision(9)<<h1->GetRMS();//小数点后9位小数
//Read data from a file & Draw points to TGraph 
{
	TCanvas *c = new TCanvas("name", "canvas", 800, 500);
	gPad->SetLogx();
	double x;
	double w;
	int i = 0;
	string line;
	TGraph *g = new TGraph();
	ifstream file_input ("ziel.dat");
	while ( getline(file_input, line) )
	{
		stringstream(line) >> x >> w;
		g->SetPoint(i,x,w);
		i++;
	}
	g->Draw();
}


//To merge/combine files from run0079 up to run0102, including every file in that range, you would use the command:
hadd sum.root run00{79..99}*.root run0{100..102}*.root
//This command correctly accounts for the change from two leading zeros to one leading zero as the run numbers move from two to three digits.
rsync -av --ignore-existing /mnt/daqtesting/pxct/stagearea/run02{25..32}*.root /mnt/analysis/e21010/pxct_data/readout_PXCT_rootfiles/
 cp /mnt/daqtesting/pxct/stagearea/run02{25..32}*.root /mnt/analysis/e21010/pxct_data/readout_PXCT_rootfiles/