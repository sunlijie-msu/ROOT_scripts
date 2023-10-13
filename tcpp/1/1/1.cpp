#include <iostream>
//#include <iomanip.h>
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
	int runstart,runstop;
	cout<<"input runstart:";
	cin>>runstart;
	cout<<"input runstop:";
	cin>>runstop;
	char resultname[80];
	sprintf(resultname,"%s%03d%s%03d%s","C:/Si24/Si22peakcali/test_",runstart,"_",runstop,".txt");
	ofstream outfile(resultname,ios::app);
	int jj=0;
	for(int jj=1;jj<7;jj++)
	{
		outfile<<'\n'<<jj;
		cout<<'\n'<<jj;
	}
}