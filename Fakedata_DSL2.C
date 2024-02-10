//1) Fill_histogram_random.C to generate Fakebackground.root

//Open a ROOT and run the following macro
//2) To generate a fake Doppler - shifted gamma peak

TFile* fin = TFile::Open("F:/out/G4_rootfiles_with_tree_Eg4156/S31_Gamma4156_Eg4155.84_Tau3.0_SP1.00_AC0.0_all.root");
TH1F* hfakepeak = new TH1F("hfakepeak", "hfakepeak", 8000, 0, 8000);
TH1F* hfakepeakfull = new TH1F("hfakepeakfull", "hfakepeakfull", 8000, 0, 8000);
tree->Draw("Clovere>>hfakepeakfull", "Clovere>4410&&Clovere<4524&&DSSD1e+DSSD2e>10000&&DSSD1e>100&&DSSD2e>100", "");
hfakepeakfull->GetXaxis()->SetRangeUser(4400, 4550);
hfakepeakfull->GetYaxis()->SetRangeUser(0, 36);
hfakepeakfull->SetLineColor(kRed);
hfakepeakfull->SetLineWidth(2);
//hfakepeakfull->Scale(1020.0 / 25318.0); // 3fs
hfakepeakfull->Scale(1020.0 / 31749.0); // 0fs
//int i = 0;
i += 1e5; cout << i << endl;
tree->Draw("Clovere", "Clovere>4410&&Clovere<4524&&DSSD1e+DSSD2e>10000&&DSSD1e>100&&DSSD2e>100", "sames", 1.7e5, i);
hfakepeakfull->GetYaxis()->SetRangeUser(0, 41);

gStyle->SetStatFormat("7.6g");
gPad->Update();

// tree->Draw("Clovere>>hfakepeak", "Clovere>0&&Clovere<8000&&DSSD1e+DSSD2e>10000&&DSSD1e>100&&DSSD2e>100", "", 1.70e5, 24e5); // manipulated for 
// tree->Draw("Clovere>>hfakepeak", "Clovere>0&&Clovere<8000&&DSSD1e+DSSD2e>10000&&DSSD1e>100&&DSSD2e>100", "", 1.70e5, 4e5);

TFile* fin = TFile::Open("F:/out/G4_rootfiles_with_tree_Eg4156/S31_Gamma4156_Eg4155.84_Tau3.0_SP1.00_AC0.0_all.root");
TH1F* hfakepeak = new TH1F("hfakepeak", "hfakepeak", 8000, 0, 8000);
tree->Draw("Clovere>>hfakepeak", "Clovere>0&&Clovere<8000&&DSSD1e+DSSD2e>10000&&DSSD1e>100&&DSSD2e>100");
hfakepeak->Scale(1020.0 / 31749.0); // 0fs
hfakepeak->Scale(510.0 / 31749.0); // 0fs
hfakepeak->Scale(1020.0 / 25318.0); // 3fs
hfakepeak->Scale(510.0 / 25318.0); // 3fs
hfakepeak->SetBinErrorOption(TH1::kPoisson);
TFile* fout = new TFile("F:/out/G4_rootfiles_with_tree_Eg4156/Fakepeak_S31_Gamma4156_Eg4155.84_Tau3.0_SP1.00_AC0.0_0.5k.root", "RECREATE");
hfakepeak->Write();
fout->Close();

//Get rid of G4 background below the peak
TFile* fin1 = TFile::Open("F:/out/G4_rootfiles_with_tree_Eg4156/Fakepeak_S31_Gamma4156_Eg4155.84_Tau3.0_SP1.00_AC0.0_0.5k.root");
TH1F* hfakepurepeak = new TH1F("hfakepurepeak", "hfakepurepeak", 8000, 0, 8000);

//TF1 *pol1=new TF1("pol1","[0]*x+[1]",0,8000);//G4 bkg
//hfakepeak->Fit("pol1","MLE","", 4280,4380);
//hfakepeak->Fit("pol1","MLE","", 4280,4380);
//hfakepeak->Fit("pol1","MLE","", 4280,4380);
//double a = pol1->GetParameter(0);
//double b = pol1->GetParameter(1);

//for (int xbin=1;xbin<=8000;xbin++){
//double x = hfakepeak->GetBinCenter(xbin);
//int Bincount = hfakepeak->GetBinContent(xbin);
//if (b+a*x>=0) hfakepurepeak->SetBinContent(xbin,Bincount-(b+a*x));
//else hfakepurepeak->SetBinContent(xbin, Bincount);
//If (hfakepeak->GetBinContent(xbin)<0) hfakepurepeak->SetBinContent(xbin, 0);
//cout<<xbin<<"   "<<Bincount<<"   "<<(b+a*x)<<"   "<<Bincount-(b+a*x)<<endl;
//}
// bin modify
for (int xbin = 4000; xbin <= 4407; xbin++) {
	hfakepurepeak->SetBinContent(xbin, hfakepeak->GetBinContent(xbin + 530));
}
for (int xbin = 4407; xbin <= 8000; xbin++) {
	hfakepurepeak->SetBinContent(xbin, hfakepeak->GetBinContent(xbin));
}
hfakepurepeak->Draw();
TFile* fout2 = new TFile("F:/out/G4_rootfiles_with_tree_Eg4156/Fakepurepeak_S31_Gamma4156_Eg4155.84_Tau3.0_SP1.00_AC0.0_0.5k.root", "RECREATE");
hfakepurepeak->Write();
fout2->Close();

//3) Combine the fake peak and fake background
TFile* fin2 = TFile::Open("F:/out/G4_rootfiles_with_tree_Eg4156/Fakebackground.root");
TFile* fin3 = TFile::Open("F:/out/G4_rootfiles_with_tree_Eg4156/Fakepurepeak_S31_Gamma4156_Eg4155.84_Tau3.0_SP1.00_AC0.0_0.5k.root");
TH1F* hfakedata = new TH1F("hfakedata", "hfakedata", 8000, 0., 8000);
fin2->cd();
hfakedata->Add(hfakedata, hbackground_linear, 1, 1);
fin3->cd();
hfakedata->Add(hfakedata, hfakepurepeak, 1, 1);
hfakedata->Draw("");
TFile* fout3 = new TFile("F:/out/G4_rootfiles_with_tree_Eg4156/Fakedata_S31_Gamma4156_Eg4155.84_Tau3.0_SP1.00_AC0.0_scaled_0.5k.root", "RECREATE");
hfakedata->Write();
fout3->Close();

