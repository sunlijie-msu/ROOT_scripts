char filename[400];
sprintf(filename, "Mg23_Gamma7333_Eg7333.20_Tau7.0_SP1.00_AD-0.9_all");

char inputFile[400];
sprintf(inputFile, "F:/out/%s.root.", filename);
TFile* _file0 = TFile::Open(inputFile);

TCanvas* canvasring = new TCanvas("canvasring", "canvasring", 890, 800);//
canvasring->cd();//
canvasring->SetTopMargin(0.03);
canvasring->SetRightMargin(0.12);
canvasring->SetLeftMargin(0.17);
canvasring->SetBottomMargin(0.15);
gStyle->SetFrameLineWidth(3);
canvasring->SetFrameLineWidth(3);
//canvasring->SetLogz();

//TH2D* h2DSSD1xy = new TH2D("h2DSSD1xy", "h2DSSD1xy", 16, 75, 125, 16, 75, 125);
TH2D* h2DSSD1xy = new TH2D("h2DSSD1xy", "h2DSSD1xy", 160, -25, 25, 160, -25, 25);

//tree->Draw("DSSD1y-100:DSSD1x-100>>h2DSSD1xy", "DSSD1e+DSSD2e>100&&DSSD1e>100&&DSSD2e>100", "surf2");
tree->Draw("DSSD1y-100:DSSD1x-100>>h2DSSD1xy", "DSSD1e+DSSD2e>100&&DSSD1e>100&&DSSD2e>100", "colz");

h2DSSD1xy->SetStats(0);
h2DSSD1xy->SetTitle("");
h2DSSD1xy->GetXaxis()->SetTitle("DSSD1 X (mm)");
h2DSSD1xy->GetYaxis()->SetTitle("DSSD1 Y (mm)");
h2DSSD1xy->GetXaxis()->CenterTitle();
h2DSSD1xy->GetYaxis()->CenterTitle();
h2DSSD1xy->GetXaxis()->SetLabelFont(132);
h2DSSD1xy->GetYaxis()->SetLabelFont(132);
h2DSSD1xy->GetXaxis()->SetLabelSize(0.05);
h2DSSD1xy->GetYaxis()->SetLabelSize(0.05);
h2DSSD1xy->GetXaxis()->SetTitleFont(132);
h2DSSD1xy->GetYaxis()->SetTitleFont(132);
h2DSSD1xy->GetXaxis()->SetTitleOffset(1.1);
h2DSSD1xy->GetYaxis()->SetTitleOffset(1.4);
h2DSSD1xy->GetXaxis()->SetTitleSize(0.06);
h2DSSD1xy->GetYaxis()->SetTitleSize(0.06);
h2DSSD1xy->GetXaxis()->SetNdivisions(505);
h2DSSD1xy->GetYaxis()->SetNdivisions(505);
h2DSSD1xy->GetXaxis()->SetRangeUser(-25, 25);
h2DSSD1xy->GetYaxis()->SetRangeUser(-25, 25);
h2DSSD1xy->SetMinimum(1);
h2DSSD1xy->SetMaximum(540);
h2DSSD1xy->SetContour(99);
TPaletteAxis* palette = (TPaletteAxis*)h2DSSD1xy->GetListOfFunctions()->FindObject("palette");
gStyle->SetPalette(kRainbow);
palette->SetTitleFont(132);
palette->SetTitleSize(0.05);
palette->SetLabelFont(132);
palette->SetLabelSize(0.05);
palette->SetLineWidth(3);
//canvasring->SaveAs("F:/out/DSL2_Fig_alpha_counts_rings_Mg23_Gamma7333_Eg7333.20_Tau7.0_SP1.00_AD-1.0_all.png");
char outputFile[400];
sprintf(outputFile, "F:/out/DSL2_Fig_alpha_counts_rings_%s.png", filename);
canvasring->SaveAs(outputFile);

tree->SetMarkerColor(1);
tree->Draw("DSSD1y-100:DSSD1x-100>>h2DSSD1xy", "DSSD1e+DSSD2e>100&&DSSD1e>100&&DSSD2e>100&&DSSD1x>96.875&&DSSD1x<103.125&&DSSD1y>96.875&&DSSD1y<103.125", "");

tree->SetMarkerColor(2);
tree->Draw("DSSD1y-100:DSSD1x-100", "DSSD1e+DSSD2e>100&&DSSD1e>100&&DSSD2e>100&&DSSD1x>96.875&&DSSD1x<103.125&&DSSD1y>103.125&&DSSD1y<106.25", "same");
tree->Draw("DSSD1y-100:DSSD1x-100", "DSSD1e+DSSD2e>100&&DSSD1e>100&&DSSD2e>100&&DSSD1x>96.875&&DSSD1x<103.125&&DSSD1y>93.75&&DSSD1y<96.875", "same");
tree->Draw("DSSD1y-100:DSSD1x-100", "DSSD1e+DSSD2e>100&&DSSD1e>100&&DSSD2e>100&&DSSD1x>103.125&&DSSD1x<106.25&&DSSD1y>96.875&&DSSD1y<103.125", "same");
tree->Draw("DSSD1y-100:DSSD1x-100", "DSSD1e+DSSD2e>100&&DSSD1e>100&&DSSD2e>100&&DSSD1x>93.75&&DSSD1x<96.875&&DSSD1y>96.875&&DSSD1y<103.125", "same");

tree->SetMarkerColor(kGreen + 1);
tree->Draw("DSSD1y-100:DSSD1x-100", "DSSD1e+DSSD2e>100&&DSSD1e>100&&DSSD2e>100&&DSSD1x>93.75&&DSSD1x<106.25&&DSSD1y>106.25&&DSSD1y<109.375", "same");
tree->Draw("DSSD1y-100:DSSD1x-100", "DSSD1e+DSSD2e>100&&DSSD1e>100&&DSSD2e>100&&DSSD1x>93.75&&DSSD1x<106.25&&DSSD1y>90.625&&DSSD1y<93.75", "same");
tree->Draw("DSSD1y-100:DSSD1x-100", "DSSD1e+DSSD2e>100&&DSSD1e>100&&DSSD2e>100&&DSSD1x>90.625&&DSSD1x<96.875&&DSSD1y>103.125&&DSSD1y<106.25", "same");
tree->Draw("DSSD1y-100:DSSD1x-100", "DSSD1e+DSSD2e>100&&DSSD1e>100&&DSSD2e>100&&DSSD1x>103.125&DSSD1x<109.375&&DSSD1y>103.125&&DSSD1y<106.25", "same");
tree->Draw("DSSD1y-100:DSSD1x-100", "DSSD1e+DSSD2e>100&&DSSD1e>100&&DSSD2e>100&&DSSD1x>90.625&&DSSD1x<93.75&&DSSD1y>96.875&&DSSD1y<103.125", "same");
tree->Draw("DSSD1y-100:DSSD1x-100", "DSSD1e+DSSD2e>100&&DSSD1e>100&&DSSD2e>100&&DSSD1x>106.25&&DSSD1x<109.375&&DSSD1y>96.875&&DSSD1y<103.125", "same");
tree->Draw("DSSD1y-100:DSSD1x-100", "DSSD1e+DSSD2e>100&&DSSD1e>100&&DSSD2e>100&&DSSD1x>90.625&&DSSD1x<96.875&&DSSD1y>93.75&&DSSD1y<96.875", "same");
tree->Draw("DSSD1y-100:DSSD1x-100", "DSSD1e+DSSD2e>100&&DSSD1e>100&&DSSD2e>100&&DSSD1x>103.125&&DSSD1x<109.375&&DSSD1y>93.75&&DSSD1y<96.875", "same");

tree->SetMarkerColor(kAzure);
tree->Draw("DSSD1y-100:DSSD1x-100", "DSSD1e+DSSD2e>100&&DSSD1e>100&&DSSD2e>100&&DSSD1x>93.75&&DSSD1x<106.25&&DSSD1y>109.375&&DSSD1y<112.5", "same");
tree->Draw("DSSD1y-100:DSSD1x-100", "DSSD1e+DSSD2e>100&&DSSD1e>100&&DSSD2e>100&&DSSD1x>93.75&&DSSD1x<106.25&&DSSD1y>87.5&&DSSD1y<90.625", "same");
tree->Draw("DSSD1y-100:DSSD1x-100", "DSSD1e+DSSD2e>100&&DSSD1e>100&&DSSD2e>100&&DSSD1x>90.625&&DSSD1x<93.75&&DSSD1y>106.25&&DSSD1y<109.375", "same");
tree->Draw("DSSD1y-100:DSSD1x-100", "DSSD1e+DSSD2e>100&&DSSD1e>100&&DSSD2e>100&&DSSD1x>106.25&DSSD1x<109.375&&DSSD1y>106.25&&DSSD1y<109.375", "same");
tree->Draw("DSSD1y-100:DSSD1x-100", "DSSD1e+DSSD2e>100&&DSSD1e>100&&DSSD2e>100&&DSSD1x>87.5&&DSSD1x<90.625&&DSSD1y>93.75&&DSSD1y<106.25", "same");
tree->Draw("DSSD1y-100:DSSD1x-100", "DSSD1e+DSSD2e>100&&DSSD1e>100&&DSSD2e>100&&DSSD1x>109.375&&DSSD1x<112.5&&DSSD1y>93.75&&DSSD1y<106.25", "same");
tree->Draw("DSSD1y-100:DSSD1x-100", "DSSD1e+DSSD2e>100&&DSSD1e>100&&DSSD2e>100&&DSSD1x>90.625&&DSSD1x<93.75&&DSSD1y>90.625&&DSSD1y<93.75", "same");
tree->Draw("DSSD1y-100:DSSD1x-100", "DSSD1e+DSSD2e>100&&DSSD1e>100&&DSSD2e>100&&DSSD1x>106.25&&DSSD1x<109.375&&DSSD1y>90.625&&DSSD1y<93.75", "same");

tree->SetMarkerColor(kViolet);
tree->Draw("DSSD1y-100:DSSD1x-100", "DSSD1e+DSSD2e>100&&DSSD1e>100&&DSSD2e>100&&DSSD1x>93.75&&DSSD1x<106.25&&DSSD1y>112.5&&DSSD1y<115.625", "same");
tree->Draw("DSSD1y-100:DSSD1x-100", "DSSD1e+DSSD2e>100&&DSSD1e>100&&DSSD2e>100&&DSSD1x>93.75&&DSSD1x<106.25&&DSSD1y>84.375&&DSSD1y<87.5", "same");
tree->Draw("DSSD1y-100:DSSD1x-100", "DSSD1e+DSSD2e>100&&DSSD1e>100&&DSSD2e>100&&DSSD1x>87.5&&DSSD1x<93.75&&DSSD1y>109.375&&DSSD1y<112.5", "same");
tree->Draw("DSSD1y-100:DSSD1x-100", "DSSD1e+DSSD2e>100&&DSSD1e>100&&DSSD2e>100&&DSSD1x>106.25&DSSD1x<112.5&&DSSD1y>109.375&&DSSD1y<112.5", "same");
tree->Draw("DSSD1y-100:DSSD1x-100", "DSSD1e+DSSD2e>100&&DSSD1e>100&&DSSD2e>100&&DSSD1x>87.5&&DSSD1x<90.625&&DSSD1y>106.25&&DSSD1y<109.375", "same");
tree->Draw("DSSD1y-100:DSSD1x-100", "DSSD1e+DSSD2e>100&&DSSD1e>100&&DSSD2e>100&&DSSD1x>109.375&&DSSD1x<112.5&&DSSD1y>106.25&&DSSD1y<109.375", "same");
tree->Draw("DSSD1y-100:DSSD1x-100", "DSSD1e+DSSD2e>100&&DSSD1e>100&&DSSD2e>100&&DSSD1x>84.375&&DSSD1x<87.5&&DSSD1y>93.75&&DSSD1y<106.25", "same");
tree->Draw("DSSD1y-100:DSSD1x-100", "DSSD1e+DSSD2e>100&&DSSD1e>100&&DSSD2e>100&&DSSD1x>112.5&&DSSD1x<115.625&&DSSD1y>93.75&&DSSD1y<106.25", "same");
tree->Draw("DSSD1y-100:DSSD1x-100", "DSSD1e+DSSD2e>100&&DSSD1e>100&&DSSD2e>100&&DSSD1x>87.5&&DSSD1x<90.625&&DSSD1y>90.625&&DSSD1y<93.75", "same");
tree->Draw("DSSD1y-100:DSSD1x-100", "DSSD1e+DSSD2e>100&&DSSD1e>100&&DSSD2e>100&&DSSD1x>109.375&&DSSD1x<112.5&&DSSD1y>90.625&&DSSD1y<93.75", "same");
tree->Draw("DSSD1y-100:DSSD1x-100", "DSSD1e+DSSD2e>100&&DSSD1e>100&&DSSD2e>100&&DSSD1x>87.5&&DSSD1x<93.75&&DSSD1y>87.5&&DSSD1y<90.625", "same");
tree->Draw("DSSD1y-100:DSSD1x-100", "DSSD1e+DSSD2e>100&&DSSD1e>100&&DSSD2e>100&&DSSD1x>106.25&&DSSD1x<112.5&&DSSD1y>87.5&&DSSD1y<90.625", "same");

TCanvas* canvaspeak = new TCanvas("canvaspeak", "canvaspeak", 1300, 800);//
canvaspeak->cd();//
canvaspeak->SetTopMargin(0.025);
canvaspeak->SetRightMargin(0.04);
canvaspeak->SetLeftMargin(0.11);
canvaspeak->SetBottomMargin(0.15);
canvaspeak->SetFrameLineWidth(3);
canvaspeak->SetLogz();

TH1D* h1DSSD_total_e_Ring1 = new TH1D("h1DSSD_total_e_Ring1", "h1DSSD_total_e_Ring1", 6000, 0, 60000);
TH1D* h1DSSD_total_e_Ring2 = new TH1D("h1DSSD_total_e_Ring2", "h1DSSD_total_e_Ring2", 6000, 0, 60000);
TH1D* h1DSSD_total_e_Ring2_sub1 = new TH1D("h1DSSD_total_e_Ring2_sub1", "h1DSSD_total_e_Ring2_sub1", 6000, 0, 60000);
TH1D* h1DSSD_total_e_Ring2_sub2 = new TH1D("h1DSSD_total_e_Ring2_sub2", "h1DSSD_total_e_Ring2_sub2", 6000, 0, 60000);
TH1D* h1DSSD_total_e_Ring2_sub3 = new TH1D("h1DSSD_total_e_Ring2_sub3", "h1DSSD_total_e_Ring2_sub3", 6000, 0, 60000);
TH1D* h1DSSD_total_e_Ring2_sub4 = new TH1D("h1DSSD_total_e_Ring2_sub4", "h1DSSD_total_e_Ring2_sub4", 6000, 0, 60000);
TH1D* h1DSSD_total_e_Ring3 = new TH1D("h1DSSD_total_e_Ring3", "h1DSSD_total_e_Ring3", 6000, 0, 60000);
TH1D* h1DSSD_total_e_Ring3_sub1 = new TH1D("h1DSSD_total_e_Ring3_sub1", "h1DSSD_total_e_Ring3_sub1", 6000, 0, 60000);
TH1D* h1DSSD_total_e_Ring3_sub2 = new TH1D("h1DSSD_total_e_Ring3_sub2", "h1DSSD_total_e_Ring3_sub2", 6000, 0, 60000);
TH1D* h1DSSD_total_e_Ring3_sub3 = new TH1D("h1DSSD_total_e_Ring3_sub3", "h1DSSD_total_e_Ring3_sub3", 6000, 0, 60000);
TH1D* h1DSSD_total_e_Ring3_sub4 = new TH1D("h1DSSD_total_e_Ring3_sub4", "h1DSSD_total_e_Ring3_sub4", 6000, 0, 60000);
TH1D* h1DSSD_total_e_Ring3_sub5 = new TH1D("h1DSSD_total_e_Ring3_sub5", "h1DSSD_total_e_Ring3_sub5", 6000, 0, 60000);
TH1D* h1DSSD_total_e_Ring3_sub6 = new TH1D("h1DSSD_total_e_Ring3_sub6", "h1DSSD_total_e_Ring3_sub6", 6000, 0, 60000);
TH1D* h1DSSD_total_e_Ring3_sub7 = new TH1D("h1DSSD_total_e_Ring3_sub7", "h1DSSD_total_e_Ring3_sub7", 6000, 0, 60000);
TH1D* h1DSSD_total_e_Ring3_sub8 = new TH1D("h1DSSD_total_e_Ring3_sub8", "h1DSSD_total_e_Ring3_sub8", 6000, 0, 60000);
TH1D* h1DSSD_total_e_Ring4 = new TH1D("h1DSSD_total_e_Ring4", "h1DSSD_total_e_Ring4", 6000, 0, 60000);
TH1D* h1DSSD_total_e_Ring4_sub1 = new TH1D("h1DSSD_total_e_Ring4_sub1", "h1DSSD_total_e_Ring4_sub1", 6000, 0, 60000);
TH1D* h1DSSD_total_e_Ring4_sub2 = new TH1D("h1DSSD_total_e_Ring4_sub2", "h1DSSD_total_e_Ring4_sub2", 6000, 0, 60000);
TH1D* h1DSSD_total_e_Ring4_sub3 = new TH1D("h1DSSD_total_e_Ring4_sub3", "h1DSSD_total_e_Ring4_sub3", 6000, 0, 60000);
TH1D* h1DSSD_total_e_Ring4_sub4 = new TH1D("h1DSSD_total_e_Ring4_sub4", "h1DSSD_total_e_Ring4_sub4", 6000, 0, 60000);
TH1D* h1DSSD_total_e_Ring4_sub5 = new TH1D("h1DSSD_total_e_Ring4_sub5", "h1DSSD_total_e_Ring4_sub5", 6000, 0, 60000);
TH1D* h1DSSD_total_e_Ring4_sub6 = new TH1D("h1DSSD_total_e_Ring4_sub6", "h1DSSD_total_e_Ring4_sub6", 6000, 0, 60000);
TH1D* h1DSSD_total_e_Ring4_sub7 = new TH1D("h1DSSD_total_e_Ring4_sub7", "h1DSSD_total_e_Ring4_sub7", 6000, 0, 60000);
TH1D* h1DSSD_total_e_Ring4_sub8 = new TH1D("h1DSSD_total_e_Ring4_sub8", "h1DSSD_total_e_Ring4_sub8", 6000, 0, 60000);
TH1D* h1DSSD_total_e_Ring5 = new TH1D("h1DSSD_total_e_Ring5", "h1DSSD_total_e_Ring5", 6000, 0, 60000);
TH1D* h1DSSD_total_e_Ring5_sub1 = new TH1D("h1DSSD_total_e_Ring5_sub1", "h1DSSD_total_e_Ring5_sub1", 6000, 0, 60000);
TH1D* h1DSSD_total_e_Ring5_sub2 = new TH1D("h1DSSD_total_e_Ring5_sub2", "h1DSSD_total_e_Ring5_sub2", 6000, 0, 60000);
TH1D* h1DSSD_total_e_Ring5_sub3 = new TH1D("h1DSSD_total_e_Ring5_sub3", "h1DSSD_total_e_Ring5_sub3", 6000, 0, 60000);
TH1D* h1DSSD_total_e_Ring5_sub4 = new TH1D("h1DSSD_total_e_Ring5_sub4", "h1DSSD_total_e_Ring5_sub4", 6000, 0, 60000);
TH1D* h1DSSD_total_e_Ring5_sub5 = new TH1D("h1DSSD_total_e_Ring5_sub5", "h1DSSD_total_e_Ring5_sub5", 6000, 0, 60000);
TH1D* h1DSSD_total_e_Ring5_sub6 = new TH1D("h1DSSD_total_e_Ring5_sub6", "h1DSSD_total_e_Ring5_sub6", 6000, 0, 60000);
TH1D* h1DSSD_total_e_Ring5_sub7 = new TH1D("h1DSSD_total_e_Ring5_sub7", "h1DSSD_total_e_Ring5_sub7", 6000, 0, 60000);
TH1D* h1DSSD_total_e_Ring5_sub8 = new TH1D("h1DSSD_total_e_Ring5_sub8", "h1DSSD_total_e_Ring5_sub8", 6000, 0, 60000);
TH1D* h1DSSD_total_e_Ring5_sub9 = new TH1D("h1DSSD_total_e_Ring5_sub9", "h1DSSD_total_e_Ring5_sub9", 6000, 0, 60000);
TH1D* h1DSSD_total_e_Ring5_sub10 = new TH1D("h1DSSD_total_e_Ring5_sub10", "h1DSSD_total_e_Ring5_sub10", 6000, 0, 60000);
TH1D* h1DSSD_total_e_Ring5_sub11 = new TH1D("h1DSSD_total_e_Ring5_sub11", "h1DSSD_total_e_Ring5_sub11", 6000, 0, 60000);
TH1D* h1DSSD_total_e_Ring5_sub12 = new TH1D("h1DSSD_total_e_Ring5_sub12", "h1DSSD_total_e_Ring5_sub12", 6000, 0, 60000);


h1DSSD_total_e_Ring1->SetLineColor(1);
tree->Draw("DSSD1e+DSSD2e>>h1DSSD_total_e_Ring1", "DSSD1e+DSSD2e>100&&DSSD1e>100&&DSSD2e>100&&DSSD1x>96.875&&DSSD1x<103.125&&DSSD1y>96.875&&DSSD1y<103.125", "");

tree->SetLineColor(2);
tree->Draw("DSSD1e+DSSD2e>>h1DSSD_total_e_Ring2_sub1", "DSSD1e+DSSD2e>100&&DSSD1e>100&&DSSD2e>100&&DSSD1x>96.875&&DSSD1x<103.125&&DSSD1y>103.125&&DSSD1y<106.25", "");
tree->Draw("DSSD1e+DSSD2e>>h1DSSD_total_e_Ring2_sub2", "DSSD1e+DSSD2e>100&&DSSD1e>100&&DSSD2e>100&&DSSD1x>96.875&&DSSD1x<103.125&&DSSD1y>93.75&&DSSD1y<96.875", "");
tree->Draw("DSSD1e+DSSD2e>>h1DSSD_total_e_Ring2_sub3", "DSSD1e+DSSD2e>100&&DSSD1e>100&&DSSD2e>100&&DSSD1x>103.125&&DSSD1x<106.25&&DSSD1y>96.875&&DSSD1y<103.125", "");
tree->Draw("DSSD1e+DSSD2e>>h1DSSD_total_e_Ring2_sub4", "DSSD1e+DSSD2e>100&&DSSD1e>100&&DSSD2e>100&&DSSD1x>93.75&&DSSD1x<96.875&&DSSD1y>96.875&&DSSD1y<103.125", "");
h1DSSD_total_e_Ring2->Add(h1DSSD_total_e_Ring2_sub1, h1DSSD_total_e_Ring2_sub2);
h1DSSD_total_e_Ring2->Add(h1DSSD_total_e_Ring2_sub3);
h1DSSD_total_e_Ring2->Add(h1DSSD_total_e_Ring2_sub4);
h1DSSD_total_e_Ring2->SetLineColor(2);
h1DSSD_total_e_Ring2->Draw("same");

tree->SetLineColor(kGreen+1);
tree->Draw("DSSD1e+DSSD2e>>h1DSSD_total_e_Ring3_sub1", "DSSD1e+DSSD2e>100&&DSSD1e>100&&DSSD2e>100&&DSSD1x>93.75&&DSSD1x<106.25&&DSSD1y>106.25&&DSSD1y<109.375", "");
tree->Draw("DSSD1e+DSSD2e>>h1DSSD_total_e_Ring3_sub2", "DSSD1e+DSSD2e>100&&DSSD1e>100&&DSSD2e>100&&DSSD1x>93.75&&DSSD1x<106.25&&DSSD1y>90.625&&DSSD1y<93.75", "");
tree->Draw("DSSD1e+DSSD2e>>h1DSSD_total_e_Ring3_sub3", "DSSD1e+DSSD2e>100&&DSSD1e>100&&DSSD2e>100&&DSSD1x>90.625&&DSSD1x<96.875&&DSSD1y>103.125&&DSSD1y<106.25", "");
tree->Draw("DSSD1e+DSSD2e>>h1DSSD_total_e_Ring3_sub4", "DSSD1e+DSSD2e>100&&DSSD1e>100&&DSSD2e>100&&DSSD1x>103.125&DSSD1x<109.375&&DSSD1y>103.125&&DSSD1y<106.25", "");
tree->Draw("DSSD1e+DSSD2e>>h1DSSD_total_e_Ring3_sub5", "DSSD1e+DSSD2e>100&&DSSD1e>100&&DSSD2e>100&&DSSD1x>90.625&&DSSD1x<93.75&&DSSD1y>96.875&&DSSD1y<103.125", "");
tree->Draw("DSSD1e+DSSD2e>>h1DSSD_total_e_Ring3_sub6", "DSSD1e+DSSD2e>100&&DSSD1e>100&&DSSD2e>100&&DSSD1x>106.25&&DSSD1x<109.375&&DSSD1y>96.875&&DSSD1y<103.125", "");
tree->Draw("DSSD1e+DSSD2e>>h1DSSD_total_e_Ring3_sub7", "DSSD1e+DSSD2e>100&&DSSD1e>100&&DSSD2e>100&&DSSD1x>90.625&&DSSD1x<96.875&&DSSD1y>93.75&&DSSD1y<96.875", "");
tree->Draw("DSSD1e+DSSD2e>>h1DSSD_total_e_Ring3_sub8", "DSSD1e+DSSD2e>100&&DSSD1e>100&&DSSD2e>100&&DSSD1x>103.125&&DSSD1x<109.375&&DSSD1y>93.75&&DSSD1y<96.875", "");
h1DSSD_total_e_Ring3->Add(h1DSSD_total_e_Ring3_sub1, h1DSSD_total_e_Ring3_sub2);
h1DSSD_total_e_Ring3->Add(h1DSSD_total_e_Ring3_sub3);
h1DSSD_total_e_Ring3->Add(h1DSSD_total_e_Ring3_sub4);
h1DSSD_total_e_Ring3->Add(h1DSSD_total_e_Ring3_sub5);
h1DSSD_total_e_Ring3->Add(h1DSSD_total_e_Ring3_sub6);
h1DSSD_total_e_Ring3->Add(h1DSSD_total_e_Ring3_sub7);
h1DSSD_total_e_Ring3->Add(h1DSSD_total_e_Ring3_sub8);
h1DSSD_total_e_Ring3->SetLineColor(kGreen+1);
h1DSSD_total_e_Ring3->Draw("same");

tree->SetLineColor(kAzure);
tree->Draw("DSSD1e+DSSD2e>>h1DSSD_total_e_Ring4_sub1", "DSSD1e+DSSD2e>100&&DSSD1e>100&&DSSD2e>100&&DSSD1x>93.75&&DSSD1x<106.25&&DSSD1y>109.375&&DSSD1y<112.5", "");
tree->Draw("DSSD1e+DSSD2e>>h1DSSD_total_e_Ring4_sub2", "DSSD1e+DSSD2e>100&&DSSD1e>100&&DSSD2e>100&&DSSD1x>93.75&&DSSD1x<106.25&&DSSD1y>87.5&&DSSD1y<90.625", "");
tree->Draw("DSSD1e+DSSD2e>>h1DSSD_total_e_Ring4_sub3", "DSSD1e+DSSD2e>100&&DSSD1e>100&&DSSD2e>100&&DSSD1x>90.625&&DSSD1x<93.75&&DSSD1y>106.25&&DSSD1y<109.375", "");
tree->Draw("DSSD1e+DSSD2e>>h1DSSD_total_e_Ring4_sub4", "DSSD1e+DSSD2e>100&&DSSD1e>100&&DSSD2e>100&&DSSD1x>106.25&&DSSD1x<109.375&&DSSD1y>106.25&&DSSD1y<109.375", "");
tree->Draw("DSSD1e+DSSD2e>>h1DSSD_total_e_Ring4_sub5", "DSSD1e+DSSD2e>100&&DSSD1e>100&&DSSD2e>100&&DSSD1x>87.5&&DSSD1x<90.625&&DSSD1y>93.75&&DSSD1y<106.25", "");
tree->Draw("DSSD1e+DSSD2e>>h1DSSD_total_e_Ring4_sub6", "DSSD1e+DSSD2e>100&&DSSD1e>100&&DSSD2e>100&&DSSD1x>109.375&&DSSD1x<112.5&&DSSD1y>93.75&&DSSD1y<106.25", "");
tree->Draw("DSSD1e+DSSD2e>>h1DSSD_total_e_Ring4_sub7", "DSSD1e+DSSD2e>100&&DSSD1e>100&&DSSD2e>100&&DSSD1x>90.625&&DSSD1x<93.75&&DSSD1y>90.625&&DSSD1y<93.75", "");
tree->Draw("DSSD1e+DSSD2e>>h1DSSD_total_e_Ring4_sub8", "DSSD1e+DSSD2e>100&&DSSD1e>100&&DSSD2e>100&&DSSD1x>106.25&&DSSD1x<109.375&&DSSD1y>90.625&&DSSD1y<93.75", "");
h1DSSD_total_e_Ring4->Add(h1DSSD_total_e_Ring4_sub1, h1DSSD_total_e_Ring4_sub2);
h1DSSD_total_e_Ring4->Add(h1DSSD_total_e_Ring4_sub3);
h1DSSD_total_e_Ring4->Add(h1DSSD_total_e_Ring4_sub4);
h1DSSD_total_e_Ring4->Add(h1DSSD_total_e_Ring4_sub5);
h1DSSD_total_e_Ring4->Add(h1DSSD_total_e_Ring4_sub6);
h1DSSD_total_e_Ring4->Add(h1DSSD_total_e_Ring4_sub7);
h1DSSD_total_e_Ring4->Add(h1DSSD_total_e_Ring4_sub8);
h1DSSD_total_e_Ring4->SetLineColor(kAzure);

tree->SetLineColor(kViolet);
tree->Draw("DSSD1e+DSSD2e>>h1DSSD_total_e_Ring5_sub1", "DSSD1e+DSSD2e>100&&DSSD1e>100&&DSSD2e>100&&DSSD1x>93.75&&DSSD1x<106.25&&DSSD1y>112.5&&DSSD1y<115.625", "");
tree->Draw("DSSD1e+DSSD2e>>h1DSSD_total_e_Ring5_sub2", "DSSD1e+DSSD2e>100&&DSSD1e>100&&DSSD2e>100&&DSSD1x>93.75&&DSSD1x<106.25&&DSSD1y>84.375&&DSSD1y<87.5", "");
tree->Draw("DSSD1e+DSSD2e>>h1DSSD_total_e_Ring5_sub3", "DSSD1e+DSSD2e>100&&DSSD1e>100&&DSSD2e>100&&DSSD1x>87.5&&DSSD1x<93.75&&DSSD1y>109.375&&DSSD1y<112.5", "");
tree->Draw("DSSD1e+DSSD2e>>h1DSSD_total_e_Ring5_sub4", "DSSD1e+DSSD2e>100&&DSSD1e>100&&DSSD2e>100&&DSSD1x>106.25&&DSSD1x<112.5&&DSSD1y>109.375&&DSSD1y<112.5", "");
tree->Draw("DSSD1e+DSSD2e>>h1DSSD_total_e_Ring5_sub5", "DSSD1e+DSSD2e>100&&DSSD1e>100&&DSSD2e>100&&DSSD1x>87.5&&DSSD1x<90.625&&DSSD1y>106.25&&DSSD1y<109.375", "");
tree->Draw("DSSD1e+DSSD2e>>h1DSSD_total_e_Ring5_sub6", "DSSD1e+DSSD2e>100&&DSSD1e>100&&DSSD2e>100&&DSSD1x>109.375&&DSSD1x<112.5&&DSSD1y>106.25&&DSSD1y<109.375", "");
tree->Draw("DSSD1e+DSSD2e>>h1DSSD_total_e_Ring5_sub7", "DSSD1e+DSSD2e>100&&DSSD1e>100&&DSSD2e>100&&DSSD1x>84.375&&DSSD1x<87.5&&DSSD1y>93.75&&DSSD1y<106.25", "");
tree->Draw("DSSD1e+DSSD2e>>h1DSSD_total_e_Ring5_sub8", "DSSD1e+DSSD2e>100&&DSSD1e>100&&DSSD2e>100&&DSSD1x>112.5&&DSSD1x<115.625&&DSSD1y>93.75&&DSSD1y<106.25", "");
tree->Draw("DSSD1e+DSSD2e>>h1DSSD_total_e_Ring5_sub9", "DSSD1e+DSSD2e>100&&DSSD1e>100&&DSSD2e>100&&DSSD1x>87.5&&DSSD1x<90.625&&DSSD1y>90.625&&DSSD1y<93.75", "");
tree->Draw("DSSD1e+DSSD2e>>h1DSSD_total_e_Ring5_sub10", "DSSD1e+DSSD2e>100&&DSSD1e>100&&DSSD2e>100&&DSSD1x>109.375&&DSSD1x<112.5&&DSSD1y>90.625&&DSSD1y<93.75", "");
tree->Draw("DSSD1e+DSSD2e>>h1DSSD_total_e_Ring5_sub11", "DSSD1e+DSSD2e>100&&DSSD1e>100&&DSSD2e>100&&DSSD1x>87.5&&DSSD1x<93.75&&DSSD1y>87.5&&DSSD1y<90.625", "");
tree->Draw("DSSD1e+DSSD2e>>h1DSSD_total_e_Ring5_sub12", "DSSD1e+DSSD2e>100&&DSSD1e>100&&DSSD2e>100&&DSSD1x>106.25&&DSSD1x<112.5&&DSSD1y>87.5&&DSSD1y<90.625", "");
h1DSSD_total_e_Ring5->Add(h1DSSD_total_e_Ring5_sub1, h1DSSD_total_e_Ring5_sub2);
h1DSSD_total_e_Ring5->Add(h1DSSD_total_e_Ring5_sub3);
h1DSSD_total_e_Ring5->Add(h1DSSD_total_e_Ring5_sub4);
h1DSSD_total_e_Ring5->Add(h1DSSD_total_e_Ring5_sub5);
h1DSSD_total_e_Ring5->Add(h1DSSD_total_e_Ring5_sub6);
h1DSSD_total_e_Ring5->Add(h1DSSD_total_e_Ring5_sub7);
h1DSSD_total_e_Ring5->Add(h1DSSD_total_e_Ring5_sub8);
h1DSSD_total_e_Ring5->Add(h1DSSD_total_e_Ring5_sub9);
h1DSSD_total_e_Ring5->Add(h1DSSD_total_e_Ring5_sub10);
h1DSSD_total_e_Ring5->Add(h1DSSD_total_e_Ring5_sub11);
h1DSSD_total_e_Ring5->Add(h1DSSD_total_e_Ring5_sub12);
h1DSSD_total_e_Ring5->SetLineColor(kViolet);



h1DSSD_total_e_Ring1->Draw();
h1DSSD_total_e_Ring2->Draw("same");
h1DSSD_total_e_Ring3->Draw("same");
h1DSSD_total_e_Ring4->Draw("same");
h1DSSD_total_e_Ring5->Draw("same");
h1DSSD_total_e_Ring1->SetStats(0);
h1DSSD_total_e_Ring1->SetTitle("");
h1DSSD_total_e_Ring1->GetXaxis()->SetTitle("DSSD1e+DSSD2e (keV)");
h1DSSD_total_e_Ring1->GetYaxis()->SetTitle("Counts per 10 keV");
h1DSSD_total_e_Ring1->GetXaxis()->CenterTitle();
h1DSSD_total_e_Ring1->GetYaxis()->CenterTitle();
h1DSSD_total_e_Ring1->GetXaxis()->SetLabelFont(132);
h1DSSD_total_e_Ring1->GetYaxis()->SetLabelFont(132);
h1DSSD_total_e_Ring1->GetXaxis()->SetLabelSize(0.05);
h1DSSD_total_e_Ring1->GetYaxis()->SetLabelSize(0.05);
h1DSSD_total_e_Ring1->GetXaxis()->SetTitleFont(132);
h1DSSD_total_e_Ring1->GetYaxis()->SetTitleFont(132);
h1DSSD_total_e_Ring1->GetXaxis()->SetTitleOffset(1.1);
h1DSSD_total_e_Ring1->GetYaxis()->SetTitleOffset(0.9);
h1DSSD_total_e_Ring1->GetXaxis()->SetTitleSize(0.06);
h1DSSD_total_e_Ring1->GetYaxis()->SetTitleSize(0.06);
h1DSSD_total_e_Ring1->GetXaxis()->SetNdivisions(505);
h1DSSD_total_e_Ring1->GetYaxis()->SetNdivisions(505);
h1DSSD_total_e_Ring1->GetYaxis()->SetTickLength(0.02);
h1DSSD_total_e_Ring1->GetXaxis()->SetRangeUser(7000, 25000);
h1DSSD_total_e_Ring1->GetYaxis()->SetRangeUser(0, 220);
gPad->RedrawAxis();
canvaspeak->SaveAs("DSL2_Fig_alpha_energy_rings.png");

TFile* _file0 = TFile::Open("Mg23_Gamma7333_Eg7333.20_Tau7.0_SP1.00_AD-0.3_all.root.");
double R1, R2, R3, R4, R5, R_total;
double Ea_min = 100, Ea_max = 30000;

Ea_min = 100; Ea_max = 30000;
R1 = tree->GetEntries(Form("DSSD1e+DSSD2e>%f&&DSSD1e+DSSD2e<%f&&DSSD1e>100&&DSSD2e>100&&DSSD1x>96.875&&DSSD1x<103.125&&DSSD1y>96.875&&DSSD1y<103.125", Ea_min, Ea_max));

Ea_min = 100; Ea_max = 30000;
R2 = tree->GetEntries(Form("DSSD1e+DSSD2e>%f&&DSSD1e+DSSD2e<%f&&DSSD1e>100&&DSSD2e>100&&DSSD1x>96.875&&DSSD1x<103.125&&DSSD1y>103.125&&DSSD1y<106.25", Ea_min, Ea_max));
R2 += tree->GetEntries(Form("DSSD1e+DSSD2e>%f&&DSSD1e+DSSD2e<%f&&DSSD1e>100&&DSSD2e>100&&DSSD1x>96.875&&DSSD1x<103.125&&DSSD1y>93.75&&DSSD1y<96.875", Ea_min, Ea_max));
R2 += tree->GetEntries(Form("DSSD1e+DSSD2e>%f&&DSSD1e+DSSD2e<%f&&DSSD1e>100&&DSSD2e>100&&DSSD1x>103.125&&DSSD1x<106.25&&DSSD1y>96.875&&DSSD1y<103.125", Ea_min, Ea_max));
R2 += tree->GetEntries(Form("DSSD1e+DSSD2e>%f&&DSSD1e+DSSD2e<%f&&DSSD1e>100&&DSSD2e>100&&DSSD1x>93.75&&DSSD1x<96.875&&DSSD1y>96.875&&DSSD1y<103.125", Ea_min, Ea_max));

Ea_min = 100; Ea_max = 30000;
R3 = tree->GetEntries(Form("DSSD1e+DSSD2e>%f&&DSSD1e+DSSD2e<%f&&DSSD1e>100&&DSSD2e>100&&DSSD1x>93.75&&DSSD1x<106.25&&DSSD1y>106.25&&DSSD1y<109.375", Ea_min, Ea_max));
R3 += tree->GetEntries(Form("DSSD1e+DSSD2e>%f&&DSSD1e+DSSD2e<%f&&DSSD1e>100&&DSSD2e>100&&DSSD1x>93.75&&DSSD1x<106.25&&DSSD1y>90.625&&DSSD1y<93.75", Ea_min, Ea_max));
R3 += tree->GetEntries(Form("DSSD1e+DSSD2e>%f&&DSSD1e+DSSD2e<%f&&DSSD1e>100&&DSSD2e>100&&DSSD1x>90.625&&DSSD1x<96.875&&DSSD1y>103.125&&DSSD1y<106.25", Ea_min, Ea_max));
R3 += tree->GetEntries(Form("DSSD1e+DSSD2e>%f&&DSSD1e+DSSD2e<%f&&DSSD1e>100&&DSSD2e>100&&DSSD1x>103.125&DSSD1x<109.375&&DSSD1y>103.125&&DSSD1y<106.25", Ea_min, Ea_max));
R3 += tree->GetEntries(Form("DSSD1e+DSSD2e>%f&&DSSD1e+DSSD2e<%f&&DSSD1e>100&&DSSD2e>100&&DSSD1x>90.625&&DSSD1x<93.75&&DSSD1y>96.875&&DSSD1y<103.125", Ea_min, Ea_max));
R3 += tree->GetEntries(Form("DSSD1e+DSSD2e>%f&&DSSD1e+DSSD2e<%f&&DSSD1e>100&&DSSD2e>100&&DSSD1x>106.25&&DSSD1x<109.375&&DSSD1y>96.875&&DSSD1y<103.125", Ea_min, Ea_max));
R3 += tree->GetEntries(Form("DSSD1e+DSSD2e>%f&&DSSD1e+DSSD2e<%f&&DSSD1e>100&&DSSD2e>100&&DSSD1x>90.625&&DSSD1x<96.875&&DSSD1y>93.75&&DSSD1y<96.875", Ea_min, Ea_max));
R3 += tree->GetEntries(Form("DSSD1e+DSSD2e>%f&&DSSD1e+DSSD2e<%f&&DSSD1e>100&&DSSD2e>100&&DSSD1x>103.125&&DSSD1x<109.375&&DSSD1y>93.75&&DSSD1y<96.875", Ea_min, Ea_max));

Ea_min = 100; Ea_max = 30000;
R4 = tree->GetEntries(Form("DSSD1e+DSSD2e>%f&&DSSD1e+DSSD2e<%f&&DSSD1e>100&&DSSD2e>100&&DSSD1x>93.75&&DSSD1x<106.25&&DSSD1y>109.375&&DSSD1y<112.5", Ea_min, Ea_max));
R4 += tree->GetEntries(Form("DSSD1e+DSSD2e>%f&&DSSD1e+DSSD2e<%f&&DSSD1e>100&&DSSD2e>100&&DSSD1x>93.75&&DSSD1x<106.25&&DSSD1y>87.5&&DSSD1y<90.625", Ea_min, Ea_max));
R4 += tree->GetEntries(Form("DSSD1e+DSSD2e>%f&&DSSD1e+DSSD2e<%f&&DSSD1e>100&&DSSD2e>100&&DSSD1x>90.625&&DSSD1x<93.75&&DSSD1y>106.25&&DSSD1y<109.375", Ea_min, Ea_max));
R4 += tree->GetEntries(Form("DSSD1e+DSSD2e>%f&&DSSD1e+DSSD2e<%f&&DSSD1e>100&&DSSD2e>100&&DSSD1x>106.25&DSSD1x<109.375&&DSSD1y>106.25&&DSSD1y<109.375", Ea_min, Ea_max));
R4 += tree->GetEntries(Form("DSSD1e+DSSD2e>%f&&DSSD1e+DSSD2e<%f&&DSSD1e>100&&DSSD2e>100&&DSSD1x>87.5&&DSSD1x<90.625&&DSSD1y>93.75&&DSSD1y<106.25", Ea_min, Ea_max));
R4 += tree->GetEntries(Form("DSSD1e+DSSD2e>%f&&DSSD1e+DSSD2e<%f&&DSSD1e>100&&DSSD2e>100&&DSSD1x>109.375&&DSSD1x<112.5&&DSSD1y>93.75&&DSSD1y<106.25", Ea_min, Ea_max));
R4 += tree->GetEntries(Form("DSSD1e+DSSD2e>%f&&DSSD1e+DSSD2e<%f&&DSSD1e>100&&DSSD2e>100&&DSSD1x>90.625&&DSSD1x<93.75&&DSSD1y>90.625&&DSSD1y<93.75", Ea_min, Ea_max));
R4 += tree->GetEntries(Form("DSSD1e+DSSD2e>%f&&DSSD1e+DSSD2e<%f&&DSSD1e>100&&DSSD2e>100&&DSSD1x>106.25&&DSSD1x<109.375&&DSSD1y>90.625&&DSSD1y<93.75", Ea_min, Ea_max));

R5 = tree->GetEntries(Form("DSSD1e+DSSD2e>%f&&DSSD1e+DSSD2e<%f&&DSSD1e>100&&DSSD2e>100&&DSSD1x>93.75&&DSSD1x<106.25&&DSSD1y>112.5&&DSSD1y<115.625", Ea_min, Ea_max));
R5 += tree->GetEntries(Form("DSSD1e+DSSD2e>%f&&DSSD1e+DSSD2e<%f&&DSSD1e>100&&DSSD2e>100&&DSSD1x>93.75&&DSSD1x<106.25&&DSSD1y>84.375&&DSSD1y<87.5", Ea_min, Ea_max));
R5 += tree->GetEntries(Form("DSSD1e+DSSD2e>%f&&DSSD1e+DSSD2e<%f&&DSSD1e>100&&DSSD2e>100&&DSSD1x>87.5&&DSSD1x<93.75&&DSSD1y>109.375&&DSSD1y<112.5", Ea_min, Ea_max));
R5 += tree->GetEntries(Form("DSSD1e+DSSD2e>%f&&DSSD1e+DSSD2e<%f&&DSSD1e>100&&DSSD2e>100&&DSSD1x>106.25&DSSD1x<112.5&&DSSD1y>109.375&&DSSD1y<112.5", Ea_min, Ea_max));
R5 += tree->GetEntries(Form("DSSD1e+DSSD2e>%f&&DSSD1e+DSSD2e<%f&&DSSD1e>100&&DSSD2e>100&&DSSD1x>87.5&&DSSD1x<90.625&&DSSD1y>106.25&&DSSD1y<109.375", Ea_min, Ea_max));
R5 += tree->GetEntries(Form("DSSD1e+DSSD2e>%f&&DSSD1e+DSSD2e<%f&&DSSD1e>100&&DSSD2e>100&&DSSD1x>109.375&&DSSD1x<112.5&&DSSD1y>106.25&&DSSD1y<109.375", Ea_min, Ea_max));
R5 += tree->GetEntries(Form("DSSD1e+DSSD2e>%f&&DSSD1e+DSSD2e<%f&&DSSD1e>100&&DSSD2e>100&&DSSD1x>84.375&&DSSD1x<87.5&&DSSD1y>93.75&&DSSD1y<106.25", Ea_min, Ea_max));
R5 += tree->GetEntries(Form("DSSD1e+DSSD2e>%f&&DSSD1e+DSSD2e<%f&&DSSD1e>100&&DSSD2e>100&&DSSD1x>112.5&&DSSD1x<115.625&&DSSD1y>93.75&&DSSD1y<106.25", Ea_min, Ea_max));
R5 += tree->GetEntries(Form("DSSD1e+DSSD2e>%f&&DSSD1e+DSSD2e<%f&&DSSD1e>100&&DSSD2e>100&&DSSD1x>87.5&&DSSD1x<90.625&&DSSD1y>90.625&&DSSD1y<93.75", Ea_min, Ea_max));
R5 += tree->GetEntries(Form("DSSD1e+DSSD2e>%f&&DSSD1e+DSSD2e<%f&&DSSD1e>100&&DSSD2e>100&&DSSD1x>109.375&&DSSD1x<112.5&&DSSD1y>90.625&&DSSD1y<93.75", Ea_min, Ea_max));
R5 += tree->GetEntries(Form("DSSD1e+DSSD2e>%f&&DSSD1e+DSSD2e<%f&&DSSD1e>100&&DSSD2e>100&&DSSD1x>87.5&&DSSD1x<93.75&&DSSD1y>87.5&&DSSD1y<90.625", Ea_min, Ea_max));
R5 += tree->GetEntries(Form("DSSD1e+DSSD2e>%f&&DSSD1e+DSSD2e<%f&&DSSD1e>100&&DSSD2e>100&&DSSD1x>106.25&&DSSD1x<112.5&&DSSD1y>87.5&&DSSD1y<90.625", Ea_min, Ea_max));

R_total = tree->GetEntries("DSSD1e+DSSD2e>100&&DSSD1e>100&&DSSD2e>100");

cout << "R1 = " << R1 << endl;
cout << "R2 = " << R2 << endl;
cout << "R3 = " << R3 << endl;
cout << "R4 = " << R4 << endl;
cout << "R5 = " << R5 << endl;
cout << "R_total = " << R_total << endl;

double R1_ratio, R2_ratio, R3_ratio, R4_ratio, R5_ratio;
R1_ratio = R1 / R_total;
R2_ratio = R2 / R_total;
R3_ratio = R3 / R_total;
R4_ratio = R4 / R_total;
R5_ratio = R5 / R_total;

cout << "R_ratios= " << R1_ratio << " " << R2_ratio << " " << R3_ratio << " " << R4_ratio << " " << R5_ratio << endl;
