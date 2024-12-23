TFile* _file0 = TFile::Open("Mg23_Gamma7333_Eg7333.20_Tau7.0_SP1.00_AC0.0_all.root.");
TCanvas* canvaspeak = new TCanvas("canvaspeak", "canvaspeak", 890, 800);//
canvaspeak->cd();//
canvaspeak->SetTopMargin(0.03);
canvaspeak->SetRightMargin(0.12);
canvaspeak->SetLeftMargin(0.17);
canvaspeak->SetBottomMargin(0.15);
gStyle->SetFrameLineWidth(3);
canvaspeak->SetFrameLineWidth(3);
canvaspeak->SetLogz();

//TH2D* h2DSSD1xy = new TH2D("h2DSSD1xy", "h2DSSD1xy", 16, 75, 125, 16, 75, 125);
TH2D* h2DSSD1xy = new TH2D("h2DSSD1xy", "h2DSSD1xy", 1600, -25, 25, 1600, -25, 25);
tree->Draw("DSSD1y-100:DSSD1x-100>>h2DSSD1xy", "DSSD1e+DSSD2e>100&&DSSD1e>100&&DSSD2e>100", "surf2");

tree->SetMarkerColor(1);
tree->Draw("DSSD1y-100:DSSD1x-100>>h2DSSD1xy", "DSSD1e+DSSD2e>100&&DSSD1e>100&&DSSD2e>100&&DSSD1x>96.875&&DSSD1x<103.125&&DSSD1y>96.875&&DSSD1y<103.125", "");

tree->SetMarkerColor(2);
tree->Draw("DSSD1y-100:DSSD1x-100", "DSSD1e+DSSD2e>100&&DSSD1e>100&&DSSD2e>100&&DSSD1x>96.875&&DSSD1x<103.125&&DSSD1y>103.125&&DSSD1y<106.25", "same");
tree->Draw("DSSD1y-100:DSSD1x-100", "DSSD1e+DSSD2e>100&&DSSD1e>100&&DSSD2e>100&&DSSD1x>96.875&&DSSD1x<103.125&&DSSD1y>93.75&&DSSD1y<96.875", "same");
tree->Draw("DSSD1y-100:DSSD1x-100", "DSSD1e+DSSD2e>100&&DSSD1e>100&&DSSD2e>100&&DSSD1x>103.125&&DSSD1x<106.25&&DSSD1y>96.875&&DSSD1y<103.125", "same");
tree->Draw("DSSD1y-100:DSSD1x-100", "DSSD1e+DSSD2e>100&&DSSD1e>100&&DSSD2e>100&&DSSD1x>93.75&&DSSD1x<96.875&&DSSD1y>96.875&&DSSD1y<103.125", "same");

tree->SetMarkerColor(kGreen+1);
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
//h2DSSD1xy->SetMaximum(100);
h2DSSD1xy->SetContour(99);
TPaletteAxis* palette = (TPaletteAxis*)h2DSSD1xy->GetListOfFunctions()->FindObject("palette");
gStyle->SetPalette(kRainbow);
palette->SetTitleFont(132);
palette->SetTitleSize(0.05);
palette->SetLabelFont(132);
palette->SetLabelSize(0.05);
palette->SetLineWidth(3);
canvaspeak->SaveAs("F:/out/G4_rootfiles_with_tree_Eg7333/DSL2_Fig_alpha_counts_rings.png");

double R1, R2, R3, R4, R5, R_total;
R1 = tree->GetEntries("DSSD1e+DSSD2e>100&&DSSD1e>100&&DSSD2e>100&&DSSD1x>96.875&&DSSD1x<103.125&&DSSD1y>96.875&&DSSD1y<103.125");

R2 = tree->GetEntries("DSSD1e+DSSD2e>100&&DSSD1e>100&&DSSD2e>100&&DSSD1x>96.875&&DSSD1x<103.125&&DSSD1y>103.125&&DSSD1y<106.25");
R2 += tree->GetEntries("DSSD1e+DSSD2e>100&&DSSD1e>100&&DSSD2e>100&&DSSD1x>96.875&&DSSD1x<103.125&&DSSD1y>93.75&&DSSD1y<96.875");
R2 += tree->GetEntries("DSSD1e+DSSD2e>100&&DSSD1e>100&&DSSD2e>100&&DSSD1x>103.125&&DSSD1x<106.25&&DSSD1y>96.875&&DSSD1y<103.125");
R2 += tree->GetEntries("DSSD1e+DSSD2e>100&&DSSD1e>100&&DSSD2e>100&&DSSD1x>93.75&&DSSD1x<96.875&&DSSD1y>96.875&&DSSD1y<103.125");

R3 = tree->GetEntries("DSSD1e+DSSD2e>100&&DSSD1e>100&&DSSD2e>100&&DSSD1x>93.75&&DSSD1x<106.25&&DSSD1y>106.25&&DSSD1y<109.375");
R3 += tree->GetEntries("DSSD1e+DSSD2e>100&&DSSD1e>100&&DSSD2e>100&&DSSD1x>93.75&&DSSD1x<106.25&&DSSD1y>90.625&&DSSD1y<93.75");
R3 += tree->GetEntries("DSSD1e+DSSD2e>100&&DSSD1e>100&&DSSD2e>100&&DSSD1x>90.625&&DSSD1x<96.875&&DSSD1y>103.125&&DSSD1y<106.25");
R3 += tree->GetEntries("DSSD1e+DSSD2e>100&&DSSD1e>100&&DSSD2e>100&&DSSD1x>103.125&DSSD1x<109.375&&DSSD1y>103.125&&DSSD1y<106.25");
R3 += tree->GetEntries("DSSD1e+DSSD2e>100&&DSSD1e>100&&DSSD2e>100&&DSSD1x>90.625&&DSSD1x<93.75&&DSSD1y>96.875&&DSSD1y<103.125");
R3 += tree->GetEntries("DSSD1e+DSSD2e>100&&DSSD1e>100&&DSSD2e>100&&DSSD1x>106.25&&DSSD1x<109.375&&DSSD1y>96.875&&DSSD1y<103.125");
R3 += tree->GetEntries("DSSD1e+DSSD2e>100&&DSSD1e>100&&DSSD2e>100&&DSSD1x>90.625&&DSSD1x<96.875&&DSSD1y>93.75&&DSSD1y<96.875");
R3 += tree->GetEntries("DSSD1e+DSSD2e>100&&DSSD1e>100&&DSSD2e>100&&DSSD1x>103.125&&DSSD1x<109.375&&DSSD1y>93.75&&DSSD1y<96.875");

R4 = tree->GetEntries("DSSD1e+DSSD2e>100&&DSSD1e>100&&DSSD2e>100&&DSSD1x>93.75&&DSSD1x<106.25&&DSSD1y>109.375&&DSSD1y<112.5");
R4 += tree->GetEntries("DSSD1e+DSSD2e>100&&DSSD1e>100&&DSSD2e>100&&DSSD1x>93.75&&DSSD1x<106.25&&DSSD1y>87.5&&DSSD1y<90.625");
R4 += tree->GetEntries("DSSD1e+DSSD2e>100&&DSSD1e>100&&DSSD2e>100&&DSSD1x>90.625&&DSSD1x<93.75&&DSSD1y>106.25&&DSSD1y<109.375");
R4 += tree->GetEntries("DSSD1e+DSSD2e>100&&DSSD1e>100&&DSSD2e>100&&DSSD1x>106.25&DSSD1x<109.375&&DSSD1y>106.25&&DSSD1y<109.375");
R4 += tree->GetEntries("DSSD1e+DSSD2e>100&&DSSD1e>100&&DSSD2e>100&&DSSD1x>87.5&&DSSD1x<90.625&&DSSD1y>93.75&&DSSD1y<106.25");
R4 += tree->GetEntries("DSSD1e+DSSD2e>100&&DSSD1e>100&&DSSD2e>100&&DSSD1x>109.375&&DSSD1x<112.5&&DSSD1y>93.75&&DSSD1y<106.25");
R4 += tree->GetEntries("DSSD1e+DSSD2e>100&&DSSD1e>100&&DSSD2e>100&&DSSD1x>90.625&&DSSD1x<93.75&&DSSD1y>90.625&&DSSD1y<93.75");
R4 += tree->GetEntries("DSSD1e+DSSD2e>100&&DSSD1e>100&&DSSD2e>100&&DSSD1x>106.25&&DSSD1x<109.375&&DSSD1y>90.625&&DSSD1y<93.75");

R5 = tree->GetEntries("DSSD1e+DSSD2e>100&&DSSD1e>100&&DSSD2e>100&&DSSD1x>93.75&&DSSD1x<106.25&&DSSD1y>112.5&&DSSD1y<115.625");
R5 += tree->GetEntries("DSSD1e+DSSD2e>100&&DSSD1e>100&&DSSD2e>100&&DSSD1x>93.75&&DSSD1x<106.25&&DSSD1y>84.375&&DSSD1y<87.5");
R5 += tree->GetEntries("DSSD1e+DSSD2e>100&&DSSD1e>100&&DSSD2e>100&&DSSD1x>87.5&&DSSD1x<93.75&&DSSD1y>109.375&&DSSD1y<112.5");
R5 += tree->GetEntries("DSSD1e+DSSD2e>100&&DSSD1e>100&&DSSD2e>100&&DSSD1x>106.25&DSSD1x<112.5&&DSSD1y>109.375&&DSSD1y<112.5");
R5 += tree->GetEntries("DSSD1e+DSSD2e>100&&DSSD1e>100&&DSSD2e>100&&DSSD1x>87.5&&DSSD1x<90.625&&DSSD1y>106.25&&DSSD1y<109.375");
R5 += tree->GetEntries("DSSD1e+DSSD2e>100&&DSSD1e>100&&DSSD2e>100&&DSSD1x>109.375&&DSSD1x<112.5&&DSSD1y>106.25&&DSSD1y<109.375");
R5 += tree->GetEntries("DSSD1e+DSSD2e>100&&DSSD1e>100&&DSSD2e>100&&DSSD1x>84.375&&DSSD1x<87.5&&DSSD1y>93.75&&DSSD1y<106.25");
R5 += tree->GetEntries("DSSD1e+DSSD2e>100&&DSSD1e>100&&DSSD2e>100&&DSSD1x>112.5&&DSSD1x<115.625&&DSSD1y>93.75&&DSSD1y<106.25");
R5 += tree->GetEntries("DSSD1e+DSSD2e>100&&DSSD1e>100&&DSSD2e>100&&DSSD1x>87.5&&DSSD1x<90.625&&DSSD1y>90.625&&DSSD1y<93.75");
R5 += tree->GetEntries("DSSD1e+DSSD2e>100&&DSSD1e>100&&DSSD2e>100&&DSSD1x>109.375&&DSSD1x<112.5&&DSSD1y>90.625&&DSSD1y<93.75");
R5 += tree->GetEntries("DSSD1e+DSSD2e>100&&DSSD1e>100&&DSSD2e>100&&DSSD1x>87.5&&DSSD1x<93.75&&DSSD1y>87.5&&DSSD1y<90.625");
R5 += tree->GetEntries("DSSD1e+DSSD2e>100&&DSSD1e>100&&DSSD2e>100&&DSSD1x>106.25&&DSSD1x<112.5&&DSSD1y>87.5&&DSSD1y<90.625");

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

cout << "R1_ratio = " << R1_ratio << endl;
cout << "R2_ratio = " << R2_ratio << endl;
cout << "R3_ratio = " << R3_ratio << endl;
cout << "R4_ratio = " << R4_ratio << endl;
cout << "R5_ratio = " << R5_ratio << endl;
