TFile* _file1 = TFile::Open("F:/e21010/pxct/run0228_0229_0230_LEGe_XtRa_MSD26_152Eu_Z2707_inChamber_vacuum_XtRa_12mm_away_window1us_XtRaCFDdelay_0.2us_for_efficiency_cal.root");

TCanvas* canvaspeak = new TCanvas("canvaspeak", "canvaspeak", 1500, 400);
canvaspeak->cd();//
canvaspeak->SetTopMargin(0.038);
canvaspeak->SetRightMargin(0.025);
canvaspeak->SetLeftMargin(0.066);
canvaspeak->SetBottomMargin(0.17);
gStyle->SetFrameLineWidth(3);
canvaspeak->SetFrameLineWidth(3);

gPad->SetLogy();
hnorth_e->SetLineWidth(3);
hnorth_e->Rebin(4);
hnorth_e->SetLineColor(1);
hnorth_e->SetStats(0);
hnorth_e->SetTitle("");
hnorth_e->GetXaxis()->SetTitle("Energy (keV)");
hnorth_e->GetYaxis()->SetTitle("Counts per 152 eV");
hnorth_e->GetXaxis()->CenterTitle();
hnorth_e->GetYaxis()->CenterTitle();
hnorth_e->GetXaxis()->SetLabelFont(132);
hnorth_e->GetYaxis()->SetLabelFont(132);
hnorth_e->GetXaxis()->SetLabelSize(0.07);
hnorth_e->GetYaxis()->SetLabelSize(0.07);
hnorth_e->GetXaxis()->SetTitleFont(132);
hnorth_e->GetYaxis()->SetTitleFont(132);
hnorth_e->GetXaxis()->SetTitleOffset(1.0);
hnorth_e->GetYaxis()->SetTitleOffset(0.41);
hnorth_e->GetXaxis()->SetTitleSize(0.08);
hnorth_e->GetYaxis()->SetTitleSize(0.08);
hnorth_e->GetYaxis()->SetTickLength(0.015);
hnorth_e->GetXaxis()->SetRangeUser(0, 2000);
hnorth_e->Scale(0.2);
hnorth_e->GetYaxis()->SetRangeUser(1, 4e5);

const Int_t numBins = hnorth_e->GetNbinsX();
const Double_t xMin = hnorth_e->GetXaxis()->GetXmin();
const Double_t xMax = hnorth_e->GetXaxis()->GetXmax();

TH1D* hnorth_e_t1bkg = new TH1D("hnorth_e_t1bkg", "hnorth_e_t1bkg", numBins, xMin, xMax);
TH1D* hnorth_e_t2bkg = new TH1D("hnorth_e_t2bkg", "hnorth_e_t2bkg", numBins, xMin, xMax);
TH1D* hnorth_e_ebkg = new TH1D("hnorth_e_ebkg", "hnorth_e_ebkg", numBins, xMin, xMax);
TH1D* hnorth_e_pureKaKb = new TH1D("hnorth_e_pureKaKb", "hnorth_e_pureKaKb", numBins, xMin, xMax);
TH1D* hnorth_e_KaKbgated = new TH1D("hnorth_e_KaKbgated", "hnorth_e_KaKbgated", numBins, xMin, xMax);

hnorth_e_KaKbgated->SetLineColor(4);
tree->Draw("north_e>>hnorth_e_KaKbgated", "north_e>0&&north_e<2000&&((lege_e>38.8&&lege_e<40.8)||(lege_e>44.8&&lege_e<47.2))&&north_t-lege_t>-600&&north_t-lege_t<350");

hnorth_e_t1bkg->SetLineColor(kGreen + 1);
tree->Draw("north_e>>hnorth_e_t1bkg", "north_e>0&&north_e<2000&&((lege_e>38.8&&lege_e<40.8)||(lege_e>44.8&&lege_e<47.2))&&north_t-lege_t>400&&north_t-lege_t<1000", "");

hnorth_e_t2bkg->SetLineColor(kOrange + 7);
tree->Draw("north_e>>hnorth_e_t2bkg", "north_e>0&&north_e<2000&&((lege_e>38.8&&lege_e<40.8)||(lege_e>44.8&&lege_e<47.2))&&north_t-lege_t>-1000&&north_t-lege_t<-650", "");

hnorth_e_ebkg->SetLineColor(kViolet);
tree->Draw("north_e>>hnorth_e_ebkg", "north_e>0&&north_e<2000&&lege_e>60&&lege_e<64.4&&north_t-lege_t>-600&&north_t-lege_t<350", "");

hnorth_e_KaKbgated->SetLineWidth(3);
hnorth_e_pureKaKb->Add(hnorth_e_KaKbgated, hnorth_e_ebkg, 1, -1);
hnorth_e_pureKaKb->Add(hnorth_e_pureKaKb, hnorth_e_t1bkg, 1, -1);
hnorth_e_pureKaKb->Add(hnorth_e_pureKaKb, hnorth_e_t2bkg, 1, -1);
hnorth_e_pureKaKb->SetLineColor(2);
hnorth_e_pureKaKb->SetLineWidth(2);

hnorth_e->Draw("hist");
//hnorth_e_KaKbgated->Draw("same");
hnorth_e_pureKaKb->Draw("same");
//hnorth_e_t1bkg->Draw("same");
//hnorth_e_t2bkg->Draw("same");
//hnorth_e_ebkg->Draw("same");



TH1D* hnorth_e_tbkg1 = new TH1D("hnorth_e_tbkg1", "hnorth_e_tbkg1", numBins, xMin, xMax);
TH1D* hnorth_e_tbkg2 = new TH1D("hnorth_e_tbkg2", "hnorth_e_tbkg2", numBins, xMin, xMax);
TH1D* hnorth_e_purebeta = new TH1D("hnorth_e_purebeta", "hnorth_e_purebeta", numBins, xMin, xMax);
TH1D* hnorth_e_betagated20 = new TH1D("hnorth_e_betagated20", "hnorth_e_betagated20", numBins, xMin, xMax);
TH1D* hnorth_e_betagated60 = new TH1D("hnorth_e_betagated60", "hnorth_e_betagated60", numBins, xMin, xMax);
TH1D* hnorth_e_betagated100 = new TH1D("hnorth_e_betagated100", "hnorth_e_betagated100", numBins, xMin, xMax);

//tree->Draw("north_e>>hnorth_e_betagated20","north_e>0&&north_e<2000&&msd26_e>20&&msd26_e<1400&&north_t-msd26_t>-500&&north_t-msd26_t<400","");

//tree->Draw("north_e>>hnorth_e_betagated60","north_e>0&&north_e<2000&&msd26_e>60&&msd26_e<1400&&north_t-msd26_t>-500&&north_t-msd26_t<400","");

tree->Draw("north_e>>hnorth_e_betagated100", "north_e>0&&north_e<2000&&msd26_e>100&&msd26_e<1400&&north_t-msd26_t>-500&&north_t-msd26_t<400", "");

tree->Draw("north_e>>hnorth_e_tbkg1", "north_e>0&&north_e<2000&&msd26_e>100&&msd26_e<1400&&north_t-msd26_t>-1000&&north_t-msd26_t<-500", "");

tree->Draw("north_e>>hnorth_e_tbkg2", "north_e>0&&north_e<2000&&msd26_e>100&&msd26_e<1400&&north_t-msd26_t>500&&north_t-msd26_t<900", "");

hnorth_e_purebeta->Add(hnorth_e_betagated100, hnorth_e_tbkg1, 1, -1);
hnorth_e_purebeta->Add(hnorth_e_purebeta, hnorth_e_tbkg2, 1, -1);
hnorth_e_purebeta->SetLineColor(kAzure);
hnorth_e_purebeta->SetLineWidth(2);

hnorth_e->Draw("hist");
//hnorth_e->Scale(0.20);
hnorth_e_purebeta->Rebin(1);
//hnorth_e_betagated20->Draw("same");
//hnorth_e_betagated60->Draw("same");
//hnorth_e_betagated100->Draw("same");
hnorth_e_purebeta->Draw("same");

hnorth_e_pureKaKb->Draw("same");

//hmsd26_e->Draw();
//tree->Draw("msd26_e", "north_e>120&&north_e<124&&msd26_e>0&&msd26_e<2000&&north_t-msd26_t>-500&&north_t-msd26_t<400", "same");
//tree->Draw("msd26_e", "north_e>243&&north_e<246&&msd26_e>0&&msd26_e<2000&&north_t-msd26_t>-500&&north_t-msd26_t<400", "same");
//tree->Draw("msd26_e", "north_e>342&&north_e<346&&msd26_e>0&&msd26_e<2000&&north_t-msd26_t>-500&&north_t-msd26_t<400", "same");

// generate figures for PRC paper
hnorth_e->GetXaxis()->SetRangeUser(100, 400);
hnorth_e->GetYaxis()->SetRangeUser(200, 3e7);
gPad->RedrawAxis();
canvaspeak->SaveAs("F:/e21010/pxct/xgamma_coincidence_152Eu_inChamber_run0228_0229_0230_100-400keV.png");

hnorth_e->GetXaxis()->SetRangeUser(400, 700);
hnorth_e->GetYaxis()->SetRangeUser(50, 1e6);
gPad->RedrawAxis();
canvaspeak->SaveAs("F:/e21010/pxct/xgamma_coincidence_152Eu_inChamber_run0228_0229_0230_400-700keV.png");

hnorth_e->GetXaxis()->SetRangeUser(700, 1000);
hnorth_e->GetYaxis()->SetRangeUser(10, 1e7);
gPad->RedrawAxis();
canvaspeak->SaveAs("F:/e21010/pxct/xgamma_coincidence_152Eu_inChamber_run0228_0229_0230_700-1000keV.png");

hnorth_e->GetXaxis()->SetRangeUser(1000, 1300);
hnorth_e->GetYaxis()->SetRangeUser(1, 4e7);
gPad->RedrawAxis();
canvaspeak->SaveAs("F:/e21010/pxct/xgamma_coincidence_152Eu_inChamber_run0228_0229_0230_1000-1300keV.png");

hnorth_e->GetXaxis()->SetRangeUser(1300, 1600);
hnorth_e->GetYaxis()->SetRangeUser(0.6, 3e7);
gPad->RedrawAxis();
canvaspeak->SaveAs("F:/e21010/pxct/xgamma_coincidence_152Eu_inChamber_run0228_0229_0230_1300-1600keV.png");


