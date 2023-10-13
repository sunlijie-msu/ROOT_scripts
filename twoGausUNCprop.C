

void twoGausUNCprop(){

  TH1D *hist = new TH1D("h","h",100,0,20);  
  TRandom3 *gen = new TRandom3();
  Double_t A,B,C;

  for(int i=0;i<1e5;i++){
    // Sample A and B
    // gen->Gaus(mean,sigma)
    A = gen->Gaus(5,1);
    B = gen->Gaus(10,2);
    C = B+A;
    hist->Fill(C);
    if(i<5) printf("A: %0.2f  B: %0.2f  C:%0.2f\n",A,B,C);
  }

  hist->Draw();
  gStyle->SetOptFit(1111);
  hist->Fit("gaus");
  // the histogram's center is the estimate of C
  // the histogram's width is the estimate of the error in C
}
