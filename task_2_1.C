void task_2_1(){

  using namespace RooFit;

  
  unsigned int cate_t = 0;
  RooRealVar m("m","",5.0,5.8);
  RooRealVar cate("cate","",-1,20);

  // B -> J\psi K MC
  TFile *fin = new TFile("/eos/uscms/store/user/cmsdas/2025/long_exercises/long-ex-bs-mumu/bupsikMc.root");
  TTree *tin = (TTree*)fin->Get("bupsikMc");

  RooDataSet* rds_mc   = new RooDataSet("rds_mc", "rds_mc", tin,  RooArgSet(cate, m));

  rds_mc = (RooDataSet*)rds_mc->reduce(m, Form("m >= 5.0 && m <= 5.8 && cate==%d", cate_t));

  // Construct Double Gaussian Model
  RooRealVar sig_mean1("sig_mean1","",5.28,5.2,5.4);
  RooRealVar sig_mean2("sig_mean2","",5.28,5.2,5.4);
  RooRealVar sig_sigma1("sig_sigma1","",0.030,0.005,0.060);
  RooRealVar sig_sigma2("sig_sigma2","",0.080,0.040,0.200);
  RooRealVar sig_frac("sig_frac","",0.9,0.5,1.0);
  RooGaussian sig_g1("sig_g1","",m,sig_mean1,sig_sigma1);
  RooGaussian sig_g2("sig_g2","",m,sig_mean2,sig_sigma2);
  RooAddPdf pdf_sig("pdf_sig","",RooArgList(sig_g1,sig_g2),RooArgList(sig_frac));

  // fit
  pdf_sig.fitTo(*rds_mc);


  // Plotting
  RooPlot *frame = m.frame(Title(" "),Bins(80));
  rds_mc->plotOn(frame, Name("t_rds_mc"));
  pdf_sig.plotOn(frame, Name("t_pdf_sig"), LineWidth(3));

  pdf_sig.plotOn(frame,Name("sig_g1"),Components(RooArgList(sig_g1)), LineColor(kRed),LineStyle(2), LineWidth(2));
  pdf_sig.plotOn(frame,Name("sig_g2"),Components(RooArgList(sig_g2)), LineColor(kGreen),LineStyle(2), LineWidth(3));



  pdf_sig.paramOn(frame,Layout(0.9,0.6,0.5)); // xmax,xmin,ymax
  frame->getAttText()->SetTextSize(0.02);

  TCanvas* canvas = new TCanvas("canvas", "", 600, 600);
  canvas->SetMargin(0.15,0.06,0.13,0.07);
        
  frame->GetYaxis()->SetTitleOffset(1.50);
  frame->GetYaxis()->SetTitle("Entries / 0.01 GeV");
  frame->GetXaxis()->SetTitleOffset(1.15);
  frame->GetXaxis()->SetLabelOffset(0.01);
  frame->GetXaxis()->SetTitle("M(#mu#muK) [GeV]");
  frame->GetXaxis()->SetTitleSize(0.043);
  frame->GetYaxis()->SetTitleSize(0.043);
  frame->Draw();
        
  TLegend* leg = new TLegend(0.58,0.77,0.93,0.92);
  leg->SetFillStyle(0);
  leg->SetLineWidth(0);
  leg->SetHeader(Form("Category %d",cate_t));
  leg->AddEntry(frame->findObject("t_rds_mc"),"Simluation","EP");
  leg->AddEntry(frame->findObject("t_pdf_sig"),"Fit","L");
  leg->Draw();
            
  canvas->Print("task_2_1.pdf");
  canvas->Print("task_2_1.png");
}

