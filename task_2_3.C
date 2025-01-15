void task_2_3(){

  unsigned int category=0;	
  using namespace RooFit;

  RooRealVar m("m","",5.0,5.8);
  RooRealVar wgt("wgt","",1.,0.,1000.);
  RooDataSet *rds_data = new RooDataSet("rds_data","",RooArgSet(m,wgt),"wgt");

  TFile *fin = new TFile("/eos/uscms/store/user/cmsdas/2025/long_exercises/long-ex-bs-mumu/bupsikData.root");
  TTree *tin = (TTree*)fin->Get("bupsikData");

  unsigned int cate_t;
  float wgt_t,m_t;
  tin->SetBranchAddress("cate",&cate_t);
  tin->SetBranchAddress("wgt",&wgt_t);
  tin->SetBranchAddress("m",&m_t);


  for(int evt=0; evt<tin->GetEntries(); evt++) {
      tin->GetEntry(evt);
      if (cate_t!=category) continue;
      if (m_t<5.0 || m_t>=5.8) continue;
      m.setVal(m_t);
      wgt.setVal(wgt_t);
      rds_data->add(RooArgSet(m,wgt),wgt_t);
     }

  delete fin;


  //category 0 MC fit params
  //sig_frac	  = 0.926022	 +/-  0.00179216	(limited)
  //sig_mean1	  = 5.27971	 +/-  5.56694e-05	(limited)
  //sig_mean2	  = 5.26646	 +/-  0.000814312	(limited)
  //sig_sigma1	  = 0.0156007	 +/-  5.44464e-05	(limited)
  //sig_sigma2	  = 0.0627436	 +/-  0.000724307	(limited)

  //category 1 MC fit params
  //sig_frac	  = 0.854807	 +/-  0.0035439	(limited)
  //sig_mean1	  = 5.27966	 +/-  7.96365e-05	(limited)
  //sig_mean2	  = 5.27309	 +/-  0.000571977	(limited)
  //sig_sigma1	  = 0.0198666	 +/-  9.6366e-05	(limited)
  //sig_sigma2	  = 0.0602973	 +/-  0.000593151	(limited)


  // corrections on the mean and resolution of the signal peak
  RooRealVar sig_shift("sig_shift","",0.,-0.02,0.02);
  RooRealVar sig_scale("sig_scale","",1.,0.8,1.2);

  RooRealVar sigmc_mean1("sigmc_mean1","",5.27971);
  RooRealVar sigmc_mean2("sigmc_mean2","",5.26646);
  RooRealVar sigmc_sigma1("sigmc_sigma1","",0.0156007);
  RooRealVar sigmc_sigma2("sigmc_sigma2","",0.0627436);


  RooAddition sig_mean1("sig_mean1","",RooArgList(sigmc_mean1,sig_shift));
  RooAddition sig_mean2("sig_mean2","",RooArgList(sigmc_mean2,sig_shift));
  RooProduct sig_sigma1("sig_sigma1","",RooArgList(sigmc_sigma1,sig_scale));
  RooProduct sig_sigma2("sig_sigma2","",RooArgList(sigmc_sigma2,sig_scale));

  RooRealVar sig_frac("sig_frac","",0.926022);
  RooGaussian sig_g1("sig_g1","",m,sig_mean1,sig_sigma1);
  RooGaussian sig_g2("sig_g2","",m,sig_mean2,sig_sigma2);
  RooAddPdf pdf_sig("pdf_sig","",RooArgList(sig_g1,sig_g2),RooArgList(sig_frac));

  RooRealVar comb_coeff("comb_coeff","",-1.2,-10.,10.);

  RooExponential pdf_comb("pdf_comb","",m,comb_coeff);


  RooRealVar jpsix_scale("jpsix_scale","",0.02,0.001,0.08);
  RooRealVar jpsix_shift("jpsix_shift","",5.13,5.12,5.16);
  RooGenericPdf pdf_jpsix("pdf_jpsix","","TMath::Erfc((@0-@1)/@2)",RooArgList(m,jpsix_shift,jpsix_scale));
  RooRealVar n_sig("n_sig","",100000,0.,1E8);
  RooRealVar n_comb("n_comb","",80000,0.,1E6);
  RooRealVar n_jpsix("n_jpsix","",20000,0.,1E5);

  RooAddPdf model("model","",RooArgList(pdf_sig,pdf_comb,pdf_jpsix),RooArgList(n_sig,n_comb,n_jpsix));

  model.fitTo(*rds_data, Extended(true), SumW2Error(true));


  //plotting
  RooPlot *frame = m.frame(Title(" "),Bins(80));
  rds_data->plotOn(frame, Name("t_rds_data"));
  model.plotOn(frame, Name("t_model"), LineWidth(3));
  model.plotOn(frame, Name("t_pdf_comb"), Components("pdf_comb"), LineWidth(3), LineStyle(2), LineColor(kGray+1));

  //model.paramOn(frame);


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
  leg->SetHeader(Form("Category %d",category));
  leg->AddEntry(frame->findObject("t_rds_data"),"Real Data","EP");
  leg->AddEntry(frame->findObject("t_pdf_sig"),"Fit","L");
  leg->AddEntry(frame->findObject("t_pdf_comb"),"Combinatorial Bkg.","L");
  leg->Draw();
									          
  canvas->Print("task_2_3.pdf");
  canvas->Print("task_2_3.png");

}
