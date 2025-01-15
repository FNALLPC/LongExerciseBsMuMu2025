void task_2_2(){

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
  tin->SetBranchAddress("m",&m_t);
  // event weight for trigger prescale correction
  tin->SetBranchAddress("wgt",&wgt_t);

  for(int evt=0; evt<tin->GetEntries(); evt++) {
      tin->GetEntry(evt);
      if (cate_t!=category) continue;
      if (m_t<5.0 || m_t>=5.8) continue;
      m.setVal(m_t);
      wgt.setVal(wgt_t);
      rds_data->add(RooArgSet(m,wgt),wgt_t);
     }

  delete fin;

  // signal PDF model
  RooRealVar sig_mean1("sig_mean1","",5.27971);
  RooRealVar sig_mean2("sig_mean2","",5.26646);
  RooRealVar sig_sigma1("sig_sigma1","",0.0156007);
  RooRealVar sig_sigma2("sig_sigma2","",0.0627436);
  RooRealVar sig_frac("sig_frac","",0.926022);
  RooGaussian sig_g1("sig_g1","",m,sig_mean1,sig_sigma1);
  RooGaussian sig_g2("sig_g2","",m,sig_mean2,sig_sigma2);
  RooAddPdf pdf_sig("pdf_sig","",RooArgList(sig_g1,sig_g2),RooArgList(sig_frac));
  RooRealVar n_sig("n_sig","",100000,0.,1E8);

  // PDF for combinatorial background (exponential)
  RooRealVar comb_coeff("comb_coeff","",-1.2,-10.,10.);
  RooExponential pdf_comb("pdf_comb","",m,comb_coeff);
  RooRealVar n_comb("n_comb","",80000,0.,1E6);

  // error function for the J\psi + X bkg
  RooRealVar jpsix_scale("jpsix_scale","",0.02,0.001,0.08);
  RooRealVar jpsix_shift("jpsix_shift","",5.13,5.12,5.16);
  RooGenericPdf pdf_jpsix("pdf_jpsix","","TMath::Erfc((@0-@1)/@2)",RooArgList(m,jpsix_shift,jpsix_scale));
  RooRealVar n_jpsix("n_jpsix","",20000,0.,1E5);

  // final PDF model
  RooAddPdf model("model","",RooArgList(pdf_sig,pdf_comb,pdf_jpsix),RooArgList(n_sig,n_comb,n_jpsix));

  model.fitTo(*rds_data, Extended(true), SumW2Error(true));

  // plotting
  RooPlot *frame = m.frame(Title(" "),Bins(80));
  rds_data->plotOn(frame, Name("t_rds_data"));
  model.plotOn(frame, Name("t_model"), LineWidth(3));
  model.plotOn(frame, Name("t_pdf_comb"), Components("pdf_comb"), LineWidth(3), LineStyle(2), LineColor(kGray+1));
  model.paramOn(frame,RooFit::Layout(0.9,0.6,0.6)); // xmax,xmin,ymax
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
  leg->SetHeader(Form("Category %d",category));
  leg->AddEntry(frame->findObject("t_rds_data"),"Real Data","EP");
  leg->AddEntry(frame->findObject("t_pdf_sig"),"Fit","L");
  leg->AddEntry(frame->findObject("t_pdf_comb"),"Combinatorial bkg.","L");
  leg->Draw();
									          
  canvas->Print("task_2_2.pdf");
  canvas->Print("task_2_2.png");

}

