//#include "common.h"

void task_2_4()
{
 
   using namespace RooFit;

   unsigned int cate = 1;
    
   RooRealVar m("m","",5.0,5.8);
   RooRealVar wgt("wgt","",1.,0.,1000.);
   RooDataSet *rds_data = new RooDataSet("rds_data","",RooArgSet(m,wgt),"wgt");
   RooDataSet *rds_mc = new RooDataSet("rds_mc","",RooArgSet(m));

   // read data
   TFile *fin = new TFile("/eos/uscms/store/user/cmsdas/2025/long_exercises/long-ex-bs-mumu/bupsikData.root");
   TTree *tin = (TTree*)fin->Get("bupsikData");
    
   unsigned int cate_t;
   float wgt_t,m_t;
   tin->SetBranchAddress("cate",&cate_t);
   tin->SetBranchAddress("wgt",&wgt_t);
   tin->SetBranchAddress("m",&m_t);
    
   for(int evt=0; evt<tin->GetEntries(); evt++) {
     tin->GetEntry(evt);
     if (cate_t!=cate) continue;
     if (m_t<5.0 || m_t>=5.8) continue;
     m.setVal(m_t);
     wgt.setVal(wgt_t);
     rds_data->add(RooArgSet(m,wgt),wgt_t);
   }
   delete fin;
   
   // read MC
   fin = new TFile("/eos/uscms/store/user/cmsdas/2025/long_exercises/long-ex-bs-mumu/bupsikMc.root");
   tin = (TTree*)fin->Get("bupsikMc");
   
   tin->SetBranchAddress("cate",&cate_t);
   tin->SetBranchAddress("m",&m_t);
    
   for(int evt=0; evt<tin->GetEntries(); evt++) {
     tin->GetEntry(evt);
     if (cate_t!=cate) continue;
     if (m_t<5.0 || m_t>=5.8) continue;
     m.setVal(m_t);
     rds_mc->add(RooArgSet(m));
   }
   delete fin;

   // MC PDF
   RooRealVar sigmc_mean1("sigmc_mean1","",5.28,5.2,5.4);
   RooRealVar sigmc_mean2("sigmc_mean2","",5.28,5.2,5.4);
   RooRealVar sigmc_sigma1("sigmc_sigma1","",0.030,0.005,0.060);
   RooRealVar sigmc_sigma2("sigmc_sigma2","",0.080,0.040,0.200);
   RooRealVar sig_frac("sig_frac","",0.9,0.5,1.0);
   RooGaussian sigmc_g1("sig_g1","",m,sigmc_mean1,sigmc_sigma1);
   RooGaussian sigmc_g2("sig_g2","",m,sigmc_mean2,sigmc_sigma2);
   RooAddPdf pdf_sigmc("pdf_sigmc","",RooArgList(sigmc_g1,sigmc_g2),RooArgList(sig_frac));
   
   pdf_sigmc.fitTo(*rds_mc);
   
   RooPlot *frame1 = m.frame(Title(" "),Bins(80));
   rds_mc->plotOn(frame1, Name("t_rds_mc"));
   pdf_sigmc.plotOn(frame1, Name("t_pdf_sigmc"), LineWidth(3));
   
   TCanvas* canvas1 = new TCanvas("canvas1", "", 600, 600);
   canvas1->SetMargin(0.15,0.06,0.13,0.07);
   
   frame1->GetYaxis()->SetTitleOffset(1.50);
   frame1->GetYaxis()->SetTitle("Entries / 0.01 GeV");
   frame1->GetXaxis()->SetTitleOffset(1.15);
   frame1->GetXaxis()->SetLabelOffset(0.01);
   frame1->GetXaxis()->SetTitle("M(#mu#muK) [GeV]");
   frame1->GetXaxis()->SetTitleSize(0.043);
   frame1->GetYaxis()->SetTitleSize(0.043);
   frame1->Draw();
   
   TLegend* leg1 = new TLegend(0.58,0.77,0.93,0.92);
   leg1->SetFillStyle(0);
   leg1->SetLineWidth(0);
   leg1->SetHeader(Form("Category %d",cate));
   leg1->AddEntry(frame1->findObject("t_rds_mc"),"Simluation","EP");
   leg1->AddEntry(frame1->findObject("t_pdf_sigmc"),"Fit","L");
   leg1->Draw();
   
   canvas1->Print("task_2_4a.pdf");
   canvas1->Print("task_2_4a.png");

   // set signal parameters constant to MC fit results
   sigmc_mean1.setConstant(true);
   sigmc_mean2.setConstant(true);
   sigmc_sigma1.setConstant(true);
   sigmc_sigma2.setConstant(true);
   sig_frac.setConstant(true);
   
   RooRealVar sig_shift("sig_shift","",0.,-0.02,0.02);
   RooRealVar sig_scale("sig_scale","",1.,0.8,1.2);
   
   RooAddition sig_mean1("sig_mean1","",RooArgList(sigmc_mean1,sig_shift));
   RooAddition sig_mean2("sig_mean2","",RooArgList(sigmc_mean2,sig_shift));
   RooProduct sig_sigma1("sig_sigma1","",RooArgList(sigmc_sigma1,sig_scale));
   RooProduct sig_sigma2("sig_sigma2","",RooArgList(sigmc_sigma2,sig_scale));
   RooGaussian sig_g1("sig_g1","",m,sig_mean1,sig_sigma1);
   RooGaussian sig_g2("sig_g2","",m,sig_mean2,sig_sigma2);
   RooAddPdf pdf_sig("pdf_sig","",RooArgList(sig_g1,sig_g2),RooArgList(sig_frac));
   
   RooRealVar comb_coeff("comb_coeff","",-1.2,-10.,10.);
   RooExponential pdf_comb("pdf_comb","",m,comb_coeff);
   
   RooRealVar jpsix_scale("jpsix_scale","",0.02,0.001,0.08);
   RooRealVar jpsix_shift("jpsix_shift","",5.13,5.12,5.16);
   RooGenericPdf pdf_jpsix("pdf_jpsix","","TMath::Erfc((@0-@1)/@2)",RooArgList(m,jpsix_shift,jpsix_scale));
   
   double n_comb_guess = rds_data->sumEntries("m>5.4")*2.;
   double n_sig_guess = rds_data->sumEntries("m>5.18&&m<5.38")-n_comb_guess/4.;
   double n_jpsix_guess = rds_data->sumEntries("m<5.18")-n_comb_guess*0.18/0.8;
   
   RooRealVar n_sig("n_sig","",n_sig_guess,0.,rds_data->sumEntries());
   RooRealVar n_comb("n_comb","",n_comb_guess,0.,rds_data->sumEntries());
   RooRealVar n_jpsix("n_jpsix","",n_jpsix_guess,0.,rds_data->sumEntries());
   RooAddPdf model("model","",RooArgList(pdf_sig,pdf_comb,pdf_jpsix),RooArgList(n_sig,n_comb,n_jpsix));
   
   model.fitTo(*rds_data, Extended(true), SumW2Error(true));
   
   RooPlot *frame2 = m.frame(Title(" "),Bins(80));
   rds_data->plotOn(frame2, Name("t_rds_data"));
   model.plotOn(frame2, Name("t_model"), LineWidth(3));
   model.plotOn(frame2, Name("t_pdf_comb"), Components("pdf_comb"), LineWidth(3), LineStyle(2), LineColor(kGray+1));
   
   TCanvas* canvas2 = new TCanvas("canvas2", "", 600, 600);
   canvas2->SetMargin(0.15,0.06,0.13,0.07);
   
   frame2->GetYaxis()->SetTitleOffset(1.50);
   frame2->GetYaxis()->SetTitle("Entries / 0.01 GeV");
   frame2->GetXaxis()->SetTitleOffset(1.15);
   frame2->GetXaxis()->SetLabelOffset(0.01);
   frame2->GetXaxis()->SetTitle("M(#mu#muK) [GeV]");
   frame2->GetXaxis()->SetTitleSize(0.043);
   frame2->GetYaxis()->SetTitleSize(0.043);
   frame2->Draw();
   
   TLegend* leg2 = new TLegend(0.58,0.77,0.93,0.92);
   leg2->SetFillStyle(0);
   leg2->SetLineWidth(0);
   leg2->SetHeader(Form("Category %d",cate));
   leg2->AddEntry(frame2->findObject("t_rds_data"),"Data","EP");
   leg2->AddEntry(frame2->findObject("t_model"),"Fit","L");
   leg2->AddEntry(frame2->findObject("t_pdf_comb"),"Combinatorial bkg.","L");
   leg2->Draw();
   
   canvas2->Print("task_2_4b.pdf");
   canvas2->Print("task_2_4b.png");
}
