#include "common.h"

void task_3_2()
{

    using namespace RooFit;	
    uint cate = 2;
    double pt_min = 20.0;
    double pt_max = 30.0;
    double eta_min = 0.;
    double eta_max = 2.5;
    
    RooRealVar m("m","",5.0,5.8);
    RooRealVar wgt("wgt","",1.,0.,1000.);
    RooDataSet *rds_data = new RooDataSet("rds_data","",RooArgSet(m,wgt),"wgt");
    RooDataSet *rds_mc = new RooDataSet("rds_mc","",RooArgSet(m));
    
    TFile *fin = new TFile("/eos/user/c/cmsdas/2024/long-ex-bph/bspsiphiData.root");
    TTree *tin = (TTree*)fin->Get("bspsiphiData");
    
    uint cate_t;
    float wgt_t,m_t,pt_t,eta_t;
    tin->SetBranchAddress("cate",&cate_t);
    tin->SetBranchAddress("wgt",&wgt_t);
    tin->SetBranchAddress("m",&m_t);
    tin->SetBranchAddress("pt",&pt_t);
    tin->SetBranchAddress("eta",&eta_t);
    
    for(int evt=0; evt<tin->GetEntries(); evt++) {
        tin->GetEntry(evt);
        if (cate_t!=cate) continue;
        if (m_t<5.0 || m_t>=5.8) continue;
        if (pt_t<pt_min || pt_t>=pt_max) continue;
        if (fabs(eta_t)<eta_min || fabs(eta_t)>=eta_max) continue;
        m.setVal(m_t);
        wgt.setVal(wgt_t);
        rds_data->add(RooArgSet(m,wgt),wgt_t);
    }
    delete fin;
    
    fin = new TFile("/eos/user/c/cmsdas/2024/long-ex-bph/bspsiphiMc.root");
    tin = (TTree*)fin->Get("bspsiphiMc");
    
    tin->SetBranchAddress("cate",&cate_t);
    tin->SetBranchAddress("m",&m_t);
    tin->SetBranchAddress("pt",&pt_t);
    tin->SetBranchAddress("eta",&eta_t);
    
    int evtcount[2] = {0,0};
    for(int evt=0; evt<tin->GetEntries(); evt++) {
        tin->GetEntry(evt);
        if (cate_t!=cate) continue;
        evtcount[0]++;
        if (m_t<5.0 || m_t>=5.8) continue;
        if (pt_t<pt_min || pt_t>=pt_max) continue;
        if (fabs(eta_t)<eta_min || fabs(eta_t)>=eta_max) continue;
        evtcount[1]++;
        m.setVal(m_t);
        rds_mc->add(RooArgSet(m));
    }
    double eff = (double)evtcount[1]/(double)evtcount[0];
    double eff_error = sqrt(eff*(1.-eff)/(double)evtcount[0]);
    delete fin;
    
    RooRealVar sigmc_mean1("sigmc_mean1","",5.37,5.2,5.5);
    RooRealVar sigmc_mean2("sigmc_mean2","",5.37,5.2,5.5);
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
    frame1->GetXaxis()->SetTitle("M(#mu#muKK) [GeV]");
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
    
    canvas1->Print("task3_2a.pdf");
    
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
    
    double n_comb_guess = rds_data->sumEntries("m>5.46||m<5.26")*0.8/0.6;
    double n_sig_guess = rds_data->sumEntries("m>5.26&&m<5.46")-n_comb_guess/4.;
    
    RooRealVar n_sig("n_sig","",n_sig_guess,0.,rds_data->sumEntries());
    RooRealVar n_comb("n_comb","",n_comb_guess,0.,rds_data->sumEntries());
    RooAddPdf model("model","",RooArgList(pdf_sig,pdf_comb),RooArgList(n_sig,n_comb));
    
    model.fitTo(*rds_data, SumW2Error(true));
    
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
    frame2->GetXaxis()->SetTitle("M(#mu#muKK) [GeV]");
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
    
    canvas2->Print("task3_2b.pdf");
    
    cout << "Category: " << cate << endl;
    cout << "pt range: " << pt_min << ", " << pt_max << endl;
    cout << "|eta| range: " << eta_min << ", " << eta_max << endl;
    cout << "Selection efficiency: " << eff << " +- " << eff_error << endl;
    cout << "Observed yield: " << n_sig.getVal() << " +- " << n_sig.getError() << endl;
}
