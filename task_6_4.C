#include "common.h"

void build_gaussian_constraint(RooWorkspace* wspace, TString varname, double mean, double sigma)
{
    double min = mean-sigma*6.;
    double max = mean+sigma*6.;
    if (min<0.) min = 0.;
    
    RooRealVar var(varname, "", mean, min, max);
    
    RooRealVar var_mean(Form("%s_mean",varname.Data()),"",mean);
    RooRealVar var_sigma(Form("%s_sigma",varname.Data()),"",sigma);
    RooGaussian gau(Form("%s_gau",varname.Data()), "", var, var_mean, var_sigma);
    wspace->import(gau);
}

void build_lognormal_constraint(RooWorkspace* wspace, TString varname, double mean, double sigma)
{
    double min = mean-sigma*6.;
    double max = mean+sigma*6.;
    if (min<0.) min = 0.;
    
    RooRealVar var(varname, "", mean, min, max);
    
    RooRealVar var_mean(Form("%s_mean",varname.Data()),"",mean);
    RooRealVar var_kappa(Form("%s_kappa",varname.Data()),"",1.+sigma/mean);
    RooLognormal lnn(Form("%s_lnn",varname.Data()), "", var, var_mean, var_kappa);
    wspace->import(lnn);
}

void build_pdf_peak(RooWorkspace *wspace, uint cate, double bdt_min)
{
    vector<TString> decay;
    vector<double> yield, yield_err;
    
    // B0 -> mu mu 
    decay.push_back("bdmmMc");
    yield.push_back(effyield[cate][N_bdmmMc]);
    yield_err.push_back(effyield[cate][dN_bdmmMc]);
    
    // Bs -> hadron hadron (Bs-> kk ; Bs-> k pi; Bs-> pi pi) 
    decay.push_back("bstohhMcBg");
    yield.push_back(effyield[cate][N_bstohhMcBg]);
    yield_err.push_back(effyield[cate][dN_bstohhMcBg]);

    // B -> hadron hadron (B-> kk ; B-> k pi; B-> pi pi)
    decay.push_back("bdtohhMcBg");
    yield.push_back(effyield[cate][N_bdtohhMcBg]);
    yield_err.push_back(effyield[cate][dN_bdtohhMcBg]);
    
    RooRealVar m("m","",4.9,5.9);
    RooRealVar wgt("wgt","",1.,0.,1000.);
    RooDataSet *rds = new RooDataSet("rds","",RooArgSet(m,wgt),"wgt");
    
    double sum_weight = 0.;
    double sum_weight_err = 0.;
    double weight_max = 0.; // get the largest weight value across samples
    for(int proc=0; proc<(int)decay.size(); proc++) {
        TFile *fin = new TFile("/eos/uscms/store/user/cmsdas/2025/long_exercises/long-ex-bs-mumu/" + decay[proc]+".root");
        TTree *tin = (TTree*)fin->Get(decay[proc]);
        double weight = yield[proc]/(double)tin->GetEntries(Form("cate==%d",cate));
        if (weight>weight_max) weight_max = weight;
        delete fin;
    }

    for(int proc=0; proc<(int)decay.size(); proc++) {
      TFile *fin = new TFile("/eos/uscms/store/user/cmsdas/2025/long_exercises/long-ex-bs-mumu/" + decay[proc]+".root");
      TTree *tin = (TTree*)fin->Get(decay[proc]);
        
        unsigned int cate_t;
        float m_t,bdt_t;
        tin->SetBranchAddress("cate",&cate_t);
        tin->SetBranchAddress("m",&m_t);
        tin->SetBranchAddress("bdt",&bdt_t);
        
        double weight = yield[proc]/(double)tin->GetEntries(Form("cate==%d",cate));
	double weight_err = yield_err[proc]/(double)tin->GetEntries(Form("cate==%d",cate))/weight_max;

	std::cout<<"weight "<<weight<<" weight_err"<<weight_err<<std::endl;
	
        for(int evt=0; evt<tin->GetEntries(); evt++) {
            tin->GetEntry(evt);
            if (cate_t!=cate) continue;
            if (bdt_t<=bdt_min) continue;
            m.setVal(m_t);
	    wgt.setVal(weight/weight_max); // rescale the event with largest weight to be 1
            rds->add(RooArgSet(m,wgt),weight/weight_max);
            
            sum_weight += weight;
            sum_weight_err += weight_err; // systematics; linear sum
        }
        delete fin;
    }

    RooKeysPdf pdf(Form("pdf_peak_%d",cate), "", m, *rds,  RooKeysPdf::NoMirror, 2.0);
    
    RooPlot *frame = m.frame(Title(" "),Bins(100));
    rds->plotOn(frame, Name("t_rds"));
    pdf.plotOn(frame, Name("t_pdf"), LineWidth(3));
    
    TCanvas* canvas = new TCanvas("canvas", "", 600, 600);
    canvas->SetMargin(0.15,0.06,0.13,0.07);
    
    frame->GetYaxis()->SetTitleOffset(1.50);
    frame->GetYaxis()->SetTitle("Entries / 0.01 GeV");
    frame->GetXaxis()->SetTitleOffset(1.15);
    frame->GetXaxis()->SetLabelOffset(0.01);
    frame->GetXaxis()->SetTitle("M(#mu#mu) [GeV]");
    frame->GetXaxis()->SetTitleSize(0.043);
    frame->GetYaxis()->SetTitleSize(0.043);
    frame->Draw();
    
    TLegend* leg = new TLegend(0.58,0.77,0.93,0.92);
    leg->SetFillStyle(0);
    leg->SetLineWidth(0);
    leg->SetHeader(Form("Category %d",cate));
    leg->AddEntry(frame->findObject("t_rds"),"Simulation","EP");
    leg->AddEntry(frame->findObject("t_pdf"),"PDF","L");
    leg->Draw();
    
    canvas->Print(Form("task6_4a_%d.pdf",cate));
    
    cout << "Category: " << cate << endl;
    cout << "BDT min: " << bdt_min << endl;
    cout << "Sum of weights: " << sum_weight << " +- " << sum_weight_err << endl;
    
    wspace->import(pdf);
    build_lognormal_constraint(wspace,Form("N_peak_%d",cate),sum_weight,sum_weight_err);
    
    delete frame;
    delete rds;
    delete leg;
    delete canvas;
}

void build_pdf_semi(RooWorkspace *wspace, uint cate, double bdt_min)
{
    vector<TString> decay;
    vector<double> yield, yield_err;
    
    decay.push_back("bskmunuMcBg");
    yield.push_back(effyield[cate][N_bskmunuMcBg]);
    yield_err.push_back(effyield[cate][dN_bskmunuMcBg]);
    
    decay.push_back("bdpimunuMcBg");
    yield.push_back(effyield[cate][N_bdpimunuMcBg]);
    yield_err.push_back(effyield[cate][dN_bdpimunuMcBg]);
    
    decay.push_back("bdpimumuMcBg");
    yield.push_back(effyield[cate][N_bdpimumuMcBg]);
    yield_err.push_back(effyield[cate][dN_bdpimumuMcBg]);
    
    decay.push_back("bupimumuMcBg");
    yield.push_back(effyield[cate][N_bupimumuMcBg]);
    yield_err.push_back(effyield[cate][dN_bupimumuMcBg]);
    
    RooRealVar m("m","",4.9,5.9);
    RooRealVar wgt("wgt","",1.,0.,1000.);
    RooDataSet *rds = new RooDataSet("rds","",RooArgSet(m,wgt),"wgt");
    
    double sum_weight = 0.;
    double sum_weight_err = 0.;
    for(int proc=0; proc<(int)decay.size(); proc++) {
        TFile *fin = new TFile("/eos/uscms/store/user/cmsdas/2025/long_exercises/long-ex-bs-mumu/"+decay[proc]+".root");
        TTree *tin = (TTree*)fin->Get(decay[proc]);
        
        uint cate_t;
        float m_t,bdt_t;
        tin->SetBranchAddress("cate",&cate_t);
        tin->SetBranchAddress("m",&m_t);
        tin->SetBranchAddress("bdt",&bdt_t);
        
        double weight = yield[proc]/(double)tin->GetEntries(Form("cate==%d",cate));
        double weight_err = yield_err[proc]/(double)tin->GetEntries(Form("cate==%d",cate));
        
        for(int evt=0; evt<tin->GetEntries(); evt++) {
            tin->GetEntry(evt);
            if (cate_t!=cate) continue;
            if (bdt_t<=bdt_min) continue;
            m.setVal(m_t);
            wgt.setVal(weight);
            rds->add(RooArgSet(m,wgt),weight);
            
            sum_weight += weight;
            sum_weight_err += weight_err; // systematics; linear sum
        }
        delete fin;
    }
    
    RooKeysPdf pdf(Form("pdf_semi_%d",cate), "", m, *rds,  RooKeysPdf::MirrorLeft, 2.0);
    
    RooPlot *frame = m.frame(Title(" "),Bins(50));
    rds->plotOn(frame, Name("t_rds"));
    pdf.plotOn(frame, Name("t_pdf"), LineWidth(3));
    
    TCanvas* canvas = new TCanvas("canvas", "", 600, 600);
    canvas->SetMargin(0.15,0.06,0.13,0.07);
    
    frame->GetYaxis()->SetTitleOffset(1.50);
    frame->GetYaxis()->SetTitle("Entries / 0.02 GeV");
    frame->GetXaxis()->SetTitleOffset(1.15);
    frame->GetXaxis()->SetLabelOffset(0.01);
    frame->GetXaxis()->SetTitle("M(#mu#mu) [GeV]");
    frame->GetXaxis()->SetTitleSize(0.043);
    frame->GetYaxis()->SetTitleSize(0.043);
    frame->Draw();
    
    TLegend* leg = new TLegend(0.58,0.77,0.93,0.92);
    leg->SetFillStyle(0);
    leg->SetLineWidth(0);
    leg->SetHeader(Form("Category %d",cate));
    leg->AddEntry(frame->findObject("t_rds"),"Simulation","EP");
    leg->AddEntry(frame->findObject("t_pdf"),"PDF","L");
    leg->Draw();
    
    canvas->Print(Form("task6_4b_%d.pdf",cate));
    
    cout << "Category: " << cate << endl;
    cout << "BDT min: " << bdt_min << endl;
    cout << "Sum of weights: " << sum_weight << " +- " << sum_weight_err << endl;
    
    wspace->import(pdf);
    build_lognormal_constraint(wspace,Form("N_semi_%d",cate),sum_weight,sum_weight_err);
    
    delete frame;
    delete rds;
    delete leg;
    delete canvas;
}

void build_pdf_bs(RooWorkspace *wspace, uint cate, double bdt_min)
{
    RooRealVar m("m","",4.9,5.9);
    RooDataSet *rds = new RooDataSet("rds","",RooArgSet(m));
    
    TFile *fin = new TFile("/eos/uscms/store/user/cmsdas/2025/long_exercises/long-ex-bs-mumu/bsmmMc.root");
    TTree *tin = (TTree*)fin->Get("bsmmMc");
        
    uint cate_t;
    float m_t,bdt_t;
    tin->SetBranchAddress("cate",&cate_t);
    tin->SetBranchAddress("m",&m_t);
    tin->SetBranchAddress("bdt",&bdt_t);
    
    int evtcount[2] = {0,0};
    for(int evt=0; evt<tin->GetEntries(); evt++) {
        tin->GetEntry(evt);
        if (cate_t!=cate) continue;
        evtcount[0]++;
        if (bdt_t<=bdt_min) continue;
        evtcount[1]++;
        m.setVal(m_t);
        rds->add(RooArgSet(m));
    }
    double eff = (double)evtcount[1]/(double)evtcount[0] * effyield[cate][Eff_bsmmMc];
    double eff_error = effyield[cate][dEff_bsmmMc]/effyield[cate][Eff_bsmmMc]*eff; // ignore MC statistics error
    delete fin;
    
    RooRealVar bs_mean1(Form("bs_%d_mean1",cate),"",5.37,5.2,5.5);
    RooRealVar bs_mean2(Form("bs_%d_mean2",cate),"",5.37,5.2,5.5);
    RooRealVar bs_sigma1(Form("bs_%d_sigma1",cate),"",0.030,0.005,0.060);
    RooRealVar bs_sigma2(Form("bs_%d_sigma2",cate),"",0.080,0.040,0.200);
    RooRealVar bs_cbalpha(Form("bs_%d_cbalpha",cate),"",1.,0.,4.);
    RooRealVar bs_cbn(Form("bs_%d_cbn",cate),"",1.,0.,4.);
    RooRealVar bs_frac(Form("bs_%d_frac",cate),"",0.7,0.,1.);
    RooGaussian bs_gaus(Form("bs_%d_gaus",cate), "", m, bs_mean1, bs_sigma1);
    RooCBShape bs_cbline(Form("bs_%d_cbline",cate), "", m, bs_mean2, bs_sigma2, bs_cbalpha, bs_cbn);
    RooAddPdf pdf(Form("pdf_bs_%d",cate),"",RooArgList(bs_gaus,bs_cbline),RooArgList(bs_frac));
    pdf.fitTo(*rds);
    
    RooPlot *frame = m.frame(Title(" "),Bins(100));
    rds->plotOn(frame, Name("t_rds"));
    pdf.plotOn(frame, Name("t_pdf"), LineWidth(3));
    
    TCanvas* canvas = new TCanvas("canvas", "", 600, 600);
    canvas->SetMargin(0.15,0.06,0.13,0.07);
    
    frame->GetYaxis()->SetTitleOffset(1.50);
    frame->GetYaxis()->SetTitle("Entries / 0.01 GeV");
    frame->GetXaxis()->SetTitleOffset(1.15);
    frame->GetXaxis()->SetLabelOffset(0.01);
    frame->GetXaxis()->SetTitle("M(#mu#mu) [GeV]");
    frame->GetXaxis()->SetTitleSize(0.043);
    frame->GetYaxis()->SetTitleSize(0.043);
    frame->Draw();
    
    TLegend* leg = new TLegend(0.58,0.77,0.93,0.92);
    leg->SetFillStyle(0);
    leg->SetLineWidth(0);
    leg->SetHeader(Form("Category %d",cate));
    leg->AddEntry(frame->findObject("t_rds"),"Simulation","EP");
    leg->AddEntry(frame->findObject("t_pdf"),"PDF","L");
    leg->Draw();
    
    canvas->Print(Form("task6_4c_%d.pdf",cate));
    
    cout << "Category: " << cate << endl;
    cout << "BDT min: " << bdt_min << endl;
    
    bs_mean1.setConstant(true);
    bs_mean2.setConstant(true);
    bs_sigma1.setConstant(true);
    bs_sigma2.setConstant(true);
    bs_cbalpha.setConstant(true);
    bs_cbn.setConstant(true);
    bs_frac.setConstant(true);
    wspace->import(pdf);
    build_gaussian_constraint(wspace,Form("Eff_bs_%d",cate),eff,eff_error);
    
    delete frame;
    delete rds;
    delete leg;
    delete canvas;
}

void fit_bupsik(RooWorkspace *wspace, uint cate, double bdt_min)
{
    RooRealVar m("m","",5.0,5.8);
    RooRealVar wgt("wgt","",1.,0.,1000.);
    RooDataSet *rds_data = new RooDataSet("rds_data","",RooArgSet(m,wgt),"wgt");
    RooDataSet *rds_mc = new RooDataSet("rds_mc","",RooArgSet(m));
    
    TFile *fin = new TFile("/eos/uscms/store/user/cmsdas/2025/long_exercises/long-ex-bs-mumu/bupsikData.root");
    TTree *tin = (TTree*)fin->Get("bupsikData");
    
    uint cate_t;
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
    
    fin = new TFile("/eos/uscms/store/user/cmsdas/2025/long_exercises/long-ex-bs-mumu/bupsikMc.root");
    tin = (TTree*)fin->Get("bupsikMc");
    
    tin->SetBranchAddress("cate",&cate_t);
    tin->SetBranchAddress("m",&m_t);
    
    int evtcount[2] = {0,0};
    for(int evt=0; evt<tin->GetEntries(); evt++) {
        tin->GetEntry(evt);
        if (cate_t!=cate) continue;
        evtcount[0]++;
        if (m_t<5.0 || m_t>=5.8) continue;
        evtcount[1]++;
        m.setVal(m_t);
        rds_mc->add(RooArgSet(m));
    }
    double eff = (double)evtcount[1]/(double)evtcount[0] * effyield[cate][Eff_bupsikMc];
    double eff_error = effyield[cate][dEff_bupsikMc]/effyield[cate][Eff_bupsikMc]*eff; // ignore MC statistics error
    delete fin;
    
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
    leg1->AddEntry(frame1->findObject("t_rds_mc"),"Simulation","EP");
    leg1->AddEntry(frame1->findObject("t_pdf_sigmc"),"Fit","L");
    leg1->Draw();
    
    canvas1->Print(Form("task6_4d_%d.pdf",cate));
    
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
    
    canvas2->Print(Form("task6_4e_%d.pdf",cate));
    
    cout << "Category: " << cate << endl;
    cout << "Selection efficiency: " << eff << " +- " << eff_error << endl;
    cout << "Observed yield: " << n_sig.getVal() << " +- " << n_sig.getError() << endl;
    
    build_gaussian_constraint(wspace,Form("N_bu_%d",cate),n_sig.getVal(),n_sig.getError());
    wspace->import(RooRealVar(Form("Eff_bu_%d",cate),"",eff));
    
    delete frame1;
    delete frame2;
    delete rds_mc;
    delete rds_data;
    delete leg1;
    delete leg2;
    delete canvas1;
    delete canvas2;
}

void build_pdf_comb(RooWorkspace *wspace, uint cate, double bdt_min)
{
    RooRealVar m("m","",4.9,5.9);
    
    RooRealVar comb_B1(Form("comb_B1_%d",cate), "", 0.5, 0. , 1);
    RooFormulaVar comb_B2(Form("comb_B2_%d",cate), "", "1.-@0", RooArgList(comb_B1));
    RooBernstein pdf(Form("pdf_comb_%d",cate), "", m, RooArgList(comb_B1, comb_B2));
    
    TFile *fin = new TFile("/eos/uscms/store/user/cmsdas/2025/long_exercises/long-ex-bs-mumu/bmmData-blind.root");
    TTree *tin = (TTree*)fin->Get("bmmData");
    
    double N_comb_guess = (double)tin->GetEntries(Form("cate==%d&&bdt>%g&&m>5.45",cate,bdt_min));
    N_comb_guess *= 1.0/0.45; // scale to full mass region
    delete fin;
 
    cout << "Category: " << cate << endl;
    wspace->import(pdf);
    wspace->import(RooRealVar(Form("N_comb_%d",cate),"",N_comb_guess,0.,N_comb_guess*10.));
}

void task_6_4()
{
  gROOT->SetBatch();
  RooWorkspace *wspace = new RooWorkspace("wspace");
    
    for(uint cate=0; cate<N_Categories; cate++) {
        double bdt_min = 0.80;
        build_pdf_peak(wspace,cate,bdt_min);
        build_pdf_semi(wspace,cate,bdt_min);
        build_pdf_bs(wspace,cate,bdt_min);
        build_pdf_comb(wspace,cate,bdt_min);
        fit_bupsik(wspace, cate, bdt_min);
    }
    
    // BF(B+ -> J/psi K+) = (1.010 +- 0.028) E-3 (PDG)
    // BF(J/psi -> mu+mu-) = (5.961 +- 0.033) E-2 (PDG)
    // Total uncertainty is 2.8%
    build_gaussian_constraint(wspace, "BF_bu", 1.010E-3*5.961E-2, 1.010E-3*5.961E-2*0.028);
    
    // fs/fu = 0.252 +- 0.012 (PDG) +- 0.015 (energy/pt dependence)
    build_gaussian_constraint(wspace, "fs_over_fu", 0.252, sqrt(pow(0.012,2)+pow(0.015,2)));
    
    wspace->writeToFile("wspace_constr.root");
    delete wspace;
}
