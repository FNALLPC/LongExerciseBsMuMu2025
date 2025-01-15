#include "common.h"

void task_6_1()
{
    double bdt_min = 0.80;
    
    TFile *fin_wspace = new TFile("wspace.root");
    RooWorkspace *wspace = (RooWorkspace*)fin_wspace->Get("wspace");
    
    RooRealVar *m = wspace->var("m");
    RooArgSet nset(*m);
    RooCategory cate("cate","");
    for(int idx=0; idx<N_Categories; idx++)
        cate.defineType(Form("c%d",idx),idx);
    
    RooSimultaneous *model = new RooSimultaneous("model", "", cate);
    
    // Main POI: Bs->mumu branching fraction
    RooRealVar *BF_bs = new RooRealVar("BF_bs","",3.57E-9,0.,2E-8);
    
    // BF(B+ -> J/psi K+) = (1.010 +- 0.028) E-3 (PDG)
    // BF(J/psi -> mu+mu-) = (5.961 +- 0.033) E-2 (PDG)
    RooRealVar *BF_bu = new RooRealVar("BF_bu","",1.010E-3 * 5.961E-2);
    
    // fs/fu = 0.252 +- 0.012 (PDG) +- 0.015 (energy/pt dependence)
    RooRealVar *fs_over_fu = new RooRealVar("fs_over_fu","",0.252);
    
    RooFormulaVar *N_bs[N_Categories];
    RooAddPdf *pdf_sum[N_Categories];
    for(int idx=0; idx<N_Categories; idx++) {
        
        RooRealVar *Eff_bs = wspace->var(Form("Eff_bs_%d",idx));
        RooRealVar *Eff_bu = wspace->var(Form("Eff_bu_%d",idx));
        RooRealVar *N_bu = wspace->var(Form("N_bu_%d",idx));
        
        N_bs[idx] = new RooFormulaVar(Form("N_bs_%d", idx), "", "@0*@1*@2*@3/@4/@5",
                                      RooArgList(*BF_bs, *N_bu, *fs_over_fu, *Eff_bs, *Eff_bu, *BF_bu));
        RooRealVar *N_peak = wspace->var(Form("N_peak_%d",idx));
        RooRealVar *N_semi = wspace->var(Form("N_semi_%d",idx));
        RooRealVar *N_comb = wspace->var(Form("N_comb_%d",idx));
        
        RooArgList pdf_list;
        pdf_list.add(*wspace->pdf(Form("pdf_bs_%d",idx)));
        pdf_list.add(*wspace->pdf(Form("pdf_peak_%d",idx)));
        pdf_list.add(*wspace->pdf(Form("pdf_semi_%d",idx)));
        pdf_list.add(*wspace->pdf(Form("pdf_comb_%d",idx)));
        
        RooArgList N_list;
        N_list.add(*N_bs[idx]);
        N_list.add(*N_peak);
        N_list.add(*N_semi);
        N_list.add(*N_comb);
        
        pdf_sum[idx] = new RooAddPdf(Form("pdf_sum_%d",idx), "", pdf_list, N_list);
        model->addPdf(*pdf_sum[idx],Form("c%d",idx));
    }

    RooDataSet *rds_data = new RooDataSet("rds_data","",RooArgSet(*m,cate));
    
    TFile *fin_data = new TFile("/eos/uscms/store/user/cmsdas/2025/long_exercises/long-ex-bs-mumu/bmmSoup10.root");
    TTree *tin = (TTree*)fin_data->Get("bmmSoup10_100");
    
    uint cate_t;
    float m_t, bdt_t;
    tin->SetBranchAddress("cate",&cate_t);
    tin->SetBranchAddress("m",&m_t);
    tin->SetBranchAddress("bdt",&bdt_t);
    
    for(int evt=0; evt<tin->GetEntries(); evt++) {
        tin->GetEntry(evt);
        if (bdt_t<=bdt_min) continue;
        cate.setIndex(cate_t);
        m->setVal(m_t);
        rds_data->add(RooArgSet(*m,cate));
    }
    
    model->fitTo(*rds_data,Extended(true),Minos(RooArgSet(*BF_bs)));
    
    TCanvas* canvas = new TCanvas("canvas", "", 1200, 600);
    canvas->Divide(4,2);
    
    RooPlot *frame[N_Categories];
    for(int idx=0; idx<N_Categories; idx++) {
        TVirtualPad* pad = canvas->cd(idx+1);
        pad->SetMargin(0.15,0.06,0.13,0.07);
        
        frame[idx] = m->frame(Title(" "),Bins(25));
        rds_data->plotOn(frame[idx], Cut(Form("cate==%d",idx)), MarkerSize(0.8));
        
        double norm = pdf_sum[idx]->expectedEvents(&nset);
	pdf_sum[idx]->plotOn(frame[idx], Normalization(norm, RooAbsReal::NumEvent), LineWidth(3));
        pdf_sum[idx]->plotOn(frame[idx], Normalization(norm, RooAbsReal::NumEvent), Components(RooArgSet(*wspace->pdf(Form("pdf_bs_%d",idx)))), DrawOption("F"), FillColor(kRed), FillStyle(3365));
        pdf_sum[idx]->plotOn(frame[idx], Normalization(norm, RooAbsReal::NumEvent), Components(RooArgSet(*wspace->pdf(Form("pdf_peak_%d",idx)))), DrawOption("F"), FillColor(kViolet-4), FillStyle(3344));
        pdf_sum[idx]->plotOn(frame[idx], Normalization(norm, RooAbsReal::NumEvent), Components(RooArgSet(*wspace->pdf(Form("pdf_semi_%d",idx)))), DrawOption("L"), LineColor(kGreen-3), LineStyle(2));        
        frame[idx]->GetYaxis()->SetTitleOffset(1.50);
        frame[idx]->GetYaxis()->SetTitle("Entries / 0.04 GeV");
        frame[idx]->GetXaxis()->SetTitleOffset(1.15);
        frame[idx]->GetXaxis()->SetLabelOffset(0.01);
        frame[idx]->GetXaxis()->SetTitle("M(#mu#mu) [GeV]");
        frame[idx]->GetXaxis()->SetTitleSize(0.043);
        frame[idx]->GetYaxis()->SetTitleSize(0.043);
        frame[idx]->Draw();
    }
    
    canvas->Print("task_6_1.pdf");
    canvas->Print("task_6_1.png");
}
