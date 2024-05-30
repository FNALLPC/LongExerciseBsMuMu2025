#include "common.h"

void do_gen_and_fit(double *res)
{
    double bdt_min = 0.8;
    
    TFile *fin_wspace = new TFile("wspace.root");
    RooWorkspace *wspace = (RooWorkspace*)fin_wspace->Get("wspace");
    
    RooRealVar *m = wspace->var("m");
    RooCategory cate("cate","");
    for(int idx=0; idx<N_Categories; idx++)
        cate.defineType(Form("c%d",idx),idx);
    
    RooSimultaneous *model = new RooSimultaneous("model", "", cate);
    
    // Main POI: Bs->mumu branching fraction
    RooRealVar *BF_bs = new RooRealVar("BF_bs","",3.57E-9,0.,3E-8);
    
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
        
        // fix the efficiencies
        Eff_bs->setConstant(true);
        Eff_bu->setConstant(true);
        
        // fix the semi/peak/bu yield
        N_bu->setConstant(true);
        N_peak->setConstant(true);
        N_semi->setConstant(true);
        
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

    RooDataSet *rds_toy = model->generate(RooArgSet(*m,cate),Extended(true));
    
    RooFitResult *res_best = model->fitTo(*rds_toy,Save(true),Extended(true),Minos(RooArgSet(*BF_bs)));
    
    res[0] = BF_bs->getVal();
    res[1] = BF_bs->getError();
    res[2] = BF_bs->getErrorHi();
    res[3] = BF_bs->getErrorLo();
    
    BF_bs->setConstant(true);
    BF_bs->setVal(0.);
    RooFitResult *res_null = model->fitTo(*rds_toy,Save(true),Extended(true));
    
    res[4] = sqrt((res_null->minNll()-res_best->minNll())*2.);
    
    delete rds_toy;
    delete model;
    for(int idx=0; idx<N_Categories; idx++) {
        delete N_bs[idx];
        delete pdf_sum[idx];
    }
    delete BF_bs;
    delete BF_bu;
    delete fs_over_fu;
    delete fin_wspace;
}

void task_6_3()
{
    TH1D *h_mean = new TH1D("h_mean","",50,0.,1E-8);
    TH1D *h_error = new TH1D("h_error","",50,0.,2E-9);
    TH1D *h_signif = new TH1D("h_signif","",50,0.,10.);
    TH1D *h_pull = new TH1D("h_pull","",50,-5.,5.);
    
    for(int iter=0; iter<1000; iter++) {
        cout << "iteration: " << iter << endl;
        double res[5];
        do_gen_and_fit(res);
        
        h_mean->Fill(res[0]);
        h_error->Fill(res[1]);
        h_signif->Fill(res[4]);
        
        double pull = (res[0]-3.57E-9);
        if (pull>0.) pull /= fabs(res[3]);
        if (pull<0.) pull /= fabs(res[2]);
        if (fabs(res[3])>0. && fabs(res[2])>0.) h_pull->Fill(pull);
    }
    
    TCanvas* canvas = new TCanvas("canvas", "", 800, 800);
    canvas->Divide(2,2);
    
    for(auto& hist : {h_mean,h_error,h_signif,h_pull}) {
        hist->GetYaxis()->SetTitleOffset(1.50);
        hist->GetYaxis()->SetTitle("# of toys");
        hist->GetXaxis()->SetTitleOffset(1.15);
        hist->GetXaxis()->SetLabelOffset(0.01);
        hist->GetXaxis()->SetTitleSize(0.043);
        hist->GetYaxis()->SetTitleSize(0.043);
        hist->SetStats(true);
        hist->SetFillColor(41);
    }
    
    canvas->cd(1)->SetMargin(0.15,0.09,0.13,0.07);
    h_mean->GetXaxis()->SetTitle("B(B_{s}#rightarrow#mu#mu) Mean");
    h_mean->Draw();
    canvas->cd(2)->SetMargin(0.15,0.09,0.13,0.07);
    h_error->GetXaxis()->SetTitle("B(B_{s}#rightarrow#mu#mu) Error");
    h_error->Draw();
    canvas->cd(3)->SetMargin(0.15,0.09,0.13,0.07);
    h_signif->GetXaxis()->SetTitle("Significance");
    h_signif->Draw();
    canvas->cd(4)->SetMargin(0.15,0.09,0.13,0.07);
    h_pull->GetXaxis()->SetTitle("Pull");
    h_pull->Fit("gaus","L");
    
    canvas->Print("task6_3.pdf");
    canvas->Print("task6_3.png");
}
