#include "common.h"

void task_4_2()
{

    using namespace RooFit;

    unsigned int cate = 0;
    double bdt_min = 0.8;
    
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
      TFile *fin = new TFile("/eos/user/c/cmsdas/2024/long-ex-bph/"+decay[proc]+".root");
      TTree *tin = (TTree*)fin->Get(decay[proc]);
      
      unsigned int cate_t;
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
    
    RooKeysPdf pdf("pdf", "", m, *rds,  RooKeysPdf::MirrorLeft, 2.0);
    
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
    leg->AddEntry(frame->findObject("t_rds"),"Simluation","EP");
    leg->AddEntry(frame->findObject("t_pdf"),"PDF","L");
    leg->Draw();
    
    canvas->Print("task_4_2.pdf");
    canvas->Print("task_4_2.png");
    
    cout << "Category: " << cate << endl;
    cout << "BDT min: " << bdt_min << endl;
    cout << "Sum of weights: " << sum_weight << " +- " << sum_weight_err << endl;
}
