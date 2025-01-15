#include "common.h"

void task_4_1()
{
  gROOT->SetBatch();
  unsigned int cate = 0;
    double bdt_min = 0.8;
    
    vector<TString> decay;
    vector<double> yield, yield_err;

    TString input_path = "/eos/uscms/store/user/cmsdas/2025/long_exercises/long-ex-bs-mumu/";
    
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
    rds->weightError(RooAbsData::SumW2);

    double weight_max = 0.; // get the largest weight value across samples
    for(int proc=0; proc<(int)decay.size(); proc++) {
        TFile *fin = new TFile(input_path + decay[proc]+".root");
        TTree *tin = (TTree*)fin->Get(decay[proc]);
        double weight = yield[proc]/(double)tin->GetEntries(Form("cate==%d",cate));
        if (weight>weight_max) weight_max = weight;
        delete fin;
    }

    for(int proc=0; proc<(int)decay.size(); proc++) {
        TFile *fin = new TFile(input_path + decay[proc]+".root");
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
    
    RooKeysPdf pdf("pdf", "", m, *rds,  RooKeysPdf::NoMirror, 2.0);
    
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
    leg->AddEntry(frame->findObject("t_rds"),"Simluation","EP");
    leg->AddEntry(frame->findObject("t_pdf"),"PDF","L");
    leg->Draw();
    
    canvas->Print("task_4_1.pdf");
    canvas->Print("task_4_1.png");
    
    cout << "Category: " << cate << endl;
    cout << "BDT min: " << bdt_min << endl;
    cout << "Sum of weights: " << sum_weight << " +- " << sum_weight_err << endl;
}
