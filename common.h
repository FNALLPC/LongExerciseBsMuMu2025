#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TAxis.h>
#include <TLegend.h>
#include <TString.h>
#include <RooRealVar.h>
#include <RooCategory.h>
#include <RooDataSet.h>
#include <RooBernstein.h>
#include <RooGaussian.h>
#include <RooLognormal.h>
#include <RooCBShape.h>
#include <RooAddPdf.h>
#include <RooProdPdf.h>
#include <RooAddition.h>
#include <RooFormulaVar.h>
#include <RooProduct.h>
#include <RooArgList.h>
#include <RooExponential.h>
#include <RooGenericPdf.h>
#include <RooKeysPdf.h>
#include <RooHistPdf.h>
#include <RooSimultaneous.h>
#include <RooPlot.h>
#include <RooFitResult.h>
#include <RooWorkspace.h>

using namespace RooFit;
using namespace std;

enum { Category, Eff_bsmmMc, dEff_bsmmMc, Eff_bdmmMc, dEff_bdmmMc, Eff_bupsikMc, dEff_bupsikMc, Eff_bspsiphiMc, dEff_bspsiphiMc, N_bupsikData, dN_bupsikData, N_bspsiphiData, dN_bspsiphiData, N_bsmmMc, dN_bsmmMc, N_bdmmMc, dN_bdmmMc, N_bstohhMcBg, dN_bstohhMcBg, N_bdtohhMcBg, dN_bdtohhMcBg, N_bskmunuMcBg, dN_bskmunuMcBg, N_bdpimunuMcBg, dN_bdpimunuMcBg, N_bdpimumuMcBg, dN_bdpimumuMcBg, N_bupimumuMcBg, dN_bupimumuMcBg, N_Columns };

const int N_Categories = 8;

const double effyield[N_Categories][N_Columns] = {
    {0,
     0.003209738, 1.419501e-05, 0.003186667, 1.453832e-05, 0.001003, 2.563e-06, 0.000435, 8e-07, 299881, 3224, 14242, 250, 13.34416, 0.5941475, 1.613997, 0.0865717, 0.2748195, 0.141854, 0.415722, 0.2112255, 2.908118, 0.9438877, 5.392937, 1.689394, 9.833555, 1.12419, 2.506161, 0.3526048, 
},
    {1,
     0.005369033, 1.836274e-05, 0.005263216, 1.869714e-05, 0.001342, 2.9844e-06, 0.0005476, 9e-07, 439299, 4762, 19345, 330, 24.43867, 1.086063, 2.918615, 0.1563351, 0.5433117, 0.2781636, 0.9562831, 0.4830151, 7.387415, 2.168292, 11.47512, 3.296641, 17.0699, 1.918935, 4.565485, 0.6345283, 
},
    {2,
     0.003793565, 1.550034e-05, 0.003774367, 1.58202e-05, 0.001267, 2.9008e-06, 0.0005476, 9e-07, 378825, 4013, 17685, 302, 15.77188, 0.7009087, 1.911719, 0.1024001, 0.3085595, 0.1605772, 0.4113678, 0.2105327, 4.772823, 1.448018, 6.575227, 1.941222, 12.33439, 1.396932, 2.812489, 0.3945325},
    {3, 0.006140739, 1.972531e-05, 0.006043221, 2.003345e-05, 0.001593, 3.2724e-06, 0.0006532, 1e-06, 476369, 5161, 21716, 365, 25.53419, 1.134107, 3.061358, 0.16391, 0.5235518, 0.2695437, 0.9531907, 0.4827444, 8.05609, 2.330973, 13.99605, 3.857029, 18.81509, 2.104402, 4.548412, 0.6313464},
    {4, 0.004318245, 1.064496e-05, 0.004324804, 1.102165e-05, 0.001446, 1.43173e-06, 0.0006311, 7e-07, 1084880, 11060, 52068, 819, 45.04996, 1.990263, 5.496646, 0.2932101, 0.7494271, 0.3841498, 1.380779, 0.6971112, 9.299132, 2.697586, 19.42981, 5.378907, 32.93601, 3.525271, 8.538031, 1.150372},
    {5, 0.006396152, 1.295791e-05, 0.006298484, 1.330974e-05, 0.00168, 1.5567e-06, 0.0006856, 7e-07, 1233267, 12693, 56125, 884, 65.28905, 2.884334, 7.832524, 0.4178027, 1.202692, 0.6133696, 2.446374, 1.231462, 17.49562, 4.878094, 34.89292, 9.325248, 48.24744, 5.136144, 12.04745, 1.618088},
    {6, 0.004408422, 8.878538e-06, 0.00437946, 8.982701e-06, 0.001539, 1.4439e-06, 0.0006823, 7e-07, 1764249, 17860, 85665, 1325, 70.27134, 3.101676, 8.504716, 0.4533649, 1.232339, 0.6297198, 2.502741, 1.260971, 18.11759, 5.160928, 33.82139, 9.185575, 51.49532, 5.479898, 12.95332, 1.738901}, 
    {7, 0.007038989, 1.122181e-05, 0.006969272, 1.134145e-05, 0.00199, 1.6615e-06, 0.000822, 8e-07, 2093494, 21319, 96148, 1484, 102.9681, 4.54431, 12.42007, 0.6620219, 2.171105, 1.104153, 4.094385, 2.058779, 26.57741, 7.387679, 63.84195, 16.81253, 75.48919, 8.000636, 19.54075, 2.616194},
};
