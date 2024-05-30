#include "common.h"

void task_3_3()
{
    using namespace RooFit;	
    uint cate = 2;
    double pt_min = 20.0;
    double pt_max = 30.0;
    double eta_min = 0.;
    double eta_max = 2.5;
    
    /*
     From task3_1.cc
     Category: 2
     pt range: 20, 30
     |eta| range: 0, 2.5
     Selection efficiency: 0.32762 +- 0.0014842
     Observed yield: 106740 +- 782.594     

     From task3_2.cc
     Category: 2
     pt range: 20, 30
     |eta| range: 0, 2.5
     Selection efficiency: 0.39079 +- 0.00154296
     Observed yield: 7593.83 +- 207.071     */

    double ret_bupsik[4], ret_bspsiphi[4];
    
    ret_bupsik[0] = 0.32762;
    ret_bupsik[1] = 0.0014842/ret_bupsik[0];
    ret_bupsik[2] = 106740;
    ret_bupsik[3] = 782.594/ret_bupsik[2];
    
    ret_bspsiphi[0] = 0.39079;
    ret_bspsiphi[1] = 0.00154296/ret_bspsiphi[0];
    ret_bspsiphi[2] = 7593.83;
    ret_bspsiphi[3] = 207.071/ret_bspsiphi[2];
    
    double baseeff_bupsik = 0.001267;
    double baseeff_bspsiphi = 0.0005476;
    
    double BF_bupsik = 1.010E-3;
    double BF_bupsik_err = 0.028/1.010;
    
    double BF_bspsiphi = 1.08E-3 * 0.492;
    double BF_bspsiphi_err = sqrt(pow(0.08/1.08,2) + pow(0.5/49.2,2));
    
    double fs = ret_bspsiphi[2]/ret_bspsiphi[0]/baseeff_bspsiphi/BF_bspsiphi;
    double fu = ret_bupsik[2]/ret_bupsik[0]/baseeff_bupsik/BF_bupsik;
    
    double fs_over_fu = fs/fu;
    double fs_over_fu_err = fs_over_fu*sqrt(pow(ret_bspsiphi[1],2)+pow(ret_bspsiphi[3],2)+pow(ret_bupsik[1],2)+pow(ret_bupsik[3],2));
    double fs_over_fu_err_common = fs_over_fu*sqrt(pow(BF_bupsik_err,2)+pow(BF_bspsiphi_err,2));
    
    cout << "Category: " << cate << endl;
    cout << "pt range: " << pt_min << ", " << pt_max << endl;
    cout << "|eta| range: " << eta_min << ", " << eta_max << endl;
    cout << "fs_over_fu: " << fs_over_fu << " +- " << fs_over_fu_err << " +- " << fs_over_fu_err_common << endl;
}
