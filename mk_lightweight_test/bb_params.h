#ifndef BB_PARAMS_H
#define BB_PARAMS_H
#include <cmath>
#include <vector>
#include <tfhe.h>
#include <tfhe_core.h>
namespace bbii {
struct BBIIParams {
    int32_t n; int32_t N; int32_t d; int32_t rho; int32_t r;
    TFheGateBootstrappingParameterSet* tfhe_params;
    BBIIParams(int32_t d_val, int32_t rho_val, int32_t N_val) : d(d_val), rho(rho_val), N(N_val) {
        // 軽量化のためパラメータ調整
        // BBII: n = 2 * d^rho
        n = 2 * std::pow(d, rho); // 強制的に小さく
        r = N / 2;
        LweParams* lp=new_LweParams(n,3.0e-5,0.5); 
        TLweParams* tlp=new_TLweParams(N,1,9.0e-9,0.5);
        TGswParams* tgp=new_TGswParams(2,10,tlp); // l=2に削減
        tfhe_params=new TFheGateBootstrappingParameterSet(10,1,lp,tgp);
    }
};
}
#endif
