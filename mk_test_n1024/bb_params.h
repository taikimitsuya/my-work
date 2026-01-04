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
        // N=1024用の標準的なパラメータ設定
        n = 500; 
        r = N / 2;
        LweParams* lp=new_LweParams(n, 2.4e-5, 0.5); 
        TLweParams* tlp=new_TLweParams(N, 1, 3.7e-9, 0.5); // N=1024 noise
        TGswParams* tgp=new_TGswParams(3, 10, tlp);        // l=3, Bg=1024
        tfhe_params=new TFheGateBootstrappingParameterSet(10, 1, lp, tgp);
    }
};
}
#endif
