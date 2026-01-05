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
    
    double lwe_alpha; double tlwe_alpha; int32_t l; int32_t Bgbit;

    BBIIParams(int32_t d_val, int32_t rho_val, int32_t N_val) : d(d_val), rho(rho_val), N(N_val) {
        n = 500; r = N / 2;
        lwe_alpha = 2.4e-5; tlwe_alpha = 3.7e-9; l = 3; Bgbit = 10;
        LweParams* lp = new_LweParams(n, lwe_alpha, 0.5); 
        TLweParams* tlp = new_TLweParams(N, 1, tlwe_alpha, 0.5); 
        TGswParams* tgp = new_TGswParams(l, Bgbit, tlp); 
        tfhe_params = new TFheGateBootstrappingParameterSet(10, 1, lp, tgp);
    }
};
}
#endif
