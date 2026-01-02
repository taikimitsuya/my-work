#ifndef BB_PARAMS_H
#define BB_PARAMS_H
#include <cmath>
#include <vector>
#include <iostream>
#include <stdexcept>
#include <complex>
// ディレクトリなしでインクルード
#include <tfhe.h>
#include <tfhe_core.h>

namespace bbii {
struct BBIIParams {
    int32_t n; int32_t N; int32_t d; int32_t rho; int32_t r;
    TFheGateBootstrappingParameterSet* tfhe_params;
    
    BBIIParams(int32_t d_val, int32_t rho_val, int32_t N_val) : d(d_val), rho(rho_val), N(N_val) {
        // n = 2 * d^rho
        // rho=20 の場合、n = 2,097,152 になります
        n = 2 * std::pow(d, rho); 
        r = N / 2; 
        
        static const int32_t k = 1;
        static const double alpha_lwe = 3.0e-5;
        static const double alpha_bk  = 9.0e-9;
        static const int32_t ks_t = 10;
        static const int32_t ks_basebit = 1;
        
        LweParams* lwe_p = new_LweParams(n, alpha_lwe, 1.0/2.0);
        TLweParams* tlwe_p = new_TLweParams(N, k, alpha_bk, 1.0/2.0);
        TGswParams* tgsw_p = new_TGswParams(3, 10, tlwe_p);
        
        tfhe_params = new TFheGateBootstrappingParameterSet(ks_t, ks_basebit, lwe_p, tgsw_p);
    }
    ~BBIIParams() { delete_gate_bootstrapping_parameters(tfhe_params); }
};

// ★ここを変更しました: (d=2, rho=20, N=1024)
inline BBIIParams* get_test_params() { return new BBIIParams(2, 20, 1024); }

}
#endif
