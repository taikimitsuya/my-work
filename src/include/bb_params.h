#ifndef BB_PARAMS_H
#define BB_PARAMS_H

#include <cmath>
#include <vector>
#include <iostream>
#include <stdexcept>
#include <complex>

// TFHEのコアヘッダ
#include <tfhe.h>
#include <tfhe_core.h>

namespace bbii {

struct BBIIParams {
    int32_t n;   // Input LWE dim: 2 * d^rho
    int32_t N;   // Ring dim (RGSW): 1024 etc.
    int32_t d;   // Base (e.g. 2)
    int32_t rho; // Recursion depth
    int32_t r;   // Batch slots

    TFheGateBootstrappingParameterSet* tfhe_params;

    // コンストラクタ
    BBIIParams(int32_t d_val, int32_t rho_val, int32_t N_val) 
        : d(d_val), rho(rho_val), N(N_val) {
        
        // BBII: n = 2 * d^rho
        n = 2 * std::pow(d, rho);
        
        // Batch容量 r: Nの約数かつ十分な大きさ
        r = N / 2; 

        // パラメータ生成 (128bitセキュリティ相当の設定例)
        // LWE: n=600くらいが必要だが、テスト用(n=16)のためノイズを小さく設定
        static const int32_t k = 1;
        static const double alpha_lwe = 3.0e-5; // LWE noise
        static const double alpha_bk  = 9.0e-9; // BK noise

        LweParams* lwe_p = new_LweParams(n, alpha_lwe, 1.0/2.0); // 1/2 interval
        TLweParams* tlwe_p = new_TLweParams(N, k, alpha_bk, 1.0/2.0);
        TGswParams* tgsw_p = new_TGswParams(3, 10, tlwe_p); // l=3, Bgbit=10

        tfhe_params = new TFheGateBootstrappingParameterSet(10, 1, lwe_p, tgsw_p);
    }

    ~BBIIParams() {
        delete_gate_bootstrapping_parameters(tfhe_params);
    }
};

// 検証用 Toy Parameters
inline BBIIParams* get_test_params() {
    // d=2, rho=3 => n=16, N=1024
    return new BBIIParams(2, 3, 1024);
}

} // namespace bbii

#endif // BB_PARAMS_H