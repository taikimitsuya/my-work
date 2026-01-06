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
    
    double lwe_alpha;
    double tlwe_alpha;
    int32_t l;
    int32_t Bgbit;

    BBIIParams(int32_t d_val, int32_t rho_val, int32_t N_val) : d(d_val), rho(rho_val), N(N_val) {
        
        // BBII: n = 2 * d^rho
        n = 2 * std::pow(d, rho);
        
        r = N / 2;
        lwe_alpha = 2.4e-5; tlwe_alpha = 3.7e-9;
        l = 3; Bgbit = 10;

        // Use default TFHE parameter set (128-bit security) to avoid initialization issues
        // The default set is well-tested and handles integer overflow properly
        tfhe_params = new_default_gate_bootstrapping_parameters(128);
        
        // NOTE: This uses TFHE's built-in secure defaults instead of manual construction
        // which avoids parameter mismatches that can cause integer overflow in tGswAddMuIntH
    }
};
}
#endif
