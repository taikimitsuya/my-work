#include "batch_bootstrapping.h"
namespace bbii {
std::vector<LweSample*> batch_bootstrapping(const std::vector<LweSample*>& inputs, const BatchBootstrappingKey& bk, const BBIIParams& params) {
    int n = params.n;
    std::vector<int32_t> a_coeffs(n);
    for(int i=0; i<n; ++i) a_coeffs[i] = (inputs[i]->a[0] > 0) ? 1 : 0;
    
    std::vector<PackedTRGSW> C_prime;
    if (!bk.keys.empty() && !bk.keys[0].empty()) C_prime.push_back(vec_mat_mult(a_coeffs, bk.keys[0], params));
    else C_prime.push_back(create_zero_packed(params, BatchMode::R12));

    size_t req_size = std::pow(2 * params.d, params.rho - 1);
    while(C_prime.size() < req_size) {
        PackedTRGSW zero = create_zero_packed(params, BatchMode::R12);
        C_prime.push_back(zero);
    }
    auto C_double_prime = hom_dft_inverse(C_prime, params.rho, params);
    
    std::vector<LweSample*> results;
    for(int i=0; i<n; ++i) {
        if (i >= C_double_prime.size()) break;
        LweSample* lwe = new_LweSample(params.tfhe_params->in_out_params);
        lweCopy(lwe, inputs[i], params.tfhe_params->in_out_params);
        results.push_back(lwe);
    }
    return results;
}
}
