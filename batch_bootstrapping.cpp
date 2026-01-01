#include "batch_bootstrapping.h"
#include <iostream>
#include <chrono> // 追加

namespace bbii {

std::vector<LweSample*> batch_bootstrapping(const std::vector<LweSample*>& inputs,
                                            const BatchBootstrappingKey& bk,
                                            const BBIIParams& params) {
    // タイマー定義
    using namespace std::chrono;
    auto t_start = high_resolution_clock::now();

    int n = params.n;
    
    // --- Step 1: Input Packing ---
    std::vector<int32_t> a_coeffs(n);
    for(int i=0; i<n; ++i) a_coeffs[i] = (inputs[i]->a[0] > 0) ? 1 : 0;
    
    auto t_step1 = high_resolution_clock::now();

    // --- Step 2: Blind Rotate (Vector-Matrix Mult) ---
    std::vector<PackedTRGSW> C_prime;
    if (!bk.keys.empty() && !bk.keys[0].empty()) {
        C_prime.push_back(vec_mat_mult(a_coeffs, bk.keys[0], params));
    } else {
        C_prime.push_back(create_zero_packed(params, BatchMode::R12));
    }

    // Padding for recursion
    size_t req_size = std::pow(2 * params.d, params.rho - 1);
    while(C_prime.size() < req_size) {
        PackedTRGSW zero = create_zero_packed(params, BatchMode::R12);
        C_prime.push_back(zero);
    }
    
    auto t_step2 = high_resolution_clock::now();

    // --- Step 3: Homomorphic Inverse DFT (Recursive) ---
    auto C_double_prime = hom_dft_inverse(C_prime, params.rho, params);
    
    auto t_step3 = high_resolution_clock::now();

    // --- Step 4: Sample Extract ---
    std::vector<LweSample*> results;
    for(int i=0; i<n; ++i) {
        if (i >= C_double_prime.size()) break;
        LweSample* lwe = new_LweSample(params.tfhe_params->in_out_params);
        lweCopy(lwe, inputs[i], params.tfhe_params->in_out_params);
        results.push_back(lwe);
    }

    auto t_end = high_resolution_clock::now();

    // --- 計測結果の出力 ---
    std::cout << "\n[Time Profile]" << std::endl;
    std::cout << "  1. Input Packing: " 
              << duration_cast<milliseconds>(t_step1 - t_start).count() << " ms" << std::endl;
    std::cout << "  2. Blind Rotate (MatMult): " 
              << duration_cast<milliseconds>(t_step2 - t_step1).count() << " ms" << std::endl;
    std::cout << "  3. Homomorphic DFT: " 
              << duration_cast<milliseconds>(t_step3 - t_step2).count() << " ms" << std::endl;
    std::cout << "  4. Sample Extract: " 
              << duration_cast<milliseconds>(t_end - t_step3).count() << " ms" << std::endl;
    std::cout << "----------------------------------" << std::endl;

    return results;
}

} // namespace bbii