#include "batch_ops.h"
#include <iostream>

namespace bbii {

// Algorithm 4.1: Vector-Matrix Mult
PackedTRGSW vec_mat_mult(const std::vector<int32_t>& a, 
                         const std::vector<PackedTRGSW>& B, 
                         const BBIIParams& params) {
    // Initialize Accumulator with 0
    PackedTRGSW acc = create_zero_packed(params, B[0].mode);

    for (size_t i = 0; i < a.size(); ++i) {
        if (i >= B.size()) break;
        if (a[i] == 1) {
            // ACC += B[i]
            trgsw_add_to(acc.cipher, B[i].cipher, params);
        }
        // if a[i] == 0, do nothing
    }
    return acc;
}

// Algorithm 5.4 Wrapper: Batch Anti-Cyclic Rotation
PackedTRGSW batch_anti_rot(const PackedTRGSW& C, int32_t delta, const BBIIParams& params) {
    // Create result container
    TRGSW* res_cipher = new_TRGSW_array(1, params.tfhe_params->tgsw_params);
    
    // 多項式回転 (X^delta)
    trgsw_mul_by_xai(res_cipher, C.cipher, delta, params);

    return PackedTRGSW(res_cipher, C.mode);
}

// Algorithm 5.5: Encrypted Vector-Matrix Mult with Exponents
std::vector<PackedTRGSW> enc_vec_mat_mult(const std::vector<std::vector<int32_t>>& M_exp, 
                                          const std::vector<PackedTRGSW>& C, 
                                          const BBIIParams& params) {
    size_t rows = M_exp.size();
    size_t cols = C.size();
    if (cols == 0) return {};

    std::vector<PackedTRGSW> result;
    result.reserve(rows);

    for (size_t i = 0; i < rows; ++i) {
        // ACC = 0
        PackedTRGSW acc = create_zero_packed(params, C[0].mode);
        
        for (size_t j = 0; j < cols; ++j) {
            // Term = C[j] * X^{M[i][j]}
            PackedTRGSW term = batch_anti_rot(C[j], M_exp[i][j], params);
            
            // ACC += Term
            trgsw_add_to(acc.cipher, term.cipher, params);
            
            // Clean up term memory (Allocated in batch_anti_rot)
            delete_TRGSW_array(1, term.cipher);
        }
        result.push_back(acc);
    }
    return result;
}

} // namespace bbii