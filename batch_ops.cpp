#include "batch_ops.h"
namespace bbii {
PackedTRGSW vec_mat_mult(const std::vector<int32_t>& a, const std::vector<PackedTRGSW>& B, const BBIIParams& params) {
    if (B.empty()) return create_zero_packed(params, BatchMode::None); // Safety
    PackedTRGSW acc = create_zero_packed(params, B[0].mode);
    for (size_t i = 0; i < a.size(); ++i) {
        if (i >= B.size()) break;
        if (a[i] == 1) trgsw_add_to(acc.cipher, B[i].cipher, params);
    }
    return acc;
}
PackedTRGSW batch_anti_rot(const PackedTRGSW& C, int32_t delta, const BBIIParams& params) {
    TGswSample* res = new_TGswSample_array(1, params.tfhe_params->tgsw_params);
    trgsw_mul_by_xai(res, C.cipher, delta, params);
    return PackedTRGSW(res, C.mode);
}
std::vector<PackedTRGSW> enc_vec_mat_mult(const std::vector<std::vector<int32_t>>& M_exp, const std::vector<PackedTRGSW>& C, const BBIIParams& params) {
    size_t rows = M_exp.size(); size_t cols = C.size();
    if (cols == 0) return {};
    std::vector<PackedTRGSW> result;
    for (size_t i = 0; i < rows; ++i) {
        PackedTRGSW acc = create_zero_packed(params, C[0].mode);
        for (size_t j = 0; j < cols; ++j) {
            PackedTRGSW term = batch_anti_rot(C[j], M_exp[i][j], params);
            trgsw_add_to(acc.cipher, term.cipher, params);
            delete_TGswSample_array(1, term.cipher);
        }
        result.push_back(acc);
    }
    return result;
}
}
