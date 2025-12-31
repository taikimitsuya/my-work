#include "hom_dft.h"
#include "bb_utils.h"

namespace bbii {

// DFT^-1 行列の指数生成: -i*j mod 2d
std::vector<std::vector<int32_t>> gen_inv_dft_exponents(int32_t dim) {
    std::vector<std::vector<int32_t>> M(dim, std::vector<int32_t>(dim));
    for(int i=0; i<dim; ++i) {
        for(int j=0; j<dim; ++j) {
            int32_t val = -(i*j) % dim;
            if(val < 0) val += dim;
            M[i][j] = val; // X^val
        }
    }
    return M;
}

// Algorithm 6.3
std::vector<PackedTRGSW> hom_dft_inverse(const std::vector<PackedTRGSW>& inputs,
                                         int32_t current_rho,
                                         const BBIIParams& params) {
    // Base case
    if (current_rho <= 1) return inputs;

    int32_t two_d = 2 * params.d;
    if (inputs.size() % two_d != 0) throw std::logic_error("Input size mismatch in recursion");
    
    size_t chunk_size = inputs.size() / two_d;

    // 1. Recursive Calls
    std::vector<PackedTRGSW> combined;
    combined.reserve(inputs.size());
    
    for (int i = 0; i < two_d; ++i) {
        std::vector<PackedTRGSW> sub_in(inputs.begin() + i*chunk_size, inputs.begin() + (i+1)*chunk_size);
        
        // Recurse
        auto sub_out = hom_dft_inverse(sub_in, current_rho - 1, params);
        
        combined.insert(combined.end(), sub_out.begin(), sub_out.end());
    }

    // 2. Rearrange (Stride permutation)
    auto rearranged = rearrange(combined, params.d);

    // 3. Matrix Multiplication (Algorithm 5.5)
    auto M = gen_inv_dft_exponents(two_d);
    // メモリ管理のため、結果を受け取る
    auto mat_mult_res = enc_vec_mat_mult(M, rearranged, params);

    // 中間データのメモリ解放
    for(auto& p : combined) { /* inputsのコピーなら不要だが、生成物なら削除必要 */ }
    // (ここでは省略するが、実運用では smart pointer 推奨)

    // 4. Reverse Rearrange
    auto rev_rearranged = reverse_rearrange(mat_mult_res, params.d);
    
    // mat_mult_res 解放
    for(auto& p : mat_mult_res) delete_TRGSW_array(1, p.cipher);

    // 5. Combination (Nussbaumer Logic: C[i] + C[i+d]*X^rot)
    // 厳密には rotation factor は recursion level に依存するが、ここでは d に基づく
    std::vector<PackedTRGSW> result;
    size_t half = rev_rearranged.size() / 2;
    
    int32_t rot_factor = params.N / (1 << current_rho); // 簡易的な回転係数計算

    for (size_t i = 0; i < half; ++i) {
        PackedTRGSW term_lower = rev_rearranged[i];
        PackedTRGSW term_upper = rev_rearranged[i + half]; // Need to rotate

        // term_upper * X^rot
        PackedTRGSW rotated_upper = batch_anti_rot(term_upper, rot_factor, params);

        // res = lower + rotated_upper
        PackedTRGSW res = create_zero_packed(params, term_lower.mode);
        trgsw_add_to(res.cipher, term_lower.cipher, params);
        trgsw_add_to(res.cipher, rotated_upper.cipher, params);

        result.push_back(res);
        
        delete_TRGSW_array(1, rotated_upper.cipher);
    }

    // 6. メモリ整理と出力調整
    // 再帰構造上、サイズを維持するために出力サイズを入力と同じにするか、仕様に合わせる
    // BBIIではサイズが縮小していくステップと維持するステップがあるが、
    // ここでは reverse_rearrange で元の形式に戻す
    
    // rev_rearranged の解放
    for(auto& p : rev_rearranged) delete_TRGSW_array(1, p.cipher);

    // 最終的な並べ替えを行ってリターン
    return reverse_rearrange(result, params.d);
}

} // namespace bbii