#include "hom_dft.h"
#include "bb_utils.h"
namespace bbii {
std::vector<std::vector<int32_t>> gen_inv_dft_exponents(int32_t dim) {
    std::vector<std::vector<int32_t>> M(dim, std::vector<int32_t>(dim));
    for(int i=0; i<dim; ++i) {
        for(int j=0; j<dim; ++j) {
            int32_t val = -(i*j) % dim;
            if(val < 0) val += dim;
            M[i][j] = val;
        }
    }
    return M;
}
std::vector<PackedTRGSW> hom_dft_inverse(const std::vector<PackedTRGSW>& inputs, int32_t current_rho, const BBIIParams& params) {
    if (current_rho <= 1) return inputs;
    int32_t two_d = 2 * params.d;
    size_t chunk_size = inputs.size() / two_d;
    std::vector<PackedTRGSW> combined;
    for (int i = 0; i < two_d; ++i) {
        std::vector<PackedTRGSW> sub_in(inputs.begin() + i*chunk_size, inputs.begin() + (i+1)*chunk_size);
        auto sub_out = hom_dft_inverse(sub_in, current_rho - 1, params);
        combined.insert(combined.end(), sub_out.begin(), sub_out.end());
    }
    auto rearranged = rearrange(combined, params.d);
    auto M = gen_inv_dft_exponents(two_d);
    auto mat_mult_res = enc_vec_mat_mult(M, rearranged, params);
    
    // 入力側(combined)はもう不要なので削除してOK
    for(auto& p : combined) delete_TGswSample_array(1, p.cipher); 

    // ★修正箇所: ここで mat_mult_res を削除してはいけない！
    // rev_rearranged は mat_mult_res のポインタを参照しているため、
    // ここで消すと rev_rearranged が空っぽのメモリを指すことになる。
    auto rev_rearranged = reverse_rearrange(mat_mult_res, params.d);
    // 削除行: for(auto& p : mat_mult_res) delete_TGswSample_array(1, p.cipher);  <-- これが原因でした

    std::vector<PackedTRGSW> result;
    size_t half = rev_rearranged.size() / 2;
    int32_t rot_factor = params.N / (1 << current_rho);
    for (size_t i = 0; i < half; ++i) {
        PackedTRGSW rotated = batch_anti_rot(rev_rearranged[i+half], rot_factor, params);
        PackedTRGSW res = create_zero_packed(params, rev_rearranged[i].mode);
        trgsw_add_to(res.cipher, rev_rearranged[i].cipher, params);
        trgsw_add_to(res.cipher, rotated.cipher, params);
        result.push_back(res);
        delete_TGswSample_array(1, rotated.cipher);
    }
    
    // 最後に rev_rearranged (＝mat_mult_resの実体) を削除
    for(auto& p : rev_rearranged) delete_TGswSample_array(1, p.cipher);
    
    return reverse_rearrange(result, params.d);
}
}
