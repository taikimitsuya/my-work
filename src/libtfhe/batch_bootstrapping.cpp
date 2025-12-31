#include "batch_bootstrapping.h"
#include <iostream>

namespace bbii {

std::vector<LweSample*> batch_bootstrapping(const std::vector<LweSample*>& inputs,
                                            const BatchBootstrappingKey& bk,
                                            const BBIIParams& params) {
    int n = params.n;
    
    // Step 1: Input Packing (LWE -> Poly Coeffs)
    // LWEのペイロード(b)は後で位相回転に使う。a項を多項式係数にする。
    std::vector<int32_t> a_coeffs(n);
    std::vector<int32_t> b_values(n);
    for(int i=0; i<n; ++i) {
        // 注: TFHEのLWE a項はTorus32。離散化してint32に丸める
        // 実際には KeySwitching が必要だが、ここでは直接利用
        a_coeffs[i] = (inputs[i]->a[0] > 0) ? 1 : 0; 
        b_values[i] = inputs[i]->b; 
    }

    // Step 2: Vec Mat Mult (Key Product)
    // 実際には Plaintext DFT が必要だが、Toy Exampleでは省略し、直接積をとる
    std::vector<PackedTRGSW> C_prime;
    if (!bk.keys.empty() && !bk.keys[0].empty()) {
        // v=1 のケース: 単一のブロックとの積
        C_prime.push_back(vec_mat_mult(a_coeffs, bk.keys[0], params));
    } else {
        C_prime.push_back(create_zero_packed(params, BatchMode::R12));
    }

    // Step 3: Homomorphic Inverse DFT
    // 入力サイズを再帰ツリーに合わせる (Padding)
    size_t req_size = std::pow(2 * params.d, params.rho - 1);
    while(C_prime.size() < req_size) {
        // Copy last or zero padding
        PackedTRGSW zero = create_zero_packed(params, BatchMode::R12);
        C_prime.push_back(zero);
    }
    
    // Run Recursive DFT
    auto C_double_prime = hom_dft_inverse(C_prime, params.rho, params);

    // Step 4: Extract LWE (Blind Rotate finalization)
    std::vector<LweSample*> results;
    for(int i=0; i<n; ++i) {
        if (i >= C_double_prime.size()) break;

        LweSample* lwe = new_LweSample(params.tfhe_params->in_out_params);
        
        // Phase shift by b_i
        // C''[i] = C''[i] * X^{b_i}
        // ここでは TRGSW 全体を回すのではなく、Extract した LWE を回す方が軽い
        // TRGSW -> LWE extraction
        // TRGSWの[0]番目のTLWEブロックを取り出す(簡易抽出)
        
        TRGSW* final_trgsw = C_double_prime[i].cipher;
        
        // TFHEの sample_extract 相当の処理を手動で
        // TLWE(s * a + e) -> LWE
        // TLweSample* tlwe = &final_trgsw->blocs[0]; // Take first block
        // tLweExtractLweSample(lwe, tlwe, params.tfhe_params->tgsw_params->tlwe_params, params.tfhe_params->in_out_params);
        
        // スタブ: 入力をコピーして返却（フロー確認用）
        // 実際には上記の抽出ロジックが入る
        lweCopy(lwe, inputs[i], params.tfhe_params->in_out_params);
        
        results.push_back(lwe);
    }

    // Clean up C_prime, C_double_prime...
    
    return results;
}

} // namespace bbii