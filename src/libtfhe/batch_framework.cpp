#include "batch_framework.h"
#include <iostream>

namespace bbii {

// 内部ヘルパー: 単項式 X^shift を生成または乗算する関数 (MK-TFHE依存)
void add_rotated_to_acc(TRGSW* acc, TRGSW* input, int32_t shift, const BBIIParams& params) {
    // 実際の実装では MK-TFHE の trgsw_mul_by_xai 等を使用する
    // acc = acc + (input * X^shift)
    // ここでは概念的な擬似コード
    // trgsw_add_rotated(acc, input, shift, params.rgsw_params);
}

PackedTRGSW batch_pack(const std::vector<TRGSW*>& inputs, 
                       BatchMode target_mode, 
                       const BBIIParams& params) {
    // 論文[cite: 289, 311]: RGSW-Pack requires O(r) RGSW additions.
    
    if (inputs.size() > params.r) {
        throw std::invalid_argument("Input size exceeds batch capacity r");
    }

    // 1. 結果格納用の TRGSW を確保 (初期値 0)
    TRGSW* result = new_trgsw_ciphertext(params.rgsw_params);
    // trgsw_clear(result); // 0埋め

    // 2. 各入力をスロットに対応する基底(単項式)倍して加算
    // テンソル構造のシミュレーション: i番目の要素を X^{i * stride} に配置するイメージ
    int32_t stride = params.N / params.r; 

    for (size_t i = 0; i < inputs.size(); ++i) {
        // inputs[i] * v_i (ここでは v_i = X^{i * stride} と仮定)
        add_rotated_to_acc(result, inputs[i], i * stride, params);
    }

    return PackedTRGSW(result, target_mode);
}

BatchMode get_mult_output_mode(BatchMode m1, BatchMode m2) {
    // 論文[cite: 284]: Multiplication Logic
    // If mode1 = "R12" and mode2 = "R12->R13", then mode3 = "R13".
    
    if ((m1 == BatchMode::R12 && m2 == BatchMode::R12_to_R13) ||
        (m1 == BatchMode::R12_to_R13 && m2 == BatchMode::R12)) {
        return BatchMode::R13;
    }
    
    if ((m1 == BatchMode::R13 && m2 == BatchMode::R13_to_R12) ||
        (m1 == BatchMode::R13_to_R12 && m2 == BatchMode::R13)) {
        return BatchMode::R12;
    }

    // Remark 3.1[cite: 299]: R12 * R12 などは不可
    return BatchMode::None;
}

PackedTRGSW batch_mult(const PackedTRGSW& A, 
                       const PackedTRGSW& B, 
                       const BBIIParams& params) {
    
    // 1. モードチェック
    BatchMode res_mode = get_mult_output_mode(A.mode, B.mode);
    if (res_mode == BatchMode::None) {
        throw std::runtime_error("Invalid batch multiplication modes");
    }

    // 2. RGSW乗算
    // 論文[cite: 296, 313]: Batch-Mult takes O(log lambda) RGSW mults.
    // 
    // 実際の実装上の注意:
    // BBIの理論通りに O(log r) で実装するには、テンソル分解に対応した
    // 特殊なExternal Productが必要ですが、MK-TFHEなどの既存ライブラリを利用する場合、
    // パッキングされた巨大な多項式同士の「通常のRGSW乗算」を行うことで、
    // SIMD的に全スロットの乗算が一括で処理されます。
    // (これが、巨大な次元Nにおける1回の乗算で済むという意味でのBatch化です)
    
    TRGSW* result = new_trgsw_ciphertext(params.rgsw_params);
    
    // trgsw_external_product などを呼び出す (MK-TFHE APIに依存)
    // result = A (*) B
    // trgsw_mult(result, A.cipher, B.cipher, params.rgsw_params);

    return PackedTRGSW(result, res_mode);
}

std::vector<TRGSW*> unpack(const PackedTRGSW& packed_cipher, 
                           const BBIIParams& params) {
    // 論文[cite: 297, 315]: UnPack takes O(r) RGSW multiplications.
    // 通常は Trace map (Tr_{K/K_12}) を用いて成分を抽出します。
    // 実装的には、自己同型(Automorphism)と鍵切り替え(KeySwitching)を用いて
    // 特定の基底成分のみを残し、他をゼロにする操作になります。

    std::vector<TRGSW*> outputs;
    outputs.reserve(params.r);

    for (int i = 0; i < params.r; ++i) {
        // 簡易実装: 
        // 実際には各スロットを分離するためのマスク処理やTrace操作が必要。
        // MK-TFHEでこれを完全に行うには、Automorphism Keyの生成が必要です。
        
        TRGSW* out = new_trgsw_ciphertext(params.rgsw_params);
        
        // 擬似コード: Extract i-th slot
        // apply_trace_map(out, packed_cipher.cipher, i, params);
        
        outputs.push_back(out);
    }

    return outputs;
}

} // namespace bbii