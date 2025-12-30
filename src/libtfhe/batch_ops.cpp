#include "batch_ops.h"
#include <iostream>

namespace bbii {

// --- Helper Functions (Internal) ---

// 単位行列やゼロ行列に相当するPacked RGSWを生成する関数 (G_0, G_1 in Paper)
PackedTRGSW get_gadget_packed(BatchMode mode, const BBIIParams& params) {
    // 実際には Gadget Matrix G をパッキングした定数暗号文を返す
    // [cite: 365] G0 = RGSW-Pack(G, ..., G, "R12") 等
    TRGSW* c = new_trgsw_ciphertext(params.rgsw_params);
    // ... Gの設定 ...
    return PackedTRGSW(c, mode); 
}

// 置換行列 C_pi を生成するヘルパー (Alg 5.2 用)
PackedTRGSW get_permutation_matrix(const std::vector<int32_t>& pi, 
                                   BatchMode mode, 
                                   const BBIIParams& params) {
    // 置換 pi に対応する「行列」を暗号化して返す。
    // [cite: 497] C_pi = sum v_i^vee * w_{pi(i)}
    TRGSW* c = new_trgsw_ciphertext(params.rgsw_params);
    // ... 実装 ...
    return PackedTRGSW(c, mode);
}


// --- Implementations ---

PackedTRGSW vec_mat_mult(const std::vector<int32_t>& a, 
                         const std::vector<PackedTRGSW>& B, 
                         const BBIIParams& params) {
    // Algorithm 4.1 [cite: 366]
    // Input: a in {0,1}^w, B (pre-processed packed ciphertexts)
    
    int w = a.size();
    if (B.size() != w) throw std::invalid_argument("Dimension mismatch in VecMatMult");

    // Line 2: ACC_0 = G_0 (Initial Accumulator)
    PackedTRGSW acc = get_gadget_packed(BatchMode::R12, params); 

    // Line 3-4: Loop
    for (int i = 0; i < w; ++i) {
        // Line 4: B_i = sum( ... )
        // a_k[i] はビットなので、これは「選択」操作 (MUX) になる。
        // ここでは B[i] が既に「その列に対応するPacked暗号文」であると仮定する簡略化を行う。
        // 厳密には、B[i] は入力ビットによって選択された暗号文である必要がある。
        
        const PackedTRGSW& b_i = B[i]; 

        // Line 5: ACC = Batch-Mult(ACC, B_i)
        acc = batch_mult(acc, b_i, params);
    }

    // Line 6: UnPack is done later or implicitly handled by return type
    return acc;
}

PackedTRGSW mult_small_ring(const std::vector<int32_t>& a, 
                            const std::vector<PackedTRGSW>& B, 
                            const BBIIParams& params) {
    // Algorithm 4.2 [cite: 404]
    
    // Line 1: a' = g^{-1}(coeffs(a))
    // 多項式の係数をビット分解(g^{-1})して、長さ w = 2d * log q のベクトルにする
    std::vector<int32_t> a_prime; 
    // ... decompose(a) ...

    // Line 2: Return VecMatMult
    return vec_mat_mult(a_prime, B, params);
}

PackedTRGSW batch_permute(const PackedTRGSW& C, 
                          const std::vector<int32_t>& permutation, 
                          const BBIIParams& params) {
    // Algorithm 5.2 [cite: 493]
    
    // Line 1: C_pi を生成 (モード変換を含む R12->R13 等)
    // モードは C.mode の双対またはペアになるものを選ぶ
    BatchMode perm_mode = (C.mode == BatchMode::R12) ? BatchMode::R12_to_R13 : BatchMode::R13_to_R12;
    PackedTRGSW c_pi = get_permutation_matrix(permutation, perm_mode, params);

    // Line 2: Batch-Mult
    return batch_mult(C, c_pi, params);
}

PackedTRGSW inv_auto(const PackedTRGSW& C, const BBIIParams& params) {
    // Algorithm 5.3 [cite: 501]
    
    // Line 2: Apply sigma to each entry (in plaintext conceptual level)
    // Line 3: RGSW-KS (Key Switching)
    // MK-TFHE の KeySwitching を呼び出して、自己同型を適用した鍵に切り替える
    TRGSW* result = new_trgsw_ciphertext(params.rgsw_params);
    
    // 擬似コード: trgsw_key_switch(result, C.cipher, params.auto_key);
    
    return PackedTRGSW(result, C.mode);
}

PackedTRGSW batch_anti_rot(const PackedTRGSW& C, 
                           int32_t delta, 
                           const BBIIParams& params) {
    // Algorithm 5.4 [cite: 511]
    
    // 1. Cyclic Rotation pi_delta
    int r = params.r;
    std::vector<int32_t> pi_delta(r);
    for(int i=0; i<r; ++i) pi_delta[i] = (i + delta) % r; // 巡回シフト

    // 2. ACC = Batch-Permute(C, pi_delta)
    PackedTRGSW acc = batch_permute(C, pi_delta, params);

    // 3. ACC' = Inv-Auto(ACC)
    PackedTRGSW acc_prime = inv_auto(acc, params);

    // 4. マスク生成 (Select Upper / Lower triangular parts)
    // ACC_+: 回転して「はみ出さなかった」部分
    // ACC_-: 回転して「はみ出して符号反転が必要な」部分
    // これらを Batch-Mult で選択・合成する
    
    // マスク用の定数暗号文を生成 (実際には事前計算すべき)
    // mask_pos: indices [delta+1, r] are 1, else 0
    // mask_neg: indices [1, delta] are 1, else 0
    PackedTRGSW mask_pos = get_gadget_packed(BatchMode::R12_to_R13, params); // 仮
    PackedTRGSW mask_neg = get_gadget_packed(BatchMode::R12_to_R13, params); // 仮

    // 5. ACC_+ = Batch-Mult(ACC, mask_pos)
    PackedTRGSW acc_plus = batch_mult(acc, mask_pos, params);

    // 6. ACC_- = Batch-Mult(ACC', mask_neg)
    PackedTRGSW acc_minus = batch_mult(acc_prime, mask_neg, params);

    // 7. Combine: C' = ACC_+ + ACC_-
    // RGSWの加算 (MK-TFHE の trgsw_add)
    TRGSW* res_cipher = new_trgsw_ciphertext(params.rgsw_params);
    // trgsw_add(res_cipher, acc_plus.cipher, acc_minus.cipher);

    return PackedTRGSW(res_cipher, C.mode);
}

std::vector<PackedTRGSW> enc_vec_mat_mult(const std::vector<std::vector<int32_t>>& M, 
                                          const std::vector<PackedTRGSW>& C, 
                                          const BBIIParams& params) {
    // Algorithm 5.5 [cite: 535]
    // M: 2d x 2d matrix (powers of xi)
    // C: Vector of Packed Ciphertexts
    
    // Line 1-2: Pre-process C_kj (Packing inputs) -> Assumed done or handled here
    
    // Line 3: ACC_0 = G (Init)
    PackedTRGSW acc = get_gadget_packed(BatchMode::R12, params);
    
    // Line 4: Loop i (rows)
    // Line 5: Loop j (sub-blocks)
    // ... (Loop structures following Alg 5.5) ...
    
    // 論文にある通り、このループ内で以下の手順を繰り返す:
    // 1. Batch-Mult: 現在の暗号文と入力 C_kj を掛ける
    // 2. Batch-Anti-Rot: 係数 M[i,k]/M[i,k+1] に基づいて回転させる [cite: 538]
    // これにより、内積 <b, x> をホーナー法のように再帰的に計算する (Section 5.3 本文参照)

    std::vector<PackedTRGSW> result;
    // ... 実装詳細 ...
    
    return result;
}

} // namespace bbii