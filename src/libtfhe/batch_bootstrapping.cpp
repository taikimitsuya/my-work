#include "batch_bootstrapping.h"
#include <iostream>

namespace bbii {

// --- Internal Helper: Plaintext DFT ---

// 平文多項式 a(X) の係数を受け取り、DFT変換後のベクトルを返す
// a_i = DFT(a) [cite: 670]
std::vector<std::vector<int32_t>> plaintext_dft(const std::vector<int32_t>& poly_coeffs, 
                                                const BBIIParams& params) {
    // 簡易実装:
    // 実際には FFT/NTT アルゴリズムを用いて、多項式 a(X) を
    // 1の根 xi 上で評価した値のベクトルに変換する。
    // ここでは、Algorithm 4.1 (VecMatMult) の入力形式に合わせて
    // ビット分解 (g^{-1}) 済みのベクトルを生成する処理まで含むとする。
    
    // 戻り値: (2d)^(rho-1) 個の要素。各要素は g^{-1} 分解されたベクトル。
    std::vector<std::vector<int32_t>> dft_result;
    
    // ... FFT implementation ...
    // ... Bit Decomposition ...

    // ダミー実装: サイズ合わせのみ
    int32_t size = std::pow(2 * params.d, params.rho - 1);
    int32_t decomp_size = 2 * params.d * 32; // log(q) approx 32
    dft_result.resize(size, std::vector<int32_t>(decomp_size, 0));

    return dft_result;
}


// --- Main Algorithm ---

std::vector<LweSample*> batch_bootstrapping(const std::vector<LweSample*>& inputs,
                                            const BatchBootstrappingKey& bk,
                                            const BBIIParams& params) {
    // Algorithm 7.1 Implementation [cite: 683, 689]

    int n = params.n;
    if (inputs.size() != n) {
        throw std::invalid_argument("Input size must match parameter n");
    }

    // Step 1: LWE Packing [cite: 653]
    // n 個の LWE 暗号文を 1つの RLWE 暗号文 (a, b) に変換する。
    // b(X) - a(X) * z(X) = message_poly(X)
    // ここでは a(X) の係数が必要。
    std::vector<int32_t> a_coeffs(n);
    std::vector<int32_t> b_values(n); // Step 7 で使用

    for (int i = 0; i < n; ++i) {
        // LWE の a パート (マスク) を集約して多項式 a(X) を構成
        // 実際のパッキング法 (Micciancio-Sorrell 18等) に依存するが、
        // 単純には i番目のLWEのマスクの一部やModSwitch後の値を使う。
        
        // 仮: LWEのa項の第1成分などを係数とする
        a_coeffs[i] = inputs[i]->a[0]; // 簡易的なマッピング
        b_values[i] = inputs[i]->b;
    }

    // Step 2: DFT(a) [cite: 689]
    // 平文多項式 a(X) を DFT 領域へ変換
    auto a_dft_vectors = plaintext_dft(a_coeffs, params);

    // Step 3 & 4: Vector-Matrix Multiplication Loop [cite: 689]
    // v' = (2d)^(rho-1) / v
    int32_t total_blocks = std::pow(2 * params.d, params.rho - 1);
    int32_t v = params.v;
    int32_t v_prime = total_blocks / v;

    std::vector<PackedTRGSW> C_prime_flat;
    C_prime_flat.reserve(total_blocks);

    for (int i = 0; i < v_prime; ++i) {
        // 部分ベクトル (a_i) の抽出 (U_i に対応)
        // サイズ v のブロック
        
        // VecMatMult の呼び出し
        // 入力: 平文ベクトル a (のDFT結果) と 暗号化された行列 B
        // Note: Algorithm 7.1 Step 4 では "VecMatMult" とあるが、
        // 実際には Algorithm 4.1 を複数回呼ぶか、バッチ処理されたバージョンを呼ぶ。
        
        // ここでは v 個の出力を得る想定
        for (int k = 0; k < v; ++k) {
            int idx = i * v + k;
            
            // a_vec: DFT結果の idx 番目の要素 (ビット分解済み)
            const auto& a_vec = a_dft_vectors[idx];
            
            // B_key: 対応するブートストラップ鍵ブロック
            // bk.keys[i] は [2d*logq] サイズのベクトルで、各要素は PackedTRGSW
            const auto& B_key = bk.keys[i]; 

            // Algorithm 4.1 [cite: 366]
            PackedTRGSW res = vec_mat_mult(a_vec, B_key, params);
            
            C_prime_flat.push_back(res);
        }
    }

    // Step 5: (Implicitly done by pushing to C_prime_flat)
    // C'_{(i-1)v + j} = C_ij

    // Step 6: Homomorphic Inverse DFT [cite: 644]
    // Recursive Nussbaumer Transform
    // 入力は (2d)^(rho-1) 個の PackedTRGSW
    std::vector<PackedTRGSW> C_double_prime = hom_dft_inverse(C_prime_flat, 
                                                              params.rho, 
                                                              params.n, // Initial n
                                                              params);

    // Step 7: Blind Rotate (Phase Shift) [cite: 689]
    // For i in [n], C''[i] = C''[i] * xi_q^{b_i}
    // C_double_prime は PackedTRGSW のリスト (論理的には n 個の RGSW)
    // ここでバッチ展開 (UnPack) して個別に回転させるか、
    // Packed のまま回転させる (Batch-Anti-Rot) 必要がある。
    
    // BBIIの文脈では、C_double_prime の各スロットが LWE サンプルに対応するため、
    // UnPack して個別に b_i を足すのが自然。
    
    std::vector<LweSample*> bootstrapped_lwes(n);

    // 全ての PackedTRGSW をアンパックしてフラットな RGSW リストにする
    int output_idx = 0;
    for (const auto& packed_c : C_double_prime) {
        // UnPack [cite: 297]
        std::vector<TRGSW*> unpacked_rgsws = unpack(packed_c, params);
        
        for (TRGSW* c : unpacked_rgsws) {
            if (output_idx >= n) break;

            // Multiply by X^{b_i} (Phase shift)
            // LWE b値 を加算
            // trgsw_mul_by_xai(c, c, b_values[output_idx], params.rgsw_params);

            // Step 8: MSB Extract 
            // RGSW から LWE を抽出
            bootstrapped_lwes[output_idx] = msb_extract(c, params);
            
            output_idx++;
        }
    }

    return bootstrapped_lwes;
}

LweSample* msb_extract(TRGSW* input_rgsw, const BBIIParams& params) {
    // MSB Extract Implementation 
    // RGSW は (RLWE, RLWE) の構成。
    // RGSW(m) * TestVector(0) -> RLWE(m) -> SampleExtract -> LWE(m)
    // または RGSW の特定の列 (通常は l-1 番目) が RLWE 暗号文そのもの。
    
    // TFHE の trgsw_sample_extract_index 相当の処理
    LweSample* output_lwe = new_lwe_sample(params.lwe_params);
    
    // 擬似コード: 
    // trgsw_extract_lwe_sample(output_lwe, input_rgsw, params.tfhe_params);
    
    return output_lwe;
}

} // namespace bbii