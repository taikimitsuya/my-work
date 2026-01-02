#include "mk_tfhe_structs.h"
#include "mk_params.h"
#include <vector>
#include <cmath>

// 依存ライブラリのインクルード
#include <tfhe.h>
#include <tfhe_core.h>
#include <polynomials.h>
#include <lweparams.h>
#include <lwekey.h>
#include <tgsw.h>

namespace bbii {

// =========================================================
// 外部関数宣言 (mk_ops.cpp で実装されているもの)
// ヘッダファイル (mk_ops.h) があればそれをincludeしてください
// =========================================================
void mk_cmux(MKRLweSample* result, 
             const TGswSampleFFT* bk_fft, 
             const MKRLweSample* in0, 
             const MKRLweSample* in1, 
             int32_t party_id,
             const TFheGateBootstrappingParameterSet* params);

void mk_rlwe_clear(MKRLweSample* result);
void mk_rlwe_copy(MKRLweSample* dst, const MKRLweSample* src);

// =========================================================
// 内部ヘルパー関数
// =========================================================

/**
 * @brief MK-RLWE 多項式回転 (X^{bar} 倍する)
 * * result = src * X^{bar}
 * * 各多項式成分に対して、単項式 X^{bar} を掛けます。
 */
void mk_rlwe_mul_by_xai(MKRLweSample* result, const MKRLweSample* src, int32_t bar, int32_t N) {
    // libtfheの関数を利用: torusPolynomialMulByXai(result, bar, src)
    // bar は 2N を法とする次数
    for (int i = 0; i <= src->k; ++i) {
        torusPolynomialMulByXai(result->parts[i], bar, src->parts[i]);
    }
}

// =========================================================
// Blind Rotate
// =========================================================

/**
 * @brief マルチキー Blind Rotate
 * * LWE暗号文 (bk_input) の位相情報を使って、アキュムレータ (acc) を回転させます。
 * * @param acc      [in/out] 初期化済みのMK-RLWEアキュムレータ (出力で上書き)
 * @param bk_input [in]     ブートストラップ対象のMK-LWE暗号文
 * @param mk_bk    [in]     全員分のブートストラップ鍵
 * @param params   [in]     パラメータ
 */
void mk_blind_rotate(MKRLweSample* acc,
                     const MKLweSample* bk_input,
                     const MKBootstrappingKey* mk_bk,
                     const TFheGateBootstrappingParameterSet* params) {
    
    int32_t k = acc->k; // パーティ数
    int32_t n = mk_bk->n_per_party; // 1人あたりのLWE次元
    int32_t N = acc->N; // RLWE次数
    int32_t _2N = 2 * N;

    // 1. 初期回転: ACC = ACC * X^{-b}
    //    bk_input->sample->b を 2N で離散化して回転量を得る
    //    符号に注意: decrypt = b - sum(as). rotate = -decrypt = -b + sum(as).
    //    TFHEの実装では、初期化で X^{-b} を掛ける。
    
    int32_t bar_b = modSwitchFromTorus32(bk_input->sample->b, _2N);
    // MKRLweSample は構造体なので、一時変数を作らず直接回転させるための一時バッファが必要
    // ここでは acc を直接操作するため、一度コピーが必要
    
    MKRLweSample* temp_acc = new MKRLweSample(k, params);
    mk_rlwe_copy(temp_acc, acc);
    
    // acc = temp_acc * X^{-bar_b}
    mk_rlwe_mul_by_xai(acc, temp_acc, -bar_b, N);

    // 2. メインループ: 各パーティの各ビットで CMUX
    //    sum(a_i * s_i) の部分を回転として追加していく
    
    for (int u = 0; u < k; ++u) { // パーティ u
        for (int i = 0; i < n; ++i) { // ビット i
            
            // LWE暗号文の a ベクトルから、該当する係数を取得
            // bk_input->sample->a はサイズ k*n の配列
            int32_t global_idx = u * n + i;
            Torus32 a_ui = bk_input->sample->a[global_idx];
            
            // 係数が0なら回転不要なのでスキップ (最適化)
            if (a_ui == 0) continue;

            // 回転量 bar_ai
            int32_t bar_ai = modSwitchFromTorus32(a_ui, _2N);
            if (bar_ai == 0) continue;

            // CMUX(BK_{u,i}, ACC * X^{bar_ai}, ACC)
            // BK_{u,i} が 1 なら X^{bar_ai} 回転、0 ならそのまま
            
            // temp_acc = acc * X^{bar_ai} (Trueの場合の入力)
            mk_rlwe_mul_by_xai(temp_acc, acc, bar_ai, N);
            
            // MK-CMUX 実行
            // result = acc (上書き)
            // bk = mk_bk->bk_fft[u][i]
            // in0 = acc (Falseの場合)
            // in1 = temp_acc (Trueの場合)
            mk_cmux(acc, mk_bk->bk_fft[u][i], acc, temp_acc, u, params);
        }
    }

    delete temp_acc;
}

// =========================================================
// Sample Extraction
// =========================================================

/**
 * @brief MK-RLWE から MK-LWE へのサンプル抽出
 * * MK-RLWEアキュムレータの定数項を取り出し、MK-LWE暗号文として整形します。
 * (LWE次元変換: N -> n_total ではなく、RLWEの係数を取り出す操作)
 */
void mk_sample_extract(MKLweSample* output, const MKRLweSample* acc, const LweParams* lwe_params) {
    int32_t k = acc->k;
    int32_t N = acc->N;

    // MK-LWE: (a'_1, ..., a'_k, b')
    // ここで a'_u は N次元ベクトル。outputは k*N 次元の巨大LWEとみなす。
    // ※注意: 入力のLWE次元 n と、抽出後のLWE次元 N (RLWEの次数) は異なります。
    // ブートストラップ後は通常、次元 N のLWEになります。
    // そのため、output の確保サイズは k*N である必要があります。
    
    // a部分の抽出
    // acc->parts[u] は多項式 (coef_0, coef_1, ..., coef_{N-1})
    // LWEのaベクトルとしては、係数の順序や符号に注意が必要だが、
    // 基本的なSample Extractでは定数項を見るため、
    // a_lwe[j] = acc->parts[u] の係数 j (符号反転等の詳細はTFHE仕様に準拠)
    
    // TFHEの標準的な抽出 (SampleExtractIndex 0):
    // LWE b = acc bodyの定数項
    // LWE a_i = acc mask_iの定数項 (または係数展開)
    
    // ここではシンプルに「定数項抽出」を実装します。
    // acc = (m_1, ..., m_k, body)
    // lwe.b = body[0]
    // lwe.a_{u, 0} = m_u[0], lwe.a_{u, j} = -m_u[N-j] ...
    // (RLWEの構造上、a_i * s = m_i * s (poly mul) の定数項は、<a_vec, s_vec> になる)

    output->sample->b = acc->parts[k]->coefsT[0];

    // 各パーティのマスク部分
    for (int u = 0; u < k; ++u) {
        TorusPolynomial* poly = acc->parts[u];
        
        // 次元 N のLWEマスクとして展開
        // output->sample->a のインデックスオフセット
        // ※ outputは 次元 k*N を想定
        int offset = u * N;
        
        // 定数項 (x^0)
        output->sample->a[offset + 0] = poly->coefsT[0];
        
        // その他の項 (x^j corresponds to index N-j with sign change due to X^N = -1)
        for (int j = 1; j < N; ++j) {
            output->sample->a[offset + j] = -poly->coefsT[N - j];
        }
    }
    
    // 分散の更新などは省略（計算困難なため）
    output->sample->current_variance = lwe_params->alpha_min * lwe_params->alpha_min; 
}


// =========================================================
// Main Bootstrapping Wrapper
// =========================================================

/**
 * @brief MK Bootstrapping エントリポイント
 * * 入力MK-LWE (次元 n*k) をブートストラップし、リフレッシュされた MK-LWE (次元 N*k) を返します。
 * * (通常、この後に KeySwitch をして次元を n*k に戻しますが、ここではBS部分のみ)
 * * @param result [out] 結果 (次元 N*k を持つMKLweSample)
 * @param input  [in]  入力 (次元 n*k)
 * @param mk_bk  [in]  全員分のブートストラップ鍵
 * @param test_vector [in] テストベクトル (多項式 v)
 * @param params [in]  パラメータ
 */
void mk_bootstrapping(MKLweSample* result,
                      const MKLweSample* input,
                      const MKBootstrappingKey* mk_bk,
                      Torus32 mu, // テストベクトルとして定数メッセージ mu を埋め込む場合
                      const TFheGateBootstrappingParameterSet* params) {

    int32_t k = input->k;
    
    // 1. アキュムレータの準備
    MKRLweSample* acc = new MKRLweSample(k, params);
    
    // テストベクトル v の作成
    // ここでは簡単のため、定数メッセージ mu を持つ多項式とします (v = mu)
    // 一般的なGate bootstrappingでは v = X^{N/2} * mu などを使います (Look-up table)
    // ここでは mu を定数項にセット
    mk_rlwe_clear(acc);
    acc->parts[k]->coefsT[0] = mu; // bodyの定数項にセット

    // 2. Blind Rotate
    mk_blind_rotate(acc, input, mk_bk, params);

    // 3. Sample Extract
    // outputのLWE次元は k*N である必要があるため、呼び出し元で適切なサイズのresultを渡すこと
    mk_sample_extract(result, acc, params->in_out_params);

    delete acc;
}

} // namespace bbii