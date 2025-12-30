#ifndef BB_PARAMS_H
#define BB_PARAMS_H

#include <cmath>
#include <stdexcept>
#include <iostream>

// MK-TFHE ライブラリのヘッダをインクルード (パスは環境に合わせて調整してください)
#include "tfhe.h"
#include "tfhe_core.h" 

namespace bbii {

/**
 * @brief Batch Bootstrapping II のパラメータ構造体
 * 論文 Section 7.2  に基づく定義
 */
struct BBIIParams {
    // --- BBII Algorithm Parameters ---

    // 入力LWEの次元。論文では n = 2 * d^rho と設定される [cite: 661, 665]。
    int32_t n;

    // 入力LWEのモジュラス。Section 7.3 では q = O(sqrt(n)) とされる 。
    // 注: TFHEの実装上は Torus32 等を使うため、これは論理的な値として保持する。
    int32_t q;

    // RGSW (Bootstrapping Key) のリング次元 [cite: 660]。
    int32_t N;

    // DFT/逆DFT のパラメータ (2d が基底となる) [cite: 663]。
    // d は2の冪乗である必要がある [cite: 664]。
    int32_t d;

    // 再帰の深さ (Recursive depth) [cite: 665]。
    int32_t rho;

    // バッチ処理可能なスロット数 r [cite: 666]。
    // 論文条件: r > 2 * d * v を満たす必要がある 。
    int32_t r;

    // Algorithm 4.2 (Multiplications over Small Rings) で使用する分割数 v 。
    int32_t v;

    // --- Underlying TFHE Parameters ---
    
    // 入出力LWE暗号文のパラメータ (次元 n)
    LweParams* lwe_params;
    
    // ブートストラップ鍵(RGSW)用のパラメータ (リング次元 N)
    TGswParams* rgsw_params;

    // キー生成や計算に必要な基本的なTFHEパラメータセット
    TFheGateBootstrappingParameterSet* tfhe_params;

    /**
     * @brief コンストラクタ
     * @param d_val パラメータ d (2の冪乗)
     * @param rho_val 再帰深度 rho
     * @param N_val RGSWの次元 N
     * @param v_val 分割数 v
     */
    BBIIParams(int32_t d_val, int32_t rho_val, int32_t N_val, int32_t v_val) 
        : d(d_val), rho(rho_val), N(N_val), v(v_val) {
        
        // n = 2 * d^rho の計算 [cite: 665]
        // pow関数の結果を整数にキャスト
        n = 2 * std::pow(d, rho);

        // q の設定 (簡易実装のため n に近い素数や2の冪乗を設定する等の調整が必要)
        // 論文では q = O(sqrt(n))  だが、実装の容易さのため調整可。
        q = 2 * n; // 仮の設定

        // r (スロット数) の計算。
        // MK-TFHEやFHEWの一般的な設定に基づき、N/2 程度がスロット数の上限となることが多いが、
        // ここでは条件 r > 2dv  を満たす最小の2の冪乗などを設定する。
        // 仮に r = N/2 とする (Nが十分大きいと仮定)。
        r = N / 2;

        if (r <= 2 * d * v) {
            std::cerr << "Warning: r <= 2dv condition violated. " 
                      << "r=" << r << ", 2dv=" << 2*d*v << std::endl;
            // 実装時は適切な r を再設定するか例外を投げる
        }

        // --- TFHEパラメータの初期化 (MK-TFHEのAPIに依存) ---
        // ここでは標準的なセキュリティパラメータ（128bit等）を生成する関数を呼ぶ想定
        // 実際には n, N に基づいてノイズ分散などを計算する必要がある。
        
        // 例: n, N を指定してパラメータを生成 (擬似コード)
        // tfhe_params = new_tfhe_gate_bootstrapping_parameters(n, N, ...);
        // lwe_params = tfhe_params->in_out_params;
        // rgsw_params = tfhe_params->tgsw_params;
    }

    ~BBIIParams() {
        // メモリ解放処理
        // delete_gate_bootstrapping_parameters(tfhe_params);
    }
};

/**
 * @brief テスト用の小さなパラメータセットを生成するファクトリ関数
 * 論文 Section 7.3 のような巨大なパラメータではなく、動作確認用。
 * 例: d=2, rho=3 -> n=16 程度の小さな系。
 */
inline BBIIParams* get_toy_bbii_params() {
    int32_t d = 2;
    int32_t rho = 3; // n = 2 * 2^3 = 16
    int32_t N = 1024; // RGSW次元 (通常はもっと大きいがテスト用)
    int32_t v = 1;    // 最小構成

    return new BBIIParams(d, rho, N, v);
}

} // namespace bbii

#endif // BB_PARAMS_H