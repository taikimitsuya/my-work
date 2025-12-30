#ifndef BATCH_BOOTSTRAPPING_H
#define BATCH_BOOTSTRAPPING_H

#include <vector>
#include <complex>

#include "tfhe_core.h"
#include "bb_params.h"
#include "batch_framework.h"
#include "batch_ops.h"
#include "hom_dft.h"
#include "bb_utils.h"

namespace bbii {

/**
 * @brief Batch Bootstrapping Key (Algorithm 7.1 Input)
 * 論文 Section 7.2 で定義される、事前計算・パッキング済みのブートストラップ鍵。
 * 秘密鍵 s の回転行列を暗号化したもの。
 */
struct BatchBootstrappingKey {
    // 構造: [v'][2d*log(q)][block_size]
    // v': ブロック分割数 (Section 7.2)
    // 2d*log(q): 行列の列数 (G^{-1} decomposition)
    // block_size: 各ブロックに含まれる行数
    std::vector<std::vector<std::vector<PackedTRGSW>>> keys;

    // キーの初期化やロード関数が必要だが、ここでは構造のみ定義
};

/**
 * @brief Algorithm 7.1: Batch Ring Bootstrapping
 * n 個の LWE 暗号文を一括でブートストラップする。
 *
 * @param inputs 入力 LWE 暗号文のリスト (サイズ n)
 * @param bk バッチブートストラップ鍵
 * @param params パラメータセット
 * @return std::vector<LweSample*> ブートストラップされた LWE 暗号文 (サイズ n)
 */
std::vector<LweSample*> batch_bootstrapping(const std::vector<LweSample*>& inputs,
                                            const BatchBootstrappingKey& bk,
                                            const BBIIParams& params);

/**
 * @brief MSB Extract (Sample Extraction) 
 * RGSW(xi^m) から LWE(m) を抽出する。
 */
LweSample* msb_extract(TRGSW* input_rgsw, const BBIIParams& params);

} // namespace bbii

#endif // BATCH_BOOTSTRAPPING_H