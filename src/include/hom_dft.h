#ifndef HOM_DFT_H
#define HOM_DFT_H

#include <vector>
#include <cmath>
#include <algorithm>

#include "bb_params.h"
#include "batch_framework.h"
#include "batch_ops.h"
#include "bb_utils.h"

namespace bbii {

/**
 * @brief 逆DFT行列の指数行列 M を生成するヘルパー関数
 * M[i][j] = -i * j (mod 2d)
 * 論文 Section 5: DFT^-1 matrix definition based on powers of xi_2d
 * * @param d パラメータ d (行列サイズは 2d x 2d)
 * @return std::vector<std::vector<int32_t>> 指数の行列
 */
std::vector<std::vector<int32_t>> gen_inverse_dft_matrix_exponents(int32_t d);

/**
 * @brief Algorithm 6.3: Homomorphic Inverse DFT (Recursive)
 * Nussbaumer Transform を用いて再帰的に逆DFTを計算する。
 *
 * @param inputs 入力暗号文ベクトル (サイズ (2d)^(rho'-1))
 * @param current_rho 現在の再帰深度 rho'
 * @param current_n 現在のリング次元 n'
 * @param params パラメータセット
 * @return std::vector<PackedTRGSW> 変換後の暗号文ベクトル
 */
std::vector<PackedTRGSW> hom_dft_inverse(const std::vector<PackedTRGSW>& inputs,
                                         int32_t current_rho,
                                         int32_t current_n,
                                         const BBIIParams& params);

} // namespace bbii

#endif // HOM_DFT_H