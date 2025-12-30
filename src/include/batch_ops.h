#ifndef BATCH_OPS_H
#define BATCH_OPS_H

#include <vector>
#include <complex>

#include "bb_params.h"
#include "batch_framework.h"
#include "bb_utils.h"

namespace bbii {

// --- Section 4: Batch Vector-Matrix Multiplication ---

/**
 * @brief Algorithm 4.1: Batch Vector-Matrix Multiplication
 * 平文バイナリベクトル a と、暗号化された行列 B の積を計算する。
 * z = M * a の計算に対応。
 *
 * @param a 平文ベクトル (ビット分解されたベクトル)
 * @param B 論文 Section 4.1 で定義される、前処理済みのPacked TRGSW行列
 * @param params パラメータ
 * @return PackedTRGSW 結果の暗号文ベクトル (連結された状態)
 * [cite: 366]
 */
PackedTRGSW vec_mat_mult(const std::vector<int32_t>& a, 
                         const std::vector<PackedTRGSW>& B, 
                         const BBIIParams& params);

/**
 * @brief Algorithm 4.2: Multiplications over Small(er) Rings
 * 小さなリング Z[xi_2d] 上での乗算を行う。
 * 内部で vec_mat_mult を呼び出す。
 *
 * @param a 平文の多項式 (係数ベクトル)
 * @param B 暗号化された多項式 (前処理済み)
 * @param params パラメータ
 * @return PackedTRGSW 乗算結果
 * [cite: 404]
 */
PackedTRGSW mult_small_ring(const std::vector<int32_t>& a, 
                            const std::vector<PackedTRGSW>& B, 
                            const BBIIParams& params);


// --- Section 5: Batch Homomorphic Anti-Cyclic Rotation & Ops ---

/**
 * @brief Algorithm 5.2: Batch-Permute
 * Packed Ciphertext 内のスロットを置換する。
 *
 * @param C 入力 Packed Ciphertext
 * @param permutation 置換マップ pi (index i -> pi(i))
 * @param params パラメータ
 * @return PackedTRGSW 置換後の暗号文
 * [cite: 493]
 */
PackedTRGSW batch_permute(const PackedTRGSW& C, 
                          const std::vector<int32_t>& permutation, 
                          const BBIIParams& params);

/**
 * @brief Algorithm 5.3: Inv-Auto (Batch Inverse Automorphism)
 * 自己同型写像 xi -> xi^{-1} を適用する (複素共役のような操作)。
 *
 * @param C 入力 Packed Ciphertext
 * @param params パラメータ (評価鍵を含む)
 * @return PackedTRGSW 変換後の暗号文
 * [cite: 501]
 */
PackedTRGSW inv_auto(const PackedTRGSW& C, const BBIIParams& params);

/**
 * @brief Algorithm 5.4: Batch-Anti-Rot
 * Packed Ciphertext に対して反巡回回転 (Anti-cyclic rotation) を行う。
 *
 * @param C 入力 (Mode: R12 or R13)
 * @param delta 回転量 (単項式 xi^delta)
 * @param params パラメータ
 * @return PackedTRGSW 回転後の暗号文
 * [cite: 511]
 */
PackedTRGSW batch_anti_rot(const PackedTRGSW& C, 
                           int32_t delta, 
                           const BBIIParams& params);

/**
 * @brief Algorithm 5.5: RGSW.EncVec-MatMult (Improved Method)
 * DFT行列のような「1の冪根の冪乗」で構成される特殊な行列 M と、
 * 暗号化されたベクトル C の積を計算する。
 *
 * @param M 平文行列 (各要素は xi の指数を表す整数)
 * @param C 入力暗号文ベクトル
 * @param params パラメータ
 * @return std::vector<PackedTRGSW> 結果の暗号文ベクトル
 * [cite: 535]
 */
std::vector<PackedTRGSW> enc_vec_mat_mult(const std::vector<std::vector<int32_t>>& M, 
                                          const std::vector<PackedTRGSW>& C, 
                                          const BBIIParams& params);

} // namespace bbii

#endif // BATCH_OPS_H