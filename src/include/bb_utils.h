#ifndef BB_UTILS_H
#define BB_UTILS_H

#include <vector>
#include <complex>
#include <cstdint>
#include <cmath>
#include <algorithm>
#include <stdexcept>

// TFHEの型定義を利用
#include "tfhe_core.h" 

namespace bbii {

/**
 * @brief 反巡回回転後のインデックスと符号を保持する構造体
 * RLWE/RGSWにおいて X^k を乗算する際の挙動を表現する
 */
struct RotatedIndex {
    int32_t index; // 新しいインデックス
    int32_t sign;  // 符号 (+1 or -1)
};

/**
 * @brief 複素数の定義 (DFT行列生成用)
 */
using Complex = std::complex<double>;

// --- Math Helpers ---

/**
 * @brief nが2の冪乗か判定
 */
bool is_power_of_two(int32_t n);

/**
 * @brief log2を計算 (整数)
 */
int32_t log2_int(int32_t n);

/**
 * @brief 1のn乗根 (Roots of Unity) のテーブルを生成
 * Nussbaumer変換やDFT行列の構築に使用 
 */
std::vector<Complex> get_roots_of_unity(int32_t n);

/**
 * @brief 反巡回回転 (Anti-cyclic Rotation) のインデックス計算
 * Z_q[X]/(X^N + 1) における X^shift * X^idx の結果
 * 
 * * @param idx 現在のインデックス (0 <= idx < N)
 * @param shift 回転量 (正の整数)
 * @param N リング次元
 * @return RotatedIndex 新しいインデックスと符号反転の有無
 */
RotatedIndex get_anti_cyclic_index(int32_t idx, int32_t shift, int32_t N);

// --- Permutation / Rearrangement Helpers (Section 6.2) ---

/**
 * @brief Section 6.2 "Rearr" アルゴリズムのインデックス置換
 * m -> 2d への表現変換。係数ベクトルのストライド変換を行う。
 * [cite: 627]
 * * @tparam T データ型 (Torus32, RGSW Ciphertextなど)
 * @param input 入力ベクトル (次元 m/2)
 * @param d ブロックサイズパラメータ
 * @return std::vector<T> 並べ替えられたベクトル
 */
template <typename T>
std::vector<T> rearrange(const std::vector<T>& input, int32_t d);

/**
 * @brief Section 6.2 "Rev-Rearr" アルゴリズム
 * 2d -> m への表現変換。Rearrの逆操作。
 * [cite: 628]
 * * @tparam T データ型
 * @param input 入力ベクトル (次元 m/(2d))
 * @param d ブロックサイズパラメータ
 * @return std::vector<T> 元の順序に戻されたベクトル
 */
template <typename T>
std::vector<T> reverse_rearrange(const std::vector<T>& input, int32_t d);

/**
 * @brief ビット反転 (FFT/NTT用)
 */
uint32_t reverse_bits(uint32_t x, int32_t bits);

} // namespace bbii

#include "bb_utils.tpp" // テンプレート実装を含む

#endif // BB_UTILS_H