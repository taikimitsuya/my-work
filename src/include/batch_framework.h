#ifndef BATCH_FRAMEWORK_H
#define BATCH_FRAMEWORK_H

#include <vector>
#include <string>
#include <stdexcept>

// MK-TFHE および BBII パラメータ
#include "tfhe_core.h"
#include "bb_params.h"

namespace bbii {

/**
 * @brief Section 3.1 で定義されるパッキングモード 
 * メッセージ空間 R1 に対して、演算用空間 R2, R3 をどう組み合わせるかを定義。
 */
enum class BatchMode {
    R12,            // x in R1 (x) R2
    R13,            // x in R1 (x) R3
    R12_to_R13,     // x in R1 (x) R2^vee (x) R3 (Key-Switching Key like)
    R13_to_R12,     // x in R1 (x) R2 (x) R3^vee
    None            // Invalid or scalar
};

/**
 * @brief バッチ化されたRGSW暗号文
 * 物理的には1つの TRGSW 暗号文だが、論理的には r 個のスロットを持つ。
 */
struct PackedTRGSW {
    TRGSW* cipher;  // MK-TFHE の TRGSW オブジェクトへのポインタ
    BatchMode mode; // 現在のパッキングモード

    PackedTRGSW(TRGSW* c, BatchMode m) : cipher(c), mode(m) {}
};

/**
 * @brief バッチパッキング (RGSW-Pack) [cite: 289]
 * 複数の RGSW 暗号文 (または平文メッセージ) を1つの PackedTRGSW にまとめる。
 * * @param inputs 入力となる r 個の TRGSW 暗号文 (各入力は R1 の要素)
 * @param target_mode 目標とするパッキングモード
 * @param params パラメータ
 * @return PackedTRGSW バッチ化された暗号文
 */
PackedTRGSW batch_pack(const std::vector<TRGSW*>& inputs, 
                       BatchMode target_mode, 
                       const BBIIParams& params);

/**
 * @brief バッチ準同型乗算 (Batch-Mult) [cite: 295]
 * 2つの PackedTRGSW の乗算を行う。
 * 対応するモードの組み合わせのみが許可される (例: R12 * R12->R13 => R13)。
 *
 * @param A 被乗数
 * @param B 乗数
 * @param params パラメータ
 * @return PackedTRGSW 乗算結果
 */
PackedTRGSW batch_mult(const PackedTRGSW& A, 
                       const PackedTRGSW& B, 
                       const BBIIParams& params);

/**
 * @brief アンパッキング (UnPack) [cite: 297]
 * PackedTRGSW から個別のスロットを取り出す。
 * * @param packed_cipher バッチ化された暗号文
 * @param params パラメータ
 * @return std::vector<TRGSW*> 個別の暗号文のリスト
 */
std::vector<TRGSW*> unpack(const PackedTRGSW& packed_cipher, 
                           const BBIIParams& params);

/**
 * @brief モード遷移の妥当性チェック [cite: 284]
 * 2つのモードの積がどのモードになるかを判定する。
 */
BatchMode get_mult_output_mode(BatchMode m1, BatchMode m2);

} // namespace bbii

#endif // BATCH_FRAMEWORK_H