// DFT/IDFT バッチ処理
void mk_homomorphic_dft_batch(std::vector<MKPackedRLWE*>& accs, const std::vector<std::vector<std::complex<double>>>& dft_matrix);
void mk_homomorphic_idft_batch(std::vector<MKPackedRLWE*>& accs, const std::vector<std::vector<std::complex<double>>>& idft_matrix);
// インクルードガード
#ifndef MK_PACKED_OPS_H
#define MK_PACKED_OPS_H

// 標準ライブラリ
#include <complex>
#include <vector>

// プロジェクトヘッダ
#include "bb_params.h"
#include "mk_ops.h"

// DFT/IDFT・Twiddle
void mk_apply_inv_twiddle(MKPackedRLWE* acc, int32_t power, int32_t N, const MKBootstrappingKey* mk_bk, const TFheGateBootstrappingParameterSet* params);
void mk_homomorphic_idft_recursive(std::vector<MKPackedRLWE*>& inputs, int32_t N, const MKBootstrappingKey* mk_bk, const TFheGateBootstrappingParameterSet* params);
// KSK生成: Automorphism後の鍵からKeySwitchingKeyを生成
void mk_fill_automorphism_ksk(BBII_KSKStruct* ksk, const std::vector<MKSecretKey*>& sks, const TFheGateBootstrappingParameterSet* params);
// 必要な関数宣言（未定義エラー回避用）
void mk_rlwe_addTo(MKPackedRLWE*, MKPackedRLWE*);
void mk_rlwe_subTo(MKPackedRLWE*, MKPackedRLWE*);

// Packed External Product
void mk_packed_external_product(MKPackedRLWE* acc, const MKPackedRGSW* rgsw, int32_t pid, const TFheGateBootstrappingParameterSet* params);

// 順列キー生成（雛形）
MKPackedRGSW* mk_create_permutation_key(const std::vector<int>& permutation, const TFheGateBootstrappingParameterSet* params);
// Automorphism用KeySwitchKey生成（雛形）
BBII_KSKStruct* mk_create_automorphism_ksk(int32_t k, int32_t N, const TFheGateBootstrappingParameterSet* params);

// DFT/IDFT定数行列生成
std::vector<std::vector<std::complex<double>>> mk_create_dft_matrix(int N, bool inverse = false);
// 多項式のインデックス反転（Automorphism）
void mk_poly_automorphism(TorusPolynomial* poly);

// 多項式の反転: P(X) → P(X^-1)
void mk_poly_inv_auto_inplace(TorusPolynomial* poly);
void mk_slice(const std::vector<MKPackedRLWE*>& input, std::vector<MKPackedRLWE*>& out_upper, std::vector<MKPackedRLWE*>& out_lower);
void mk_butterfly(MKPackedRLWE* u, MKPackedRLWE* v, const TFheGateBootstrappingParameterSet* params);
void mk_apply_twiddle(MKPackedRLWE* acc, int32_t power, int32_t N, const MKBootstrappingKey* mk_bk, const TFheGateBootstrappingParameterSet* params);
void mk_homomorphic_dft_recursive(std::vector<MKPackedRLWE*>& inputs, int32_t N, const MKBootstrappingKey* mk_bk, const TFheGateBootstrappingParameterSet* params);
void mk_inv_auto(MKPackedRLWE* acc, const BBII_KSKStruct* ksk, const TFheGateBootstrappingParameterSet* params);
void mk_homomorphic_dft(MKPackedRLWE* acc, const std::vector<std::vector<std::complex<double>>>& dft_matrix);
void mk_homomorphic_idft(MKPackedRLWE* acc, const std::vector<std::vector<std::complex<double>>>& idft_matrix);
void mk_batch_anti_rot(MKPackedRLWE* acc, const MKPackedRGSW* perm_key, const BBII_KSKStruct* ksk, const TFheGateBootstrappingParameterSet* params);
void mk_vec_mat_mult(
    MKPackedRLWE* acc,
    const std::vector<MKPackedRGSW*>& bk_list,
    const std::vector<int32_t>& coeffs,
    const TFheGateBootstrappingParameterSet* params
);

// Batch-Permute: permutation keyによる並べ替え（外部積ラッパー）

void mk_batch_permute(MKPackedRLWE* acc, const MKPackedRGSW* perm_key, const TFheGateBootstrappingParameterSet* params);

#endif
