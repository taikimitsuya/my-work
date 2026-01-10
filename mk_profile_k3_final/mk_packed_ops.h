#ifndef MK_PACKED_OPS_H
#define MK_PACKED_OPS_H

#include <complex>
#include <vector>
#include "bb_params.h"
#include "mk_ops.h"

namespace bbii {
void mk_apply_inv_twiddle(MKPackedRLWE* acc, int32_t power, int32_t N, const MKBootstrappingKey* mk_bk, const TFheGateBootstrappingParameterSet* params);
void mk_homomorphic_idft_recursive(
	std::vector<MKPackedRLWE*>& inputs,
	int32_t N,
	const MKBootstrappingKey* mk_bk,
	const TFheGateBootstrappingParameterSet* params
);
// KSK生成: Automorphism後の鍵からKeySwitchingKeyを生成
void mk_fill_automorphism_ksk(MKKeySwitchKey* ksk, const std::vector<MKSecretKey*>& sks, const TFheGateBootstrappingParameterSet* params);
struct MKPackedRLWE;
struct MKPackedRGSW;
struct MKKeySwitchKey;

// Packed External Product
void mk_packed_external_product(MKPackedRLWE* acc, const MKPackedRGSW* rgsw, int32_t pid, const TFheGateBootstrappingParameterSet* params);

// 順列キー生成（雛形）
MKPackedRGSW* mk_create_permutation_key(const std::vector<int>& permutation, const TFheGateBootstrappingParameterSet* params);
// Automorphism用KeySwitchKey生成（雛形）
MKKeySwitchKey* mk_create_automorphism_ksk(int32_t k, int32_t N, const TFheGateBootstrappingParameterSet* params);

// DFT/IDFT定数行列生成
std::vector<std::vector<std::complex<double>>> mk_create_dft_matrix(int N, bool inverse = false);
// 多項式のインデックス反転（Automorphism）
void mk_poly_automorphism(TorusPolynomial* poly);

// 多項式の反転: P(X) → P(X^-1)
void mk_poly_inv_auto_inplace(TorusPolynomial* poly);
void mk_slice(const std::vector<bbii::MKPackedRLWE*>& input, std::vector<bbii::MKPackedRLWE*>& out_upper, std::vector<bbii::MKPackedRLWE*>& out_lower);
void mk_butterfly(bbii::MKPackedRLWE* u, bbii::MKPackedRLWE* v, const TFheGateBootstrappingParameterSet* params);
void mk_apply_twiddle(bbii::MKPackedRLWE* acc, int32_t power, int32_t N, const bbii::MKBootstrappingKey* mk_bk, const TFheGateBootstrappingParameterSet* params);
void mk_homomorphic_dft_recursive(std::vector<bbii::MKPackedRLWE*>& inputs, int32_t N, const bbii::MKBootstrappingKey* mk_bk, const TFheGateBootstrappingParameterSet* params);
void mk_inv_auto(MKPackedRLWE* acc, const MKKeySwitchKey* ksk, const TFheGateBootstrappingParameterSet* params);
void mk_homomorphic_dft(MKPackedRLWE* acc, const std::vector<std::vector<std::complex<double>>>& dft_matrix);
void mk_homomorphic_idft(MKPackedRLWE* acc, const std::vector<std::vector<std::complex<double>>>& idft_matrix);
void mk_batch_anti_rot(MKPackedRLWE* acc, const MKPackedRGSW* perm_key, const MKKeySwitchKey* ksk, const TFheGateBootstrappingParameterSet* params);
void mk_vec_mat_mult(
 MKPackedRLWE* acc,
 const std::vector<MKPackedRGSW*>& bk_list,
 const std::vector<int32_t>& coeffs,
 const TFheGateBootstrappingParameterSet* params
);

// Batch-Permute: permutation keyによる並べ替え（外部積ラッパー）
void mk_batch_permute(MKPackedRLWE* acc, const MKPackedRGSW* perm_key, const TFheGateBootstrappingParameterSet* params);

// Inv-Auto: インデックス反転（雛形）
void mk_inv_auto(MKPackedRLWE* acc, const MKKeySwitchKey* ksk);

// 再帰的DFT（雛形）


// ヘルパー関数
void mk_slice(const std::vector<MKPackedRLWE*>& input, std::vector<MKPackedRLWE*>& out_upper, std::vector<MKPackedRLWE*>& out_lower);
void mk_butterfly(MKPackedRLWE* u, MKPackedRLWE* v, const TFheGateBootstrappingParameterSet* params);
void mk_apply_twiddle(MKPackedRLWE* acc, int32_t power, int32_t N, const MKBootstrappingKey* mk_bk, const TFheGateBootstrappingParameterSet* params);

} // namespace bbii



#include <complex>
#include <vector>
#include "bb_params.h"
#include "mk_ops.h"

namespace bbii {
struct MKPackedRLWE;
struct MKPackedRGSW;
struct MKKeySwitchKey;

// Packed External Product
void mk_packed_external_product(MKPackedRLWE* acc, const MKPackedRGSW* rgsw, int32_t pid, const TFheGateBootstrappingParameterSet* params);

// 順列キー生成（雛形）
MKPackedRGSW* mk_create_permutation_key(const std::vector<int>& permutation, const TFheGateBootstrappingParameterSet* params);
// Automorphism用KeySwitchKey生成（雛形）
MKKeySwitchKey* mk_create_automorphism_ksk(int32_t k, int32_t N, const TFheGateBootstrappingParameterSet* params);

// DFT/IDFT定数行列生成
std::vector<std::vector<std::complex<double>>> mk_create_dft_matrix(int N, bool inverse);
// 多項式のインデックス反転（Automorphism）
void mk_poly_automorphism(TorusPolynomial* poly);

// 多項式の反転: P(X) → P(X^-1)
void mk_poly_inv_auto_inplace(TorusPolynomial* poly);

// Inv-Auto: Automorphism+KeySwitching
void mk_inv_auto(MKPackedRLWE* acc, const MKKeySwitchKey* ksk, const TFheGateBootstrappingParameterSet* params);
// Homomorphic DFT/IDFT（雛形）
void mk_homomorphic_dft(MKPackedRLWE* acc, const std::vector<std::vector<std::complex<double>>>& dft_matrix);
void mk_homomorphic_idft(MKPackedRLWE* acc, const std::vector<std::vector<std::complex<double>>>& idft_matrix);
// Batch-Anti-Rot: 反巡回シフト（本体）
void mk_batch_anti_rot(MKPackedRLWE* acc, const MKPackedRGSW* perm_key, const MKKeySwitchKey* ksk, const TFheGateBootstrappingParameterSet* params);
void mk_vec_mat_mult(
 MKPackedRLWE* acc,
 const std::vector<MKPackedRGSW*>& bk_list,
 const std::vector<int32_t>& coeffs,
 const TFheGateBootstrappingParameterSet* params
);

// Batch-Permute: permutation keyによる並べ替え（外部積ラッパー）
void mk_batch_permute(MKPackedRLWE* acc, const MKPackedRGSW* perm_key, const TFheGateBootstrappingParameterSet* params);

// Inv-Auto: インデックス反転（雛形）
void mk_inv_auto(MKPackedRLWE* acc, const MKKeySwitchKey* ksk);

// 再帰的DFT（雛形）
constexpr int BASE_SIZE = 8;
void mk_homomorphic_dft_recursive(MKPackedRLWE* acc, int current_N, const TFheGateBootstrappingParameterSet* params);

} // namespace bbii

#endif
