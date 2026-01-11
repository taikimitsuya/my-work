#ifndef MK_PACKED_OPS_H
#define MK_PACKED_OPS_H

#include <complex>
#undef complex
#include <vector>
#include "bb_params.h"
#include "mk_ops.h"

// 再帰型DFT（ダミー: 実際は分割統治）
void mk_homomorphic_dft_recursive(std::vector<MKPackedRLWE*>& inputs, int32_t N, const MKBootstrappingKey* mk_bk, const TFheGateBootstrappingParameterSet* params);


// DFT/IDFT バッチ処理
void mk_homomorphic_dft_batch(std::vector<MKPackedRLWE*>& accs, const std::vector<std::vector<std::complex<double>>>& dft_matrix);
void mk_homomorphic_idft_batch(std::vector<MKPackedRLWE*>& accs, const std::vector<std::vector<std::complex<double>>>& idft_matrix);

// DFT/IDFT定数行列生成
std::vector<std::vector<std::complex<double>>> mk_create_dft_matrix(int N, bool inverse);

// Homomorphic DFT/IDFT
void mk_homomorphic_dft(MKPackedRLWE* acc, const std::vector<std::vector<std::complex<double>>>& dft_matrix);
void mk_homomorphic_idft(MKPackedRLWE* acc, const std::vector<std::vector<std::complex<double>>>& idft_matrix);

// Batch-Anti-Rot
void mk_batch_anti_rot(MKPackedRLWE* acc, const MKPackedRGSW* perm_key, const BBII_KSKStruct* ksk, const TFheGateBootstrappingParameterSet* params);

// VecMatMult
void mk_vec_mat_mult(MKPackedRLWE* acc, const std::vector<MKPackedRGSW*>& bk_list, const std::vector<int32_t>& coeffs, const TFheGateBootstrappingParameterSet* params);

// Perm/Key cache
MKPackedRGSW* get_perm_key_cached(MKBootstrappingKey* mk_bk, const std::vector<int>& permutation, const TFheGateBootstrappingParameterSet* params);
BBII_KSKStruct* get_ksk_cached(MKBootstrappingKey* mk_bk, int delta, int k, int N, const TFheGateBootstrappingParameterSet* params);

// DFT/IDFT・Twiddle
void mk_apply_inv_twiddle(MKPackedRLWE* acc, int32_t power, int32_t N, const MKBootstrappingKey* mk_bk, const TFheGateBootstrappingParameterSet* params);
void mk_homomorphic_idft_recursive(std::vector<MKPackedRLWE*>& inputs, int32_t N, const MKBootstrappingKey* mk_bk, const TFheGateBootstrappingParameterSet* params);

// Automorphism Key Switching
void mk_fill_automorphism_ksk(BBII_KSKStruct* ksk, const std::vector<MKSecretKey*>& sks, const TFheGateBootstrappingParameterSet* params);
BBII_KSKStruct* mk_create_automorphism_ksk(int32_t k, int32_t N, const TFheGateBootstrappingParameterSet* params);

// RLWE演算
void mk_rlwe_addTo(MKPackedRLWE*, MKPackedRLWE*);
void mk_rlwe_subTo(MKPackedRLWE*, MKPackedRLWE*);

// 外部積
void mk_packed_external_product(MKPackedRLWE* acc, const MKPackedRGSW* rgsw, int32_t pid, const TFheGateBootstrappingParameterSet* params);

// パーミュテーションキー生成
MKPackedRGSW* mk_create_permutation_key(const std::vector<int>& permutation, const TFheGateBootstrappingParameterSet* params);

// Batch-Permute: permutation keyによる並べ替え（外部積ラッパー）
void mk_batch_permute(MKPackedRLWE* acc, const MKPackedRGSW* perm_key, const TFheGateBootstrappingParameterSet* params);

#endif // MK_PACKED_OPS_H
