#include "mk_methods.h"
#include "mk_profiler.h"
#include <cmath>
#include <iostream>
#include <chrono>
#include <tfhe.h>
#include <tfhe_core.h>

#include "mk_packed_ops.h"
#include <polynomials.h>

namespace bbii {
// DFTベースBlind Rotate（雛形・流れのみ）
void mk_blind_rotate_dft(bbii::MKRLweSample* acc, const bbii::MKRLweSample* bk_input, const bbii::MKBootstrappingKey* mk_bk, const TFheGateBootstrappingParameterSet* params) {
    int32_t k = acc->k, N = acc->N;
    // 1. accをパック型にラップ
    MKPackedRLWE* acc_packed = new MKPackedRLWE(k, params, BBIIMode::R12);
    mk_rlwe_copy(acc_packed->sample, acc);

    // 2. DFT行列生成
    auto dft_mat = mk_create_dft_matrix(N, false);
    auto idft_mat = mk_create_dft_matrix(N, true);

    // 3. Homomorphic DFT
    mk_homomorphic_dft(acc_packed, dft_mat);


    // 4. Batch-Anti-Rot（テスト用ダミー: perm_key, kskをその場で生成）
    // 実際はdeltaやbk_input等から適切なperm_key, kskを選択
    std::vector<int> dummy_perm(N); for(int i=0;i<N;++i) dummy_perm[i]=i; // 恒等順列
    MKPackedRGSW* perm_key = new MKPackedRGSW(params, BBIIMode::R12_TO_R13); // ダミー
    MKKeySwitchKey* ksk = new MKKeySwitchKey(k, N, params); // ダミー
    mk_batch_anti_rot(acc_packed, perm_key, ksk, params);
    delete perm_key;
    delete ksk;

    // 5. Homomorphic IDFT
    mk_homomorphic_idft(acc_packed, idft_mat);

    // 6. 結果を書き戻す
    mk_rlwe_copy(acc, acc_packed->sample);
    delete acc_packed;
}
void mk_cmux(MKRLweSample*, const TGswSampleFFT*, const MKRLweSample*, const MKRLweSample*, int32_t, const TFheGateBootstrappingParameterSet*);
void mk_rlwe_clear(MKRLweSample*); 
void mk_rlwe_copy(MKRLweSample*, const MKRLweSample*);
void mk_mul_xai(MKRLweSample* r, const MKRLweSample* s, int32_t b, int32_t N){
    for(int i=0;i<=s->k;++i) {
        if (r->parts[i]->N != s->parts[i]->N) {
            TorusPolynomial* old = r->parts[i];
            r->parts[i] = new_TorusPolynomial(s->parts[i]->N);
            torusPolynomialClear(r->parts[i]);
            delete_TorusPolynomial(old);
        }
        torusPolynomialMulByXai(r->parts[i], b, s->parts[i]);
    }
    r->N = s->N;
}

void mk_blind_rotate(MKRLweSample* acc, const MKRLweSample* bk_input, const MKBootstrappingKey* mk_bk, const TFheGateBootstrappingParameterSet* params) {
    auto br_start = std::chrono::high_resolution_clock::now();
    double extprod_start = global_profiler.time_external_product;

    int32_t k=acc->k, n=mk_bk->n_per_party, N=acc->N, _2N=2*N;
    int32_t bar_b = modSwitchFromTorus32(bk_input->parts[k]->coefsT[0], _2N);

    // accをBBII用にラップ
    MKPackedRLWE* acc_packed = new MKPackedRLWE(k, params, BBIIMode::R12);
    mk_rlwe_copy(acc_packed->sample, acc);

    int32_t a0 = ((_2N - bar_b) % _2N);
    mk_mul_xai(acc_packed->sample, acc_packed->sample, a0, N);

    // BBII型: iごとに全uのbk_packedをまとめて一括外部積
    std::cout << "  [Progress] Blind Rotate (BBII) loop start..." << std::endl;
    for (int i = 0; i < n; ++i) {
        std::vector<MKPackedRGSW*> bk_vec;
        std::vector<int32_t> coeff_vec;
        for (int u = 0; u < k; ++u) {
            int32_t bar_ai = modSwitchFromTorus32(bk_input->parts[u]->coefsT[i], _2N);
            if (bar_ai == 0) continue;
            bk_vec.push_back(mk_bk->bk_packed[u][i]);
            coeff_vec.push_back(bar_ai);
        }
        if (!bk_vec.empty()) {
            mk_vec_mat_mult(acc_packed, bk_vec, coeff_vec, params);
        }
    }

    // 結果を書き戻す
    mk_rlwe_copy(acc, acc_packed->sample);
    delete acc_packed;

    auto br_end = std::chrono::high_resolution_clock::now();
    double br_total = std::chrono::duration<double, std::milli>(br_end - br_start).count();
    double extprod_end = global_profiler.time_external_product;
    double extprod_diff = extprod_end - extprod_start;
    global_profiler.time_blind_rotate_control = br_total - extprod_diff;
}

void mk_sample_extract(MKRLweSample* output, const MKRLweSample* acc, const MKBootstrappingKey* mk_bk, const TFheGateBootstrappingParameterSet* params, const LweParams* lwe_params) {
    auto start = std::chrono::high_resolution_clock::now();
    double extract0 = global_profiler.time_sample_extract;

    int32_t k = acc->k, N = acc->N;
    MKPackedRLWE* acc_packed = new MKPackedRLWE(k, params, BBIIMode::R12);
    mk_rlwe_copy(acc_packed->sample, acc);

    // DFT/IDFT行列
    auto dft_mat = mk_create_dft_matrix(N, false);
    auto idft_mat = mk_create_dft_matrix(N, true);

    // DFT
    mk_homomorphic_dft(acc_packed, dft_mat);

    // delta, perm_key, kskの取得（例: delta=1の回転）
    int delta = 1;
    std::vector<int> permutation(N);
    for(int i=0;i<N;++i) permutation[i] = (i+delta)%N;
    // perm_key/kskをキャッシュ経由で取得
    MKPackedRGSW* perm_key = bbii::get_perm_key_cached(const_cast<bbii::MKBootstrappingKey*>(mk_bk), permutation, params);
    MKKeySwitchKey* ksk = bbii::get_ksk_cached(const_cast<bbii::MKBootstrappingKey*>(mk_bk), delta, k, N, params);

    // Batch-Anti-Rot
    mk_batch_anti_rot(acc_packed, perm_key, ksk, params);

    // IDFT
    mk_homomorphic_idft(acc_packed, idft_mat);

    mk_rlwe_copy(const_cast<MKRLweSample*>(acc), acc_packed->sample);

    delete acc_packed;
}
} 
