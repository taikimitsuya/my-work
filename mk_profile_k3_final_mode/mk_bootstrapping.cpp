#include "mk_methods.h"
#include "mk_profiler.h"
#include <cmath>
#include <iostream>
#include <chrono>
#include <tfhe.h>
#include <tfhe_core.h>
#include <polynomials.h>
#include "mk_tfhe_structs.h"
#include "mk_packed_ops.h"

namespace bbii {
void mk_cmux(MKRLweSample*, const MKPackedRGSW*, const MKRLweSample*, const MKRLweSample*, int32_t, const TFheGateBootstrappingParameterSet*);
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

void mk_blind_rotate(MKRLweSample* acc, const MKLweSample* bk_input, const MKBootstrappingKey* mk_bk, const TFheGateBootstrappingParameterSet* params) {
    auto br_start = std::chrono::high_resolution_clock::now();
    // Blind Rotate内でのexternal_product合計のみを差分で計測
    double extprod_start = global_profiler.time_external_product;

    int32_t k=acc->k, n=mk_bk->n_per_party, N=acc->N, _2N=2*N;
    int32_t bar_b = modSwitchFromTorus32(bk_input->sample->b, _2N);
    
    MKRLweSample* temp_acc = new MKRLweSample(k, params);
    mk_rlwe_copy(temp_acc, acc);
    
    int32_t a0 = ((_2N - bar_b) % _2N);
    mk_mul_xai(acc, temp_acc, a0, N);

    // --- BBII論文 Algorithm 4.1 の3-4行目に相当する桁ごとの分解・行列積ループ ---
    // 例: L = log_B(q) 桁分解数, B = 基数, n = パーティごとのビット長
    const int L = 3; // 例: 桁数（パラメータに応じて調整）
    const int B = 1 << 10; // 例: 基数（パラメータに応じて調整）
    for (int ell = 0; ell < L; ++ell) { // 各桁ごと
        for (int u = 0; u < k; ++u) { // 各参加者
            for (int j = 0; j < n; ++j) { // 各ビット
                // 桁ごとにbar_aを分解
                int32_t aij = bk_input->sample->a[u*n+j];
                int32_t digit = (aij >> (ell * 10)) & (B - 1); // 10bitごとに分解
                if (digit == 0) continue;
                // BBII: 各桁ごとにBKを適用
                mk_vec_mat_mult(acc, mk_bk->bk_packed[u][j], digit, params);
            }
        }
    }
    
    delete temp_acc; 

    auto br_end = std::chrono::high_resolution_clock::now();
    double br_total = std::chrono::duration<double, std::milli>(br_end - br_start).count();
    double extprod_end = global_profiler.time_external_product;
    double extprod_diff = extprod_end - extprod_start;
    // Blind Rotate (Control) = Blind Rotate全体 - Blind Rotate内でのexternal_product合計
    global_profiler.time_blind_rotate_control = br_total - extprod_diff;
}

void mk_sample_extract(MKLweSample* output, const MKRLweSample* acc, const LweParams* lwe_params) {
    auto start = std::chrono::high_resolution_clock::now();
    double extract0 = global_profiler.time_sample_extract;

    int32_t k = acc->k; int32_t N = acc->N;
    output->sample->b = acc->parts[k]->coefsT[0];
    for (int u = 0; u < k; ++u) {
        TorusPolynomial* poly = acc->parts[u];
        for (int j = 0; j < N; ++j) { output->sample->a[u*N+j] = (j==0)? poly->coefsT[0] : -poly->coefsT[N-j]; }
    }

    auto end = std::chrono::high_resolution_clock::now();
    double extract1 = global_profiler.time_sample_extract;
    global_profiler.time_sample_extract = (extract1 - extract0) + std::chrono::duration<double, std::milli>(end - start).count();
}

void mk_bootstrapping(MKLweSample* res, const MKLweSample* in, const MKBootstrappingKey* bk, Torus32 mu, const TFheGateBootstrappingParameterSet* p) {
    auto pack_start = std::chrono::high_resolution_clock::now();
    double pack0 = global_profiler.time_input_packing;
    int32_t k=in->k;
    MKRLweSample* acc=new MKRLweSample(k,p); 
    mk_rlwe_clear(acc); 
    acc->parts[k]->coefsT[0]=mu;
    auto pack_end = std::chrono::high_resolution_clock::now();
    double pack1 = global_profiler.time_input_packing;
    global_profiler.time_input_packing = (pack1 - pack0) + std::chrono::duration<double, std::milli>(pack_end - pack_start).count();

    mk_blind_rotate(acc, in, bk, p);

    mk_sample_extract(res, acc, p->in_out_params);
    delete acc;
}

// --- VecMatMultの雛形（再掲）---
void mk_vec_mat_mult(MKRLweSample* acc, MKPackedRGSW* packed_bk, int32_t digit, const TFheGateBootstrappingParameterSet* params) {
    // digitが0ならガジェット行列Gを使う
    MKPackedRGSW* op = nullptr;
    if (digit == 0) {
        // マルチキー用ガジェット行列を生成
        op = mk_generate_gadget_rgsw(params, acc->k, packed_bk->mode);
    } else {
        // digit回加算（またはスカラー乗算）
        // ここでは単純にdigit回加算する例
        op = new MKPackedRGSW(*packed_bk); // コピー
        // 実際はdigit回加算やスカラー乗算に最適化可能
    }
    // accとopのBatch-Mult（RGSW×RGSW→RGSW）
    bbii::mk_rgsw_batch_mult(reinterpret_cast<MKPackedRGSW*>(acc), reinterpret_cast<const MKPackedRGSW*>(acc), op, params);
    if (digit == 0 || op != packed_bk) delete op;
}

// --- BBII: RGSW×RGSW→RGSWの行列積（Batch-Mult）具体例 ---
void mk_batch_mult(MKPackedRGSW* acc, MKPackedRGSW* op, const TFheGateBootstrappingParameterSet* params) {
    // acc, opともにRGSW（行列）として扱う
    // ここでは各行（RLWE）ごとにExternal Productを適用し、新しいRGSWを構築する例
    // ※本来はTGswSampleFFT*の多次元配列を使うが、雛形として簡略化
    // 例: acc->sample->all_sample[i] = ExternalProduct(acc->sample->all_sample[i], op->sample, ...)
    // 必要に応じてmode遷移もここで管理
    // --- 擬似コード ---
    // for (int row = 0; row < 行数; ++row) {
    //     tGswFFTExternMulToTLwe(acc->sample->all_sample[row], op->sample, ...);
    // }
}

// --- BBII: マルチキー用ガジェット行列生成関数の具体例 ---
MKPackedRGSW* mk_generate_gadget_rgsw(const TFheGateBootstrappingParameterSet* params, int k, BBIIMode mode) {
    // (k+1)次元のガジェット行列Gを生成し、MKPackedRGSWとして返す
    MKPackedRGSW* gadget = new MKPackedRGSW(params, mode);
    // 例: 各行にガジェットベクトル（[0, ..., 0, g]）をセット
    // for (int i = 0; i <= k; ++i) {
    //     set_gadget_vector(gadget->sample->all_sample[i], ...);
    // }
    return gadget;
}
}