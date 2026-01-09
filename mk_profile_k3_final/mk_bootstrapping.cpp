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

void mk_sample_extract(MKRLweSample* output, const MKRLweSample* acc, const LweParams* lwe_params) {
    auto start = std::chrono::high_resolution_clock::now();
    double extract0 = global_profiler.time_sample_extract;

    int32_t k = acc->k; int32_t N = acc->N;
    output->parts[k]->coefsT[0] = acc->parts[k]->coefsT[0];
    for (int u = 0; u < k; ++u) {
        TorusPolynomial* poly = acc->parts[u];
        for (int j = 0; j < N; ++j) { output->parts[u]->coefsT[j] = (j==0)? poly->coefsT[0] : -poly->coefsT[N-j]; }
    }

    auto end = std::chrono::high_resolution_clock::now();
    double extract1 = global_profiler.time_sample_extract;
    global_profiler.time_sample_extract = (extract1 - extract0) + std::chrono::duration<double, std::milli>(end - start).count();
}

void mk_bootstrapping(MKRLweSample* res, const MKRLweSample* in, const MKBootstrappingKey* bk, Torus32 mu, const TFheGateBootstrappingParameterSet* p) {
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
} 
