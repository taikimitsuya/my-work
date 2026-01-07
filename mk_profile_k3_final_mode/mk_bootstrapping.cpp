#include "mk_methods.h"
#include "mk_profiler.h"
#include <cmath>
#include <iostream>
#include <chrono>
#include <tfhe.h>
#include <tfhe_core.h>
#include <polynomials.h>
#include "mk_tfhe_structs.h"

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

    for (int u = 0; u < k; ++u) { 
        if(u==0) std::cout << "  [Progress] Blind Rotate loop start..." << std::endl;
        for (int i = 0; i < n; ++i) { 
            int32_t bar_ai = modSwitchFromTorus32(bk_input->sample->a[u*n+i], _2N);
            if (bar_ai == 0) continue;
            
            mk_mul_xai(temp_acc, acc, bar_ai, N);
            MKPackedRGSW* packed_ct = mk_bk->bk_packed[u][i];
            // 必要なら packed_ct->mode で分岐可能
            mk_cmux(acc, packed_ct, acc, temp_acc, u, params);
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
}