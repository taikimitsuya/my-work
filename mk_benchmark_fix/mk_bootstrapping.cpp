#include "mk_tfhe_structs.h"
#include "mk_params.h"
#include <cmath>
#include <chrono> 
#include <iostream>
#include <tfhe.h>
#include <tfhe_core.h>
#include <polynomials.h>

namespace bbii {
void mk_cmux(MKRLweSample*, const TGswSampleFFT*, const MKRLweSample*, const MKRLweSample*, int32_t, const TFheGateBootstrappingParameterSet*);
void mk_rlwe_clear(MKRLweSample*); 
void mk_rlwe_copy(MKRLweSample*, const MKRLweSample*);
void mk_mul_xai(MKRLweSample* r, const MKRLweSample* s, int32_t b, int32_t N){ 
    for(int i=0;i<=s->k;++i) torusPolynomialMulByXai(r->parts[i], b, s->parts[i]); 
}

void mk_blind_rotate(MKRLweSample* acc, const MKLweSample* bk_input, const MKBootstrappingKey* mk_bk, const TFheGateBootstrappingParameterSet* params) {
    int32_t k=acc->k, n=mk_bk->n_per_party, N=acc->N, _2N=2*N;
    int32_t bar_b = modSwitchFromTorus32(bk_input->sample->b, _2N);
    
    MKRLweSample* temp_acc = new MKRLweSample(k, params);
    mk_rlwe_copy(temp_acc, acc);
    mk_mul_xai(acc, temp_acc, -bar_b, N);

    auto start_loop = std::chrono::high_resolution_clock::now();
    long long op_count = 0;

    for (int u = 0; u < k; ++u) { 
        for (int i = 0; i < n; ++i) { 
            int32_t global_idx = u * n + i;
            if (bk_input->sample->a[global_idx] == 0) continue;
            int32_t bar_ai = modSwitchFromTorus32(bk_input->sample->a[global_idx], _2N);
            if (bar_ai == 0) continue;

            mk_mul_xai(temp_acc, acc, bar_ai, N);
            mk_cmux(acc, mk_bk->bk_fft[u][i], acc, temp_acc, u, params);
            op_count++;
        }
    }
    auto end_loop = std::chrono::high_resolution_clock::now();
    double total_ms = std::chrono::duration<double, std::milli>(end_loop - start_loop).count();
    double avg_per_op = (op_count > 0) ? total_ms / op_count : 0.0;

    std::cout << "    [Granularity 2] BlindRotate Loop Total: " << total_ms << " ms" << std::endl;
    std::cout << "    [Granularity 3] Avg per Atomic Op (CMUX): " << avg_per_op << " ms (" 
              << op_count << " ops executed)" << std::endl;

    delete temp_acc;
}

void mk_sample_extract(MKLweSample* output, const MKRLweSample* acc, const LweParams* lwe_params) {
    int32_t k = acc->k; int32_t N = acc->N;
    output->sample->b = acc->parts[k]->coefsT[0];
    for (int u = 0; u < k; ++u) {
        TorusPolynomial* poly = acc->parts[u];
        for (int j = 0; j < N; ++j) {
            output->sample->a[u*N+j] = (j==0)? poly->coefsT[0] : -poly->coefsT[N-j];
        }
    }
    output->sample->current_variance = lwe_params->alpha_min * lwe_params->alpha_min; 
}

void mk_bootstrapping(MKLweSample* res, const MKLweSample* in, const MKBootstrappingKey* bk, Torus32 mu, const TFheGateBootstrappingParameterSet* p) {
    int32_t k=in->k;
    MKRLweSample* acc=new MKRLweSample(k,p); mk_rlwe_clear(acc); acc->parts[k]->coefsT[0]=mu;
    mk_blind_rotate(acc, in, bk, p);
    mk_sample_extract(res, acc, p->in_out_params);
    delete acc;
}
} 
