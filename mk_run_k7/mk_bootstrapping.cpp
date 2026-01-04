#include "mk_methods.h"
#include <cmath>
#include <iostream>
#include <tfhe/tfhe.h>
#include <tfhe/tfhe_core.h>
#include <tfhe/polynomials.h>
namespace bbii {
void mk_cmux(MKRLweSample*, const TGswSampleFFT*, const MKRLweSample*, const MKRLweSample*, int32_t, const TFheGateBootstrappingParameterSet*);
void mk_rlwe_clear(MKRLweSample*); 
void mk_rlwe_copy(MKRLweSample*, const MKRLweSample*);
void mk_mul_xai(MKRLweSample* r, const MKRLweSample* s, int32_t b, int32_t N){ for(int i=0;i<=s->k;++i) torusPolynomialMulByXai(r->parts[i], b, s->parts[i]); }
void mk_blind_rotate(MKRLweSample* acc, const MKLweSample* bk_input, const MKBootstrappingKey* mk_bk, const TFheGateBootstrappingParameterSet* params) {
    int32_t k=acc->k, n=mk_bk->n_per_party, N=acc->N, _2N=2*N;
    int32_t bar_b = modSwitchFromTorus32(bk_input->sample->b, _2N);
    MKRLweSample* temp_acc = new MKRLweSample(k, params);
    mk_rlwe_copy(temp_acc, acc);
    mk_mul_xai(acc, temp_acc, -bar_b, N);
    
    // パーティ数ループ (k=7)
    for (int u = 0; u < k; ++u) { 
        std::cout << "  [Blind Rotate] Processing Party " << u+1 << "/" << k << "..." << std::endl;
        for (int i = 0; i < n; ++i) { 
            int32_t bar_ai = modSwitchFromTorus32(bk_input->sample->a[u*n+i], _2N);
            if (bar_ai == 0) continue;
            mk_mul_xai(temp_acc, acc, bar_ai, N);
            mk_cmux(acc, mk_bk->bk_fft[u][i], acc, temp_acc, u, params);
        }
    }
}
void mk_sample_extract(MKLweSample* output, const MKRLweSample* acc, const LweParams* lwe_params) {
    int32_t k = acc->k; int32_t N = acc->N;
    output->sample->b = acc->parts[k]->coefsT[0];
    for (int u = 0; u < k; ++u) {
        TorusPolynomial* poly = acc->parts[u];
        for (int j = 0; j < N; ++j) { output->sample->a[u*N+j] = (j==0)? poly->coefsT[0] : -poly->coefsT[N-j]; }
    }
}
void mk_bootstrapping(MKLweSample* res, const MKLweSample* in, const MKBootstrappingKey* bk, Torus32 mu, const TFheGateBootstrappingParameterSet* p) {
    int32_t k=in->k;
    MKRLweSample* acc=new MKRLweSample(k,p); mk_rlwe_clear(acc); acc->parts[k]->coefsT[0]=mu;
    mk_blind_rotate(acc, in, bk, p);
    mk_sample_extract(res, acc, p->in_out_params);
}
} 
