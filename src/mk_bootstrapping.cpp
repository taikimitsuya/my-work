#include "mk_tfhe_structs.h"
#include "mk_params.h"
#include <vector>
#include <cmath>
#include <tfhe.h>
#include <tfhe_core.h>
#include <polynomials.h>
#include <lweparams.h>
#include <lwekey.h>
#include <tgsw.h>

namespace bbii {

void mk_cmux(MKRLweSample* result, const TGswSampleFFT* bk_fft, const MKRLweSample* in0, const MKRLweSample* in1, int32_t party_id, const TFheGateBootstrappingParameterSet* params);
void mk_rlwe_clear(MKRLweSample* result);
void mk_rlwe_copy(MKRLweSample* dst, const MKRLweSample* src);

void mk_rlwe_mul_by_xai(MKRLweSample* result, const MKRLweSample* src, int32_t bar, int32_t N) {
    for (int i = 0; i <= src->k; ++i) torusPolynomialMulByXai(result->parts[i], bar, src->parts[i]);
}

void mk_blind_rotate(MKRLweSample* acc,
                     const MKLweSample* bk_input,
                     const MKBootstrappingKey* mk_bk,
                     const TFheGateBootstrappingParameterSet* params) {
    int32_t k = acc->k;
    int32_t n = mk_bk->n_per_party; 
    int32_t N = acc->N;
    int32_t _2N = 2 * N;

    int32_t bar_b = modSwitchFromTorus32(bk_input->sample->b, _2N);
    MKRLweSample* temp_acc = new MKRLweSample(k, params);
    mk_rlwe_copy(temp_acc, acc);
    mk_rlwe_mul_by_xai(acc, temp_acc, -bar_b, N);

    for (int u = 0; u < k; ++u) { 
        for (int i = 0; i < n; ++i) { 
            int32_t global_idx = u * n + i;
            Torus32 a_ui = bk_input->sample->a[global_idx];
            if (a_ui == 0) continue;
            int32_t bar_ai = modSwitchFromTorus32(a_ui, _2N);
            if (bar_ai == 0) continue;

            mk_rlwe_mul_by_xai(temp_acc, acc, bar_ai, N);
            mk_cmux(acc, mk_bk->bk_fft[u][i], acc, temp_acc, u, params);
        }
    }
    delete temp_acc;
}

void mk_sample_extract(MKLweSample* output, const MKRLweSample* acc, const LweParams* lwe_params) {
    int32_t k = acc->k;
    int32_t N = acc->N;
    output->sample->b = acc->parts[k]->coefsT[0];
    for (int u = 0; u < k; ++u) {
        TorusPolynomial* poly = acc->parts[u];
        int offset = u * N;
        output->sample->a[offset + 0] = poly->coefsT[0];
        for (int j = 1; j < N; ++j) {
            output->sample->a[offset + j] = -poly->coefsT[N - j];
        }
    }
    output->sample->current_variance = lwe_params->alpha_min * lwe_params->alpha_min; 
}

void mk_bootstrapping(MKLweSample* result,
                      const MKLweSample* input,
                      const MKBootstrappingKey* mk_bk,
                      Torus32 mu, 
                      const TFheGateBootstrappingParameterSet* params) {
    int32_t k = input->k;
    MKRLweSample* acc = new MKRLweSample(k, params);
    mk_rlwe_clear(acc);
    acc->parts[k]->coefsT[0] = mu; 

    mk_blind_rotate(acc, input, mk_bk, params);
    mk_sample_extract(result, acc, params->in_out_params);
    delete acc;
}

} 
