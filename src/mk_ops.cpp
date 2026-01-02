#include "mk_tfhe_structs.h"
#include <vector>
#include <iostream>
#include <tfhe.h>
#include <tfhe_core.h>
#include <polynomials.h>
#include <tlwe.h>
#include <tgsw.h>

namespace bbii {

void mk_rlwe_clear(MKRLweSample* result) {
    for (int i = 0; i <= result->k; ++i) torusPolynomialClear(result->parts[i]);
}

void mk_rlwe_copy(MKRLweSample* dst, const MKRLweSample* src) {
    for (int i = 0; i <= dst->k; ++i) torusPolynomialCopy(dst->parts[i], src->parts[i]);
}

void mk_rlwe_addTo(MKRLweSample* result, const MKRLweSample* sample) {
    for (int i = 0; i <= result->k; ++i) torusPolynomialAddTo(result->parts[i], sample->parts[i]);
}

void mk_rlwe_subTo(MKRLweSample* result, const MKRLweSample* sample) {
    for (int i = 0; i <= result->k; ++i) torusPolynomialSubTo(result->parts[i], sample->parts[i]);
}

void mk_external_product(MKRLweSample* result, 
                         const TGswSampleFFT* bk_fft, 
                         const MKRLweSample* in_acc, 
                         int32_t party_id, 
                         const TFheGateBootstrappingParameterSet* params) {
    int32_t k = in_acc->k;
    TLweParams* tlwe_params = params->tgsw_params->tlwe_params;
    TLweSample* temp_tlwe = new_TLweSample(tlwe_params);

    mk_rlwe_clear(result);

    for (int i = 0; i <= k; ++i) {
        const TorusPolynomial* current_poly = in_acc->parts[i];
        torusPolynomialClear(&temp_tlwe->a[0]); 
        torusPolynomialCopy(temp_tlwe->b, current_poly); 
        temp_tlwe->current_variance = 0; 

        tGswFFTExternMulToTLwe(temp_tlwe, bk_fft, params->tgsw_params);

        torusPolynomialAddTo(result->parts[i], temp_tlwe->b);
        torusPolynomialAddTo(result->parts[party_id], &temp_tlwe->a[0]);
    }
    delete_TLweSample(temp_tlwe);
}

void mk_cmux(MKRLweSample* result, 
             const TGswSampleFFT* bk_fft, 
             const MKRLweSample* in0, 
             const MKRLweSample* in1, 
             int32_t party_id,
             const TFheGateBootstrappingParameterSet* params) {
    MKRLweSample temp_diff(in0->k, params);
    mk_rlwe_copy(&temp_diff, in1);
    mk_rlwe_subTo(&temp_diff, in0); 

    MKRLweSample temp_prod(in0->k, params);
    mk_external_product(&temp_prod, bk_fft, &temp_diff, party_id, params);

    mk_rlwe_copy(result, in0);
    mk_rlwe_addTo(result, &temp_prod);
}

} 
