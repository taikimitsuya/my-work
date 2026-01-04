#include "mk_methods.h"
#include "mk_profiler.h" // プロファイラー読み込み
#include <iostream>
#include <chrono>
#include <tfhe/tfhe.h>
#include <tfhe/tfhe_core.h>
#include <tfhe/polynomials.h>
#include <tfhe/tlwe.h>
#include <tfhe/tgsw.h>

namespace bbii {
void mk_rlwe_clear(MKRLweSample* r){ for(int i=0;i<=r->k;++i) torusPolynomialClear(r->parts[i]); }
void mk_rlwe_copy(MKRLweSample* d, const MKRLweSample* s){ for(int i=0;i<=d->k;++i) torusPolynomialCopy(d->parts[i], s->parts[i]); }
void mk_rlwe_addTo(MKRLweSample* r, const MKRLweSample* s){ for(int i=0;i<=r->k;++i) torusPolynomialAddTo(r->parts[i], s->parts[i]); }
void mk_rlwe_subTo(MKRLweSample* r, const MKRLweSample* s){ for(int i=0;i<=r->k;++i) torusPolynomialSubTo(r->parts[i], s->parts[i]); }

void mk_external_product(MKRLweSample* res, const TGswSampleFFT* bk, const MKRLweSample* acc, int32_t pid, const TFheGateBootstrappingParameterSet* p) {
    // 【重要】ここの計算が最も重い (Homomorphic DFT / FFT部分)
    auto start = std::chrono::high_resolution_clock::now();

    TLweParams* tlp = const_cast<TLweParams*>(p->tgsw_params->tlwe_params);
    TLweSample* tmp = new_TLweSample(tlp);
    mk_rlwe_clear(res);
    for(int i=0;i<=acc->k;++i){
        torusPolynomialClear(&tmp->a[0]); torusPolynomialCopy(tmp->b, acc->parts[i]); tmp->current_variance=0;
        tGswFFTExternMulToTLwe(tmp, bk, const_cast<TGswParams*>(p->tgsw_params));
        torusPolynomialAddTo(res->parts[i], tmp->b); torusPolynomialAddTo(res->parts[pid], &tmp->a[0]);
    }
    delete_TLweSample(tmp);

    auto end = std::chrono::high_resolution_clock::now();
    global_profiler.time_external_product += std::chrono::duration<double, std::milli>(end - start).count();
}

void mk_cmux(MKRLweSample* res, const TGswSampleFFT* bk, const MKRLweSample* in0, const MKRLweSample* in1, int32_t pid, const TFheGateBootstrappingParameterSet* p) {
    MKRLweSample diff(in0->k, p); mk_rlwe_copy(&diff, in1); mk_rlwe_subTo(&diff, in0);
    MKRLweSample prod(in0->k, p); mk_external_product(&prod, bk, &diff, pid, p);
    mk_rlwe_copy(res, in0); mk_rlwe_addTo(res, &prod);
}
} 
