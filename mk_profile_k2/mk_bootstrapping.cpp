#include "mk_methods.h"
#include "mk_profiler.h"
#include <cmath>
#include <iostream>
#include <chrono>
#include <tfhe/tfhe.h>
#include <tfhe/tfhe_core.h>
#include <tfhe/polynomials.h>

namespace bbii {
void mk_cmux(MKRLweSample*, const TGswSampleFFT*, const MKRLweSample*, const MKRLweSample*, int32_t, const TFheGateBootstrappingParameterSet*);
void mk_rlwe_clear(MKRLweSample*); 
void mk_rlwe_copy(MKRLweSample*, const MKRLweSample*);
void mk_mul_xai(MKRLweSample* r, const MKRLweSample* s, int32_t b, int32_t N){ for(int i=0;i<=s->k;++i) torusPolynomialMulByXai(r->parts[i], b, s->parts[i]); }

void mk_blind_rotate(MKRLweSample* acc, const MKLweSample* bk_input, const MKBootstrappingKey* mk_bk, const TFheGateBootstrappingParameterSet* params) {
    // Blind Rotate 全体の時間計測開始
    auto br_start = std::chrono::high_resolution_clock::now();

    int32_t k=acc->k, n=mk_bk->n_per_party, N=acc->N, _2N=2*N;
    int32_t bar_b = modSwitchFromTorus32(bk_input->sample->b, _2N);
    MKRLweSample* temp_acc = new MKRLweSample(k, params);
    mk_rlwe_copy(temp_acc, acc);
    mk_mul_xai(acc, temp_acc, -bar_b, N);

    for (int u = 0; u < k; ++u) { 
        if(u==0 && k>5) std::cout << "  [Progress] Party processing started..." << std::endl;
        for (int i = 0; i < n; ++i) { 
            int32_t bar_ai = modSwitchFromTorus32(bk_input->sample->a[u*n+i], _2N);
            if (bar_ai == 0) continue;
            mk_mul_xai(temp_acc, acc, bar_ai, N);
            mk_cmux(acc, mk_bk->bk_fft[u][i], acc, temp_acc, u, params);
        }
    }
    
    auto br_end = std::chrono::high_resolution_clock::now();
    double total_br_time = std::chrono::duration<double, std::milli>(br_end - br_start).count();
    
    // Control部分の時間 = 全体 - (重いExternalProductの時間)
    // ExternalProductの時間は mk_ops.cpp で別途加算されているため
    // ここでは「純粋な制御時間」を計算するために引き算はせず、
    // 表示の際に調整するか、Blind Rotate全体として記録する。
    // 今回は「MatMult以外」を強調したいため、ここでは何も足さず、表示で工夫します。
}

void mk_sample_extract(MKLweSample* output, const MKRLweSample* acc, const LweParams* lwe_params) {
    auto start = std::chrono::high_resolution_clock::now();

    int32_t k = acc->k; int32_t N = acc->N;
    output->sample->b = acc->parts[k]->coefsT[0];
    for (int u = 0; u < k; ++u) {
        TorusPolynomial* poly = acc->parts[u];
        for (int j = 0; j < N; ++j) { output->sample->a[u*N+j] = (j==0)? poly->coefsT[0] : -poly->coefsT[N-j]; }
    }

    auto end = std::chrono::high_resolution_clock::now();
    global_profiler.time_sample_extract += std::chrono::duration<double, std::milli>(end - start).count();
}

void mk_bootstrapping(MKLweSample* res, const MKLweSample* in, const MKBootstrappingKey* bk, Torus32 mu, const TFheGateBootstrappingParameterSet* p) {
    // 1. Input Packing (Initialization)
    auto pack_start = std::chrono::high_resolution_clock::now();
    int32_t k=in->k;
    MKRLweSample* acc=new MKRLweSample(k,p); mk_rlwe_clear(acc); acc->parts[k]->coefsT[0]=mu;
    auto pack_end = std::chrono::high_resolution_clock::now();
    global_profiler.time_input_packing += std::chrono::duration<double, std::milli>(pack_end - pack_start).count();

    // 2. Blind Rotate (Control + MatMult)
    auto br_start = std::chrono::high_resolution_clock::now();
    mk_blind_rotate(acc, in, bk, p);
    auto br_end = std::chrono::high_resolution_clock::now();
    
    // Blind Rotate全体の時間から、ExternalProductの時間を引いたものを「制御コスト」とする
    double br_total = std::chrono::duration<double, std::milli>(br_end - br_start).count();
    global_profiler.time_blind_rotate_control = br_total - global_profiler.time_external_product;

    // 3. Sample Extract
    mk_sample_extract(res, acc, p->in_out_params);
}
} 
