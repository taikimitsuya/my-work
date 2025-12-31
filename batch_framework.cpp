#include "batch_framework.h"

// 内部ヘルパー: TGswSample内のTLweSampleへアクセス
static inline TLweSample* get_tlwe_sample(TGswSample* sample, int index) {
    return &sample->all_sample[index];
}

namespace bbii {

PackedTRGSW create_zero_packed(const BBIIParams& params, BatchMode mode) {
    TGswSample* c = new_TGswSample_array(1, params.tfhe_params->tgsw_params);
    tGswClear(c, params.tfhe_params->tgsw_params);
    return PackedTRGSW(c, mode);
}

// 多項式加算: TorusPolynomial を使用
void torus_poly_add_to(TorusPolynomial* res, const TorusPolynomial* src, int32_t N) {
    for (int i = 0; i < N; ++i) {
        res->coefsT[i] += src->coefsT[i];
    }
}

void trgsw_add_to(TGswSample* res, const TGswSample* A, const BBIIParams& params) {
    const TGswParams* tgsw_p = params.tfhe_params->tgsw_params;
    const TLweParams* tlwe_p = tgsw_p->tlwe_params;
    int k = tlwe_p->k;
    int l = tgsw_p->l;
    int block_count = (k + 1) * l;
    int N = params.N;

    for (int i = 0; i < block_count; ++i) {
        TLweSample* res_tlwe = &res->all_sample[i];
        const TLweSample* A_tlwe = &A->all_sample[i];

        // a part (array of k TorusPolynomial)
        for (int j = 0; j < k; ++j) {
            torus_poly_add_to(&res_tlwe->a[j], &A_tlwe->a[j], N);
        }
        // b part (pointer to TorusPolynomial)
        torus_poly_add_to(res_tlwe->b, A_tlwe->b, N);

        res_tlwe->current_variance += A_tlwe->current_variance;
    }
}

// 多項式回転: TorusPolynomial を使用
void torus_poly_mul_by_xai(TorusPolynomial* res, const TorusPolynomial* src, int32_t delta, int32_t N) {
    delta %= (2 * N);
    if (delta < 0) delta += 2 * N;
    
    for (int i = 0; i < N; ++i) {
        int deg = i + delta;
        if (deg < N) res->coefsT[deg] = src->coefsT[i];
        else if (deg < 2 * N) res->coefsT[deg - N] = -src->coefsT[i];
        else res->coefsT[deg - 2*N] = src->coefsT[i];
    }
}

void trgsw_mul_by_xai(TGswSample* res, const TGswSample* input, int32_t delta, const BBIIParams& params) {
    const TGswParams* tgsw_p = params.tfhe_params->tgsw_params;
    const TLweParams* tlwe_p = tgsw_p->tlwe_params;
    int k = tlwe_p->k;
    int l = tgsw_p->l;
    int block_count = (k + 1) * l;
    int N = params.N;

    // Temporary buffer: コンストラクタで初期化 (Nを渡す)
    // クラス内部でメモリ確保されるため、手動new/deleteは不要
    TorusPolynomial temp_poly(N);

    for (int i = 0; i < block_count; ++i) {
        TLweSample* res_tlwe = &res->all_sample[i];
        const TLweSample* in_tlwe = &input->all_sample[i];

        // Process 'a' polynomials
        for (int j = 0; j < k; ++j) {
            torus_poly_mul_by_xai(&temp_poly, &in_tlwe->a[j], delta, N);
            for(int p=0; p<N; ++p) res_tlwe->a[j].coefsT[p] = temp_poly.coefsT[p];
        }
        
        // Process 'b' polynomial
        torus_poly_mul_by_xai(&temp_poly, in_tlwe->b, delta, N);
        for(int p=0; p<N; ++p) res_tlwe->b->coefsT[p] = temp_poly.coefsT[p];

        res_tlwe->current_variance = in_tlwe->current_variance;
    }
    // temp_poly のデストラクタが自動的にメモリを解放します
}

}
