// batch_framework.cpp
#include "batch_framework.h"

namespace bbii {

PackedTRGSW create_zero_packed(const BBIIParams& params, BatchMode mode) {
    TRGSW* c = new_TRGSW_array(1, params.tfhe_params->tgsw_params);
    trgswClear(c, params.tfhe_params->tgsw_params);
    return PackedTRGSW(c, mode);
}

// TRGSWの加算実装 (行列の各要素を加算)
void trgsw_add_to(TRGSW* res, const TRGSW* A, const BBIIParams& params) {
    const TGswParams* tgsw_p = params.tfhe_params->tgsw_params;
    const TLweParams* tlwe_p = tgsw_p->tlwe_params;
    
    // TRGSWは (k+1)*l 個の TLWE からなる
    int block_count = (tlwe_p->k + 1) * tgsw_p->l;

    for (int i = 0; i < block_count; ++i) {
        // TLwe同士の加算: res->blocs[i] += A->blocs[i]
        // 各係数を足す
        for (int j = 0; j <= tlwe_p->k; ++j) { // a_0...a_k, b
             TorusPoly* res_poly = &res->blocs[i].a[j]; // TFHEの構造に依存(a[k]がbかも)
             // 注: libtfheでは blocs[i] は TLweSample
             // TLweSample { TorusPoly* a; TorusPoly* b; ... }
             // a はサイズ k, b は単体
             
             // 簡易アクセスのため、TFHEの関数を使うのが安全だが、
             // ここでは内部構造を操作する:
             // tLweAddTo(&res->blocs[i], &A->blocs[i], tlwe_p); // もしあれば
             
             // ない場合は手動:
             const TorusPoly* a_poly_src = (j < tlwe_p->k) ? &A->blocs[i].a[j] : &A->blocs[i].b;
             TorusPoly* res_poly_dst = (j < tlwe_p->k) ? &res->blocs[i].a[j] : &res->blocs[i].b;
             
             for (int p = 0; p < params.N; ++p) {
                 res_poly_dst->coefsT[p] += a_poly_src->coefsT[p];
             }
        }
        // Current variance update (approx)
        res->blocs[i].current_variance += A->blocs[i].current_variance;
    }
}

// 多項式回転の実装 (Negacyclic Shift)
// X^N = -1 mod (X^N + 1)
void torus_poly_mul_by_xai(TorusPoly* res, const TorusPoly* src, int32_t delta, int32_t N) {
    // delta を [0, 2N) に正規化
    delta %= (2 * N);
    if (delta < 0) delta += 2 * N;

    for (int i = 0; i < N; ++i) {
        // term: src[i] * X^i * X^delta = src[i] * X^{i+delta}
        int deg = i + delta;
        
        // Handle X^N = -1
        if (deg < N) {
            res->coefsT[deg] = src->coefsT[i];
        } else if (deg < 2 * N) {
            res->coefsT[deg - N] = -src->coefsT[i];
        } else { // deg >= 2N
             res->coefsT[deg - 2*N] = src->coefsT[i];
        }
    }
}

void trgsw_mul_by_xai(TRGSW* res, const TRGSW* input, int32_t delta, const BBIIParams& params) {
    const TGswParams* tgsw_p = params.tfhe_params->tgsw_params;
    const TLweParams* tlwe_p = tgsw_p->tlwe_params;
    int block_count = (tlwe_p->k + 1) * tgsw_p->l;
    int N = params.N;

    for (int i = 0; i < block_count; ++i) {
        // a vector
        for (int k = 0; k < tlwe_p->k; ++k) {
            // 一時バッファを使って In-place 対応
            TorusPoly temp_poly;
            // new ではなく stack alloc or pre-alloc推奨だが簡易的に
            temp_poly.coefsT = new Torus32[N]; 
            
            torus_poly_mul_by_xai(&temp_poly, &input->blocs[i].a[k], delta, N);
            
            // Copy back
            for(int p=0; p<N; ++p) res->blocs[i].a[k].coefsT[p] = temp_poly.coefsT[p];
            delete[] temp_poly.coefsT;
        }
        // b polynomial
        TorusPoly temp_poly_b;
        temp_poly_b.coefsT = new Torus32[N];
        torus_poly_mul_by_xai(&temp_poly_b, &input->blocs[i].b, delta, N);
        for(int p=0; p<N; ++p) res->blocs[i].b.coefsT[p] = temp_poly_b.coefsT[p];
        delete[] temp_poly_b.coefsT;
        
        // Variance copy
        res->blocs[i].current_variance = input->blocs[i].current_variance;
    }
}

} // namespace bbii