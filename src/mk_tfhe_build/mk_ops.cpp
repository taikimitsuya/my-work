#include "mk_tfhe_structs.h"
#include <vector>
#include <iostream>

// libtfhe の内部関数 (FFT, 多項式演算) を利用
#include <tfhe.h>
#include <tfhe_core.h>
#include <polynomials.h>
#include <tlwe.h>
#include <tgsw.h>

namespace bbii {

// ==========================================
// ヘルパー関数: 多項式ベクトルの基本演算
// ==========================================

// MK-RLWE をゼロクリアする
void mk_rlwe_clear(MKRLweSample* result) {
    for (int i = 0; i <= result->k; ++i) {
        torusPolynomialClear(result->parts[i]);
    }
}

// MK-RLWE のコピー (dst = src)
void mk_rlwe_copy(MKRLweSample* dst, const MKRLweSample* src) {
    if (dst->k != src->k) {
        throw std::runtime_error("MKRLweSample dimension mismatch in copy");
    }
    for (int i = 0; i <= dst->k; ++i) {
        torusPolynomialCopy(dst->parts[i], src->parts[i]);
    }
}

// MK-RLWE の加算 (result += sample)
void mk_rlwe_addTo(MKRLweSample* result, const MKRLweSample* sample) {
    for (int i = 0; i <= result->k; ++i) {
        torusPolynomialAddTo(result->parts[i], sample->parts[i]);
    }
}

// MK-RLWE の減算 (result -= sample)
void mk_rlwe_subTo(MKRLweSample* result, const MKRLweSample* sample) {
    for (int i = 0; i <= result->k; ++i) {
        torusPolynomialSubTo(result->parts[i], sample->parts[i]);
    }
}

// ==========================================
// コア演算: External Product
// ==========================================

/**
 * @brief マルチキー External Product
 * * 指定されたパーティのTRGSW鍵 (bk_fft) を、MK-RLWEアキュムレータ (in_acc) に掛け合わせます。
 * MK-RLWEのすべての多項式成分に対して、TRGSWとの積を計算し、適切に再配置します。
 * * @param result     [out] 結果を格納するMKRLweSample (初期化済みであること)
 * @param bk_fft     [in]  特定のパーティのTRGSWサンプル (FFT形式)
 * @param in_acc     [in]  入力となるMKRLweSample
 * @param party_id   [in]  bk_fft の持ち主のパーティID (0 ~ k-1)
 * @param params     [in]  パラメータ
 */
void mk_external_product(MKRLweSample* result, 
                         const TGswSampleFFT* bk_fft, 
                         const MKRLweSample* in_acc, 
                         int32_t party_id, 
                         const TFheGateBootstrappingParameterSet* params) {
    
    int32_t k = in_acc->k;
    int32_t N = in_acc->N;
    
    // 計算用の一時的なTLWEサンプル (シングルキー用、k=1)
    // 構造: (a, b) のペア。N個の係数を持つ多項式2つ。
    TLweParams* tlwe_params = params->tgsw_params->tlwe_params;
    TLweSample* temp_tlwe = new_TLweSample(tlwe_params);

    // 結果蓄積用バッファは、呼び出し元でクリアされていることを期待しない場合、
    // ここでゼロクリアするか、呼び出し元に任せるか。
    // 通常は result = 0 から始めて加算していくため、ここではクリアしてから計算します。
    mk_rlwe_clear(result);

    // MK-RLWEの各成分 (a_0, ..., a_{k-1}, b) に対してループ
    // 合計 k+1 回の Single-Key External Product を実行
    for (int i = 0; i <= k; ++i) {
        const TorusPolynomial* current_poly = in_acc->parts[i];

        // 1. temp_tlwe を (0, current_poly) にセットする
        //    mask (a) は 0, body (b) は current_poly
        torusPolynomialClear(&temp_tlwe->a[0]); // a=0
        torusPolynomialCopy(temp_tlwe->b, current_poly); // b=current_poly
        temp_tlwe->current_variance = 0; // ノイズ管理は簡易化

        // 2. Single-Key External Product を実行
        //    temp_tlwe = bk_fft * temp_tlwe
        //    計算結果は temp_tlwe に上書きされます (出力: u, v)
        //    u = temp_tlwe->a[0], v = temp_tlwe->b
        tGswFFTExternMulToTLwe(temp_tlwe, bk_fft, params->tgsw_params);

        // 3. 結果を MK-RLWE に配置
        //    理論: current_poly * bk = (u, v)
        //    これは、u * s_{party_id} + v = current_poly * message という関係
        //    したがって、
        //      v (body部分) は、元の位置 i に加算
        //      u (mask部分) は、party_id に対応する位置に加算
        
        // v を result->parts[i] に加算
        torusPolynomialAddTo(result->parts[i], temp_tlwe->b);

        // u を result->parts[party_id] に加算
        // ※ party_id は 0 ~ k-1 なので、parts配列のインデックスと一致します
        torusPolynomialAddTo(result->parts[party_id], &temp_tlwe->a[0]);
    }

    // メモリ解放
    delete_TLweSample(temp_tlwe);
}


// ==========================================
// CMUX (Controlled Multiplexer)
// ==========================================

/**
 * @brief MK-RLWE CMUX
 * * result = in0 + bk_fft * (in1 - in0)
 * もし bk_fft が 1 (True) を暗号化していれば in1、0 (False) なら in0 を選択
 * * @param result     [out] 結果
 * @param bk_fft     [in]  選択ビットのTRGSW (FFT)
 * @param in0        [in]  Falseの場合の入力
 * @param in1        [in]  Trueの場合の入力
 * @param party_id   [in]  bk_fft の持ち主ID
 * @param params     [in]  パラメータ
 */
void mk_cmux(MKRLweSample* result, 
             const TGswSampleFFT* bk_fft, 
             const MKRLweSample* in0, 
             const MKRLweSample* in1, 
             int32_t party_id,
             const TFheGateBootstrappingParameterSet* params) {
    
    // 1. 差分計算: temp = in1 - in0
    //    MKRLweSample は動的に確保する必要があるため、スタックではなくヒープ推奨だが
    //    ここではラッパーとして管理
    MKRLweSample temp_diff(in0->k, params);
    mk_rlwe_copy(&temp_diff, in1);
    mk_rlwe_subTo(&temp_diff, in0); // temp = in1 - in0

    // 2. 積計算: temp_prod = bk_fft * temp_diff
    MKRLweSample temp_prod(in0->k, params);
    mk_external_product(&temp_prod, bk_fft, &temp_diff, party_id, params);

    // 3. 最終加算: result = in0 + temp_prod
    mk_rlwe_copy(result, in0);
    mk_rlwe_addTo(result, &temp_prod);
    
    // (temp_diff, temp_prod はデストラクタで解放される)
}

} // namespace bbii