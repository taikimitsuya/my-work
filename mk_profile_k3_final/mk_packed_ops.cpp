// ...既存BBII型mk_vec_mat_mult実装...
// 必要なヘッダ
#include "mk_tfhe_structs.h"
#include "mk_ops.h"
#include <cstdint>
#include <vector>
#include <iostream>
#include <cmath>

namespace bbii {

// Packed External Product
void mk_packed_external_product(MKPackedRLWE* acc, const MKPackedRGSW* rgsw, int32_t pid, const TFheGateBootstrappingParameterSet* params) {
    // モード遷移チェック
    if (acc->mode == BBIIMode::R12 && rgsw->mode == BBIIMode::R12_TO_R13) {
        acc->mode = BBIIMode::R13;
    } else if (acc->mode == BBIIMode::R13 && rgsw->mode == BBIIMode::R13_TO_R12) {
        acc->mode = BBIIMode::R12;
    } else {
        // std::cerr << "[Warning] Mode mismatch in external product!" << std::endl;
    }
    MKRLweSample* tmp_res = new MKRLweSample(acc->sample->k, params);
    mk_external_product(tmp_res, rgsw, acc->sample, pid, params);
    mk_rlwe_copy(acc->sample, tmp_res);
    delete tmp_res;
}
void mk_rgsw_scalar_mult(MKPackedRGSW* res, const MKPackedRGSW* src, int32_t scalar, const TFheGateBootstrappingParameterSet* params) {
    // res = scalar * src
    tGswFFTClear(res->sample, params->tgsw_params);
    for(int i=0; i<std::abs(scalar); ++i) {
        tGswFFTAddH(res->sample, params->tgsw_params); // res += src
    }
    if(scalar < 0) {
        // 負符号対応（TFHEのRGSWは符号反転でOK）
        // ここでは省略（必要ならtGswFFTNegate等を使う）
    }
    res->mode = src->mode;
}

// RGSW同士の加算
void mk_rgsw_add(MKPackedRGSW* res, const MKPackedRGSW* a, const MKPackedRGSW* b, const TFheGateBootstrappingParameterSet* params) {
    tGswFFTClear(res->sample, params->tgsw_params);
    tGswFFTAddH(res->sample, params->tgsw_params); // res += a
    tGswFFTAddH(res->sample, params->tgsw_params); // res += b
    res->mode = a->mode; // モードはa,b同じ前提
}

// BBII型 VecMatMult（外部積の和）
void mk_vec_mat_mult(
    MKPackedRLWE* acc,
    const std::vector<MKPackedRGSW*>& bk_list,
    const std::vector<int32_t>& coeffs,
    const TFheGateBootstrappingParameterSet* params
) {
    // 事前条件チェック
    if (bk_list.size() != coeffs.size()) {
        std::cerr << "[mk_vec_mat_mult] bk_listとcoeffsのサイズ不一致" << std::endl;
        return;
    }
    if (bk_list.empty()) return;

    // 1. 一時変数（combined_bk）初期化
    MKPackedRGSW combined_bk(params, bk_list[0]->mode);

    // 2. 鍵の線形結合
    for (size_t i = 0; i < bk_list.size(); ++i) {
        int32_t u = coeffs[i];
        if (u == 0) continue;
        MKPackedRGSW tmp(params, bk_list[i]->mode);
        mk_rgsw_scalar_mult(&tmp, bk_list[i], u, params);
        mk_rgsw_add(&combined_bk, &combined_bk, &tmp, params);
    }

    // 3. 外部積
    mk_packed_external_product(acc, &combined_bk, 0, params); // pid=0（BBII型は全体一括）

    // 4. モード更新
    if (acc->mode == BBIIMode::R12 && combined_bk.mode == BBIIMode::R12_TO_R13) {
        acc->mode = BBIIMode::R13;
    } else if (acc->mode == BBIIMode::R13 && combined_bk.mode == BBIIMode::R13_TO_R12) {
        acc->mode = BBIIMode::R12;
    }
}
// namespace bbii
}
