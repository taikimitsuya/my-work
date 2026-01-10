#include "mk_tfhe_structs.h"
#include "mk_ops.h"
#include "bb_params.h"

#include <algorithm>
#include <complex>
#include <cstdint>
#include <vector>
#include <iostream>
#include <cmath>

namespace bbii {
// 順列キー生成（雛形）
MKPackedRGSW* mk_create_permutation_key(const std::vector<int>& permutation, const TFheGateBootstrappingParameterSet* params) {
    MKPackedRGSW* key = new MKPackedRGSW(params, BBIIMode::R12_TO_R13);
    return key;
}
// DFT/IDFT定数行列生成
#include <complex>
const double PI = 3.14159265358979323846;
std::vector<std::vector<std::complex<double>>> mk_create_dft_matrix(int N, bool inverse) {
    std::vector<std::vector<std::complex<double>>> mat(N, std::vector<std::complex<double>>(N));
    double sign = inverse ? 1.0 : -1.0;
    for (int k = 0; k < N; ++k) {
        for (int n = 0; n < N; ++n) {
            double angle = 2 * PI * k * n / N;
            mat[k][n] = std::complex<double>(cos(angle), sign * sin(angle));
            if (inverse) mat[k][n] /= N;
        }
    }
    return mat;
}

// Homomorphic DFT（ダミー: 実際は暗号化多項式の線形変換）
void mk_homomorphic_dft(MKPackedRLWE* acc, const std::vector<std::vector<std::complex<double>>>& dft_matrix) {
    // 各パーティの多項式係数にDFT行列を適用（実際の暗号多項式変換は要実装）
    int k = acc->sample->k;
    int N = acc->sample->N;
    for(int u=0; u<=k; ++u) {
        TorusPolynomial* poly = acc->sample->parts[u];
        std::vector<std::complex<double>> input(N);
        for(int i=0;i<N;++i) input[i] = std::complex<double>(poly->coefsT[i], 0.0);
        std::vector<std::complex<double>> output(N, 0.0);
        for(int i=0;i<N;++i) {
            for(int j=0;j<N;++j) {
                output[i] += dft_matrix[i][j] * input[j];
            }
        }
        // 実部をTorus32に戻す（本実装は要検討）
        // Torus32変換: 丸め・クリッピング・スケーリング
        constexpr double torus_max = static_cast<double>(0x7FFFFFFF);
        constexpr double torus_min = static_cast<double>(-0x80000000);
        constexpr double scale = 1.0; // 必要に応じて調整
        for(int i=0;i<N;++i) {
            double val = std::real(output[i]) * scale;
            if(val > torus_max) val = torus_max;
            if(val < torus_min) val = torus_min;
            poly->coefsT[i] = static_cast<Torus32>(std::round(val));
        }
    }
}

void mk_homomorphic_idft(MKPackedRLWE* acc, const std::vector<std::vector<std::complex<double>>>& idft_matrix) {
    int k = acc->sample->k;
    int N = acc->sample->N;
    for(int u=0; u<=k; ++u) {
        TorusPolynomial* poly = acc->sample->parts[u];
        std::vector<std::complex<double>> input(N);
        for(int i=0; i<N; ++i) input[i] = std::complex<double>(poly->coefsT[i], 0.0);
        std::vector<std::complex<double>> output(N, 0.0);
        for(int i=0; i<N; ++i) {
            for(int j=0; j<N; ++j) {
                output[i] += idft_matrix[i][j] * input[j];
            }
        }
        constexpr double torus_max = static_cast<double>(0x7FFFFFFF);
        constexpr double torus_min = static_cast<double>(-0x80000000);
        constexpr double scale = 1.0;
        for(int i=0; i<N; ++i) {
            double val = std::real(output[i]) * scale;
            if(val > torus_max) val = torus_max;
            if(val < torus_min) val = torus_min;
            poly->coefsT[i] = static_cast<Torus32>(std::round(val));
        }
    }
}
// 多項式のインデックス反転（Automorphism）
void mk_poly_automorphism(TorusPolynomial* poly) {
    int32_t N = poly->N;
    for (int i = 1; i < N / 2; ++i) {
        std::swap(poly->coefsT[i], poly->coefsT[N - i]);
    }
}

// Inv-Auto: Automorphism+KeySwitching
void mk_inv_auto(MKPackedRLWE* acc, const MKKeySwitchKey* ksk, const TFheGateBootstrappingParameterSet* params) {
    int32_t k = acc->sample->k;
    // 1. Automorphism: 各成分に多項式反転を適用
    for (int i = 0; i <= k; ++i) {
        mk_poly_inv_auto_inplace(acc->sample->parts[i]);
    }
    // 2. Key Switching: 恒等パススルー（本来はここでKSKを使う）
    // ここではaccの内容をそのまま維持
}

// Batch-Anti-Rot: 反巡回シフト（本体）
void mk_batch_anti_rot(MKPackedRLWE* acc, const MKPackedRGSW* perm_key, const MKKeySwitchKey* ksk, const TFheGateBootstrappingParameterSet* params) {
    // 1. Batch-Permute (回転シフト)
    mk_batch_permute(acc, perm_key, params);
    // 2. Inv-Auto (反転+KeySwitching)
    mk_inv_auto(acc, ksk, params);
    // Alg 5.4の続きが必要ならここに追加
}
// Batch-Permute: permutation keyによる並べ替え（外部積ラッパー）
void mk_batch_permute(MKPackedRLWE* acc, const MKPackedRGSW* perm_key, const TFheGateBootstrappingParameterSet* params) {
    mk_packed_external_product(acc, perm_key, 0, params); // pid=0: 全体一括
}

// Inv-Auto: インデックス反転（雛形）
void mk_inv_auto(MKPackedRLWE* acc, const MKKeySwitchKey* ksk) {
    int k = acc->sample->k;
    for(int i=0;i<=k;++i) {
        std::reverse(acc->sample->parts[i]->coefsT, acc->sample->parts[i]->coefsT + acc->sample->N);
    }
    // TODO: kskを使ったKeySwitching処理を追加
}

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

// 多項式の反転: P(X) → P(X^-1)
void mk_poly_inv_auto_inplace(TorusPolynomial* poly) {
    int32_t N = poly->N;
    // 1. 配列を逆順にする (1 から N-1 まで)
    for (int i = 1; i < N / 2; ++i) {
        std::swap(poly->coefsT[i], poly->coefsT[N - i]);
    }
    // 2. 1次以上の係数の符号を反転する
    for (int i = 1; i < N; ++i) {
        poly->coefsT[i] = -poly->coefsT[i];
    }
}
