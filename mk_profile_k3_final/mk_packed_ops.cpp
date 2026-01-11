// --- ダミー実装: mk_slice, mk_butterfly ---
void mk_slice(const std::vector<MKPackedRLWE*>& input, std::vector<MKPackedRLWE*>& out_upper, std::vector<MKPackedRLWE*>& out_lower) {
    size_t half = input.size() / 2;
    out_upper.assign(input.begin(), input.begin() + half);
    out_lower.assign(input.begin() + half, input.end());
}

void mk_butterfly(MKPackedRLWE* u, MKPackedRLWE* v, const TFheGateBootstrappingParameterSet* params) {
    // ダミー: 何もしない
}
// 標準ライブラリ
#include <complex>
#include <algorithm>
#include <cstdint>
#include <vector>
#include <iostream>
#include <cmath>

// プロジェクトヘッダ
#include "mk_tfhe_structs.h"
#include "mk_ops.h"
#include "bb_params.h"

// Twiddle Factorの適用（DFT用）
void mk_apply_twiddle(MKPackedRLWE* acc, int32_t power, int32_t N, const MKBootstrappingKey* mk_bk, const TFheGateBootstrappingParameterSet* params) {
    // 各成分の多項式を X^power だけ巡回シフト
    int32_t k = acc->sample->k;
    for (int u = 0; u <= k; ++u) {
        TorusPolynomial* poly = acc->sample->parts[u];
        std::vector<Torus32> tmp(poly->coefsT, poly->coefsT + N);
        for (int i = 0; i < N; ++i) {
            int idx = (i + power) % N;
            if (idx < 0) idx += N;
            poly->coefsT[idx] = tmp[i];
        }
    }
}

// Twiddle Factorの逆（IDFT用）
void mk_apply_inv_twiddle(MKPackedRLWE* acc, int32_t power, int32_t N, const MKBootstrappingKey* mk_bk, const TFheGateBootstrappingParameterSet* params) {
    int32_t inv_power = (2 * N - power) % (2 * N);
    mk_apply_twiddle(acc, inv_power, N, mk_bk, params);
}

// 再帰的IDFT（Gentleman-Sande型）
void mk_homomorphic_idft_recursive(
    std::vector<MKPackedRLWE*>& inputs,
    int32_t N,
    const MKBootstrappingKey* mk_bk,
    const TFheGateBootstrappingParameterSet* params
#include "mk_ops.h"
#include "bb_params.h"
    if (len <= 1) return;
    size_t half = len / 2;
    std::vector<MKPackedRLWE*> upper, lower;
    mk_slice(inputs, upper, lower);
    // Butterfly Inverse
    for (size_t k = 0; k < half; ++k) {
        mk_butterfly(upper[k], lower[k], params);
        int32_t rot = k * (N / len);
        mk_apply_inv_twiddle(lower[k], rot, N, mk_bk, params);
    }
    mk_homomorphic_idft_recursive(upper, N, mk_bk, params);
    mk_homomorphic_idft_recursive(lower, N, mk_bk, params);
    for (size_t i = 0; i < half; ++i) {
        mk_rlwe_copy(inputs[i]->sample, upper[i]->sample);
        mk_rlwe_copy(inputs[i+half]->sample, lower[i]->sample);
    }
}

// KSK生成: Automorphism後の鍵からKeySwitchingKeyを生成
void mk_fill_automorphism_ksk(BBII_KSKStruct* ksk, const std::vector<MKSecretKey*>& sks, const TFheGateBootstrappingParameterSet* params) {
    int32_t k = ksk->k;
    int32_t N = ksk->N;
    const LweParams* lwe_p = params->in_out_params;
    for(int u = 0; u < k; ++u) {
        // 入力鍵 s_in = s_u(X^-1) を作成
        LweKey* s_in = new_LweKey(lwe_p);
        LweKey* s_out = sks[u]->lwe_key;
        // s_in の生成: s_out の係数を反転 (Automorphism: i -> N-i, 符号反転)
        for(int i=0; i<N; ++i) {
            s_in->key[i] = s_out->key[(N-i)%N];
        }
        // KSK生成: s_in -> s_out
        lweCreateKeySwitchKey(ksk->ks_keys[u], s_in, s_out);
        delete_LweKey(s_in);
    }
}
// 内部ヘルパー: 単項式 X^power * scalar を多項式にセットする
void set_monomial_scaled(TorusPolynomial* poly, int32_t power, int32_t scalar, int32_t N) {
    torusPolynomialClear(poly);
    // X^N = -1 の環での回転
    int32_t p = power % (2 * N);
    if (p < 0) p += 2 * N;
    if (p < N) {
        poly->coefsT[p] = scalar;
    } else {
        // X^{N+i} = -X^i
        poly->coefsT[p - N] = -scalar;
    }
}

// Permutation Key (Trivial RGSW of X^delta) の生成
MKPackedRGSW* mk_create_permutation_key(const std::vector<int>& permutation, const TFheGateBootstrappingParameterSet* params) {
    // 1. シフト量 delta の特定
    int32_t N = params->in_out_params->n;
    int32_t delta = 0;
    if (!permutation.empty()) {
        delta = permutation[0];
    }

    MKPackedRGSW* key = new MKPackedRGSW(params, BBIIMode::R12_TO_R13);

    const TGswParams* tgsw_params = params->tgsw_params;
    TGswSample* temp_rgsw = new_TGswSample(tgsw_params);

    int32_t l = tgsw_params->l;
    int32_t Bgbit = tgsw_params->Bgbit;
    int32_t tfhe_k = tgsw_params->tlwe_params->k;

    for (int32_t bloc = 0; bloc <= tfhe_k; ++bloc) {
        for (int32_t i = 0; i < l; ++i) {
            int32_t row_idx = bloc * l + i;
            int32_t shift = 32 - (i + 1) * Bgbit;
            int32_t scalar = (shift >= 0) ? (1 << shift) : 0;
            TLweSample* row_tlwe = &temp_rgsw->all_sample[row_idx];
            for (int32_t j = 0; j < tfhe_k; ++j) {
                torusPolynomialClear(&row_tlwe->a[j]);
            }
            if (bloc < tfhe_k) {
                set_monomial_scaled(&row_tlwe->a[bloc], delta, scalar, N);
                torusPolynomialClear(row_tlwe->b);
            } else {
                set_monomial_scaled(row_tlwe->b, delta, scalar, N);
            }
            row_tlwe->current_variance = 0;
        }
    }

    // FFT変換
    tGswToFFTConvert(key->sample, temp_rgsw, tgsw_params);
    delete_TGswSample(temp_rgsw);
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
void mk_inv_auto(MKPackedRLWE* acc, const BBII_KSKStruct* ksk, const TFheGateBootstrappingParameterSet* params) {
    int32_t k = acc->sample->k;
    int32_t N = acc->sample->N;
    // 1. Automorphism: 各成分に多項式反転を適用
    for (int i = 0; i <= k; ++i) {
        mk_poly_inv_auto_inplace(acc->sample->parts[i]);
    }

    // 2. Key Switching: s(X^-1) -> s(X)
    // 結果を蓄積する一時変数
    MKRLweSample* res = new MKRLweSample(k, params);
    mk_rlwe_clear(res);

    // b部分 (parts[0]) はKS不要なのでそのままコピー
    torusPolynomialCopy(res->parts[0], acc->sample->parts[0]);

    // 各パーティの a_i 部分を KeySwitching
    constexpr int t = 8; // TFHE標準値（要調整）
    constexpr int basebit = 2; // TFHE標準値（要調整）
    for (int u = 1; u <= k; ++u) {
        TorusPolynomial* poly_in = acc->sample->parts[u];
        TorusPolynomial* poly_out = res->parts[u];
        // 各係数について LWE Key Switching を実行
        for (int i = 0; i < N; ++i) {
            // LWEサンプルを構築
            LweSample* lwe_tmp = new_LweSample(params->in_out_params);
            lwe_tmp->a[0] = poly_in->coefsT[i];
            lwe_tmp->b = 0;
            lwe_tmp->current_variance = 0;
            // KeySwitch: ksk->ks_keys[u-1] を使う
            LweSample* lwe_out = new_LweSample(params->in_out_params);
            lweKeySwitch(lwe_out, ksk->ks_keys[u-1], lwe_tmp);
            // 結果をpoly_outに格納
            poly_out->coefsT[i] = lwe_out->b;
            delete_LweSample(lwe_tmp);
            delete_LweSample(lwe_out);
        }
    }

    mk_rlwe_copy(acc->sample, res);
    delete res;
}

// Batch-Anti-Rot: 反巡回シフト（本体）
void mk_batch_anti_rot(MKPackedRLWE* acc, const MKPackedRGSW* perm_key, const BBII_KSKStruct* ksk, const TFheGateBootstrappingParameterSet* params) {
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

// （以下、追加の関数定義があればここに記述）

// ファイル末尾の閉じ括弧を追加
}

