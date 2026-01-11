#ifndef BB_PARAMS_H
#define BB_PARAMS_H
#include <cmath>
#include <vector>
#include <tfhe/tfhe.h>
#include <tfhe/tfhe_core.h>

// [追加] BBIIモード管理用enum
enum class BBIIMode {
    R12,
    R13,
    R12_TO_R13,
    R13_TO_R12
};

struct BBIIParams {
    int32_t n; int32_t N; int32_t d; int32_t rho; int32_t r;
    TFheGateBootstrappingParameterSet* tfhe_params;
    double lwe_alpha;
    double tlwe_alpha;
    int32_t l;
    int32_t Bgbit;

    // 切り替え用フラグ
    static constexpr bool use_custom_params = true; // trueで手動生成, falseで公式API

    // BBII論文推奨値（例: CGGI19, Table 2）
    // lwe_alpha = 2^{-15} ≈ 3.05e-5
    // tlwe_alpha = 2^{-25} ≈ 2.98e-8
    // l = 2, Bgbit = 10 など（論文値に応じて調整）
    BBIIParams(int32_t d_val, int32_t rho_val, int32_t N_val,
              double lwe_alpha_val = 3.05e-5,
              double tlwe_alpha_val = 2.98e-8,
              int32_t l_val = 2,
              int32_t Bgbit_val = 10) :
        d(d_val), rho(rho_val), N(N_val),
        lwe_alpha(lwe_alpha_val), tlwe_alpha(tlwe_alpha_val), l(l_val), Bgbit(Bgbit_val) {
        // BBII: n = 2 * d^rho
        n = 2 * std::pow(d, rho);
        r = N / 2;

        if constexpr (use_custom_params) {
            // 手動パラメータ生成（BBII論文値で生成）
            LweParams* lp = new_LweParams(n, lwe_alpha, 0.5);
            TLweParams* tlp = new_TLweParams(N, 1, tlwe_alpha, 0.5);
            TGswParams* tgp = new_TGswParams(l, Bgbit, tlp);
            tfhe_params = new TFheGateBootstrappingParameterSet(10, 1, lp, tgp);
        } else {
            // 公式API（推奨パラメータのみ）
            tfhe_params = new_default_gate_bootstrapping_parameters(n);
        }
    }
};
#endif // BB_PARAMS_H
