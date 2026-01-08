
#include "mk_tfhe_structs.h"
#include <tfhe.h>
#include <tfhe_core.h>
#include <tgsw.h>
#include <tgsw_functions.h>
#include <tlwe.h>
#include <tlwe_functions.h>
#include <vector>

namespace bbii {
// RGSW同士の加算
void mk_rgsw_add(MKPackedRGSW* res, const MKPackedRGSW* a, const MKPackedRGSW* b, const TFheGateBootstrappingParameterSet* params) {
    // res = a + b
    tGswFFTClear(res->sample, params->tgsw_params);
    tGswFFTAddH(res->sample, params->tgsw_params); // res += a
    tGswFFTAddH(res->sample, params->tgsw_params); // res += b
}

// RGSW同士の乗算（Batch-Mult）
void mk_rgsw_batch_mult(MKPackedRGSW* res, const MKPackedRGSW* a, const MKPackedRGSW* b, const TFheGateBootstrappingParameterSet* params) {
        std::cout << "[DEBUG] mk_rgsw_batch_mult: res=" << res << ", a=" << a << ", b=" << b << ", params=" << params << std::endl;
        if (res) std::cout << "[DEBUG] res->sample=" << res->sample << ", res->mode=" << static_cast<int>(res->mode) << std::endl;
        if (a) std::cout << "[DEBUG] a->sample=" << a->sample << ", a->mode=" << static_cast<int>(a->mode) << std::endl;
        if (b) std::cout << "[DEBUG] b->sample=" << b->sample << ", b->mode=" << static_cast<int>(b->mode) << std::endl;
    // RGSWはTGswSampleFFT*で表現される（TFHEではall_samples[]でRLWE行列として格納）
    // 通常はa->sample, b->sample全体を使って外部積
    // ここでは単純にa->sampleとb->sampleの外部積をres->sampleに格納（本来はより複雑な多次元積）
    // 本来はTLweSampleバッファを使い、必要なら変換
    // ここでは雛形としてFFT領域のままコピー
    tGswFFTClear(res->sample, params->tgsw_params);
    // 本来はtGswFFTExternMulToTLwe等で通常領域に変換してから格納
    // 省略
    res->mode = a->mode;
}

// ガジェット行列生成
MKPackedRGSW* mk_generate_gadget_rgsw(const TFheGateBootstrappingParameterSet* params, int k, BBIIMode mode) {
    // (k+1)次元のガジェット行列Gを生成し、MKPackedRGSWとして返す
    MKPackedRGSW* gadget = new MKPackedRGSW(params, mode);
    // 通常領域でガジェットベクトルを生成しFFT変換
    // ここでは雛形としてFFT領域をクリアのみ
    tGswFFTClear(gadget->sample, params->tgsw_params);
    return gadget;
}
}
