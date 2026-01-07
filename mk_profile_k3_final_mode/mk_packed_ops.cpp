#include "mk_tfhe_structs.h"
#include <tfhe.h>
#include <tfhe_core.h>
#include <vector>

namespace bbii {
// RGSW同士の加算
void mk_rgsw_add(MKPackedRGSW* res, const MKPackedRGSW* a, const MKPackedRGSW* b, const TFheGateBootstrappingParameterSet* params) {
    // TGswSampleFFT*の加算（要素ごと）
    // 例: tGswFFTSampleAdd(res->sample, a->sample, b->sample, params->tgsw_params);
    // 実際は各行ごとに加算する必要あり
}

// RGSW同士の乗算（Batch-Mult）
void mk_rgsw_batch_mult(MKPackedRGSW* res, const MKPackedRGSW* a, const MKPackedRGSW* b, const TFheGateBootstrappingParameterSet* params) {
    // 各行ごとにmk_external_productを適用し、RGSW×RGSW→RGSWの行列積を構築
    // 例: for (int row = 0; row < ...; ++row) { ... }
    // 擬似コード: tGswFFTExternMulToTLwe(res->sample->all_sample[row], b->sample, ...);
}
}
