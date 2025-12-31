// batch_framework.h
#ifndef BATCH_FRAMEWORK_H
#define BATCH_FRAMEWORK_H

#include "bb_params.h"

namespace bbii {

enum class BatchMode { R12, R13, R12_to_R13, R13_to_R12 };

struct PackedTRGSW {
    TRGSW* cipher;
    BatchMode mode;

    PackedTRGSW(TRGSW* c, BatchMode m) : cipher(c), mode(m) {}
    // デストラクタ等は省略（ポインタ管理は上位で行う想定）
};

// ゼロのTRGSWを生成
PackedTRGSW create_zero_packed(const BBIIParams& params, BatchMode mode);

// TRGSW同士の加算: res += A
void trgsw_add_to(TRGSW* res, const TRGSW* A, const BBIIParams& params);

// 多項式回転: res = input * X^delta (反巡回)
void trgsw_mul_by_xai(TRGSW* res, const TRGSW* input, int32_t delta, const BBIIParams& params);

} // namespace bbii
#endif