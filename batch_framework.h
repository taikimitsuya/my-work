#ifndef BATCH_FRAMEWORK_H
#define BATCH_FRAMEWORK_H
#include "bb_params.h"

// TFHEの型を確実に認識させる
struct TGswSample;
struct TorusPolynomial; // TorusPoly -> TorusPolynomial

namespace bbii {
enum class BatchMode { R12, R13, R12_to_R13, R13_to_R12, None };
struct PackedTRGSW {
    TGswSample* cipher; BatchMode mode;
    PackedTRGSW() : cipher(nullptr), mode(BatchMode::None) {}
    PackedTRGSW(TGswSample* c, BatchMode m) : cipher(c), mode(m) {}
};
PackedTRGSW create_zero_packed(const BBIIParams& params, BatchMode mode);
void trgsw_add_to(TGswSample* res, const TGswSample* A, const BBIIParams& params);
void trgsw_mul_by_xai(TGswSample* res, const TGswSample* input, int32_t delta, const BBIIParams& params);
}
#endif
