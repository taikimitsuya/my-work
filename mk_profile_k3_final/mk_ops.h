#ifndef MK_OPS_H
#define MK_OPS_H

namespace bbii {
struct MKRLweSample;
struct MKPackedRGSW;
}

#include "mk_tfhe_structs.h"



namespace bbii {
// RLWEクリア
void mk_rlwe_clear(bbii::MKRLweSample* r);
// 外部積
void mk_external_product(bbii::MKRLweSample* res, const bbii::MKPackedRGSW* rgsw, const bbii::MKRLweSample* acc, int32_t pid, const TFheGateBootstrappingParameterSet* params);
// RLWEコピー
void mk_rlwe_copy(bbii::MKRLweSample* dest, const bbii::MKRLweSample* src);
}

#endif
