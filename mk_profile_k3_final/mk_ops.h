#ifndef MK_OPS_H
#define MK_OPS_H
#include "mk_tfhe_structs.h"

namespace bbii {
// 外部積
void mk_external_product(MKRLweSample* res, const MKPackedRGSW* rgsw, const MKRLweSample* acc, int32_t pid, const TFheGateBootstrappingParameterSet* params);
// RLWEコピー
void mk_rlwe_copy(MKRLweSample* dest, const MKRLweSample* src);
}

#endif
