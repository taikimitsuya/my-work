#ifndef MK_OPS_H
#define MK_OPS_H
#include "mk_tfhe_structs.h"
#include <tfhe.h>
#include <tfhe_core.h>
#include <polynomials.h>
#include <tlwe.h>
#include <tgsw.h>

namespace bbii {
void mk_external_product(MKRLweSample* res, const MKPackedRGSW* bk, const MKRLweSample* acc, int32_t pid, const TFheGateBootstrappingParameterSet* p);
}

#endif // MK_OPS_H
