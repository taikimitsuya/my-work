#ifndef MK_PACKED_OPS_H
#define MK_PACKED_OPS_H
#include "mk_tfhe_structs.h"
#include <tfhe.h>
#include <tfhe_core.h>

namespace bbii {
void mk_rgsw_add(MKPackedRGSW* res, const MKPackedRGSW* a, const MKPackedRGSW* b, const TFheGateBootstrappingParameterSet* params);
void mk_rgsw_batch_mult(MKPackedRGSW* res, const MKPackedRGSW* a, const MKPackedRGSW* b, const TFheGateBootstrappingParameterSet* params);
MKPackedRGSW* mk_generate_gadget_rgsw(const TFheGateBootstrappingParameterSet* params, int k, BBIIMode mode);
void mk_vec_mat_mult(MKRLweSample* acc, MKPackedRGSW* packed_bk, int32_t digit, const TFheGateBootstrappingParameterSet* params);
}
#endif
