#ifndef MK_PACKED_OPS_H
#define MK_PACKED_OPS_H
#include "mk_tfhe_structs.h"
#include "mk_ops.h"

namespace bbii {
void mk_vec_mat_mult(
	MKPackedRLWE* acc,
	const std::vector<MKPackedRGSW*>& bk_list,
	const std::vector<int32_t>& coeffs,
	const TFheGateBootstrappingParameterSet* params
);
}

#endif
