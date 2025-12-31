#ifndef BATCH_OPS_H
#define BATCH_OPS_H
#include "batch_framework.h"
#include <vector>
namespace bbii {
PackedTRGSW vec_mat_mult(const std::vector<int32_t>& a, const std::vector<PackedTRGSW>& B, const BBIIParams& params);
PackedTRGSW batch_anti_rot(const PackedTRGSW& C, int32_t delta, const BBIIParams& params);
std::vector<PackedTRGSW> enc_vec_mat_mult(const std::vector<std::vector<int32_t>>& M_exp, const std::vector<PackedTRGSW>& C, const BBIIParams& params);
}
#endif
