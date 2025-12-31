#ifndef BATCH_OPS_H
#define BATCH_OPS_H

#include "batch_framework.h"
#include <vector>

namespace bbii {

// Algorithm 4.1: 平文ビットベクトル a と 暗号化行列 B の積
PackedTRGSW vec_mat_mult(const std::vector<int32_t>& a, 
                         const std::vector<PackedTRGSW>& B, 
                         const BBIIParams& params);

// Algorithm 5.4: 反巡回回転 (X^delta 倍)
PackedTRGSW batch_anti_rot(const PackedTRGSW& C, int32_t delta, const BBIIParams& params);

// Algorithm 5.5: 特殊行列(DFT)との積
std::vector<PackedTRGSW> enc_vec_mat_mult(const std::vector<std::vector<int32_t>>& M_exponents, 
                                          const std::vector<PackedTRGSW>& C, 
                                          const BBIIParams& params);
}
#endif