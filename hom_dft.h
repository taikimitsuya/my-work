#ifndef HOM_DFT_H
#define HOM_DFT_H
#include "batch_ops.h"
namespace bbii {
std::vector<PackedTRGSW> hom_dft_inverse(const std::vector<PackedTRGSW>& inputs, int32_t current_rho, const BBIIParams& params);
}
#endif
