#ifndef BATCH_BOOTSTRAPPING_H
#define BATCH_BOOTSTRAPPING_H
#include "hom_dft.h"
namespace bbii {
struct BatchBootstrappingKey { std::vector<std::vector<PackedTRGSW>> keys; };
std::vector<LweSample*> batch_bootstrapping(const std::vector<LweSample*>& inputs, const BatchBootstrappingKey& bk, const BBIIParams& params);
}
#endif
