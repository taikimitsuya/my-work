#ifndef MK_PARAMS_H
#define MK_PARAMS_H
#include "bb_params.h"
#include <vector>
#include <stdexcept>
namespace bbii {
struct MKParams {
    int32_t k; int32_t n_per_party; int32_t total_n; int32_t N; BBIIParams* sk_params;
    MKParams(int32_t parties, BBIIParams* base_params) : k(parties), sk_params(base_params) {
        this->n_per_party = base_params->n; this->N = base_params->N;
        this->total_n = this->k * this->n_per_party;
    }
    ~MKParams() { }
    TFheGateBootstrappingParameterSet* get_tfhe_params() const { return sk_params->tfhe_params; }
};
inline MKParams* get_mk_test_params(int32_t k, int32_t d, int32_t rho, int32_t N) {
    return new MKParams(k, new BBIIParams(d, rho, N));
}
} 
#endif
