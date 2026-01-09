#ifndef MK_METHODS_H
#define MK_METHODS_H
#include "mk_tfhe_structs.h"
#include "mk_params.h"
namespace bbii {
    void mk_lwe_sym_encrypt(MKRLweSample* result, Torus32 message, const MKSecretKey* sk, int32_t party_id, const TFheGateBootstrappingParameterSet* params, int32_t n_per_party);
    Torus32 mk_lwe_decrypt(const MKRLweSample* ciphertext, const std::vector<MKSecretKey*>& all_keys, const TFheGateBootstrappingParameterSet* params);
    void mk_bootstrapping(MKRLweSample* result, const MKRLweSample* input, const MKBootstrappingKey* mk_bk, Torus32 mu, const TFheGateBootstrappingParameterSet* params);
}
#endif
