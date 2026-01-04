#ifndef MK_METHODS_H
#define MK_METHODS_H
#include "mk_tfhe_structs.h"
#include "mk_params.h"
namespace bbii {
    void mk_lwe_sym_encrypt(MKLweSample* result, Torus32 message, const MKSecretKey* sk, int32_t party_id, const TFheGateBootstrappingParameterSet* params);
    Torus32 mk_lwe_decrypt(const MKLweSample* ciphertext, const std::vector<MKSecretKey*>& all_keys, const TFheGateBootstrappingParameterSet* params);
    void mk_bootstrapping(MKLweSample* result, const MKLweSample* input, const MKBootstrappingKey* mk_bk, Torus32 mu, const TFheGateBootstrappingParameterSet* params);
}
#endif
