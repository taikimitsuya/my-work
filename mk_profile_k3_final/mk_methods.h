#ifndef MK_METHODS_H
#define MK_METHODS_H
#include "mk_tfhe_structs.h"
#include "mk_params.h"
#include "mk_tfhe_structs.h"


// DFTベースBlind Rotate（BBII本体）
void mk_blind_rotate_dft(MKRLweSample* acc, const MKRLweSample* bk_input, const MKBootstrappingKey* mk_bk, const TFheGateBootstrappingParameterSet* params);
void mk_lwe_sym_encrypt(MKRLweSample* result, Torus32 message, const MKSecretKey* sk, int32_t party_id, const TFheGateBootstrappingParameterSet* params, int32_t n_per_party);
Torus32 mk_lwe_decrypt(const MKRLweSample* ciphertext, const std::vector<MKSecretKey*>& all_keys, const TFheGateBootstrappingParameterSet* params);
void mk_batch_bootstrapping(
	std::vector<MKRLweSample*>& results,
	const std::vector<MKRLweSample*>& inputs,
	const MKBootstrappingKey* mk_bk,
	const std::vector<Torus32>& mus,
	const TFheGateBootstrappingParameterSet* params
);
#endif
