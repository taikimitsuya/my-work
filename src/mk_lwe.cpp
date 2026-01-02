#include "mk_tfhe_structs.h"
#include <vector>
#include <cmath>
#include <random>
#include <limits>
#include <tfhe.h>
#include <tfhe_core.h>

namespace bbii {

Torus32 get_uniform_random_torus32() {
    static std::random_device rd;
    static std::mt19937 gen(rd());
    static std::uniform_int_distribution<int32_t> dist(
        std::numeric_limits<int32_t>::min(), 
        std::numeric_limits<int32_t>::max()
    );
    return dist(gen);
}

void mk_lwe_sym_encrypt(MKLweSample* result, 
                        Torus32 message, 
                        const MKSecretKey* sk, 
                        int32_t party_id, 
                        const TFheGateBootstrappingParameterSet* params) {
    const LweParams* lwe_params = params->in_out_params;
    int32_t n = result->n_per_party; 
    int32_t k = result->k;           
    int32_t total_n = n * k;         

    for (int i = 0; i < total_n; ++i) result->sample->a[i] = 0;

    int32_t offset = party_id * n;
    Torus32 multi_key_product = 0;

    for (int i = 0; i < n; ++i) {
        Torus32 a_val = get_uniform_random_torus32();
        result->sample->a[offset + i] = a_val;
        if (sk->lwe_key->key[i] == 1) {
            multi_key_product += a_val;
        }
    }

    double alpha = lwe_params->alpha_min;
    Torus32 error = gaussian32(0, alpha); 
    result->sample->b = multi_key_product + message + error;
    result->sample->current_variance = alpha * alpha;
}

Torus32 mk_lwe_decrypt(const MKLweSample* ciphertext, 
                       const std::vector<MKSecretKey*>& all_keys, 
                       const TFheGateBootstrappingParameterSet* params) {
    int32_t n = ciphertext->n_per_party;
    int32_t k = ciphertext->k;
    if (all_keys.size() != (size_t)k) throw std::runtime_error("Key mismatch.");

    Torus32 phase = ciphertext->sample->b;
    for (int j = 0; j < k; ++j) {
        int32_t offset = j * n;
        const LweKey* current_sk = all_keys[j]->lwe_key;
        for (int i = 0; i < n; ++i) {
            Torus32 a_val = ciphertext->sample->a[offset + i];
            if (current_sk->key[i] == 1) phase -= a_val;
        }
    }
    return phase;
}

void mk_lwe_noiseless_trivial(MKLweSample* result, Torus32 message) {
    int32_t total_n = result->k * result->n_per_party;
    for (int i = 0; i < total_n; ++i) result->sample->a[i] = 0;
    result->sample->b = message;
    result->sample->current_variance = 0.0;
}

} 
