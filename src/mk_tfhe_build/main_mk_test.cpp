#include <iostream>
#include <vector>
#include <cmath>
#include <cassert>

#include "mk_params.h"
#include "mk_tfhe_structs.h"

namespace bbii {
    void mk_lwe_sym_encrypt(MKLweSample* result, Torus32 message, const MKSecretKey* sk, int32_t party_id, const TFheGateBootstrappingParameterSet* params);
    Torus32 mk_lwe_decrypt(const MKLweSample* ciphertext, const std::vector<MKSecretKey*>& all_keys, const TFheGateBootstrappingParameterSet* params);
    void mk_bootstrapping(MKLweSample* result, const MKLweSample* input, const MKBootstrappingKey* mk_bk, Torus32 mu, const TFheGateBootstrappingParameterSet* params);
}

using namespace bbii;
using namespace std;

double verify_extracted_sample(const bbii::MKLweSample* extracted, 
                               const vector<bbii::MKSecretKey*>& keys, 
                               const bbii::MKParams* mk_p) {
    int32_t N = mk_p->N;
    int32_t k = mk_p->k;
    
    Torus32 phase = extracted->sample->b;

    for (int u = 0; u < k; ++u) {
        const IntPolynomial* s_poly = keys[u]->rlwe_key->key; // Fixed: removed &
        int offset = u * N;
        for (int i = 0; i < N; ++i) {
            Torus32 a_val = extracted->sample->a[offset + i];
            int32_t s_val = s_poly->coefs[i]; 
            phase -= a_val * s_val;
        }
    }
    return t32tod(phase);
}

int main() {
    cout << "=== Multi-Key FHE (MK-TFHE) Test Start ===" << endl;

    int32_t k = 2;
    int32_t d = 2;
    int32_t rho = 5; 
    int32_t N = 1024;
    
    cout << "Generating Parameters (k=" << k << ", n=" << 2*pow(d,rho) << ", N=" << N << ")..." << endl;
    MKParams* mk_params = get_mk_test_params(k, d, rho, N);
    const TFheGateBootstrappingParameterSet* tfhe_params = mk_params->get_tfhe_params();

    cout << "Generating Keys..." << endl;
    vector<MKSecretKey*> secret_keys(k);
    MKBootstrappingKey* mk_bk = new MKBootstrappingKey(k, mk_params->n_per_party, tfhe_params);

    for(int i=0; i<k; ++i) {
        secret_keys[i] = new MKSecretKey(tfhe_params);
        mk_bk->generateKeyForParty(i, secret_keys[i], tfhe_params);
    }

    cout << "Encrypting message..." << endl;
    MKLweSample* input_ct = new MKLweSample(k, mk_params->n_per_party, tfhe_params);
    
    double message_double = 0.25;
    Torus32 message = dtot32(message_double);
    mk_lwe_sym_encrypt(input_ct, message, secret_keys[0], 0, tfhe_params);

    Torus32 decrypted_raw = mk_lwe_decrypt(input_ct, secret_keys, tfhe_params);
    cout << "  Input Decrypted Phase: " << t32tod(decrypted_raw) 
         << " (Expected: " << message_double << ")" << endl;

    cout << "Running MK Bootstrapping (This may take time)..." << endl;
    MKLweSample* output_ct = new MKLweSample(k, N, tfhe_params);
    
    mk_bootstrapping(output_ct, input_ct, mk_bk, dtot32(0.5), tfhe_params);

    double result_phase = verify_extracted_sample(output_ct, secret_keys, mk_params);
    
    cout << "Result Phase: " << result_phase << endl;
    cout << "Test Finished Successfully." << endl;

    delete input_ct;
    delete output_ct;
    delete mk_bk;
    delete mk_params;
    for(auto sk : secret_keys) delete sk;

    return 0;
}
