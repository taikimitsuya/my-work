#include <iostream>
#include <vector>
#include <cmath>
#include <chrono>
#include "mk_params.h"
#include "mk_tfhe_structs.h"

// Explicitly declare functions inside the bbii namespace so main() can find them
namespace bbii {
    void mk_lwe_sym_encrypt(MKLweSample* result, Torus32 message, const MKSecretKey* sk, int32_t party_id, const TFheGateBootstrappingParameterSet* params);
    Torus32 mk_lwe_decrypt(const MKLweSample* ciphertext, const std::vector<MKSecretKey*>& all_keys, const TFheGateBootstrappingParameterSet* params);
    void mk_bootstrapping(MKLweSample* result, const MKLweSample* input, const MKBootstrappingKey* mk_bk, Torus32 mu, const TFheGateBootstrappingParameterSet* params);
}

using namespace std;

int main() {
    cout << "=== Multi-Key TFHE Benchmark ===" << endl;

    int32_t k = 2; 
    int32_t d = 3; 
    int32_t rho = 4; 
    int32_t N = 1024;
    
    cout << "Params: k=" << k << ", N=" << N << ", d=" << d << ", rho=" << rho << endl;

    cout << "Generating Keys..." << endl;
    auto start_setup = chrono::high_resolution_clock::now();

    bbii::MKParams* mp = bbii::get_mk_test_params(k, d, rho, N);
    
    vector<bbii::MKSecretKey*> sks(k); 
    bbii::MKBootstrappingKey* bk = new bbii::MKBootstrappingKey(k, mp->n_per_party, mp->get_tfhe_params());
    
    for(int i=0; i<k; ++i){ 
        sks[i] = new bbii::MKSecretKey(mp->get_tfhe_params()); 
        bk->generateKeyForParty(i, sks[i], mp->get_tfhe_params()); 
    }

    auto end_setup = chrono::high_resolution_clock::now();
    cout << "  [Time] KeyGen (Total): " << chrono::duration<double, milli>(end_setup - start_setup).count() << " ms" << endl;

    bbii::MKLweSample* in = new bbii::MKLweSample(k, mp->n_per_party, mp->get_tfhe_params());
    
    // Call functions explicitly with bbii::
    bbii::mk_lwe_sym_encrypt(in, dtot32(0.25), sks[0], 0, mp->get_tfhe_params());

    cout << "Running MK Bootstrapping..." << endl;
    
    bbii::MKLweSample* out = new bbii::MKLweSample(k, N, mp->get_tfhe_params());
    
    auto start_boot = chrono::high_resolution_clock::now();
    
    bbii::mk_bootstrapping(out, in, bk, dtot32(0.5), mp->get_tfhe_params());
    
    auto end_boot = chrono::high_resolution_clock::now();

    cout << "  [Granularity 1] Bootstrapping (End-to-End): " << chrono::duration<double, milli>(end_boot - start_boot).count() << " ms" << endl;

    return 0;
}
