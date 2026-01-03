#include <iostream>
#include <vector>
#include <cmath>
#include "mk_methods.h" 

using namespace std;

double verify(const bbii::MKLweSample* extracted, const vector<bbii::MKSecretKey*>& keys, const bbii::MKParams* mk_p) {
    int32_t N = mk_p->N; int32_t k = mk_p->k;
    Torus32 phase = extracted->sample->b;
    for (int u = 0; u < k; ++u) {
        const IntPolynomial* s_poly = keys[u]->rlwe_key->key; 
        for (int i = 0; i < N; ++i) phase -= extracted->sample->a[u*N+i] * s_poly->coefs[i];
    }
    return t32tod(phase);
}

int main() {
    cout << "=== Multi-Key TFHE Test ===" << endl;
    int32_t k = 2; int32_t d = 3; int32_t rho = 4; int32_t N = 1024;
    
    bbii::MKParams* mp = bbii::get_mk_test_params(k, d, rho, N);
    vector<bbii::MKSecretKey*> sks(k); 
    bbii::MKBootstrappingKey* bk = new bbii::MKBootstrappingKey(k, mp->n_per_party, mp->get_tfhe_params());
    
    for(int i=0; i<k; ++i){ 
        sks[i] = new bbii::MKSecretKey(mp->get_tfhe_params()); 
        bk->generateKeyForParty(i, sks[i], mp->get_tfhe_params()); 
    }

    // Explicit bbii::MKLweSample
    bbii::MKLweSample* in = new bbii::MKLweSample(k, mp->n_per_party, mp->get_tfhe_params());
    bbii::mk_lwe_sym_encrypt(in, dtot32(0.25), sks[0], 0, mp->get_tfhe_params());

    cout << "Running MK Bootstrapping..." << endl;
    bbii::MKLweSample* out = new bbii::MKLweSample(k, N, mp->get_tfhe_params());
    
    bbii::mk_bootstrapping(out, in, bk, dtot32(0.5), mp->get_tfhe_params());
    
    cout << "Result: " << verify(out, sks, mp) << endl;
    return 0;
}
