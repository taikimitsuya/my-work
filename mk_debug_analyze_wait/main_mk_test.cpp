#include <iostream>
#include <vector>
#include "mk_methods.h" 
using namespace std;
int main() {
    cout << "DEBUG: Start Main" << endl;
    int32_t k = 2; int32_t d = 3; int32_t rho = 4; int32_t N = 1024;
    
    bbii::MKParams* mp = bbii::get_mk_test_params(k, d, rho, N);
    vector<bbii::MKSecretKey*> sks(k); 
    bbii::MKBootstrappingKey* bk = new bbii::MKBootstrappingKey(k, mp->n_per_party, mp->get_tfhe_params());
    
    cout << "DEBUG: KeyGen" << endl;
    for(int i=0; i<k; ++i){ 
        sks[i] = new bbii::MKSecretKey(mp->get_tfhe_params()); 
        bk->generateKeyForParty(i, sks[i], mp->get_tfhe_params()); 
    }

    cout << "DEBUG: Encrypt" << endl;
    bbii::MKLweSample* in = new bbii::MKLweSample(k, mp->n_per_party, mp->get_tfhe_params());
    bbii::mk_lwe_sym_encrypt(in, 1000, sks[0], 0, mp->get_tfhe_params());

    bbii::MKLweSample* out = new bbii::MKLweSample(k, N, mp->get_tfhe_params());
    
    cout << "DEBUG: Calling mk_bootstrapping..." << endl;
    bbii::mk_bootstrapping(out, in, bk, 1000, mp->get_tfhe_params());
    
    cout << "DEBUG: Finished Main" << endl;
    return 0;
}
