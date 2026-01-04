#include <iostream>
#include <vector>
#include "mk_methods.h" 
using namespace std;
int main() {
    cout << "=== Start Minimal Memory Test ===" << endl;
    
    // パラメータを極小にする (N=1024だと落ちるため 500 に)
    // k=1 (シングルキー)でまず動くか試す
    int32_t k = 1; 
    int32_t d = 2; 
    int32_t rho = 2; 
    int32_t N = 500; 
    
    cout << "Params: N=" << N << ", k=" << k << endl;

    bbii::MKParams* mp = bbii::get_mk_test_params(k, d, rho, N);
    vector<bbii::MKSecretKey*> sks(k); 
    bbii::MKBootstrappingKey* bk = new bbii::MKBootstrappingKey(k, mp->n_per_party, mp->get_tfhe_params());
    
    cout << "LOG: Generating Keys..." << endl;
    for(int i=0; i<k; ++i){ 
        sks[i] = new bbii::MKSecretKey(mp->get_tfhe_params()); 
        bk->generateKeyForParty(i, sks[i], mp->get_tfhe_params()); 
    }

    cout << "LOG: Encrypting..." << endl;
    bbii::MKLweSample* in = new bbii::MKLweSample(k, mp->n_per_party, mp->get_tfhe_params());
    bbii::mk_lwe_sym_encrypt(in, 1000, sks[0], 0, mp->get_tfhe_params());

    bbii::MKLweSample* out = new bbii::MKLweSample(k, N, mp->get_tfhe_params());
    
    cout << "LOG: Bootstrapping..." << endl;
    bbii::mk_bootstrapping(out, in, bk, 1000, mp->get_tfhe_params());
    
    cout << "=== SUCCESS: Finished Main ===" << endl;
    return 0;
}
