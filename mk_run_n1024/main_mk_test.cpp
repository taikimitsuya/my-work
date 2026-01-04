#include <iostream>
#include <vector>
#include "mk_methods.h" 
using namespace std;
int main() {
    cout << "=== MK-TFHE N=1024 Local Test ===" << endl;
    int32_t k = 2; int32_t d = 3; int32_t rho = 4; int32_t N = 1024;
    cout << "Params: N=" << N << ", k=" << k << endl;

    bbii::MKParams* mp = bbii::get_mk_test_params(k, d, rho, N);
    vector<bbii::MKSecretKey*> sks(k); 
    bbii::MKBootstrappingKey* bk = new bbii::MKBootstrappingKey(k, mp->n_per_party, mp->get_tfhe_params());
    
    cout << "[Step 1] Generating Keys..." << endl;
    for(int i=0; i<k; ++i){ 
        sks[i] = new bbii::MKSecretKey(mp->get_tfhe_params()); 
        bk->generateKeyForParty(i, sks[i], mp->get_tfhe_params()); 
    }

    cout << "[Step 2] Encrypting..." << endl;
    bbii::MKLweSample* in = new bbii::MKLweSample(k, mp->n_per_party, mp->get_tfhe_params());
    bbii::mk_lwe_sym_encrypt(in, 1000, sks[0], 0, mp->get_tfhe_params());

    bbii::MKLweSample* out = new bbii::MKLweSample(k, N, mp->get_tfhe_params());
    
    cout << "[Step 3] Bootstrapping (Heavy Calculation)..." << endl;
    bbii::mk_bootstrapping(out, in, bk, 1000, mp->get_tfhe_params());
    
    cout << "=== SUCCESS: Finished! ===" << endl;
    return 0;
}
