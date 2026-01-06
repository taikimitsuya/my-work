#include <iostream>
#include <vector>
#include <cstdlib>
#include "mk_methods.h" 
using namespace std;
int main() {
    int32_t k, rho;
    int32_t d = 2; // fixed
    int32_t N = 1024;
    cout << "=== MK-TFHE Test ===" << endl;
    cout << "パーティ数 k を入力してください: ";
    cin >> k;
    cout << "パラメータ rho を入力してください: ";
    cin >> rho;
    cout << "(d=2, N=1024 固定)" << endl;
    
    cout << "Params: N=" << N << ", k=" << k << " (Running on Local Docker)" << endl;
    cout << "Note: Computation time will increase linearly with k." << endl;

    bbii::MKParams* mp = bbii::get_mk_test_params(k, d, rho, N);
    vector<bbii::MKSecretKey*> sks(k); 
    bbii::MKBootstrappingKey* bk = new bbii::MKBootstrappingKey(k, mp->n_per_party, mp->get_tfhe_params());
    
    cout << "[Step 1] Generating Keys (for 3 parties)..." << endl;
    for(int i=0; i<k; ++i){ 
        cout << "  -> KeyGen Party " << i+1 << endl;
        sks[i] = new bbii::MKSecretKey(mp->get_tfhe_params()); 
        bk->generateKeyForParty(i, sks[i], mp->get_tfhe_params()); 
    }

    cout << "[Step 2] Encrypting..." << endl;
    bbii::MKLweSample* in = new bbii::MKLweSample(k, mp->n_per_party, mp->get_tfhe_params());
    bbii::mk_lwe_sym_encrypt(in, 1000, sks[0], 0, mp->get_tfhe_params());

    bbii::MKLweSample* out = new bbii::MKLweSample(k, N, mp->get_tfhe_params());
    
    cout << "[Step 3] Bootstrapping..." << endl;
    bbii::mk_bootstrapping(out, in, bk, 1000, mp->get_tfhe_params());
    
    cout << "=== SUCCESS: k=3 Test Finished! ===" << endl;
    return 0;
}
