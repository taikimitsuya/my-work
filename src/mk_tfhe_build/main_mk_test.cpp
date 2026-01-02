#include <iostream>
#include <vector>
#include <cmath>
#include <limits>
#include "mk_params.h"
#include "mk_tfhe_structs.h"

namespace bbii {
    void mk_lwe_sym_encrypt(MKLweSample*, Torus32, const MKSecretKey*, int32_t, const TFheGateBootstrappingParameterSet*);
    Torus32 mk_lwe_decrypt(const MKLweSample*, const std::vector<MKSecretKey*>&, const TFheGateBootstrappingParameterSet*);
    void mk_bootstrapping(MKLweSample*, const MKLweSample*, const MKBootstrappingKey*, Torus32, const TFheGateBootstrappingParameterSet*);
}

using namespace bbii; 
using namespace std;

// 結果検証用関数
double verify(const bbii::MKLweSample* ex, const vector<bbii::MKSecretKey*>& k, const bbii::MKParams* p) {
    Torus32 ph = ex->sample->b; 
    int32_t N = p->N;
    for(int u=0; u<p->k; ++u) {
        for(int i=0; i<N; ++i) {
            ph -= ex->sample->a[u*N+i] * k[u]->rlwe_key->key->coefs[i];
        }
    }
    return t32tod(ph);
}

int main() {
    int32_t k, d, rho, N;

    cout << "========================================" << endl;
    cout << "   Multi-Key FHE (MK-TFHE) Parameter Setup" << endl;
    cout << "========================================" << endl;

    // 1. パーティ数 (k)
    cout << "Enter number of parties (k) [e.g., 2]: ";
    while(!(cin >> k) || k < 1) {
        cout << "Invalid input. k must be >= 1: ";
        cin.clear(); cin.ignore(numeric_limits<streamsize>::max(), '\n');
    }

    // 2. 分解の底 (d)
    cout << "Enter decomposition base (d) [e.g., 2]: ";
    while(!(cin >> d) || d < 2) {
        cout << "Invalid input. d must be >= 2: ";
        cin.clear(); cin.ignore(numeric_limits<streamsize>::max(), '\n');
    }

    // 3. 分解の深さ (rho)
    cout << "Enter decomposition depth (rho) [e.g., 5]: ";
    while(!(cin >> rho) || rho < 1) {
        cout << "Invalid input. rho must be >= 1: ";
        cin.clear(); cin.ignore(numeric_limits<streamsize>::max(), '\n');
    }

    // 4. 多項式の次数 (N)
    cout << "Enter polynomial degree (N) [e.g., 1024]: ";
    while(!(cin >> N) || N < 1) {
        cout << "Invalid input. N must be >= 1: ";
        cin.clear(); cin.ignore(numeric_limits<streamsize>::max(), '\n');
    }

    // 計算されるLWE次元の確認
    int32_t n = 2 * pow(d, rho);
    cout << "\n----------------------------------------" << endl;
    cout << " Configuration Summary:" << endl;
    cout << "  Parties (k)      : " << k << endl;
    cout << "  LWE Dimension (n): " << n << " (calculated as 2 * d^rho)" << endl;
    cout << "  RLWE Degree (N)  : " << N << endl;
    cout << "----------------------------------------" << endl;
    cout << "Generating Parameters and Keys..." << endl;

    // パラメータ生成
    MKParams* mp = get_mk_test_params(k, d, rho, N);
    
    // 鍵生成
    vector<MKSecretKey*> sks(k); 
    MKBootstrappingKey* bk = new MKBootstrappingKey(k, mp->n_per_party, mp->get_tfhe_params());
    
    for(int i=0; i<k; ++i){ 
        sks[i] = new MKSecretKey(mp->get_tfhe_params()); 
        bk->generateKeyForParty(i, sks[i], mp->get_tfhe_params()); 
    }
    
    // 暗号化 (Party 0)
    cout << "Encrypting message (m=0.25) by Party 0..." << endl;
    bbii::MKLweSample* in = new bbii::MKLweSample(k, mp->n_per_party, mp->get_tfhe_params());
    
    mk_lwe_sym_encrypt(in, dtot32(0.25), sks[0], 0, mp->get_tfhe_params());
    cout << "  -> Input Decrypted Check: " << t32tod(mk_lwe_decrypt(in, sks, mp->get_tfhe_params())) << endl;
    
    // Bootstrapping
    cout << "Running MK Bootstrapping (with test vector m=0.5)..." << endl;
    bbii::MKLweSample* out = new bbii::MKLweSample(k, N, mp->get_tfhe_params());
    
    mk_bootstrapping(out, in, bk, dtot32(0.5), mp->get_tfhe_params());
    
    // 結果検証
    cout << "  -> Result Phase: " << verify(out, sks, mp) << endl;
    cout << "Test Finished." << endl;
    
    return 0;
}
