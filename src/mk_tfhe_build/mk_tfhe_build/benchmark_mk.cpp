#include <iostream>
#include <vector>
#include <cmath>
#include <chrono>
#include "mk_params.h"
#include "mk_tfhe_structs.h"

namespace bbii {
    void mk_lwe_sym_encrypt(MKLweSample*, Torus32, const MKSecretKey*, int32_t, const TFheGateBootstrappingParameterSet*);
    Torus32 mk_lwe_decrypt(const MKLweSample*, const std::vector<MKSecretKey*>&, const TFheGateBootstrappingParameterSet*);
    void mk_bootstrapping(MKLweSample*, const MKLweSample*, const MKBootstrappingKey*, Torus32, const TFheGateBootstrappingParameterSet*);
}

using namespace bbii; using namespace std;

int main() {
    cout << "=== Multi-Key TFHE Benchmark (Your Implementation) ===" << endl;

    // 計測用固定パラメータ (比較のため固定)
    int32_t k = 2;    // 2パーティ
    int32_t d = 3;    // 分解の底 (Single Key defaultに近い設定)
    int32_t rho = 4;  // 分解の深さ
    int32_t N = 1024; // RLWE次数

    cout << "Params: k=" << k << ", N=" << N << ", d=" << d << ", rho=" << rho << endl;

    // セットアップ計測
    cout << "Generating Keys..." << endl;
    auto start_setup = chrono::high_resolution_clock::now();

    MKParams* mp = get_mk_test_params(k, d, rho, N);
    vector<MKSecretKey*> sks(k); 
    MKBootstrappingKey* bk = new MKBootstrappingKey(k, mp->n_per_party, mp->get_tfhe_params());
    for(int i=0; i<k; ++i){ 
        sks[i] = new MKSecretKey(mp->get_tfhe_params()); 
        bk->generateKeyForParty(i, sks[i], mp->get_tfhe_params()); 
    }

    auto end_setup = chrono::high_resolution_clock::now();
    cout << "  [Time] KeyGen (Total): " << chrono::duration<double, milli>(end_setup - start_setup).count() << " ms" << endl;

    // 暗号化
    bbii::MKLweSample* in = new bbii::MKLweSample(k, mp->n_per_party, mp->get_tfhe_params());
    mk_lwe_sym_encrypt(in, dtot32(0.25), sks[0], 0, mp->get_tfhe_params());

    // ブートストラップ全体計測
    cout << "Running MK Bootstrapping..." << endl;
    bbii::MKLweSample* out = new bbii::MKLweSample(k, N, mp->get_tfhe_params());
    
    auto start_boot = chrono::high_resolution_clock::now();
    
    mk_bootstrapping(out, in, bk, dtot32(0.5), mp->get_tfhe_params());
    
    auto end_boot = chrono::high_resolution_clock::now();
    double time_boot = chrono::duration<double, milli>(end_boot - start_boot).count();

    cout << "  [Granularity 1] Bootstrapping (End-to-End): " << time_boot << " ms" << endl;

    return 0;
}
