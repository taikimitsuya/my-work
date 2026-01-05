#include <iostream>
#include <vector>
#include <chrono>
#include <iomanip>
#include "mk_methods.h"
#include "mk_profiler.h"

using namespace std;

MKProfiler global_profiler;

int main() {
    cout << "=== MK-TFHE Time Profiling (k=2) FIXED ===" << endl;
    
    int32_t k = 2; 
    int32_t d = 3; 
    int32_t rho = 4; 
    int32_t N = 1024; 
    
    cout << "Initializing Params..." << endl;
    bbii::MKParams* mp = bbii::get_mk_test_params(k, d, rho, N);
    bbii::BBIIParams* bp = mp->sk_params;

    cout << "==================================" << endl;
    cout << "  Parties (k)        : " << k << endl;
    cout << "  TLWE Dimension (N) : " << bp->N << endl;
    cout << "==================================" << endl;

    cout << "Generating Keys..." << endl;
    vector<bbii::MKSecretKey*> sks(k); 
    bbii::MKBootstrappingKey* bk = new bbii::MKBootstrappingKey(k, mp->n_per_party, mp->get_tfhe_params());
    
    auto kg_start = std::chrono::high_resolution_clock::now();
    for(int i=0; i<k; ++i){ 
        sks[i] = new bbii::MKSecretKey(mp->get_tfhe_params()); 
        bk->generateKeyForParty(i, sks[i], mp->get_tfhe_params()); 
    }
    auto kg_end = std::chrono::high_resolution_clock::now();
    global_profiler.time_keygen = std::chrono::duration<double, std::milli>(kg_end - kg_start).count();

    cout << "Encrypting..." << endl;
    bbii::MKLweSample* in = new bbii::MKLweSample(k, mp->n_per_party, mp->get_tfhe_params());
    bbii::mk_lwe_sym_encrypt(in, 1000, sks[0], 0, mp->get_tfhe_params());

    bbii::MKLweSample* out = new bbii::MKLweSample(k, N, mp->get_tfhe_params());
    
    cout << "Running Bootstrapping (This will take time)..." << endl;
    auto bs_start = std::chrono::high_resolution_clock::now();
    bbii::mk_bootstrapping(out, in, bk, 1000, mp->get_tfhe_params());
    auto bs_end = std::chrono::high_resolution_clock::now();
    double total_bs_time = std::chrono::duration<double, std::milli>(bs_end - bs_start).count();

    cout << endl;
    cout << "==================================" << endl;
    cout << "         [Time Profile]           " << endl;
    cout << "==================================" << endl;
    cout << fixed << setprecision(2);
    
    cout << "  0. Key Generation        : " << global_profiler.time_keygen << " ms" << endl;
    cout << "  --------------------------------" << endl;
    cout << "  1. Input Packing         : " << global_profiler.time_input_packing << " ms" << endl;
    cout << "  2. Blind Rotate (Control): " << global_profiler.time_blind_rotate_control << " ms" << endl;
    cout << "  3. Homomorphic DFT/Mult  : " << global_profiler.time_external_product << " ms" << endl;
    cout << "  4. Sample Extract        : " << global_profiler.time_sample_extract << " ms" << endl;
    cout << "----------------------------------" << endl;
    cout << "Total Execution Time       : " << total_bs_time << " ms" << endl;
    cout << "==================================" << endl;

    return 0;
}
