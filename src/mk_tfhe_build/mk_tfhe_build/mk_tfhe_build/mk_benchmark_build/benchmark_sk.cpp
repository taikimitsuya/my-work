#include <iostream>
#include <chrono>
#include <vector>
#include <tfhe.h>
#include <tfhe_io.h>

using namespace std;

int main() {
    cout << "=== Single-Key TFHE Benchmark (Baseline) ===" << endl;
    
    const int minimum_lambda = 100; 
    TFheGateBootstrappingParameterSet* params = new_default_gate_bootstrapping_parameters(minimum_lambda);

    cout << "Params: N=" << params->tgsw_params->tlwe_params->N 
         << ", n=" << params->in_out_params->n << endl;

    cout << "Generating Keys..." << endl;
    auto start_keygen = chrono::high_resolution_clock::now();
    TFheGateBootstrappingSecretKeySet* key = new_random_gate_bootstrapping_secret_keyset(params);
    auto end_keygen = chrono::high_resolution_clock::now();
    cout << "  [Time] KeyGen: " << chrono::duration<double, milli>(end_keygen - start_keygen).count() << " ms" << endl;

    LweSample* ciphertext = new_gate_bootstrapping_ciphertext(params);
    bootsSymEncrypt(ciphertext, 1, key);

    LweSample* result = new_gate_bootstrapping_ciphertext(params);
    
    cout << "Running Bootstrapping (NAND gate)..." << endl;
    auto start_boot = chrono::high_resolution_clock::now();

    bootsNAND(result, ciphertext, ciphertext, key);

    auto end_boot = chrono::high_resolution_clock::now();
    cout << "  [Time] Single-Key Bootstrapping: " << chrono::duration<double, milli>(end_boot - start_boot).count() << " ms" << endl;

    delete_gate_bootstrapping_ciphertext(ciphertext);
    delete_gate_bootstrapping_ciphertext(result);
    delete_gate_bootstrapping_secret_keyset(key);
    delete_gate_bootstrapping_parameters(params);

    return 0;
}
