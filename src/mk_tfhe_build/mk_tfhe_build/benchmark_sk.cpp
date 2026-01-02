#include <iostream>
#include <chrono>
#include <vector>
#include <tfhe/tfhe.h>
#include <tfhe/tfhe_io.h>

using namespace std;

int main() {
    cout << "=== Single-Key TFHE Benchmark (Baseline) ===" << endl;
    
    // パラメータ設定 (MKと条件を合わせるため N=1024 近辺の標準設定を使用)
    // 128bit security settings
    const int minimum_lambda = 100; 
    TFheGateBootstrappingParameterSet* params = new_default_gate_bootstrapping_parameters(minimum_lambda);

    cout << "Params: N=" << params->tgsw_params->tlwe_params->N 
         << ", n=" << params->in_out_params->n << endl;

    // 1. 鍵生成 (KeyGen)
    cout << "Generating Keys..." << endl;
    auto start_keygen = chrono::high_resolution_clock::now();
    
    TFheGateBootstrappingSecretKeySet* key = new_random_gate_bootstrapping_secret_keyset(params);
    
    auto end_keygen = chrono::high_resolution_clock::now();
    double time_keygen = chrono::duration<double, milli>(end_keygen - start_keygen).count();
    cout << "  [Time] KeyGen: " << time_keygen << " ms" << endl;

    // 2. 暗号化 (Encrypt)
    LweSample* ciphertext = new_gate_bootstrapping_ciphertext(params);
    bootsSymEncrypt(ciphertext, 1, key); // encrypt '1'

    // 3. ブートストラップ (Bootstrapping)
    // NANDゲートなどを通すことでブートストラップを実行
    LweSample* result = new_gate_bootstrapping_ciphertext(params);
    
    cout << "Running Bootstrapping (NAND gate)..." << endl;
    auto start_boot = chrono::high_resolution_clock::now();

    // bootsNAND は内部で BlindRotate + KeySwitch を行います
    bootsNAND(result, ciphertext, ciphertext, key);

    auto end_boot = chrono::high_resolution_clock::now();
    double time_boot = chrono::duration<double, milli>(end_boot - start_boot).count();

    cout << "  [Time] Single-Key Bootstrapping: " << time_boot << " ms" << endl;

    // 後始末
    delete_gate_bootstrapping_ciphertext(ciphertext);
    delete_gate_bootstrapping_ciphertext(result);
    delete_gate_bootstrapping_secret_keyset(key);
    delete_gate_bootstrapping_parameters(params);

    return 0;
}
