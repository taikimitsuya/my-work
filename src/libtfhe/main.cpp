#include "batch_bootstrapping.h"
#include <iostream>
#include <vector>

using namespace bbii;

int main() {
    // 1. パラメータ (n=16, d=2, rho=3)
    BBIIParams* params = get_test_params();
    std::cout << "Params: n=" << params->n << ", N=" << params->N << std::endl;

    // 2. 鍵生成 (TFHE鍵)
    TFheGateBootstrappingSecretKeySet* key = new_random_gate_bootstrapping_secret_keyset(params->tfhe_params);

    // 3. Batch Bootstrapping Key 生成 (Stub)
    BatchBootstrappingKey bk;
    // ダミーの鍵を1つ生成
    int vec_len = params->n; // 簡易的な長さ
    std::vector<PackedTRGSW> block_key;
    for(int i=0; i<vec_len; ++i) {
        TRGSW* c = new_TRGSW_array(1, params->tfhe_params->tgsw_params);
        // 0を暗号化
        trgswSymEncryptInt(c, 0, 0, key->tgsw_key); 
        block_key.push_back(PackedTRGSW(c, BatchMode::R12));
    }
    bk.keys.push_back(block_key);

    // 4. 入力暗号文生成
    std::vector<LweSample*> inputs(params->n);
    for(int i=0; i<params->n; ++i) {
        inputs[i] = new_LweSample(params->tfhe_params->in_out_params);
        // Message: 0 or 1
        int msg = (i % 2); 
        // LWE Encrypt: mu = msg * (1/8) などのスケーリング
        bootsSymEncrypt(inputs[i], msg, key);
    }

    std::cout << "Running Batch Bootstrapping..." << std::endl;
    
    // 5. 実行
    auto outputs = batch_bootstrapping(inputs, bk, *params);

    std::cout << "Done. Verifying..." << std::endl;

    // 6. 検証
    int correct = 0;
    for(int i=0; i<params->n; ++i) {
        int dec = bootsSymDecrypt(outputs[i], key);
        if (dec == (i%2)) correct++;
    }

    std::cout << "Accuracy (Stub): " << correct << "/" << params->n << std::endl;

    // Cleanup
    delete_gate_bootstrapping_secret_keyset(key);
    delete params;
    
    return 0;
}