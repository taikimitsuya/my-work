#include "batch_bootstrapping.h"
#include <iostream>
using namespace bbii;

int main() {
    BBIIParams* params = get_test_params();
    std::cout << "Params: n=" << params->n << ", N=" << params->N << std::endl;
    TFheGateBootstrappingSecretKeySet* key = new_random_gate_bootstrapping_secret_keyset(params->tfhe_params);
    
    // 生成される鍵 (Toy Example)
    BatchBootstrappingKey bk;
    std::vector<PackedTRGSW> block_key;
    for(int i=0; i<params->n; ++i) {
        TGswSample* c = new_TGswSample_array(1, params->tfhe_params->tgsw_params);
        tGswSymEncryptInt(c, 0, 0, key->tgsw_key); 
        block_key.push_back(PackedTRGSW(c, BatchMode::R12));
    }
    bk.keys.push_back(block_key);

    // 入力データ生成
    std::vector<LweSample*> inputs(params->n);
    for(int i=0; i<params->n; ++i) {
        inputs[i] = new_LweSample(params->tfhe_params->in_out_params);
        bootsSymEncrypt(inputs[i], (i%2), key);
    }

    std::cout << "Running Batch Bootstrapping..." << std::endl;
    auto outputs = batch_bootstrapping(inputs, bk, *params);
    
    std::cout << "Done. Checking accuracy..." << std::endl;

    int correct = 0;
    // ★修正: params->n ではなく、実際の出力サイズに合わせてループする
    for(size_t i=0; i<outputs.size(); ++i) {
        int dec = bootsSymDecrypt(outputs[i], key);
        // 入力の i 番目と比較 (i < params->n の範囲内のみ)
        if (i < (size_t)params->n && dec == (i%2)) correct++;
        
        // メモリ解放 (これを忘れるとリークするが、今回のクラッシュ原因ではない)
        delete_LweSample(outputs[i]);
    }
    std::cout << "Accuracy: " << correct << "/" << outputs.size() << std::endl;

    // Cleanup
    for(auto p : inputs) delete_LweSample(p);
    // block_key内のTGswSample解放は省略しているが、実運用では必要
    delete_gate_bootstrapping_secret_keyset(key);
    delete params;
    
    return 0;
}
