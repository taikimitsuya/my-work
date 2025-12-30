#include <iostream>
#include <vector>
#include <cassert>
#include <chrono>
#include <random>

// MK-TFHE ライブラリ (パスは環境に合わせて変更)
#include "tfhe.h"
#include "tfhe_io.h"

// BBII プロジェクトヘッダー
#include "bb_params.h"
#include "bb_utils.h"
#include "batch_framework.h"
#include "batch_ops.h"
#include "batch_bootstrapping.h"

using namespace bbii;
using namespace std;

// --- Helper: テスト用の鍵生成 (簡易版) ---
// 本来は論文 Section 7.2 に従い、秘密鍵 s のDFT変換結果を回転行列として暗号化する必要がある。
// ここでは構造体の初期化と、ダミーデータの投入のみを行うスタブ関数とする。
void generate_test_keys(LweKey*& lwe_key, BatchBootstrappingKey& bk, const BBIIParams& params) {
    cout << "Generating keys..." << endl;

    // 1. 標準的な TFHE 秘密鍵の生成
    lwe_key = new_lwe_key(params.lwe_params);
    lweKeyGen(lwe_key);

    // 2. Batch Bootstrapping Key (BK) の構造確保
    // 論文: v' = (2d)^(rho-1) / v
    int32_t total_blocks = std::pow(2 * params.d, params.rho - 1);
    int32_t v_prime = total_blocks / params.v;
    int32_t decomp_size = 2 * params.d * 32; // 仮: log(q) approx 32 (実際にはパラメータ依存)
    
    bk.keys.resize(v_prime);
    for (int i = 0; i < v_prime; ++i) {
        bk.keys[i].resize(decomp_size);
        for (int j = 0; j < decomp_size; ++j) {
            // 各エントリは PackedTRGSW
            // 実際にはここで「秘密鍵の回転」を RGSW 暗号化したものをセットする。
            // (実装負荷が高いため、ここでは 0 の暗号化で代用する)
            TRGSW* dummy_cipher = new_trgsw_ciphertext(params.rgsw_params);
            trgswSymEncryptInt(dummy_cipher, 0, 0, lwe_key); // 0 を暗号化
            
            // モードは R12->R13 等、Algorithm 4.1 の要件に合わせる
            bk.keys[i][j] = PackedTRGSW(dummy_cipher, BatchMode::R12_to_R13);
        }
    }
    
    cout << "Key generation (stub) completed." << endl;
}

// --- Helper: ノイズ計測 ---
double check_noise(LweSample* ciphertext, int32_t expected_message, LweKey* key, int32_t mod_size) {
    // 復号 (Phase計算)
    Torus32 phase = lwePhase(ciphertext, key);
    // メッセージの期待値 (Torus)
    Torus32 expected_torus = dtot32(double(expected_message) / mod_size);
    // 差分
    Torus32 diff = phase - expected_torus;
    if (diff > 1073741824) diff = -diff; // Abs check for Torus32
    
    // おおよそのノイズ量を返す
    return t32tod(diff);
}

int main() {
    // 1. パラメータ設定 (テスト用の軽量パラメータ)
    // d=2, rho=3 -> n=16, r=512 (if N=1024)
    BBIIParams* params = get_toy_bbii_params();
    
    cout << "=== Batch Bootstrapping II Implementation Test ===" << endl;
    cout << "Parameters: n=" << params->n 
         << ", d=" << params->d 
         << ", rho=" << params->rho 
         << ", r=" << params->r << endl;

    // 2. 鍵生成
    LweKey* secret_key;
    BatchBootstrappingKey bk;
    generate_test_keys(secret_key, bk, *params);

    // 3. 入力データの作成 (n 個の LWE 暗号文)
    int n = params->n;
    std::vector<LweSample*> inputs(n);
    std::vector<int32_t> original_messages(n);

    // 乱数生成器
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(0, 1); // 1bit message

    cout << "Encrypting " << n << " input ciphertexts..." << endl;
    for (int i = 0; i < n; ++i) {
        inputs[i] = new_lwe_sample(params->lwe_params);
        original_messages[i] = dis(gen);
        
        // メッセージを Torus に変換 (1/4 or -1/4 などをスケーリング)
        // ここでは単純に 0 -> 0, 1 -> 1/2 とする (Bootstrapping Gateによる)
        Torus32 mu = original_messages[i] ? 1073741824 : -1073741824; // Modulo switch trick
        
        // 高めのノイズを入れて暗号化 (ブートストラップの効果を確認するため)
        double alpha = 0.01; // Large noise
        lweSymEncrypt(inputs[i], mu, alpha, secret_key);
    }

    // 4. Batch Bootstrapping 実行
    cout << "Running Batch Bootstrapping..." << endl;
    auto start_time = std::chrono::high_resolution_clock::now();

    std::vector<LweSample*> outputs;
    try {
        outputs = batch_bootstrapping(inputs, bk, *params);
    } catch (const std::exception& e) {
        cerr << "Error during bootstrapping: " << e.what() << endl;
        return 1;
    }

    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end_time - start_time;
    cout << "Finished in " << elapsed.count() << " seconds." << endl;
    cout << "Amortized time: " << elapsed.count() / n << " seconds per LWE." << endl;

    // 5. 結果検証
    int correct_count = 0;
    cout << "\n--- Verification ---" << endl;
    for (int i = 0; i < n; ++i) {
        // 復号 (本来はBootstrapping後のキーで復号するが、
        // KeySwitchingを省略している場合は元のキーで復号できる場合もある。
        // Algorithm 7.1 は "Refresh" なので、同じ秘密鍵で復号可能と仮定)
        
        int32_t decrypted = lweSymDecrypt(outputs[i], original_messages[i] ? 1073741824 : -1073741824, secret_key);
        // Note: lweSymDecrypt は位相から近いメッセージを返す

        // 簡易チェック: 復号結果が 0 か 1 か
        // 実際には位相を直接見て判定する
        Torus32 phase = lwePhase(outputs[i], secret_key);
        int32_t dec_msg = (phase > 0) ? 1 : 0; // Simple sign check

        if (dec_msg == original_messages[i]) {
            correct_count++;
        } else {
            // Debug info
            // cout << "Fail at " << i << ": Expected " << original_messages[i] << ", Got " << dec_msg << endl;
        }
    }

    cout << "Accuracy: " << correct_count << " / " << n << " (" << (double)correct_count/n * 100 << "%)" << endl;

    // 6. メモリ解放
    delete_lwe_key(secret_key);
    for (auto c : inputs) delete_lwe_sample(c);
    for (auto c : outputs) delete_lwe_sample(c);
    delete params;

    return (correct_count == n) ? 0 : 1;
}