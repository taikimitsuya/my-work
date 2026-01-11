#undef MKKeySwitchKey
#include <iostream>
#include <cstdlib>
#include <vector>
#include <chrono>
#include <iomanip>
#include "mk_methods.h"
#include "mk_tfhe_structs.h"
#include "mk_packed_ops.h"
#include "mk_profiler.h"




int main() {
    // DFTベースBlind Rotateのテスト呼び出し（デバッグ用）
    // extern void mk_test_blind_rotate_dft(const TFheGateBootstrappingParameterSet* params);
    // TFheGateBootstrappingParameterSet* test_params = new_default_gate_bootstrapping_parameters(110);
    // mk_test_blind_rotate_dft(test_params);
    // delete_gate_bootstrapping_parameters(test_params);
    int32_t k, rho;
    int32_t d = 2; // fixed
    int32_t N = 1024;
    int32_t REPEAT;
    std::cout << "=== MK-TFHE Time Profiling ===" << std::endl;
    std::cout << "パーティ数 k を入力してください: ";
    std::cin >> k;
    std::cout << "パラメータ rho を入力してください: ";
    std::cin >> rho;
    std::cout << "試行回数を入力してください: ";
    std::cin >> REPEAT;
    std::cout << "(d=2, N=1024 固定)" << std::endl;

    // 平均用バッファ
    double sum_keygen = 0, sum_inputpack = 0, sum_blind = 0, sum_dft = 0, sum_extract = 0, sum_total = 0;
        for(int rep=1; rep<=REPEAT; ++rep) {
        // 各ループの最初でプロファイラ値をリセット
        global_profiler.time_keygen = 0;
        global_profiler.time_input_packing = 0;
        global_profiler.time_blind_rotate_control = 0;
        global_profiler.time_external_product = 0;
        global_profiler.time_sample_extract = 0;
        std::cout << "\n[Run " << rep << "/" << REPEAT << "]" << std::endl;
        std::cout << "Initializing Params..." << std::endl;
        MKParams* mp = get_mk_test_params(k, d, rho, N);
        BBIIParams* bp = mp->sk_params;

        std::cout << "==================================" << std::endl;
        std::cout << "         [Parameters]            " << std::endl;
        std::cout << "==================================" << std::endl;
        std::cout << "  Parties (k)              : " << k << std::endl;
        std::cout << "  Recursion Depth (d)      : " << d << std::endl;
        std::cout << "  Partition Size (rho)     : " << rho << std::endl;
        std::cout << "  TLWE Dimension (N)       : " << N << std::endl;
        std::cout << "  n per party              : " << mp->n_per_party << std::endl;
        std::cout << "==================================" << std::endl;
        std::cout << std::endl;

        std::cout << "Generating Keys..." << std::endl;
        std::vector<MKSecretKey*> sks(k); 
        auto* bk = new MKBootstrappingKey(k, mp->n_per_party, mp->get_tfhe_params());
        // KeyGen: 差分計測
        auto kg_start = std::chrono::high_resolution_clock::now();
        double keygen0 = global_profiler.time_keygen;
        for(int i=0; i<k; ++i){ 
            sks[i] = new MKSecretKey(mp->get_tfhe_params()); 
            bk->generateKeyForParty(i, sks[i], mp->get_tfhe_params()); 
        }
        // KSK生成（Automorphism用KeySwitchingKeyを正しくセット）
        // 例: delta=1用KSKを生成してキャッシュに登録
        int delta = 1; // 必要に応じて変更
        auto* ksk = new BBII_KSKStruct(k, mp->n_per_party, mp->get_tfhe_params());
        mk_fill_automorphism_ksk(ksk, sks, mp->get_tfhe_params());
        bk->ksk_cache[delta] = ksk;
        auto kg_end = std::chrono::high_resolution_clock::now();
        double keygen1 = global_profiler.time_keygen;
        global_profiler.time_keygen = (keygen1 - keygen0) + std::chrono::duration<double, std::milli>(kg_end - kg_start).count();

        std::cout << "Encrypting..." << std::endl;
        // バッチ用入力生成
        int batch_size = 8;
        std::vector<Torus32> plains(batch_size, 1000);
        std::vector<MKRLweSample*> ins(batch_size);
        for(int i=0;i<batch_size;++i) {
            ins[i] = new MKRLweSample(k, mp->get_tfhe_params());
            mk_lwe_sym_encrypt(ins[i], plains[i], sks[0], 0, mp->get_tfhe_params(), mp->n_per_party);
        }

        std::vector<MKRLweSample*> outs(batch_size);
        for(int i=0;i<batch_size;++i) outs[i] = new MKRLweSample(k, mp->get_tfhe_params());

        std::cout << "Running Batch Bootstrapping (This will take time)..." << std::endl;
        auto bs_start = std::chrono::high_resolution_clock::now();
        mk_batch_bootstrapping(outs, ins, bk, plains, mp->get_tfhe_params());
        auto bs_end = std::chrono::high_resolution_clock::now();
        double total_bs_time = std::chrono::duration<double, std::milli>(bs_end - bs_start).count();

        for(int b=0;b<batch_size;++b) {
            std::cout << "[DEBUG] After Bootstrapping. outs[" << b << "]->k=" << outs[b]->k << ", outs[" << b << "]->N=" << outs[b]->N << std::endl;
            std::cout << "[DEBUG] outs[" << b << "]->parts[0]->coefsT[0..7]: ";
            for(int i=0;i<8 && i<outs[b]->N;++i) std::cout << outs[b]->parts[0]->coefsT[i] << " ";
            std::cout << std::endl;

            // 復号（デバッグ用出力追加）
            std::cout << "[DEBUG] Before Decryption: plain=" << plains[b] << std::endl;
            Torus32 decrypted = mk_lwe_decrypt(outs[b], sks, mp->get_tfhe_params());
            std::cout << "[DEBUG] After Decryption: decrypted=" << decrypted << std::endl;
            std::cout << std::endl;
            std::cout << "[Decryption Check]" << std::endl;
            std::cout << "  Original Plaintext : " << plains[b] << std::endl;
            std::cout << "  Decrypted Value    : " << decrypted << std::endl;
            if (plains[b] == decrypted) {
                std::cout << "  [OK] Decryption matches original plaintext." << std::endl;
            } else {
                std::cout << "  [NG] Decryption does NOT match original plaintext!" << std::endl;
            }
            std::cout << std::endl;
        }
        std::cout << "==================================" << std::endl;
        std::cout << "         [Time Profile]           " << std::endl;
        std::cout << "==================================" << std::endl;
        std::cout << std::fixed << std::setprecision(16);
        std::cout << "  0. Key Generation        : " << global_profiler.time_keygen << " ms" << std::endl;
        std::cout << "  --------------------------------" << std::endl;
        std::cout << "  1. Input Packing         : " << global_profiler.time_input_packing << " ms" << std::endl;
        std::cout << "  2. Blind Rotate (Control): " << global_profiler.time_blind_rotate_control << " ms" << std::endl;
        std::cout << "  3. Homomorphic DFT/Mult  : " << global_profiler.time_external_product << " ms" << std::endl;
        std::cout << "  4. Sample Extract        : " << global_profiler.time_sample_extract << " ms" << std::endl;
        std::cout << "----------------------------------" << std::endl;
        std::cout << "Total Execution Time       : " << total_bs_time << " ms" << std::endl;
        std::cout << "==================================" << std::endl;

        sum_keygen   += global_profiler.time_keygen;
        sum_inputpack+= global_profiler.time_input_packing;
        sum_blind    += global_profiler.time_blind_rotate_control;
        sum_dft      += global_profiler.time_external_product;
        sum_extract  += global_profiler.time_sample_extract;
        sum_total    += total_bs_time;

        // メモリ解放（必要に応じて）
        for(int i=0;i<batch_size;++i) {
            delete ins[i];
            delete outs[i];
        }
        for(int i=0;i<k;++i) delete sks[i];
        delete bk;
        delete mp;
    }
    std::cout << std::endl;
    std::cout << "===== [平均 Time Profile] =====" << std::endl;
    std::cout << std::fixed << std::setprecision(16);
    std::cout << "  0. Key Generation        : " << (sum_keygen/REPEAT) << " ms" << std::endl;
    std::cout << "  1. Input Packing         : " << (sum_inputpack/REPEAT) << " ms" << std::endl;
    std::cout << "  2. Blind Rotate (Control): " << (sum_blind/REPEAT) << " ms" << std::endl;
    std::cout << "  3. Homomorphic DFT/Mult  : " << (sum_dft/REPEAT) << " ms" << std::endl;
    std::cout << "  4. Sample Extract        : " << (sum_extract/REPEAT) << " ms" << std::endl;
    std::cout << "----------------------------------" << std::endl;
    std::cout << "Total Execution Time       : " << (sum_total/REPEAT) << " ms" << std::endl;
    std::cout << "==================================" << std::endl;


    return 0;
}
