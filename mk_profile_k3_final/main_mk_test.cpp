#include <iostream>
#include <cstdlib>
#include <vector>
#include <chrono>
#include <iomanip>
#include "mk_methods.h"
#include "mk_profiler.h"

using namespace std;

MKProfiler global_profiler;

int main() {
        // DFTベースBlind Rotateのテスト呼び出し
        extern void mk_test_blind_rotate_dft(const TFheGateBootstrappingParameterSet* params);
        // テスト用パラメータセット（N=1024固定）
        TFheGateBootstrappingParameterSet* test_params = new_default_gate_bootstrapping_parameters(110);
        mk_test_blind_rotate_dft(test_params);
        delete_gate_bootstrapping_parameters(test_params);
    int32_t k, rho;
    int32_t d = 2; // fixed
    int32_t N = 1024;
    int32_t REPEAT;
    cout << "=== MK-TFHE Time Profiling ===" << endl;
    cout << "パーティ数 k を入力してください: ";
    cin >> k;
    cout << "パラメータ rho を入力してください: ";
    cin >> rho;
    cout << "試行回数を入力してください: ";
    cin >> REPEAT;
    cout << "(d=2, N=1024 固定)" << endl;

    // 平均用バッファ
    double sum_keygen = 0, sum_inputpack = 0, sum_blind = 0, sum_dft = 0, sum_extract = 0, sum_total = 0;
        for(int rep=1; rep<=REPEAT; ++rep) {
        // 各ループの最初でプロファイラ値をリセット
        global_profiler.time_keygen = 0;
        global_profiler.time_input_packing = 0;
        global_profiler.time_blind_rotate_control = 0;
        global_profiler.time_external_product = 0;
        global_profiler.time_sample_extract = 0;
        cout << "\n[Run " << rep << "/" << REPEAT << "]" << endl;
        cout << "Initializing Params..." << endl;
        bbii::MKParams* mp = bbii::get_mk_test_params(k, d, rho, N);
        bbii::BBIIParams* bp = mp->sk_params;

        cout << "==================================" << endl;
        cout << "         [Parameters]            " << endl;
        cout << "==================================" << endl;
        cout << "  Parties (k)              : " << k << endl;
        cout << "  Recursion Depth (d)      : " << d << endl;
        cout << "  Partition Size (rho)     : " << rho << endl;
        cout << "  TLWE Dimension (N)       : " << N << endl;
        cout << "  n per party              : " << mp->n_per_party << endl;
        cout << "==================================" << endl;
        cout << endl;

        cout << "Generating Keys..." << endl;
        vector<bbii::MKSecretKey*> sks(k); 
        bbii::MKBootstrappingKey* bk = new bbii::MKBootstrappingKey(k, mp->n_per_party, mp->get_tfhe_params());
        // KeyGen: 差分計測
        auto kg_start = std::chrono::high_resolution_clock::now();
        double keygen0 = global_profiler.time_keygen;
        for(int i=0; i<k; ++i){ 
            sks[i] = new bbii::MKSecretKey(mp->get_tfhe_params()); 
            bk->generateKeyForParty(i, sks[i], mp->get_tfhe_params()); 
        }
        auto kg_end = std::chrono::high_resolution_clock::now();
        double keygen1 = global_profiler.time_keygen;
        global_profiler.time_keygen = (keygen1 - keygen0) + std::chrono::duration<double, std::milli>(kg_end - kg_start).count();

        cout << "Encrypting..." << endl;
        bbii::MKRLweSample* in = new bbii::MKRLweSample(k, mp->get_tfhe_params());
        bbii::mk_lwe_sym_encrypt(in, 1000, sks[0], 0, mp->get_tfhe_params(), mp->n_per_party);

        bbii::MKRLweSample* out = new bbii::MKRLweSample(k, mp->get_tfhe_params());
        cout << "Running Bootstrapping (This will take time)..." << endl;
        auto bs_start = std::chrono::high_resolution_clock::now();
        bbii::mk_bootstrapping(out, in, bk, 1000, mp->get_tfhe_params());
        auto bs_end = std::chrono::high_resolution_clock::now();
        double total_bs_time = std::chrono::duration<double, std::milli>(bs_end - bs_start).count();

        cout << endl;
        cout << "==================================" << endl;
        cout << "         [Time Profile]           " << endl;
        cout << "==================================" << endl;
        cout << fixed << setprecision(16);
        cout << "  0. Key Generation        : " << global_profiler.time_keygen << " ms" << endl;
        cout << "  --------------------------------" << endl;
        cout << "  1. Input Packing         : " << global_profiler.time_input_packing << " ms" << endl;
        cout << "  2. Blind Rotate (Control): " << global_profiler.time_blind_rotate_control << " ms" << endl;
        cout << "  3. Homomorphic DFT/Mult  : " << global_profiler.time_external_product << " ms" << endl;
        cout << "  4. Sample Extract        : " << global_profiler.time_sample_extract << " ms" << endl;
        cout << "----------------------------------" << endl;
        cout << "Total Execution Time       : " << total_bs_time << " ms" << endl;
        cout << "==================================" << endl;

        sum_keygen   += global_profiler.time_keygen;
        sum_inputpack+= global_profiler.time_input_packing;
        sum_blind    += global_profiler.time_blind_rotate_control;
        sum_dft      += global_profiler.time_external_product;
        sum_extract  += global_profiler.time_sample_extract;
        sum_total    += total_bs_time;

        // メモリ解放（必要に応じて）
        delete in;
        delete out;
        for(int i=0;i<k;++i) delete sks[i];
        delete bk;
        delete mp;
    }
    cout << endl;
    cout << "===== [平均 Time Profile] =====" << endl;
    cout << fixed << setprecision(16);
    cout << "  0. Key Generation        : " << (sum_keygen/REPEAT) << " ms" << endl;
    cout << "  1. Input Packing         : " << (sum_inputpack/REPEAT) << " ms" << endl;
    cout << "  2. Blind Rotate (Control): " << (sum_blind/REPEAT) << " ms" << endl;
    cout << "  3. Homomorphic DFT/Mult  : " << (sum_dft/REPEAT) << " ms" << endl;
    cout << "  4. Sample Extract        : " << (sum_extract/REPEAT) << " ms" << endl;
    cout << "----------------------------------" << endl;
    cout << "Total Execution Time       : " << (sum_total/REPEAT) << " ms" << endl;
    cout << "==================================" << endl;
    return 0;
}
