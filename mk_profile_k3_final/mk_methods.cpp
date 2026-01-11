#include "mk_methods.h"
#include "mk_profiler.h"
#include <iostream>
#include <chrono>

void mk_batch_bootstrapping(
    std::vector<MKRLweSample*>& results,
    const std::vector<MKRLweSample*>& inputs,
    const MKBootstrappingKey* mk_bk,
    const std::vector<Torus32>& mus,
    const TFheGateBootstrappingParameterSet* params
) {
    int batch_size = 8; // paramsから計算する場合は d^rho
    int32_t k = inputs[0]->k;
    // --- Input Packing ---
    auto t_pack_start = std::chrono::high_resolution_clock::now();
    std::vector<MKPackedRLWE*> packed_inputs;
    for(size_t i=0; i<batch_size; ++i) {
        MKPackedRLWE* acc = new MKPackedRLWE(k, params, BBIIMode::R12);
        mk_rlwe_clear(acc->sample);
        if (i < inputs.size()) {
            acc->sample->parts[k]->coefsT[0] = mus[i];
        }
        packed_inputs.push_back(acc);
    }
    auto t_pack_end = std::chrono::high_resolution_clock::now();
    global_profiler.time_input_packing = std::chrono::duration<double, std::milli>(t_pack_end - t_pack_start).count();

    // --- DFT ---
    int N = params->in_out_params->n;
    auto dft_matrix = mk_create_dft_matrix(N, false);
    auto t_dft_start = std::chrono::high_resolution_clock::now();
    mk_homomorphic_dft_batch(packed_inputs, dft_matrix);
    auto t_dft_end = std::chrono::high_resolution_clock::now();
    global_profiler.time_external_product = std::chrono::duration<double, std::milli>(t_dft_end - t_dft_start).count();

    // --- IDFT ---
    auto idft_matrix = mk_create_dft_matrix(N, true);
    auto t_idft_start = std::chrono::high_resolution_clock::now();
    mk_homomorphic_idft_batch(packed_inputs, idft_matrix);
    auto t_idft_end = std::chrono::high_resolution_clock::now();
    // DFT/IDFT合算で記録したい場合は加算、分けたい場合は別途項目を用意
    global_profiler.time_external_product += std::chrono::duration<double, std::milli>(t_idft_end - t_idft_start).count();

    // --- Sample Extract（仮: RLWE→LWEコピーのみ）---
    auto t_extract_start = std::chrono::high_resolution_clock::now();
    for(size_t i=0; i<results.size() && i<packed_inputs.size(); ++i) {
        mk_rlwe_copy(results[i], packed_inputs[i]->sample);
    }
    auto t_extract_end = std::chrono::high_resolution_clock::now();
    global_profiler.time_sample_extract = std::chrono::duration<double, std::milli>(t_extract_end - t_extract_start).count();

    // メモリ解放
    for(auto* acc : packed_inputs) delete acc;
}
