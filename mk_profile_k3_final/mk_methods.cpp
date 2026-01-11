#include "mk_methods.h"
#include "mk_profiler.h"
#include <iostream>
#include <chrono>

namespace bbii {
void mk_batch_bootstrapping(
    std::vector<MKRLweSample*>& results,
    const std::vector<MKRLweSample*>& inputs,
    const MKBootstrappingKey* mk_bk,
    const std::vector<Torus32>& mus,
    const TFheGateBootstrappingParameterSet* params
) {
    int batch_size = 8; // paramsから計算する場合は d^rho
    int32_t k = inputs[0]->k;
    std::vector<MKPackedRLWE*> packed_inputs;
    for(size_t i=0; i<batch_size; ++i) {
        MKPackedRLWE* acc = new MKPackedRLWE(k, params, BBIIMode::R12);
        mk_rlwe_clear(acc->sample);
        if (i < inputs.size()) {
            acc->sample->parts[k]->coefsT[0] = mus[i];
        }
        packed_inputs.push_back(acc);
    }
    // DFT/IDFT バッチ処理例（本来はBlind Rotate/DFT/IDFT/Extractを組み合わせる）
    int N = params->in_out_params->n;
    auto dft_matrix = mk_create_dft_matrix(N, false);
    auto idft_matrix = mk_create_dft_matrix(N, true);
    mk_homomorphic_dft_batch(packed_inputs, dft_matrix);
    mk_homomorphic_idft_batch(packed_inputs, idft_matrix);

    // 結果の抽出（仮: Sample Extract未実装）
    for(size_t i=0; i<results.size() && i<packed_inputs.size(); ++i) {
        mk_rlwe_copy(results[i], packed_inputs[i]->sample);
    }

    // メモリ解放
    for(auto* acc : packed_inputs) delete acc;
}
}
