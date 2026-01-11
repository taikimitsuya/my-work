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
    // ここで一括Blind Rotate (DFT) を呼び出す（仮: 未実装）
    // mk_blind_rotate_dft_batch(packed_inputs, ...);

    // 結果の抽出（仮: 未実装）
    // for (auto* acc : packed_inputs) { ... }

    // メモリ解放
    for(auto* acc : packed_inputs) delete acc;
}
}
