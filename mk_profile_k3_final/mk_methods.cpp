#include "mk_methods.h"
#include "mk_profiler.h"
#include <iostream>
#include <chrono>

namespace bbii {
void mk_bootstrapping(MKRLweSample* result, const MKRLweSample* input, const MKBootstrappingKey* mk_bk, Torus32 mu, const TFheGateBootstrappingParameterSet* params) {
    auto pack_start = std::chrono::high_resolution_clock::now();
    double pack0 = global_profiler.time_input_packing;
    int32_t k = input->k;
    MKRLweSample* acc = new MKRLweSample(k, params);
    mk_rlwe_clear(acc);
    acc->parts[k]->coefsT[0] = mu;
    auto pack_end = std::chrono::high_resolution_clock::now();
    double pack1 = global_profiler.time_input_packing;
    global_profiler.time_input_packing = (pack1 - pack0) + std::chrono::duration<double, std::milli>(pack_end - pack_start).count();

    // DFTベースBlind Rotate
    mk_blind_rotate_dft(acc, input, mk_bk, params);

    // 結果をresultにコピー
    mk_rlwe_copy(result, acc);
    delete acc;
}
}
