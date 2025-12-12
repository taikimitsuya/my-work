#ifndef BATCH_PARAMS_H
#define BATCH_PARAMS_H

#include "TFHEParams.h" // MK-TFHEの既存パラメータを想定

struct BatchParams {
    // 論文のパラメータ設定 [cite: 662-667]
    int32_t n;      // Input LWE dimension
    int32_t N;      // Ring dimension (RGSW)
    int32_t q;      // Modulus
    int32_t d;      // Parameter for DFT/Inverse-DFT (2d | m)
    int32_t rho;    // Recursive depth (n = 2d^rho)
    int32_t r;      // Max slots for batching (r = O(sqrt(N/q)))
    int32_t v;      // Number of inputs for vector-matrix mult

    // Tensor Ring modes [cite: 278]
    enum Mode { MODE_R12, MODE_R13, MODE_R12_TO_R13, MODE_R13_TO_R12 };
};

#endif
