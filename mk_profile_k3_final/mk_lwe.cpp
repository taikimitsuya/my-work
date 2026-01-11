#include "mk_methods.h"
#include "mk_profiler.h"
#include <iostream>
#include <cstdlib>
#include <vector>
#include <chrono>
#include <tfhe/tfhe.h>


Torus32 get_rnd() { 
    uint32_t r1 = (uint32_t)rand();
    uint32_t r2 = (uint32_t)rand();
    return (Torus32)((r1 << 16) | (r2 & 0xFFFF));
}

void mk_lwe_sym_encrypt(MKRLweSample* res, Torus32 msg, const MKSecretKey* sk, int32_t pid, const TFheGateBootstrappingParameterSet* p, int32_t n_per_party) {
    auto start = std::chrono::high_resolution_clock::now();
    static bool seeded = false;
    if(!seeded) { srand(12345); seeded = true; }
    // 各パーティの0次係数にLWE的な値をセット
    for(int u=0; u<=res->k; ++u) {
        for(int i=0; i<n_per_party; ++i) {
            res->parts[u]->coefsT[i] = 0;
        }
    }
    // pid番目のパーティの0次係数に平文+ノイズをセット
    Torus32 prod = 0;
    for(int i=0; i<n_per_party; ++i) {
        Torus32 a = get_rnd();
        res->parts[pid]->coefsT[i] = a;
        if(sk->lwe_key->key[i]) prod += a;
    }
    res->parts[pid]->coefsT[0] = prod + msg + gaussian32(0, p->in_out_params->alpha_min);
    auto end = std::chrono::high_resolution_clock::now();
    global_profiler.time_encrypt += std::chrono::duration<double, std::milli>(end - start).count();
}
Torus32 mk_lwe_decrypt(const MKRLweSample* c, const std::vector<MKSecretKey*>& keys, const TFheGateBootstrappingParameterSet* p) {
    int32_t k = c->k;
    int32_t N = p->in_out_params->n; // n_per_party
    Torus32 b = c->parts[k]->coefsT[0];
    Torus32 sum_as = 0;
    for(int u=0; u<k; ++u) {
        const LweKey* key = keys[u]->lwe_key;
        const TorusPolynomial* poly = c->parts[u];
        for(int i=0; i<N; ++i) {
            sum_as += poly->coefsT[i] * key->key[i];
        }
    }
    Torus32 phase = b - sum_as;
    return phase;
}
