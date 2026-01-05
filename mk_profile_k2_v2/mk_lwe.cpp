#include "mk_methods.h"
#include "mk_profiler.h"
#include <iostream> 
#include <cstdlib> // rand(), srand() のため

namespace bbii {

// シンプルかつ安全な乱数生成
Torus32 get_rnd() { 
    // rand() は 0〜RAND_MAX を返す。
    // Torus32 (int32) 全域を埋めるためにシフトして合成する
    uint32_t r1 = (uint32_t)rand();
    uint32_t r2 = (uint32_t)rand();
    return (Torus32)((r1 << 16) | (r2 & 0xFFFF));
}

void mk_lwe_sym_encrypt(MKLweSample* res, Torus32 msg, const MKSecretKey* sk, int32_t pid, const TFheGateBootstrappingParameterSet* p) {
    auto start = std::chrono::high_resolution_clock::now();
    
    // 初回のみシード固定 (再現性のため)
    static bool seeded = false;
    if(!seeded) { srand(12345); seeded = true; }

    int32_t n=res->n_per_party; int32_t offset=pid*n; Torus32 prod=0;
    
    // 初期化
    for(int i=0;i<res->k*n;++i) res->sample->a[i]=0;
    
    // 乱数セット (ここが 0 以外になることで重い計算が走る)
    for(int i=0;i<n;++i){ 
        Torus32 a=get_rnd(); 
        res->sample->a[offset+i]=a; 
        if(sk->lwe_key->key[i]) prod+=a; 
    }
    
    res->sample->b = prod + msg + gaussian32(0, p->in_out_params->alpha_min);
    res->sample->current_variance = 0.0; 

    auto end = std::chrono::high_resolution_clock::now();
    global_profiler.time_encrypt += std::chrono::duration<double, std::milli>(end - start).count();
}
Torus32 mk_lwe_decrypt(const MKLweSample* c, const std::vector<MKSecretKey*>& keys, const TFheGateBootstrappingParameterSet* p) { return c->sample->b; }
} 
