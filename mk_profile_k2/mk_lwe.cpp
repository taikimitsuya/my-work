#include "mk_methods.h"
#include "mk_profiler.h"
#include <iostream> 
namespace bbii {
Torus32 get_rnd() { return 12345; }
void mk_lwe_sym_encrypt(MKLweSample* res, Torus32 msg, const MKSecretKey* sk, int32_t pid, const TFheGateBootstrappingParameterSet* p) {
    // 暗号化の計測
    auto start = std::chrono::high_resolution_clock::now();
    
    int32_t n=res->n_per_party; int32_t offset=pid*n; Torus32 prod=0;
    for(int i=0;i<res->k*n;++i) res->sample->a[i]=0;
    for(int i=0;i<n;++i){ Torus32 a=get_rnd(); res->sample->a[offset+i]=a; if(sk->lwe_key->key[i]) prod+=a; }
    res->sample->b = prod + msg + gaussian32(0, p->in_out_params->alpha_min);
    res->sample->current_variance = 0.0; 

    auto end = std::chrono::high_resolution_clock::now();
    global_profiler.time_encrypt += std::chrono::duration<double, std::milli>(end - start).count();
}
Torus32 mk_lwe_decrypt(const MKLweSample* c, const std::vector<MKSecretKey*>& keys, const TFheGateBootstrappingParameterSet* p) { return c->sample->b; }
} 
