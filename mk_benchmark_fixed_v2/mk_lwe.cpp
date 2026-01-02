#include "mk_methods.h"
#include <vector>
#include <cmath>
#include <random>
#include <limits>
namespace bbii {
Torus32 get_rnd() {
    static std::random_device rd; static std::mt19937 gen(rd());
    static std::uniform_int_distribution<int32_t> dist(std::numeric_limits<int32_t>::min(), std::numeric_limits<int32_t>::max());
    return dist(gen);
}
void mk_lwe_sym_encrypt(MKLweSample* res, Torus32 msg, const MKSecretKey* sk, int32_t pid, const TFheGateBootstrappingParameterSet* p) {
    int32_t n=res->n_per_party; int32_t offset=pid*n; Torus32 prod=0;
    for(int i=0;i<res->k*n;++i) res->sample->a[i]=0;
    for(int i=0;i<n;++i){ Torus32 a=get_rnd(); res->sample->a[offset+i]=a; if(sk->lwe_key->key[i]) prod+=a; }
    res->sample->b = prod + msg + gaussian32(0, p->in_out_params->alpha_min);
    res->sample->current_variance = pow(p->in_out_params->alpha_min,2);
}
Torus32 mk_lwe_decrypt(const MKLweSample* c, const std::vector<MKSecretKey*>& keys, const TFheGateBootstrappingParameterSet* p) {
    Torus32 ph=c->sample->b; int32_t n=c->n_per_party;
    for(int j=0;j<c->k;++j) for(int i=0;i<n;++i) if(keys[j]->lwe_key->key[i]) ph-=c->sample->a[j*n+i];
    return ph;
}
} 
