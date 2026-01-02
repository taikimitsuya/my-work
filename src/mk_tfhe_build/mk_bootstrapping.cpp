#include "mk_tfhe_structs.h"
#include "mk_params.h"
#include <cmath>
#include <tfhe.h>
#include <tfhe_core.h>
#include <polynomials.h>
namespace bbii {
void mk_cmux(MKRLweSample*, const TGswSampleFFT*, const MKRLweSample*, const MKRLweSample*, int32_t, const TFheGateBootstrappingParameterSet*);
void mk_rlwe_clear(MKRLweSample*); void mk_rlwe_copy(MKRLweSample*, const MKRLweSample*);
void mk_mul_xai(MKRLweSample* r, const MKRLweSample* s, int32_t b, int32_t N){ for(int i=0;i<=s->k;++i) torusPolynomialMulByXai(r->parts[i], b, s->parts[i]); }
void mk_bootstrapping(MKLweSample* res, const MKLweSample* in, const MKBootstrappingKey* bk, Torus32 mu, const TFheGateBootstrappingParameterSet* p) {
    int32_t k=in->k, n=in->n_per_party, N=p->tgsw_params->tlwe_params->N, _2N=2*N;
    MKRLweSample* acc=new MKRLweSample(k,p); mk_rlwe_clear(acc); acc->parts[k]->coefsT[0]=mu;
    MKRLweSample* t=new MKRLweSample(k,p); mk_rlwe_copy(t,acc);
    mk_mul_xai(acc,t,-modSwitchFromTorus32(in->sample->b,_2N),N);
    for(int u=0;u<k;++u) for(int i=0;i<n;++i){
        int32_t bar=modSwitchFromTorus32(in->sample->a[u*n+i],_2N);
        if(bar){ mk_mul_xai(t,acc,bar,N); mk_cmux(acc,bk->bk_fft[u][i],acc,t,u,p); }
    }
    res->sample->b=acc->parts[k]->coefsT[0];
    for(int u=0;u<k;++u){
        for(int j=0;j<N;++j) res->sample->a[u*N+j] = (j==0)? acc->parts[u]->coefsT[0] : -acc->parts[u]->coefsT[N-j];
    }
    delete acc; delete t;
}
} 
