#ifndef MK_TFHE_STRUCTS_H
#define MK_TFHE_STRUCTS_H
#include <vector>
#include <iostream>
#include "bb_params.h"
#include <tfhe/tfhe.h>
#include <tfhe/tfhe_core.h>
namespace bbii {
struct MKRLweSample {
    int32_t k; int32_t N; TorusPolynomial** parts;
    MKRLweSample(int32_t parties, const TFheGateBootstrappingParameterSet* params) {
        this->k = parties; this->N = params->tgsw_params->tlwe_params->N;
        this->parts = new TorusPolynomial*[k + 1];
        for (int i = 0; i <= k; ++i) { this->parts[i] = new_TorusPolynomial(N); torusPolynomialClear(this->parts[i]); }
    }
    ~MKRLweSample() { /* Leak safe */ }
};
struct MKLweSample {
    LweSample* sample; int32_t k; int32_t n_per_party; int32_t* my_array; 
    MKLweSample(int32_t parties, int32_t n, const TFheGateBootstrappingParameterSet* params) 
        : k(parties), n_per_party(n), my_array(nullptr) {
        sample = new_LweSample(params->in_out_params);
        int32_t total_n = k * n;
        my_array = new int32_t[total_n]; 
        sample->a = my_array;
        for(int i=0; i<total_n; ++i) sample->a[i] = 0;
        sample->b = 0; sample->current_variance = 0.0;
    }
    ~MKLweSample() { if(sample) { sample->a = nullptr; delete sample; } }
};
struct MKSecretKey {
    LweKey* lwe_key; TGswKey* rlwe_key; 
    MKSecretKey(const TFheGateBootstrappingParameterSet* params) {
        lwe_key = new_LweKey(params->in_out_params); lweKeyGen(lwe_key);
        rlwe_key = new_TGswKey(params->tgsw_params); tGswKeyGen(rlwe_key);
    }
    ~MKSecretKey() { /* Leak safe */ }
};
struct MKBootstrappingKey {
    int32_t k; int32_t n_per_party; std::vector<std::vector<TGswSampleFFT*>> bk_fft;
    MKBootstrappingKey(int32_t parties, int32_t n, const TFheGateBootstrappingParameterSet* params) : k(parties), n_per_party(n) {
        bk_fft.resize(k);
        for(int i=0;i<k;++i){ 
            bk_fft[i].resize(n); 
            for(int j=0;j<n;++j){ bk_fft[i][j]=new_TGswSampleFFT(params->tgsw_params); tGswFFTClear(bk_fft[i][j], params->tgsw_params); } 
        }
    }
    ~MKBootstrappingKey() { /* Leak safe */ }
    void generateKeyForParty(int pid, const MKSecretKey* sk, const TFheGateBootstrappingParameterSet* params) {
        for(int j=0; j<n_per_party; ++j) {
            TGswSample* t = new_TGswSample(params->tgsw_params);
            tGswSymEncryptInt(t, sk->lwe_key->key[j], params->tgsw_params->tlwe_params->alpha_min, sk->rlwe_key);
            tGswToFFTConvert(bk_fft[pid][j], t, params->tgsw_params);
            delete_TGswSample(t);
        }
    }
};
} 
#endif
