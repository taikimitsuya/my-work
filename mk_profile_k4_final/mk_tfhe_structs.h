#ifndef MK_TFHE_STRUCTS_H
#define MK_TFHE_STRUCTS_H
#include <vector>
#include <iostream>
#include "bb_params.h"
#include "mk_utils.h" // 自作アロケータ読み込み
#include <tfhe.h>
#include <tfhe_core.h>

namespace bbii {
struct MKRLweSample {
    int32_t k; int32_t N; TorusPolynomial** parts;
    MKRLweSample(int32_t parties, const TFheGateBootstrappingParameterSet* params) {
        this->k = parties; this->N = params->tgsw_params->tlwe_params->N;
        this->parts = new TorusPolynomial*[k + 1];
        for (int i = 0; i <= k; ++i) { 
            this->parts[i] = new_TorusPolynomial(N); 
            torusPolynomialClear(this->parts[i]); 
        }
    }
    ~MKRLweSample() { 
        if (this->parts) {
            for (int i = 0; i <= k; ++i) { 
                if (this->parts[i]) {
                    delete_TorusPolynomial(this->parts[i]);
                    this->parts[i] = nullptr;
                }
            }
            delete[] this->parts;
            this->parts = nullptr;
        }
    }
};

struct MKLweSample {
    LweSample* sample; int32_t k; int32_t n_per_party; LweParams* alloc_params;
    
    MKLweSample(int32_t parties, int32_t n, const TFheGateBootstrappingParameterSet* params) 
        : k(parties), n_per_party(n), alloc_params(nullptr) {
        int32_t total_n = k * n;
        // TFHE 側で期待される "n" を拡張した LweParams を作成し、それで LweSample を生成する
        alloc_params = new_LweParams(total_n, params->in_out_params->alpha_min, 0.5);
        sample = new_LweSample(alloc_params);
        for (int i = 0; i < total_n; ++i) sample->a[i] = 0;
        sample->b = 0; sample->current_variance = 0.0;
    }
    ~MKLweSample() { 
        if (sample) {
            delete_LweSample(sample);
            sample = nullptr;
        }
        if (alloc_params) {
            delete_LweParams(alloc_params);
            alloc_params = nullptr;
        }
    }
};

struct MKSecretKey {
    LweKey* lwe_key; TGswKey* rlwe_key; 
    MKSecretKey(const TFheGateBootstrappingParameterSet* params) {
        lwe_key = new_LweKey(params->in_out_params); lweKeyGen(lwe_key);
        rlwe_key = new_TGswKey(params->tgsw_params); tGswKeyGen(rlwe_key);
    }
};

struct MKBootstrappingKey {
    int32_t k; int32_t n_per_party; std::vector<std::vector<TGswSampleFFT*>> bk_fft;
    MKBootstrappingKey(int32_t parties, int32_t n, const TFheGateBootstrappingParameterSet* params) : k(parties), n_per_party(n) {
        bk_fft.resize(k);
        for(int i=0;i<k;++i){ 
            bk_fft[i].resize(n); 
            for(int j=0;j<n;++j){ 
                bk_fft[i][j]=new_TGswSampleFFT(params->tgsw_params); 
                tGswFFTClear(bk_fft[i][j], params->tgsw_params); 
            } 
        }
    }
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
