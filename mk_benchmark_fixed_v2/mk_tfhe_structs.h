#ifndef MK_TFHE_STRUCTS_H
#define MK_TFHE_STRUCTS_H
#include <vector>
#include <iostream>
#include <stdexcept>
#include <cstdlib>
#include "bb_params.h"
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
        if(parts){ for(int i=0;i<=k;++i) delete_TorusPolynomial(parts[i]); delete[] parts; } 
    }
};

struct MKLweSample {
    LweSample* sample; 
    int32_t k; int32_t n_per_party;
    
    MKLweSample(int32_t parties, int32_t n, const TFheGateBootstrappingParameterSet* params) 
        : k(parties), n_per_party(n) {
        
        // 【修正】C++のnewを使い、正しいコンストラクタ引数(LweParams*)を渡す
        sample = new LweSample(params->in_out_params);
        
        // 初期確保された配列を削除し、マルチキーサイズで再確保
        // (LweSampleのコンストラクタが内部で new int32_t[] していると仮定)
        if(sample->a) delete[] sample->a;

        int32_t total_n = k * n;
        sample->a = new int32_t[total_n]; 
        
        // ゼロ初期化
        for(int i=0; i<total_n; ++i) sample->a[i] = 0;
        sample->b = 0;
        sample->current_variance = 0.0;
    }

    ~MKLweSample() { 
        if(sample){ 
            if(sample->a) delete[] sample->a; 
            // sampleポインタ自体を解放 (C++クラスとして扱うためdelete)
            // ※ポインタの中身を先に消したので、sample->a = nullptr にしておくと安全だが
            // ここではTFHEの実装に依存しないよう、sampleのデストラクタが呼ばれても
            // 二重解放にならないように注意が必要。
            // 自分で確保した 'a' は今消したので、nullptrを入れておく
            sample->a = nullptr; 
            delete sample;
        } 
    }
};

struct MKSecretKey {
    LweKey* lwe_key; TGswKey* rlwe_key; 
    MKSecretKey(const TFheGateBootstrappingParameterSet* params) {
        lwe_key = new_LweKey(params->in_out_params); lweKeyGen(lwe_key);
        rlwe_key = new_TGswKey(params->tgsw_params); tGswKeyGen(rlwe_key);
    }
    ~MKSecretKey() { delete_LweKey(lwe_key); delete_TGswKey(rlwe_key); }
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
    ~MKBootstrappingKey() { for(auto& v:bk_fft) for(auto* p:v) delete_TGswSampleFFT(p); }
    
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
