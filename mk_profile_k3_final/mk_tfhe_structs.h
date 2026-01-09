#ifndef MK_TFHE_STRUCTS_H
#define MK_TFHE_STRUCTS_H
#include <vector>
#include <iostream>
#include "bb_params.h"
#include "mk_utils.h" // 自作アロケータ読み込み

#include <tfhe.h>
#include <tgsw.h>
#include <tfhe_core.h>

namespace bbii {

// --- 1. BBIIモード ---
// enum class BBIIMode は bb_params.h で定義済み

// --- 2. RLWEサンプル（BBII型） ---
struct MKRLweSample {
    int32_t k; // パーティ数
    int32_t N; // 多項式次数
    std::vector<TorusPolynomial*> parts; // 各パーティのRLWE部分
    MKRLweSample(int32_t parties, const TFheGateBootstrappingParameterSet* params) : k(parties), N(params->in_out_params->n) {
        parts.resize(k+1);
        for(int i=0;i<=k;++i) {
            parts[i] = new_TorusPolynomial(N);
            torusPolynomialClear(parts[i]);
        }
    }
    ~MKRLweSample() {
        for(int i=0;i<=k;++i) {
            if(parts[i]) delete_TorusPolynomial(parts[i]);
        }
    }
};

// --- 3. Packed RGSW（BBII用） ---
struct MKPackedRGSW {
    TGswSampleFFT* sample; // 実体の暗号文
    BBIIMode mode;         // BBIIモード情報
    MKPackedRGSW(const TFheGateBootstrappingParameterSet* params, BBIIMode m) : mode(m) {
        sample = new_TGswSampleFFT(params->tgsw_params);
        tGswFFTClear(sample, params->tgsw_params);
    }
    ~MKPackedRGSW() {
        if (sample) delete_TGswSampleFFT(sample);
    }
};

// --- 4. 秘密鍵構造体（参考） ---
struct MKSecretKey {
    LweKey* lwe_key;
    TGswKey* rlwe_key;
    MKSecretKey(const TFheGateBootstrappingParameterSet* params) {
        lwe_key = new_LweKey(params->in_out_params); lweKeyGen(lwe_key);
        rlwe_key = new_TGswKey(params->tgsw_params); tGswKeyGen(rlwe_key);
    }
    ~MKSecretKey() {
        if(lwe_key) delete_LweKey(lwe_key);
        if(rlwe_key) delete_TGswKey(rlwe_key);
    }
};

// --- 5. Packed RLWE（BBII用アキュムレータ） ---
struct MKPackedRLWE {
    MKRLweSample* sample;
    BBIIMode mode;
    MKPackedRLWE(int32_t parties, const TFheGateBootstrappingParameterSet* params, BBIIMode m) : mode(m) {
        sample = new MKRLweSample(parties, params);
    }
    ~MKPackedRLWE() {
        if (sample) delete sample;
    }
};

// --- 6. BBII用パックド鍵 ---
struct MKBootstrappingKey {
    int32_t k;
    int32_t n_per_party;
    std::vector<std::vector<TGswSampleFFT*>> bk_fft;
    std::vector<std::vector<MKPackedRGSW*>> bk_packed;
    MKBootstrappingKey(int32_t parties, int32_t n, const TFheGateBootstrappingParameterSet* params) : k(parties), n_per_party(n) {
        bk_fft.resize(k);
        bk_packed.resize(k);
        for(int i=0;i<k;++i){
            bk_fft[i].resize(n);
            bk_packed[i].resize(n);
            for(int j=0;j<n;++j){
                bk_fft[i][j]=new_TGswSampleFFT(params->tgsw_params);
                tGswFFTClear(bk_fft[i][j], params->tgsw_params);
                bk_packed[i][j]=new MKPackedRGSW(params, BBIIMode::R12_TO_R13);
            }
        }
    }
    ~MKBootstrappingKey() {
        for(int i=0;i<k;++i){
            for(int j=0;j<n_per_party;++j){
                if(bk_fft[i][j]) delete_TGswSampleFFT(bk_fft[i][j]);
                if(bk_packed[i][j]) delete bk_packed[i][j];
            }
        }
    }
    void generateKeyForParty(int pid, const MKSecretKey* sk, const TFheGateBootstrappingParameterSet* params) {
        for(int j=0; j<this->n_per_party; ++j) {
            TGswSample* t = new_TGswSample(params->tgsw_params);
            tGswSymEncryptInt(t, sk->lwe_key->key[j], params->tgsw_params->tlwe_params->alpha_min, sk->rlwe_key);
            tGswToFFTConvert(this->bk_fft[pid][j], t, params->tgsw_params);
            delete_TGswSample(t);
        }
    }
};
} 
#endif
