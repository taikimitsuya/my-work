#ifndef MK_TFHE_STRUCTS_H
#define MK_TFHE_STRUCTS_H

#include <vector>
#include <iostream>
#include <stdexcept>
#include <cstdlib> // ★追加: malloc/freeのため

// 既存のBBIIパラメータとTFHEライブラリをインクルード
#include "bb_params.h"
#include <tfhe.h>
#include <tfhe_core.h>

namespace bbii {

/**
 * @brief マルチキーRLWEサンプル (MK-GLWE)
 */
struct MKRLweSample {
    int32_t k;              
    int32_t N;              
    TorusPolynomial** parts;

    MKRLweSample(int32_t parties, const TFheGateBootstrappingParameterSet* params) {
        if (!params || !params->tgsw_params || !params->tgsw_params->tlwe_params) {
            throw std::runtime_error("Invalid parameters passed to MKRLweSample");
        }

        this->k = parties;
        this->N = params->tgsw_params->tlwe_params->N;
        
        this->parts = new TorusPolynomial*[k + 1];
        
        for (int i = 0; i <= k; ++i) {
            this->parts[i] = new_TorusPolynomial(N);
            torusPolynomialClear(this->parts[i]);
        }
    }

    ~MKRLweSample() {
        if (parts) {
            for (int i = 0; i <= k; ++i) {
                delete_TorusPolynomial(parts[i]);
            }
            delete[] parts;
        }
    }

    MKRLweSample(const MKRLweSample&) = delete;
    MKRLweSample& operator=(const MKRLweSample&) = delete;

    TorusPolynomial* getPart(int index) const {
        if (index < 0 || index > k) {
            throw std::out_of_range("MKRLweSample index out of bounds");
        }
        return parts[index];
    }
};

/**
 * @brief マルチキーLWEサンプル
 */
struct MKLweSample {
    LweSample* sample; 
    int32_t k;
    int32_t n_per_party;

    MKLweSample(int32_t parties, int32_t n, const TFheGateBootstrappingParameterSet* params) 
        : k(parties), n_per_party(n) {
        
        int32_t total_n = k * n;
        
        // ★修正箇所: new LweSample ではなく malloc を使用
        // LweSample構造体自体のメモリを確保
        sample = (LweSample*)std::malloc(sizeof(LweSample));

        if (!sample) {
             throw std::runtime_error("Failed to allocate memory for MKLweSample");
        }

        // 内部の配列は new[] で確保 (デストラクタで delete[] するため整合性は取れる)
        sample->a = new int32_t[total_n]; // Torus32 は int32_t のエイリアス
        sample->b = 0;
        sample->current_variance = 0.0;
        
        // ゼロ初期化
        for(int i=0; i<total_n; ++i) sample->a[i] = 0;
    }

    ~MKLweSample() {
        if (sample) {
            // 配列は new[] で確保したので delete[] で解放
            if (sample->a) delete[] sample->a;
            
            // ★修正箇所: 構造体自体は malloc で確保したので free で解放
            std::free(sample);
        }
    }
};

/**
 * @brief 個々のパーティの秘密鍵
 */
struct MKSecretKey {
    LweKey* lwe_key;   
    TGswKey* rlwe_key; 

    MKSecretKey(const TFheGateBootstrappingParameterSet* params) {
        const LweParams* lwe_p = params->in_out_params;
        const TGswParams* tgsw_p = params->tgsw_params;

        lwe_key = new_LweKey(lwe_p);
        lweKeyGen(lwe_key);

        rlwe_key = new_TGswKey(tgsw_p);
        tGswKeyGen(rlwe_key);
    }

    ~MKSecretKey() {
        delete_LweKey(lwe_key);
        delete_TGswKey(rlwe_key);
    }
};

/**
 * @brief マルチキー・ブートストラップ鍵 (MK-BK)
 */
struct MKBootstrappingKey {
    int32_t k;            
    int32_t n_per_party;  
    
    std::vector<std::vector<TGswSampleFFT*>> bk_fft;

    MKBootstrappingKey(int32_t parties, int32_t n, const TFheGateBootstrappingParameterSet* params) 
        : k(parties), n_per_party(n) {
        
        bk_fft.resize(k);
        for (int i = 0; i < k; ++i) {
            bk_fft[i].resize(n_per_party);
            for (int j = 0; j < n_per_party; ++j) {
                bk_fft[i][j] = new_TGswSampleFFT(params->tgsw_params);
                tGswFFTClear(bk_fft[i][j], params->tgsw_params);
            }
        }
    }

    ~MKBootstrappingKey() {
        for (int i = 0; i < k; ++i) {
            for (int j = 0; j < n_per_party; ++j) {
                delete_TGswSampleFFT(bk_fft[i][j]);
            }
        }
    }

    void generateKeyForParty(int party_id, const MKSecretKey* sk, const TFheGateBootstrappingParameterSet* params) {
        if (party_id < 0 || party_id >= k) return;

        const LweKey* lwe_k = sk->lwe_key;
        const TGswKey* rlwe_k = sk->rlwe_key;
        const TGswParams* tgsw_p = params->tgsw_params;

        for (int j = 0; j < n_per_party; ++j) {
            int32_t s_j = lwe_k->key[j];

            TGswSample* temp_bk = new_TGswSample(tgsw_p);
            tGswSymEncryptInt(temp_bk, s_j, params->tgsw_params->tlwe_params->alpha_min, rlwe_k);

            tGswToFFTConvert(bk_fft[party_id][j], temp_bk, tgsw_p);

            delete_TGswSample(temp_bk);
        }
    }
};

} // namespace bbii

#endif // MK_TFHE_STRUCTS_H