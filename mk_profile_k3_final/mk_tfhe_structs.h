#ifndef MK_TFHE_STRUCTS_H
#define MK_TFHE_STRUCTS_H

#include <vector>
#include <iostream>

#include "bb_params.h"
#include "mk_utils.h" // 自作アロケータ読み込み
#include "mk_packed_ops.h"

#include <tfhe.h>
#include <tgsw.h>
#include <map>
#include <tfhe_core.h>

namespace bbii {
// 前方宣言
struct MKPackedRGSW;
struct MKKeySwitchKey;
struct MKBootstrappingKey;

// --- 追加: Automorphism用 Key Switching Key ---
struct MKKeySwitchKey {
    int32_t k;
    int32_t N;
    std::vector<LweKeySwitchKey*> ks_keys; // 各パーティごとのKSK
    MKKeySwitchKey(int32_t parties, int32_t N_, const TFheGateBootstrappingParameterSet* params) : k(parties), N(N_) {
        ks_keys.resize(k);
        constexpr int basebit = 2; // TFHE標準値（要調整）
        constexpr int t = 8;       // TFHE標準値（要調整）
        for(int i=0;i<k;++i) {
            ks_keys[i] = new_LweKeySwitchKey(params->in_out_params->n, t, basebit, params->in_out_params);
            // 実際には正しい鍵で初期化が必要
        }
    }
    ~MKKeySwitchKey() {
        for(int i=0;i<k;++i) {
            if(ks_keys[i]) delete_LweKeySwitchKey(ks_keys[i]);
        }
    }
};

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
    // perm_key/kskキャッシュ（順列→鍵, delta→鍵）
    std::map<std::vector<int>, MKPackedRGSW*> perm_key_cache;
    std::map<int, MKKeySwitchKey*> ksk_cache;
    int32_t k;
    int32_t n_per_party;
    std::vector<std::vector<TGswSampleFFT*>> bk_fft;
    std::vector<std::vector<MKPackedRGSW*>> bk_packed;
    // 自己同型用KeySwitchingKey（各パーティ・各自己同型写像ごと）
    std::vector<std::vector<MKKeySwitchKey*>> auto_ksk;
    MKBootstrappingKey(int32_t parties, int32_t n, const TFheGateBootstrappingParameterSet* params) : k(parties), n_per_party(n) {
        bk_fft.resize(k);
        bk_packed.resize(k);
        auto_ksk.resize(k); // 各パーティ分
        for(int i=0;i<k;++i){
            bk_fft[i].resize(n);
            bk_packed[i].resize(n);
            // 自己同型用KSKは未初期化（必要に応じて生成）
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
            for(auto* ksk : auto_ksk[i]) {
                if(ksk) delete ksk;
            }
        }
        for(auto& p : perm_key_cache) if(p.second) delete p.second;
        for(auto& p : ksk_cache) if(p.second) delete p.second;
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


// ---- インライン関数群 ----

// perm_keyキャッシュ取得
inline MKPackedRGSW* get_perm_key_cached(MKBootstrappingKey* bk, const std::vector<int>& permutation, const TFheGateBootstrappingParameterSet* params) {
    auto it = bk->perm_key_cache.find(permutation);
    if(it != bk->perm_key_cache.end()) return it->second;
    MKPackedRGSW* key = mk_create_permutation_key(permutation, params);
    bk->perm_key_cache[permutation] = key;
    return key;
}

// kskキャッシュ取得（delta等で区別）
inline MKKeySwitchKey* get_ksk_cached(MKBootstrappingKey* bk, int delta, int32_t k, int32_t N, const TFheGateBootstrappingParameterSet* params) {
    auto it = bk->ksk_cache.find(delta);
    if(it != bk->ksk_cache.end()) return it->second;
    MKKeySwitchKey* ksk = mk_create_automorphism_ksk(k, N, params);
    bk->ksk_cache[delta] = ksk;
    return ksk;
}

// Automorphism用KSK生成（雛形）
inline MKKeySwitchKey* mk_create_automorphism_ksk(int32_t k, int32_t N, const TFheGateBootstrappingParameterSet* params) {
    // 実際はs(X^-1)→s(X)変換用のKSKを生成
    return new MKKeySwitchKey(k, N, params);
}
} // namespace bbii

#endif
