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
 * * シングルキーのRLWEは (a, b) のペアですが、
 * k人のパーティがいる場合、(a_1, a_2, ..., a_k, b) の k+1 個の多項式になります。
 * * 主にBlind Rotateのアキュムレータとして使用されます。
 */
struct MKRLweSample {
    int32_t k;              // パーティ数 (参加者数)
    int32_t N;              // 多項式の次数
    TorusPolynomial** parts; // 多項式の配列 (サイズ k+1)

    // コンストラクタ
    MKRLweSample(int32_t parties, const TFheGateBootstrappingParameterSet* params) {
        if (!params || !params->tgsw_params || !params->tgsw_params->tlwe_params) {
            throw std::runtime_error("Invalid parameters passed to MKRLweSample");
        }

        this->k = parties;
        this->N = params->tgsw_params->tlwe_params->N;
        
        // k+1 個の多項式ポインタを確保 (a_1...a_k, b)
        this->parts = new TorusPolynomial*[k + 1];
        
        for (int i = 0; i <= k; ++i) {
            // libtfheの関数を使ってアライメントされたメモリを確保
            this->parts[i] = new_TorusPolynomial(N);
            
            // 初期化（ゼロ埋め）
            torusPolynomialClear(this->parts[i]);
        }
    }

    // デストラクタ
    ~MKRLweSample() {
        if (parts) {
            for (int i = 0; i <= k; ++i) {
                delete_TorusPolynomial(parts[i]);
            }
            delete[] parts;
        }
    }

    // コピー禁止 (ポインタ管理が複雑になるため、必要ならDeep Copy関数を別途実装)
    MKRLweSample(const MKRLweSample&) = delete;
    MKRLweSample& operator=(const MKRLweSample&) = delete;

    // 特定の多項式へのアクセス (0 <= index <= k)
    // index < k  : マスク部 (a_{index+1})
    // index == k : ボディ部 (b)
    TorusPolynomial* getPart(int index) const {
        if (index < 0 || index > k) {
            throw std::out_of_range("MKRLweSample index out of bounds");
        }
        return parts[index];
    }
};

/**
 * @brief マルチキーLWEサンプル
 * * ブートストラップ後の出力や、回路の入出力として使用します。
 * 構造: (a_1, ..., a_k, b) ここで各要素は torus32 の係数ではなくベクトルですが、
 * 簡略化のため、長さ n*k の結合されたマスク a と、ボディ b を持つ標準LWEとして扱えるようにラップします。
 * または、拡張性を考慮して論理的に分離します。ここでは結合型として定義します。
 */
struct MKLweSample {
    LweSample* sample; // サイズ N_total = n * k のLWEサンプル
    int32_t k;
    int32_t n_per_party;

    MKLweSample(int32_t parties, int32_t n, const TFheGateBootstrappingParameterSet* params) 
        : k(parties), n_per_party(n) {
        
        // 全体の次元 = パーティ数 * 1人あたりのLWE次元
        int32_t total_n = k * n;
        
        // 標準LWE関数を使って確保（パラメータはダミーで次元だけ合わせるか、カスタム生成）
        // ここでは raw 生成を用います

        // ★修正箇所: new LweSample ではなく malloc を使用
        // LweSample構造体自体のメモリを確保
        sample = (LweSample*)std::malloc(sizeof(LweSample));

        //sample = new_LweSample(&params->in_out_params->lwe_params[0]); 
        // ※注意: 正確には total_n サイズの LweParams が必要ですが、
        // 簡易実装として、既存の構造体を流用しつつ、サイズ管理は外部で行う想定です。
        // もし厳密に行うなら、以下のように手動確保します：
        
        if (sample) delete_LweSample(sample); // params由来のものを破棄
        
        // total_n サイズで手動確保
        //sample = new LweSample;
        //sample->a = new int32_t[total_n];
        //sample->b = 0;
        //sample->current_variance = 0.0;
        
        // 内部の配列は new[] で確保 (デストラクタで delete[] するため整合性は取れる)
        sample->a = new int32_t[total_n]; // Torus32 は int32_t のエイリアス
        sample->b = 0;
        sample->current_variance = 0.0;

        // ゼロ初期化
        for(int i=0; i<total_n; ++i) sample->a[i] = 0;
    }

    /*~MKLweSample() {
        if (sample) {
            delete[] sample->a;
            delete sample;
        }
    }*/

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
 * * シミュレーション実行用です。
 * LWEの秘密鍵と、RLWEの秘密鍵を持ちます。
 */
struct MKSecretKey {
    LweKey* lwe_key;   // 1人分のLWE鍵 (サイズ n)
    TGswKey* rlwe_key; // 1人分のRLWE鍵 (BK生成用)

    MKSecretKey(const TFheGateBootstrappingParameterSet* params) {
        // パラメータ抽出
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
 * * 全パーティのブートストラップ鍵の集合体です。
 * 構造: vector< vector<TGswSampleFFT*> >
 * user_keys[i][j] : パーティ i の LWE秘密鍵の j 番目のビットを暗号化したもの
 * * ※計算高速化のため、FFT領域 (TGswSampleFFT) で保持します。
 */
struct MKBootstrappingKey {
    int32_t k;            // パーティ数
    int32_t n_per_party;  // 1人あたりのLWE次元 (rho次第, 例: n=64, 512...)
    
    // [パーティID][LWEビットID] -> TRGSWサンプル(FFT)
    std::vector<std::vector<TGswSampleFFT*>> bk_fft;

    // コンストラクタ: 枠だけ確保
    MKBootstrappingKey(int32_t parties, int32_t n, const TFheGateBootstrappingParameterSet* params) 
        : k(parties), n_per_party(n) {
        
        bk_fft.resize(k);
        for (int i = 0; i < k; ++i) {
            bk_fft[i].resize(n_per_party);
            for (int j = 0; j < n_per_party; ++j) {
                // FFT形式のTRGSWサンプルを確保
                bk_fft[i][j] = new_TGswSampleFFT(params->tgsw_params);
                // 初期化 (0)
                tGswFFTClear(bk_fft[i][j], params->tgsw_params);
            }
        }
    }

    // デストラクタ
    ~MKBootstrappingKey() {
        for (int i = 0; i < k; ++i) {
            for (int j = 0; j < n_per_party; ++j) {
                delete_TGswSampleFFT(bk_fft[i][j]);
            }
        }
    }

    // 鍵生成メソッド (シミュレーション用)
    // 指定したパーティIDの秘密鍵 sk を使って、BKを生成・格納する
    void generateKeyForParty(int party_id, const MKSecretKey* sk, const TFheGateBootstrappingParameterSet* params) {
        if (party_id < 0 || party_id >= k) return;

        const LweKey* lwe_k = sk->lwe_key;
        const TGswKey* rlwe_k = sk->rlwe_key;
        const TGswParams* tgsw_p = params->tgsw_params;

        for (int j = 0; j < n_per_party; ++j) {
            // LWE秘密鍵のビット s_j (0 or 1) を取得
            int32_t s_j = lwe_k->key[j];

            // s_j を TRGSW で暗号化 (一時変数は標準ドメイン)
            TGswSample* temp_bk = new_TGswSample(tgsw_p);
            tGswSymEncryptInt(temp_bk, s_j, params->tgsw_params->tlwe_params->alpha_min, rlwe_k);

            // FFTドメインに変換して格納
            tGswToFFTConvert(bk_fft[party_id][j], temp_bk, tgsw_p);

            delete_TGswSample(temp_bk);
        }
    }
};

} // namespace bbii

#endif // MK_TFHE_STRUCTS_H