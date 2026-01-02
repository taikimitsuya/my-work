#include "mk_tfhe_structs.h"
#include <vector>
#include <cmath>
#include <random>

// libtfheの内部関数を利用するためのinclude
// (乱数生成などのため)
#include <tfhe.h>
#include <tfhe_core.h>

namespace bbii {

/**
 * @brief MK-LWE 対称鍵暗号化
 * * 指定されたパーティ(party_id)の秘密鍵を使ってメッセージを暗号化します。
 * 結果のMK-LWE暗号文において、party_idに対応するブロック以外のa部分は0になります。
 * * @param result     [out] 暗号化結果を格納するMKLweSample
 * @param message    [in]  暗号化する平文 (Torus32形式)
 * @param sk         [in]  暗号化を行うパーティの秘密鍵
 * @param party_id   [in]  暗号化を行うパーティのID (0 ~ k-1)
 * @param params     [in]  パラメータセット
 */
void mk_lwe_sym_encrypt(MKLweSample* result, 
                        Torus32 message, 
                        const MKSecretKey* sk, 
                        int32_t party_id, 
                        const TFheGateBootstrappingParameterSet* params) {
    
    // パラメータ取得
    const LweParams* lwe_params = params->in_out_params;
    int32_t n = result->n_per_party; // 1人あたりの次元 (例: 512)
    int32_t k = result->k;           // 全パーティ数
    int32_t total_n = n * k;         // 全体の次元

    // 1. 暗号文全体(aベクトル)をゼロリセット
    for (int i = 0; i < total_n; ++i) {
        result->sample->a[i] = 0;
    }

    // 2. 指定されたパーティの領域だけ乱数(a_i)をセットし、内積(a_i * s_i)を計算
    // offset: 配列全体の中での、このパーティの開始位置
    int32_t offset = party_id * n;
    Torus32 multi_key_product = 0;

    // a_i を一様乱数で生成しつつ、s_i との内積をとる
    for (int i = 0; i < n; ++i) {
        // 一様乱数生成 (libtfheの関数を使用)
        // a_i の要素
        Torus32 a_val = split_random_generator::lwe_random(); // もしくは uniform_random()
        
        // 結果の該当位置にセット
        result->sample->a[offset + i] = a_val;

        // 内積計算: a_val * s_i[i]
        // sk->lwe_key->key[i] は 0 か 1 なので、1のときだけ加算
        if (sk->lwe_key->key[i] == 1) {
            multi_key_product += a_val;
        }
    }

    // 3. エラー(Gaussian noise)を生成
    double alpha = lwe_params->alpha_min;
    Torus32 error = gaussian32(0, alpha); // libtfheの関数

    // 4. b = a*s + m + e を計算
    // MK-LWEの構成上、他のパーティのaは0なので、計算した multi_key_product だけでよい
    result->sample->b = multi_key_product + message + error;

    // 現在の分散値を記録（TFHEの仕様に準拠）
    result->sample->current_variance = alpha * alpha;
}

/**
 * @brief MK-LWE 復号
 * * 参加者全員の秘密鍵を使用して、MK-LWE暗号文を復号します。
 * (シミュレーション用、または最終的な結合復号処理用)
 * * @param ciphertext [in] 復号対象のMKLweSample
 * @param all_keys   [in] 全パーティのMKSecretKeyポインタのリスト
 * @param params     [in] パラメータセット
 * @return Torus32        復号された位相(Phase)
 */
Torus32 mk_lwe_decrypt(const MKLweSample* ciphertext, 
                       const std::vector<MKSecretKey*>& all_keys, 
                       const TFheGateBootstrappingParameterSet* params) {
    
    int32_t n = ciphertext->n_per_party;
    int32_t k = ciphertext->k;

    // エラーチェック
    if (all_keys.size() != (size_t)k) {
        throw std::runtime_error("Number of keys provided does not match the ciphertext parameters.");
    }

    // 1. b で初期化
    Torus32 phase = ciphertext->sample->b;

    // 2. 各パーティごとの成分 (a_j * s_j) を減算する
    //    Phase = b - sum(a_j * s_j)
    for (int j = 0; j < k; ++j) {
        int32_t offset = j * n;
        const LweKey* current_sk = all_keys[j]->lwe_key;

        for (int i = 0; i < n; ++i) {
            // 暗号文の a ベクトルの該当部分
            Torus32 a_val = ciphertext->sample->a[offset + i];
            
            // 秘密鍵との積を引く
            if (current_sk->key[i] == 1) {
                phase -= a_val;
            }
        }
    }

    // 3. 復号結果（Phase）を返す
    // 必要に応じて、呼び出し元で approxPhase(phase) などで丸めを行う
    return phase;
}

/**
 * @brief (テスト用) 自明なMK-LWE暗号文の作成
 * ノイズなし、a=0 でメッセージのみを持つ暗号文を作成します。
 */
void mk_lwe_noiseless_trivial(MKLweSample* result, Torus32 message) {
    int32_t total_n = result->k * result->n_per_party;
    
    for (int i = 0; i < total_n; ++i) {
        result->sample->a[i] = 0;
    }
    result->sample->b = message;
    result->sample->current_variance = 0.0;
}

} // namespace bbii