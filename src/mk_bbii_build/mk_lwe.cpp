#include "mk_tfhe_structs.h"
#include <vector>
#include <cmath>
#include <random> // ★追加: 標準乱数ライブラリ
#include <limits> // ★追加: 型の最大値・最小値取得のため

// libtfheの内部関数を利用するためのinclude
#include <tfhe.h>
#include <tfhe_core.h>

namespace bbii {

// ★ヘルパー関数: 安全な乱数生成器
// 関数内でstatic宣言することで、呼び出し毎の再初期化を防ぎます
Torus32 get_uniform_random_torus32() {
    static std::random_device rd;
    static std::mt19937 gen(rd());
    // int32_t (Torus32) の全範囲で一様乱数を生成
    static std::uniform_int_distribution<int32_t> dist(
        std::numeric_limits<int32_t>::min(), 
        std::numeric_limits<int32_t>::max()
    );
    return dist(gen);
}

void mk_lwe_sym_encrypt(MKLweSample* result, 
                        Torus32 message, 
                        const MKSecretKey* sk, 
                        int32_t party_id, 
                        const TFheGateBootstrappingParameterSet* params) {
    
    // パラメータ取得
    const LweParams* lwe_params = params->in_out_params;
    int32_t n = result->n_per_party; 
    int32_t k = result->k;           
    int32_t total_n = n * k;         

    // 1. 暗号文全体(aベクトル)をゼロリセット
    for (int i = 0; i < total_n; ++i) {
        result->sample->a[i] = 0;
    }

    // 2. 指定されたパーティの領域だけ乱数(a_i)をセットし、内積(a_i * s_i)を計算
    int32_t offset = party_id * n;
    Torus32 multi_key_product = 0;

    for (int i = 0; i < n; ++i) {
        // ★修正箇所: C++標準ライブラリで乱数生成
        Torus32 a_val = get_uniform_random_torus32();
        
        // 結果の該当位置にセット
        result->sample->a[offset + i] = a_val;

        // 内積計算
        if (sk->lwe_key->key[i] == 1) {
            multi_key_product += a_val;
        }
    }

    // 3. エラー(Gaussian noise)を生成
    // gaussian32 は tfhe_core.h に含まれているはずですが、
    // もしここでもエラーが出る場合は、同様に標準ライブラリで実装可能です。
    double alpha = lwe_params->alpha_min;
    Torus32 error = gaussian32(0, alpha); 

    // 4. b = a*s + m + e を計算
    result->sample->b = multi_key_product + message + error;

    result->sample->current_variance = alpha * alpha;
}

Torus32 mk_lwe_decrypt(const MKLweSample* ciphertext, 
                       const std::vector<MKSecretKey*>& all_keys, 
                       const TFheGateBootstrappingParameterSet* params) {
    
    int32_t n = ciphertext->n_per_party;
    int32_t k = ciphertext->k;

    if (all_keys.size() != (size_t)k) {
        throw std::runtime_error("Number of keys provided does not match the ciphertext parameters.");
    }

    Torus32 phase = ciphertext->sample->b;

    for (int j = 0; j < k; ++j) {
        int32_t offset = j * n;
        const LweKey* current_sk = all_keys[j]->lwe_key;

        for (int i = 0; i < n; ++i) {
            Torus32 a_val = ciphertext->sample->a[offset + i];
            
            if (current_sk->key[i] == 1) {
                phase -= a_val;
            }
        }
    }

    return phase;
}

void mk_lwe_noiseless_trivial(MKLweSample* result, Torus32 message) {
    int32_t total_n = result->k * result->n_per_party;
    
    for (int i = 0; i < total_n; ++i) {
        result->sample->a[i] = 0;
    }
    result->sample->b = message;
    result->sample->current_variance = 0.0;
}

} // namespace bbii