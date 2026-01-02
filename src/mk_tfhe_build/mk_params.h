#ifndef MK_PARAMS_H
#define MK_PARAMS_H

#include "bb_params.h"
#include <vector>
#include <stdexcept>

namespace bbii {

/**
 * @brief マルチキーFHEのパラメータ構造体
 * * 既存の BBIIParams (シングルキー設定) を保持しつつ、
 * 参加者数 (k) に基づく拡張されたパラメータを管理します。
 */
struct MKParams {
    // 参加者数 (Number of parties)
    int32_t k;

    // 1人あたりのLWE次元 (Single-key n)
    // rho によって決まる値 (例: rho=5 -> n=64)
    int32_t n_per_party;

    // 全体のLWE次元 (Total n = k * n_per_party)
    int32_t total_n;

    // RLWEの多項式次数 (N)
    int32_t N;

    // 基礎となるシングルキーパラメータへのポインタ
    // (LWEパラメータ、BKパラメータ、分解パラメータd, rhoなどを含む)
    BBIIParams* sk_params;

    /**
     * @brief コンストラクタ
     * @param parties 参加者数
     * @param base_params 基礎となるBBIIパラメータ (d, rho, N 等の設定済み)
     */
    MKParams(int32_t parties, BBIIParams* base_params) : k(parties), sk_params(base_params) {
        if (parties < 1) {
            throw std::runtime_error("Number of parties (k) must be at least 1.");
        }
        if (!base_params) {
            throw std::runtime_error("Base BBIIParams cannot be null.");
        }

        // シングルキーパラメータから値をコピー
        this->n_per_party = base_params->n;
        this->N = base_params->N;

        // マルチキー全体の次元を計算
        this->total_n = this->k * this->n_per_party;
    }

    /**
     * @brief デストラクタ
     * sk_params の所有権がこのクラスにあると仮定して解放します。
     * (外部で管理する場合はコメントアウトしてください)
     */
    ~MKParams() {
        if (sk_params) {
            delete sk_params;
            sk_params = nullptr;
        }
    }

    // 便利なアクセサ: TFHEの低レベルパラメータセットを取得
    TFheGateBootstrappingParameterSet* get_tfhe_params() const {
        return sk_params->tfhe_params;
    }
};

/**
 * @brief テスト用のマルチキーパラメータ生成関数
 * * 指定された構成で MKParams を生成して返します。
 * * @param k   パーティ数 (例: 2, 3...)
 * @param d   分解の底 (例: 2)
 * @param rho 分解の回数/深さ (例: 5) -> n = 2 * d^rho
 * @param N   RLWE環の次数 (例: 1024, 2048)
 * @return MKParams* 生成されたパラメータへのポインタ
 */
inline MKParams* get_mk_test_params(int32_t k, int32_t d, int32_t rho, int32_t N) {
    // 1. まずシングルキー用のパラメータを生成
    //    (bb_params.h の構造体を使用)
    BBIIParams* base = new BBIIParams(d, rho, N);

    // 2. マルチキー用ラッパーで包んで返す
    return new MKParams(k, base);
}

} // namespace bbii

#endif // MK_PARAMS_H