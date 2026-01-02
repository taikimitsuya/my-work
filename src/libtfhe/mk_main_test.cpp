#include <iostream>
#include <vector>
#include <cmath>
#include <cassert>

// これまで作成したヘッダ
#include "mk_params.h"
#include "mk_tfhe_structs.h"

// 外部関数の宣言 (ヘッダファイル化していないため)
namespace bbii {
    void mk_lwe_sym_encrypt(MKLweSample* result, Torus32 message, const MKSecretKey* sk, int32_t party_id, const TFheGateBootstrappingParameterSet* params);
    Torus32 mk_lwe_decrypt(const MKLweSample* ciphertext, const std::vector<MKSecretKey*>& all_keys, const TFheGateBootstrappingParameterSet* params);
    void mk_bootstrapping(MKLweSample* result, const MKLweSample* input, const MKBootstrappingKey* mk_bk, Torus32 mu, const TFheGateBootstrappingParameterSet* params);
}

using namespace bbii;
using namespace std;

// Sample Extract後の巨大LWE (次元 k*N) を復号して検証するヘルパー関数
// (mk_lwe_decrypt は 次元 k*n 用なので、ここではRLWE鍵を使って復号します)
double verify_extracted_sample(const bbii::MKLweSample* extracted, 
                               const vector<MKSecretKey*>& keys, 
                               const MKParams* mk_p) {
    int32_t N = mk_p->N;
    int32_t k = mk_p->k;
    
    // b' - sum(a'_i * s'_i)
    Torus32 phase = extracted->sample->b;

    for (int u = 0; u < k; ++u) {
        // ユーザーuのRLWE秘密鍵 (多項式) を取得
        const IntPolynomial* s_poly = /*&*/keys[u]->rlwe_key->key; // TGswKey -> IntPoly

        // LWEのaベクトル (サイズ k*N) のうち、このユーザーのブロック
        int offset = u * N;
        
        for (int i = 0; i < N; ++i) {
            Torus32 a_val = extracted->sample->a[offset + i];
            int32_t s_val = s_poly->coefs[i]; // 多項式の係数がそのままLWEの鍵になる
            
            phase -= a_val * s_val;
        }
    }
    
    return t32tod(phase);
}

int main() {
    cout << "=== Multi-Key FHE (MK-TFHE) Test Start ===" << endl;

    // 1. パラメータ設定
    // k=2 (2人), d=2, rho=3 (n=16), N=1024
    // ※ rhoが小さいと計算は早いが、安全性や精度は低い。テスト用設定。
    // 本格稼働時は rho=5以上推奨。
    int32_t k = 2;
    int32_t d = 2;
    int32_t rho = 5; 
    int32_t N = 1024;
    
    cout << "Generating Parameters (k=" << k << ", n=" << 2*pow(d,rho) << ", N=" << N << ")..." << endl;
    MKParams* mk_params = get_mk_test_params(k, d, rho, N);
    const TFheGateBootstrappingParameterSet* tfhe_params = mk_params->get_tfhe_params();

    // 2. 鍵生成 (Secret Key & Bootstrapping Key)
    cout << "Generating Keys..." << endl;
    vector<MKSecretKey*> secret_keys(k);
    MKBootstrappingKey* mk_bk = new MKBootstrappingKey(k, mk_params->n_per_party, tfhe_params);

    for(int i=0; i<k; ++i) {
        // 秘密鍵生成
        secret_keys[i] = new MKSecretKey(tfhe_params);
        // BK生成
        mk_bk->generateKeyForParty(i, secret_keys[i], tfhe_params);
    }

    // 3. 暗号化 (Encryption)
    // Party 0 が メッセージ 0.25 (1/4) を暗号化する
    cout << "Encrypting message..." << endl;
    bbii::MKLweSample* input_ct = new bbii::MKLweSample(k, mk_params->n_per_party, tfhe_params);
    
    double message_double = 0.25;
    Torus32 message = dtot32(message_double);
    
    // Party 0 が暗号化 (他のパーティ部分は0になる)
    mk_lwe_sym_encrypt(input_ct, message, secret_keys[0], 0, tfhe_params);

    // 暗号化直後の復号確認 (n次元)
    Torus32 decrypted_raw = mk_lwe_decrypt(input_ct, secret_keys, tfhe_params);
    cout << "  Input Decrypted Phase: " << t32tod(decrypted_raw) 
         << " (Expected: " << message_double << ")" << endl;


    // 4. MK Bootstrapping (Blind Rotate + Extract)
    // テストベクトルとして、メッセージ v = 0.25 をセットして回転させる
    // 成功すれば、出力の位相は 0.25 付近になるはず (エラーは蓄積する)
    cout << "Running MK Bootstrapping (This may take time)..." << endl;
    
    // 出力用: 次元は k * N になる (Sample Extractの結果)
    bbii::MKLweSample* output_ct = new bbii::MKLweSample(k, N, tfhe_params); // 注意: 第2引数は N
    
    // テストベクトル v = 0.25 (回転の基準)
    // 実際には入力暗号文の位相分だけ回転するため、結果は v * X^{-message} の定数項などになるが
    // ここでは簡易テストとして「入力の位相」を正しく反映して回転できたかを確認する。
    // Bootstrappingのテスト: v = 0 (test vector) に入力 0.25 を加えると...
    // 通常のGate Bootstrappingの構成とは異なり、ここでは
    // 「入力暗号文によってアキュムレータを回転させる機能」を確認する。
    
    // テスト: アキュムレータの初期値を v=0.5 とし、入力暗号文(0.25)で回転させる
    // 結果の位相が変化することを確認する。
    mk_bootstrapping(output_ct, input_ct, mk_bk, dtot32(0.5), tfhe_params);

    // 5. 結果検証
    // Sample ExtractされたLWE (次元 k*N) を復号
    double result_phase = verify_extracted_sample(output_ct, secret_keys, mk_params);
    
    cout << "Result Phase: " << result_phase << endl;
    cout << "  (Note: Value depends on rotation amount and test vector logic)" << endl;
    cout << "  Rotation worked if result is not random noise." << endl;

    // クリーンアップ
    delete input_ct;
    delete output_ct;
    delete mk_bk;
    delete mk_params; // base paramsも削除される
    for(auto sk : secret_keys) delete sk;

    cout << "Test Finished Successfully." << endl;
    return 0;
}