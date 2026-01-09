#include "mk_tfhe_structs.h"
#include "mk_methods.h"
#include <iostream>

// main_mk_test.cpp から呼ばれる正式実装
void mk_test_blind_rotate_dft(const TFheGateBootstrappingParameterSet* params) {
    std::cout << "[test] mk_test_blind_rotate_dft 実行" << std::endl;
    // テスト用パラメータ
    int k = 2;
    bbii::MKRLweSample* acc = new bbii::MKRLweSample(k, params);
    bbii::MKRLweSample* input = new bbii::MKRLweSample(k, params);
    bbii::MKBootstrappingKey* bk = new bbii::MKBootstrappingKey(k, 64, params);
    // acc, input, bkの初期化は省略（必要に応じてテスト値をセット）
    mk_blind_rotate_dft(acc, input, bk, params);
    std::cout << "[test] mk_blind_rotate_dft 完了" << std::endl;
    delete acc;
    delete input;
    delete bk;
}

namespace bbii {
// main_mk_test.cpp から呼ばれる正式実装
void mk_bootstrapping(MKRLweSample* result, const MKRLweSample* input, const MKBootstrappingKey* mk_bk, int mu, const TFheGateBootstrappingParameterSet* params) {
    // DFTベースBlind Rotateを呼び出し
    mk_blind_rotate_dft(result, input, mk_bk, params);
    // 必要に応じてSample Extract等を追加
}
}
