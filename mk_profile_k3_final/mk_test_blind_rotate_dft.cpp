#include "mk_tfhe_structs.h"
#include "mk_packed_ops.h"
#include "mk_bootstrapping.cpp"
#include <iostream>

using namespace bbii;

void mk_test_blind_rotate_dft(const TFheGateBootstrappingParameterSet* params) {
    int k = 2, N = params->in_out_params->n;
    std::cout << "[TEST] DFT-based Blind Rotate" << std::endl;
    MKRLweSample acc(k, params);
    // テスト用: accに単純な値をセット
    for(int u=0; u<=k; ++u) for(int i=0; i<N; ++i) acc.parts[u]->coefsT[i] = (i==0 ? 1 : 0);
    MKRLweSample bk_input(k, params);
    for(int u=0; u<=k; ++u) for(int i=0; i<N; ++i) bk_input.parts[u]->coefsT[i] = 0;
    MKBootstrappingKey bk(k, N, params);
    mk_blind_rotate_dft(&acc, &bk_input, &bk, params);
    std::cout << "[TEST] acc after DFT-based Blind Rotate: ";
    for(int i=0; i<8 && i<N; ++i) std::cout << acc.parts[0]->coefsT[i] << " ";
    std::cout << "..." << std::endl;
}
