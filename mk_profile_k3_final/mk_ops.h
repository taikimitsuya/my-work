

#ifndef MK_OPS_H
#define MK_OPS_H

// 標準ライブラリ
#include <cstdint>

// プロジェクトヘッダ
#include "mk_tfhe_structs.h"


// RLWEクリア
void mk_rlwe_clear(MKRLweSample* r);
// 外部積
void mk_external_product(MKRLweSample* res, const MKPackedRGSW* rgsw, const MKRLweSample* acc, int32_t pid, const TFheGateBootstrappingParameterSet* params);
// RLWEコピー
void mk_rlwe_copy(MKRLweSample* dest, const MKRLweSample* src);
void mk_rlwe_addTo(MKRLweSample* r, const MKRLweSample* s);
void mk_rlwe_subTo(MKRLweSample* r, const MKRLweSample* s);
void mk_external_product(MKRLweSample* res, const TGswSampleFFT* bk, const MKRLweSample* acc, int32_t pid, const TFheGateBootstrappingParameterSet* p);
void mk_cmux(MKRLweSample* res, const TGswSampleFFT* bk, const MKRLweSample* in0, const MKRLweSample* in1, int32_t pid, const TFheGateBootstrappingParameterSet* p);


#endif
