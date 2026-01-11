#include "mk_params.h"
#include "bb_params.h"

MKParams* get_mk_test_params(int32_t k, int32_t d, int32_t rho, int32_t N) {
    return new MKParams(k, new BBIIParams(d, rho, N));
}
