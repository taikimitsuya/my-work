#include "bb_utils.h"

namespace bbii {

bool is_power_of_two(int32_t n) {
    return n > 0 && (n & (n - 1)) == 0;
}

int32_t log2_int(int32_t n) {
    if (!is_power_of_two(n)) {
        throw std::invalid_argument("Input must be power of two");
    }
    int32_t res = 0;
    while (n >>= 1) ++res;
    return res;
}

std::vector<Complex> get_roots_of_unity(int32_t n) {
    std::vector<Complex> roots(n);
    const double pi = std::acos(-1);
    for (int i = 0; i < n; ++i) {
        double angle = 2.0 * pi * i / n;
        roots[i] = Complex(std::cos(angle), std::sin(angle));
    }
    return roots;
}

RotatedIndex get_anti_cyclic_index(int32_t idx, int32_t shift, int32_t N) {
    // 論文 Section 2.5: Anti-Rot(a, z)
    // 多項式環 X^N + 1 = 0 (mod q) 上での回転
    // X^i * X^shift = X^{i+shift}
    
    int32_t raw_idx = idx + shift;
    int32_t div = raw_idx / N;
    int32_t mod = raw_idx % N;
    
    // Nを超えるごとに符号が反転する (X^N = -1)
    // div が偶数なら符号そのまま(+1), 奇数なら反転(-1)
    int32_t sign = (div % 2 == 0) ? 1 : -1;
    
    return {mod, sign};
}

uint32_t reverse_bits(uint32_t x, int32_t bits) {
    uint32_t reversed = 0;
    for (int i = 0; i < bits; ++i) {
        if ((x >> i) & 1) {
            reversed |= (1 << (bits - 1 - i));
        }
    }
    return reversed;
}

} // namespace bbii