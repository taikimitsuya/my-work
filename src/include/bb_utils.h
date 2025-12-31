// bb_utils.h
#ifndef BB_UTILS_H
#define BB_UTILS_H

#include <vector>
#include <cstdint>
#include <stdexcept>

namespace bbii {

// 1次元配列を d x (n/d) の行列とみなして転置する (m -> 2d)
template <typename T>
std::vector<T> rearrange(const std::vector<T>& input, int32_t d) {
    if (input.empty()) return {};
    size_t n = input.size();
    if (n % d != 0) throw std::invalid_argument("Input size must be multiple of d");

    size_t rows = d;
    size_t cols = n / d;
    std::vector<T> output(n);
    
    // Output[j*rows + i] = Input[i*cols + j]
    for (size_t i = 0; i < rows; ++i) {
        for (size_t j = 0; j < cols; ++j) {
            output[j * rows + i] = input[i * cols + j];
        }
    }
    return output;
}

// 逆転置 (2d -> m)
template <typename T>
std::vector<T> reverse_rearrange(const std::vector<T>& input, int32_t d) {
    if (input.empty()) return {};
    size_t n = input.size();
    if (n % d != 0) throw std::invalid_argument("Input size must be multiple of d");

    size_t rows = d; 
    size_t cols = n / d;
    std::vector<T> output(n);
    
    // Output[i*cols + j] = Input[j*rows + i]
    for (size_t j = 0; j < cols; ++j) {
        for (size_t i = 0; i < rows; ++i) {
            output[i * cols + j] = input[j * rows + i];
        }
    }
    return output;
}

} // namespace bbii
#endif