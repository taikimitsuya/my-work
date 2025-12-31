#ifndef BB_UTILS_H
#define BB_UTILS_H
#include <vector>
#include <cstdint>
#include <stdexcept>
namespace bbii {
template <typename T>
std::vector<T> rearrange(const std::vector<T>& input, int32_t d) {
    if (input.empty()) return {};
    size_t n = input.size();
    if (n % d != 0) throw std::invalid_argument("Input size must be multiple of d");
    size_t rows = d; size_t cols = n / d;
    std::vector<T> output(n);
    for (size_t i = 0; i < rows; ++i) {
        for (size_t j = 0; j < cols; ++j) output[j * rows + i] = input[i * cols + j];
    }
    return output;
}
template <typename T>
std::vector<T> reverse_rearrange(const std::vector<T>& input, int32_t d) {
    if (input.empty()) return {};
    size_t n = input.size();
    if (n % d != 0) throw std::invalid_argument("Input size must be multiple of d");
    size_t rows = d; size_t cols = n / d;
    std::vector<T> output(n);
    for (size_t j = 0; j < cols; ++j) {
        for (size_t i = 0; i < rows; ++i) output[i * cols + j] = input[j * rows + i];
    }
    return output;
}
}
#endif
