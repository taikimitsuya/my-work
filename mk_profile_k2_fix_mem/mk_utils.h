#ifndef MK_UTILS_H
#define MK_UTILS_H
#include <cstdlib>
#include <iostream>

// 32バイトアライメントでメモリを確保する関数
inline int32_t* alloc_aligned_array(size_t size) {
    void* ptr = nullptr;
    // posix_memalign(ポインタのアドレス, アライメント, サイズ)
    // TFHEは通常32バイト(AVX)または16バイト(SSE)のアライメントを好む
    if (posix_memalign(&ptr, 32, size * sizeof(int32_t)) != 0) {
        std::cerr << "Memory allocation failed!" << std::endl;
        exit(1);
    }
    return static_cast<int32_t*>(ptr);
}

// 解放用
inline void free_aligned_array(int32_t* ptr) {
    free(ptr);
}
#endif
