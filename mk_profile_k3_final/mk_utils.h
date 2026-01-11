#ifndef MK_UTILS_H
#define MK_UTILS_H
#include <cstdlib>
#include <iostream>
#include <tfhe/tfhe.h>

// 32バイト(256bit)アライメントでメモリを確保する
// これにより AVX 命令でのクラッシュを防ぐ
inline int32_t* my_new_Torus32_array(int32_t size) {
    void* ptr = nullptr;
    // posix_memalign(ポインタのアドレス, アライメント, サイズ)
    if (posix_memalign(&ptr, 32, size * sizeof(int32_t)) != 0) {
        std::cerr << "Memory allocation failed!" << std::endl;
        exit(1);
    }
    return static_cast<int32_t*>(ptr);
}

// 解放関数 (posix_memalign で確保したメモリは free で解放可能)
inline void my_delete_Torus32_array(int32_t* ptr) {
    free(ptr);
}
#endif
