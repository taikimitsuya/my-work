#ifndef MK_UTILS_H
#define MK_UTILS_H
#include <cstdlib>
#include <iostream>
#include <tfhe.h>

// 32バイトアライメント確保
inline int32_t* my_new_Torus32_array(int32_t size) {
    void* ptr = nullptr;
    if (posix_memalign(&ptr, 32, size * sizeof(int32_t)) != 0) {
        std::cerr << "Memory allocation failed!" << std::endl;
        exit(1);
    }
    return static_cast<int32_t*>(ptr);
}
inline void my_delete_Torus32_array(int32_t* ptr) {
    free(ptr);
}
#endif
