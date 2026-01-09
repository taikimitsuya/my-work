cmake_minimum_required(VERSION 3.10)
project(mk_profile_k3_final CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(SOURCES
    main_mk_test.cpp
    mk_bootstrapping.cpp
    mk_lwe.cpp
    mk_ops.cpp
    mk_packed_ops.cpp
)
add_executable(main_mk_test ${SOURCES})

# TFHE include
include_directories(${CMAKE_SOURCE_DIR}/../src/include)

# TFHE static lib

# TFHE static libを絶対パスで直接リンク（libtfheサブディレクトリ）
target_link_libraries(main_mk_test PRIVATE /workspaces/my-work/tfhe/build/libtfhe/libtfhe_static.a -lm -lpthread)

# 依存する他のヘッダも必要ならinclude_directoriesで追加
