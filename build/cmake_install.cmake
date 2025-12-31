# Install script for directory: /workspaces/my-work/src

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "optim")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "1")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

# Set default install directory permissions.
if(NOT DEFINED CMAKE_OBJDUMP)
  set(CMAKE_OBJDUMP "/usr/bin/objdump")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/tfhe" TYPE FILE PERMISSIONS OWNER_READ OWNER_WRITE GROUP_READ WORLD_READ FILES
    "/workspaces/my-work/src/include/BatchOperations.h"
    "/workspaces/my-work/src/include/BatchParams.h"
    "/workspaces/my-work/src/include/NussbaumerDFT.h"
    "/workspaces/my-work/src/include/TensorRGSW.h"
    "/workspaces/my-work/src/include/batch_bootstrapping.h"
    "/workspaces/my-work/src/include/batch_framework.h"
    "/workspaces/my-work/src/include/batch_ops.h"
    "/workspaces/my-work/src/include/bb_params.h"
    "/workspaces/my-work/src/include/bb_utils.h"
    "/workspaces/my-work/src/include/hom_dft.h"
    "/workspaces/my-work/src/include/lagrangehalfc_arithmetic.h"
    "/workspaces/my-work/src/include/lwe-functions.h"
    "/workspaces/my-work/src/include/lwebootstrappingkey.h"
    "/workspaces/my-work/src/include/lwekey.h"
    "/workspaces/my-work/src/include/lwekeyswitch.h"
    "/workspaces/my-work/src/include/lweparams.h"
    "/workspaces/my-work/src/include/lwesamples.h"
    "/workspaces/my-work/src/include/mkTFHEfunctions.h"
    "/workspaces/my-work/src/include/mkTFHEkeygen.h"
    "/workspaces/my-work/src/include/mkTFHEkeys.h"
    "/workspaces/my-work/src/include/mkTFHEparams.h"
    "/workspaces/my-work/src/include/mkTFHEsamples.h"
    "/workspaces/my-work/src/include/numeric_functions.h"
    "/workspaces/my-work/src/include/polynomials.h"
    "/workspaces/my-work/src/include/polynomials_arithmetic.h"
    "/workspaces/my-work/src/include/tfhe.h"
    "/workspaces/my-work/src/include/tfhe_core.h"
    "/workspaces/my-work/src/include/tfhe_garbage_collector.h"
    "/workspaces/my-work/src/include/tfhe_gate_bootstrapping_functions.h"
    "/workspaces/my-work/src/include/tfhe_gate_bootstrapping_structures.h"
    "/workspaces/my-work/src/include/tfhe_generic_streams.h"
    "/workspaces/my-work/src/include/tfhe_generic_templates.h"
    "/workspaces/my-work/src/include/tfhe_io.h"
    "/workspaces/my-work/src/include/tgsw.h"
    "/workspaces/my-work/src/include/tgsw_functions.h"
    "/workspaces/my-work/src/include/tlwe.h"
    "/workspaces/my-work/src/include/tlwe_functions.h"
    )
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("/workspaces/my-work/build/libtfhe/cmake_install.cmake")

endif()

if(CMAKE_INSTALL_COMPONENT)
  set(CMAKE_INSTALL_MANIFEST "install_manifest_${CMAKE_INSTALL_COMPONENT}.txt")
else()
  set(CMAKE_INSTALL_MANIFEST "install_manifest.txt")
endif()

string(REPLACE ";" "\n" CMAKE_INSTALL_MANIFEST_CONTENT
       "${CMAKE_INSTALL_MANIFEST_FILES}")
file(WRITE "/workspaces/my-work/build/${CMAKE_INSTALL_MANIFEST}"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
