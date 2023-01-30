# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT


cmake_minimum_required(VERSION 3.17)
cmake_policy(VERSION 3.18...3.25)
project(strawC)
set(CMAKE_CXX_STANDARD 17)

include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/CompilerWarnings.cmake)
include(cmake/TestStdlibFeatures.cmake)

option(STRAW_BUILD_EXAMPLES "Build straw example(s)" OFF)
option(STRAW_ENABLE_TESTING "Build unit tests" ON)
option(STRAW_USE_ZLIBNG "Use zlib-ng instead of zlib for decompression" OFF)
option(STRAW_USE_CURL "Use cURL to read remote .hic files" OFF)

add_library(straw_project_options INTERFACE)
target_compile_features(straw_project_options INTERFACE cxx_std_${CMAKE_CXX_STANDARD})

add_library(straw_project_warnings INTERFACE)
target_compile_options(
        straw_project_warnings
        INTERFACE # C++ warnings
        $<$<COMPILE_LANGUAGE:CXX>:${STRAW_PROJECT_WARNINGS_CXX}>
        # C warnings
        $<$<COMPILE_LANGUAGE:C>:${STRAW_PROJECT_WARNINGS_C}>
        # Cuda warnings
        $<$<COMPILE_LANGUAGE:CUDA>:${STRAW_PROJECT_WARNINGS_CUDA}>)

add_subdirectory(src)

if(STRAW_BUILD_EXAMPLES)
    add_subdirectory(examples)
endif()

if(STRAW_ENABLE_TESTING)
    enable_testing()
    add_subdirectory(test)
endif()
