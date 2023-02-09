# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT


cmake_minimum_required(VERSION 3.18)
cmake_policy(VERSION 3.18...3.25)

project(
        hicxx
        LANGUAGES CXX
        VERSION 0.0.1
        HOMEPAGE_URL https://github.com/robomics/hicxx
        DESCRIPTION "hicxx: C++17 library and CLI tool to read .hic files.")

set(CMAKE_CXX_STANDARD 17)

include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/CompilerWarnings.cmake)
include(cmake/TestStdlibFeatures.cmake)

option(HICXX_BUILD_EXAMPLES "Build hicxx example(s)" OFF)
option(HICXX_ENABLE_TESTING "Build unit tests" ON)
option(HICXX_USE_ZLIBNG "Use zlib-ng instead of zlib for decompression" OFF)
option(HICXX_USE_CURL "Use cURL to read remote .hic files" OFF)

add_library(hicxx_project_options INTERFACE)
target_compile_features(hicxx_project_options INTERFACE cxx_std_${CMAKE_CXX_STANDARD})

add_library(hicxx_project_warnings INTERFACE)
target_compile_options(
        hicxx_project_warnings
        INTERFACE # C++ warnings
        $<$<COMPILE_LANGUAGE:CXX>:${HICXX_PROJECT_WARNINGS_CXX}>
        # C warnings
        $<$<COMPILE_LANGUAGE:C>:${HICXX_PROJECT_WARNINGS_C}>
        # Cuda warnings
        $<$<COMPILE_LANGUAGE:CUDA>:${HICXX_PROJECT_WARNINGS_CUDA}>)

add_subdirectory(src)

if(HICXX_BUILD_EXAMPLES)
    add_subdirectory(examples)
endif()

if(HICXX_ENABLE_TESTING)
    include(cmake/FetchTestDatasets.cmake)
    enable_testing()
    add_subdirectory(test)
endif()
