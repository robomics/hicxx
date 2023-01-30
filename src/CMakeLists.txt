# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

add_library(hicxx INTERFACE)

include(FetchContent)

FetchContent_Declare(
    fmt
    URL https://github.com/fmtlib/fmt/archive/refs/tags/9.1.0.tar.gz
    URL_HASH SHA256=5dea48d1fcddc3ec571ce2058e13910a0d4a6bab4cc09a809d8b1dd1c88ae6f2)

FetchContent_Declare(
    zlib-ng
    URL https://github.com/zlib-ng/zlib-ng/archive/refs/tags/2.0.6.tar.gz
    URL_HASH SHA256=8258b75a72303b661a238047cb348203d88d9dddf85d480ed885f375916fcab6)

if(HICXX_USE_CURL)
    find_package(CURL REQUIRED)
endif()

if(HICXX_USE_ZLIBNG)
    set(ZLIB_ENABLE_TESTS OFF)
    set(WITH_GTEST OFF)
    FetchContent_MakeAvailable(fmt zlib-ng)
else()
    find_package(ZLIB REQUIRED)
    FetchContent_MakeAvailable(fmt)
endif()

target_sources(
    hicxx
    INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/hic_file_impl.hpp
              ${CMAKE_CURRENT_SOURCE_DIR}/hic_file_stream_impl.hpp
              ${CMAKE_CURRENT_SOURCE_DIR}/hic_footer_impl.hpp
              ${CMAKE_CURRENT_SOURCE_DIR}/hic_header_impl.hpp
              ${CMAKE_CURRENT_SOURCE_DIR}/hic_matrix_selector_impl.hpp
              ${CMAKE_CURRENT_SOURCE_DIR}/include/hicxx/hicxx.hpp
              ${CMAKE_CURRENT_SOURCE_DIR}/include/hicxx/internal/common.hpp
              ${CMAKE_CURRENT_SOURCE_DIR}/include/hicxx/internal/hic_footer.hpp
              ${CMAKE_CURRENT_SOURCE_DIR}/include/hicxx/internal/hic_matrix_selector.hpp
              ${CMAKE_CURRENT_SOURCE_DIR}/include/hicxx/internal/hic_file_stream.hpp
              ${CMAKE_CURRENT_SOURCE_DIR}/include/hicxx/internal/hic_header.hpp)

target_include_directories(hicxx INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(
    hicxx
    INTERFACE hicxx_project_options
              hicxx_project_warnings
              fmt::fmt-header-only
              $<$<BOOL:${HICXX_USE_CURL}>:CURL::libcurl>
              $<IF:$<BOOL:${HICXX_USE_ZLIBNG}>,zlib,ZLIB::ZLIB>)

target_compile_definitions(hicxx INTERFACE FMT_ENFORCE_COMPILE_STRING $<$<BOOL:${HICXX_USE_CURL}>:HICXX_USE_CURL>
                                           $<$<BOOL:${HICXX_USE_ZLIBNG}>:HICXX_USE_ZLIBNG>)
