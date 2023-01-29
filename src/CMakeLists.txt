# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

add_library(straw INTERFACE)

include(FetchContent)

FetchContent_Declare(
    filestream
    URL ${PROJECT_SOURCE_DIR}/external/filestream-fa7b6e1.tar.xz
    URL_HASH SHA256=6a4927277b551403b357e68b61e56b7d2d6a7e9eb8e65d4d768df72e8951f543)

FetchContent_Declare(
    fmt
    URL https://github.com/fmtlib/fmt/archive/refs/tags/9.1.0.tar.gz
    URL_HASH SHA256=5dea48d1fcddc3ec571ce2058e13910a0d4a6bab4cc09a809d8b1dd1c88ae6f2)

FetchContent_Declare(
    zlib-ng
    URL https://github.com/zlib-ng/zlib-ng/archive/refs/tags/2.0.6.tar.gz
    URL_HASH SHA256=8258b75a72303b661a238047cb348203d88d9dddf85d480ed885f375916fcab6)

set(ENABLE_TESTING OFF)
if(STRAW_USE_ZLIBNG)
    set(ZLIB_ENABLE_TESTS OFF)
    set(WITH_GTEST OFF)
    FetchContent_MakeAvailable(filestream fmt zlib-ng)
else()
    find_package(ZLIB REQUIRED)
    FetchContent_MakeAvailable(filestream fmt)
endif()

target_sources(
    straw
    INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/hic_file_impl.hpp
              ${CMAKE_CURRENT_SOURCE_DIR}/hic_file_stream_impl.hpp
              ${CMAKE_CURRENT_SOURCE_DIR}/matrix_zoom_data_impl.hpp
              ${CMAKE_CURRENT_SOURCE_DIR}/include/straw/straw.h
              ${CMAKE_CURRENT_SOURCE_DIR}/include/straw/internal/common.h)

target_include_directories(straw INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(
    straw
    INTERFACE straw_project_options
              filestream
              fmt::fmt-header-only
              $<IF:$<BOOL:${STRAW_USE_ZLIBNG}>,zlib,ZLIB::ZLIB>)

target_compile_definitions(straw INTERFACE FMT_ENFORCE_COMPILE_STRING $<$<BOOL:${STRAW_USE_ZLIBNG}>:STRAW_USE_ZLIBNG>)

