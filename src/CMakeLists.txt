# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

find_package(ZLIB REQUIRED)

add_library(straw INTERFACE)

include(FetchContent)

FetchContent_Declare(
        filestream
        URL ${PROJECT_SOURCE_DIR}/external/filestream-fa7b6e1.tar.xz
        URL_HASH SHA256=6a4927277b551403b357e68b61e56b7d2d6a7e9eb8e65d4d768df72e8951f543)

set(ENABLE_TESTING OFF)
FetchContent_MakeAvailable(filestream)

target_sources(
    straw
    INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/straw_impl.hpp
              ${CMAKE_CURRENT_SOURCE_DIR}/hic_file_impl.hpp
              ${CMAKE_CURRENT_SOURCE_DIR}/hic_file_stream_impl.hpp
              ${CMAKE_CURRENT_SOURCE_DIR}/matrix_zoom_data_impl.hpp
              ${CMAKE_CURRENT_SOURCE_DIR}/include/straw/straw.h
              ${CMAKE_CURRENT_SOURCE_DIR}/include/straw/internal/common.h)

target_include_directories(straw INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(straw INTERFACE straw_project_options filestream ZLIB::ZLIB)
