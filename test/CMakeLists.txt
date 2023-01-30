# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

include(FetchContent)

FetchContent_Declare(
    Catch2
    GIT_REPOSITORY https://github.com/catchorg/Catch2.git
    GIT_TAG v3.3.0)

FetchContent_MakeAvailable(Catch2)

include(CTest)
include(Catch)

add_executable(straw_test_main)

target_sources(
    straw_test_main
    PRIVATE filestream_local_test.cpp
            filestream_remote_test.cpp
            hic_file_test.cpp
            hic_file_stream_test.cpp
            matrix_zoom_data_test.cpp)

target_link_libraries(
    straw_test_main
    PRIVATE straw_project_options
            straw_project_warnings
            straw
            Catch2::Catch2WithMain)

catch_discover_tests(straw_test_main WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../)
