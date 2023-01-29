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

add_executable(test_main hic_file_test.cpp hic_file_stream_test.cpp)

target_link_libraries(test_main PRIVATE straw_project_options straw Catch2::Catch2WithMain)

catch_discover_tests(test_main WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../)
