# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT


set(SOURCE_FILES examples/main.cpp src/hicxx.cpp)
add_executable(hicxx_example example.cpp)
target_link_libraries(
    hicxx_example
    PRIVATE hicxx_project_options
    PUBLIC hicxx)
