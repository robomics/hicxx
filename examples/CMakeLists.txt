# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

add_executable(hicxx_dump hicxx_dump.cpp)
target_link_libraries(
    hicxx_dump
    PRIVATE hicxx_project_options
    PUBLIC hicxx)
