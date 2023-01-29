cmake_minimum_required(VERSION 3.14) # CMake version check
cmake_policy(VERSION 3.14...3.25)
project(strawC) # Create project "simple_example"
set(CMAKE_CXX_STANDARD 14) # Enable c++14 standard

option(STRAW_BUILD_EXAMPLES "Build straw example(s)" OFF)
option(STRAW_ENABLE_TESTING "Build unit tests" ON)
option(STRAW_USE_ZLIBNG "Use zlib-ng instead of zlib for decompression" OFF)

add_library(straw_project_options INTERFACE)
target_compile_features(straw_project_options INTERFACE cxx_std_${CMAKE_CXX_STANDARD})

add_subdirectory(src)

if(STRAW_BUILD_EXAMPLES)
    add_subdirectory(examples)
endif()

if(STRAW_ENABLE_TESTING)
    enable_testing()
    add_subdirectory(test)
endif()
