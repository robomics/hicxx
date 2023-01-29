# Add main.cpp file of project root directory as source file
set(SOURCE_FILES examples/main.cpp src/straw.cpp)
add_executable(straw_example example.cpp)
target_link_libraries(
    straw_example
    PRIVATE straw_project_options
    PUBLIC straw)
