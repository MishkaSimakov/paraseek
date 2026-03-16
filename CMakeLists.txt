cmake_minimum_required(VERSION 3.31)
project(paraseek)

add_compile_definitions(PROBLEMS_PATH="/Users/mihailsimakov/Documents/Programs/misc/miplib")
set(CMAKE_CXX_STANDARD 23)

add_subdirectory(src)
add_subdirectory(tests)
add_subdirectory(benchmark)
