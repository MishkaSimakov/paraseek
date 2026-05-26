cmake_minimum_required(VERSION 3.31)
project(paraseek)

if(NOT DEFINED PROBLEMS_PATH)
    message(FATAL_ERROR "You must specify path to MIPLIB problems in PROBLEMS_PATH variable.")
endif()

add_compile_definitions(PROBLEMS_PATH=${PROBLEMS_PATH})
add_compile_definitions(LOG_PATH="${CMAKE_CURRENT_LIST_DIR}/log")
set(CMAKE_CXX_STANDARD 23)

include(cmake/load_boost.cmake)
include(cmake/load_google_benchmark.cmake)
include(cmake/load_google_test.cmake)
include(cmake/load_highs.cmake)

add_subdirectory(src)
add_subdirectory(tests)
add_subdirectory(benchmark)
add_subdirectory(playground)
