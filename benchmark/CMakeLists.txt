include(load_google_benchmark.cmake)

add_executable(paraseek_benchmark benchmark.cpp)

target_link_libraries(paraseek_benchmark benchmark::benchmark paraseek)