add_library(paraseek)

target_sources(paraseek PUBLIC
        similarity/Hamming.cpp
        similarity/Jaccard.cpp
        similarity/Union.cpp
        ProblemsNames.cpp
)
target_include_directories(paraseek PUBLIC .)

add_executable(main main.cpp)
target_link_libraries(main paraseek)