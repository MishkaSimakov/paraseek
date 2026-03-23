add_library(paraseek)

target_sources(paraseek PUBLIC
        problems/ProblemsNames.cpp
        utils/Hamming.cpp
        utils/Printing.cpp
        seekers/Result.cpp
)
target_include_directories(paraseek PUBLIC .)
