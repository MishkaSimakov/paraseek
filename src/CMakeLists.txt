add_library(paraseek)

target_sources(paraseek PUBLIC
        problems/ProblemsNames.cpp
        utils/Printing.cpp
        vars/Result.cpp
)
target_link_libraries(paraseek PUBLIC highs)
target_include_directories(paraseek PUBLIC .)