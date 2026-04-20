find_package(SCIP REQUIRED)

add_library(paraseek)

target_sources(paraseek PUBLIC
        problems/ProblemsNames.cpp
        utils/Printing.cpp
        seekers/Result.cpp
)
target_include_directories(paraseek PUBLIC . ${SCIP_INCLUDE_DIRS})
target_link_libraries(paraseek PUBLIC ${SCIP_LIBRARIES})