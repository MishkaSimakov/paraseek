enable_testing()

add_executable(tests)

target_sources(tests PUBLIC
        unit/HammingTests.cpp
        unit/HighsTests.cpp
        unit/MPSReaderTests.cpp

        unit/nonzeros/ReduceNonzerosTests.cpp

        unit/vars/ReducerTests.cpp
        unit/vars/SmallRowsTests.cpp
        unit/vars/BigRowsTests.cpp
        unit/vars/AllRowsTests.cpp
)
target_link_libraries(tests PRIVATE GTest::gtest_main paraseek highs Boost::boost)
target_include_directories(tests PRIVATE .)
set_target_properties(tests PROPERTIES CMAKE_CXX_FLAGS "-fsanitize=address")

include(GoogleTest)
gtest_discover_tests(tests)