enable_testing()

add_executable(tests)

target_sources(tests PUBLIC
        unit/HammingTests.cpp
        unit/TablesTests.cpp
        unit/ReducerTests.cpp
        unit/HighsTests.cpp
)
target_link_libraries(tests PRIVATE GTest::gtest_main paraseek highs Boost::boost)
target_include_directories(tests PRIVATE .)
set_target_properties(tests PROPERTIES CMAKE_CXX_FLAGS "-fsanitize=address")

include(GoogleTest)
gtest_discover_tests(tests)