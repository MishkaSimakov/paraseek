include(load_google_test.cmake)

enable_testing()

add_executable(tests)

target_sources(tests PUBLIC
        SimilarityTests.cpp
        FastHammingTests.cpp
        TablesTests.cpp
)
target_link_libraries(tests GTest::gtest_main paraseek)
set_target_properties(tests PROPERTIES CMAKE_CXX_FLAGS "-fsanitize=address")

include(GoogleTest)
gtest_discover_tests(tests)