function(add_playground_target name filepath)
    add_executable(paraseek_playground_${name} ${filepath})
    target_link_libraries(paraseek_playground_${name} PUBLIC paraseek)
endfunction()

add_playground_target(solve_all solve_all.cpp)
add_playground_target(solve_one solve_one.cpp)
