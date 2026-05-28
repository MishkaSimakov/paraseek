function(add_playground_target name filepath)
    add_executable(paraseek_playground_${name} ${filepath})
    target_link_libraries(paraseek_playground_${name} PUBLIC paraseek)
endfunction()

add_playground_target(splitters splitters.cpp)
add_playground_target(reduce_all reduce_all.cpp)
add_playground_target(reduce_nonzeros reduce_nonzeros.cpp)

add_playground_target(highs_presolve highs.cpp)
target_link_libraries(paraseek_playground_highs_presolve PRIVATE highs)

add_playground_target(plot_entries_growth plot/plot_entries_growth.cpp)
add_playground_target(plot_tables_considered_count plot/plot_tables_considered_count.cpp)
add_playground_target(plot_brute_force_vs_tables_small plot/plot_brute_force_vs_tables_small.cpp)
add_playground_target(plot_brute_force_time plot/plot_brute_force_time.cpp)
add_playground_target(plot_tables_time plot/plot_tables_time.cpp)
add_playground_target(plot_problems_reduction plot/plot_problems_reduction.cpp)
add_playground_target(plot_nonzeros_reduction plot/plot_nonzeros_reduction.cpp)

