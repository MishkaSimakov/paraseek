function(add_playground_target name filepath)
    add_executable(paraseek_playground_${name} ${filepath})
    target_link_libraries(paraseek_playground_${name} PUBLIC paraseek)
endfunction()

add_playground_target(solve_all solve_all.cpp)
add_playground_target(solve_one solve_one.cpp)
add_playground_target(brute_force_vs_tables brute_force_vs_tables.cpp)
add_playground_target(splitters splitters.cpp)
add_playground_target(reduce_all reduce_all.cpp)
add_playground_target(lsh lsh.cpp)
add_playground_target(scip scip.cpp)

add_playground_target(plot_entries_growth plot/plot_entries_growth.cpp)
add_playground_target(plot_tables_considered_count plot/plot_tables_considered_count.cpp)
add_playground_target(plot_brute_force_vs_tables_small plot/plot_brute_force_vs_tables_small.cpp)
add_playground_target(plot_brute_force_time plot/plot_brute_force_time.cpp)
add_playground_target(plot_tables_time plot/plot_tables_time.cpp)

