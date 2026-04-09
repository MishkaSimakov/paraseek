#include <cassert>
#include <chrono>
#include <map>
#include <print>

#include "ExpressionDisjointSet.h"
#include "problems/ProblemMatrix.h"
#include "problems/ProblemsNames.h"
#include "seekers/BruteForce.h"
#include "seekers/Tables.h"
#include "utils/Printing.h"

int main() {
  std::ofstream os(paths::log("runtimes.csv"));
  std::println(os, "problem_name,time,small_rows_time,big_rows_time");

  for (size_t problem_index = 0; problem_index < problems_names.size();
       ++problem_index) {
    const auto problem =
        get_problem_matrices(problems_names[problem_index], true);
    const auto [n, d] = problem.A.shape();

    std::println("  size: {} x {} (nz = {})", n, d, problem.A.nonzero_count());

    seekers::TablesParameters params{
        .groups_count = 4,
        .max_small_row_size = 8,
    };

    auto seeker = seekers::Tables(2, params);
    auto result = seeker.seek(problem.A);

    std::println("  done!");

    auto normalized_result = seekers::normalize_result(result);
    std::println("  size = {}", normalized_result.singular.size());


  }
}
