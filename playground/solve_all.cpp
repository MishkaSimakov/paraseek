#include <cassert>
#include <chrono>
#include <map>
#include <print>

#include "ExpressionDisjointSet.h"
#include "Reducer.h"
#include "problems/ProblemMatrix.h"
#include "problems/ProblemsNames.h"
#include "seekers/BruteForce.h"
#include "seekers/Tables.h"
#include "utils/Printing.h"

using namespace std::chrono_literals;

int main() {
  std::ofstream bf_os(paths::log("bf_runtime.csv"));
  std::println(bf_os, "problem_name,time");

  std::ofstream tables_os(paths::log("tables_runtime.csv"));
  std::println(tables_os, "problem_name,time");

  for (size_t problem_index = 0; problem_index < benchmark_set.size();
       ++problem_index) {
    const auto& problem_name = benchmark_set[problem_index];

    std::println("{}/{}: {}", problem_index + 1, benchmark_set.size(),
                 problem_name);

    auto problem = get_problem(problem_name, true);

    const auto [n, d] = problem.A.shape();

    std::println("  size: {} x {} (nz = {})", n, d, problem.A.nonzero_count());

    // solve using tables
    seekers::TablesParameters tables_params{
        .groups_count = 4,
        .max_small_row_size = 8,
    };

    const auto tables_duration = timing::timeit([&] {
      seekers::Tables<double, seekers::DoubleHasher>(2, tables_params)
          .seek(problem.A);
    });

    // solve using brute force
    seekers::BruteForceParameters bf_params{
        .size_limit = 1'000'000, .deadline = timing::Deadline::after(600s)};

    const auto bf_duration = timing::timeit(
        [&] { seekers::BruteForce<double>(2, bf_params).seek(problem.A); });

    std::println(bf_os, "{},{}", problem_name, bf_duration.count());
    bf_os.flush();

    std::println(tables_os, "{},{}", problem_name, tables_duration.count());
    tables_os.flush();
  }
}
