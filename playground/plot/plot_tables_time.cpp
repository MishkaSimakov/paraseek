#include <cassert>
#include <chrono>
#include <map>
#include <print>

#include "ExpressionDisjointSet.h"
#include "Reducer.h"
#include "problems/ProblemMatrix.h"
#include "problems/ProblemStatistics.h"
#include "problems/ProblemsNames.h"
#include "seekers/BruteForce.h"
#include "seekers/Tables.h"
#include "utils/Printing.h"

using namespace std::chrono_literals;

int main() {
  constexpr size_t max_diff = 2;

  const auto filename = std::format("tables_time_{}.csv", max_diff);
  std::ofstream os(paths::log(filename));
  std::println(os,
               "problem,rows_count,cols_count,nonzeros_count,groups_squared,"
               "tables_time");

  for (size_t problem_index = 0; problem_index < problems_names.size();
       ++problem_index) {
    const auto& problem_name = problems_names[problem_index];

    std::println("{}/{}: {}", problem_index + 1, problems_names.size(),
                 problem_name);
    auto problem = get_problem(problem_name, true);

    const auto [n, d] = problem.A.shape();
    std::println("  size: {} x {} (nz = {})", n, d, problem.A.nonzero_count());

    seekers::TablesParameters params{
        .groups_count = 4,
        .max_small_row_size = 8,
        .entries_reduction = true,
        .log_prefix = "",
        .log_entries_growth = false,
    };

    auto tables_time = timing::timeit([&] {
      seekers::Tables<double, seekers::DoubleHasher>(max_diff, params)
          .seek(problem.A);
    });

    std::println(os, "{},{},{},{},{},{}", problem_name, n, d,
                 problem.A.nonzero_count(), groups_squared(problem.A, max_diff),
                 tables_time.count());
    os.flush();
  }
}
