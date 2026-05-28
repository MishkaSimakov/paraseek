#include <chrono>
#include <fstream>
#include <print>

#include "problems/ArchivedProblem.h"
#include "problems/ProblemStatistics.h"
#include "problems/ProblemsNames.h"
#include "problems/ReplaceInequalities.h"
#include "utils/Paths.h"
#include "vars/AllRows.h"

using namespace std::chrono_literals;

int main() {
  std::ofstream os(paths::log("tables_considered_count.csv"));
  std::println(os, "problem,groups_squared,groups_count,considered_count");

  for (size_t problem_index = 0; problem_index < benchmark_set.size();
       ++problem_index) {
    const auto& problem_name = benchmark_set[problem_index];

    std::println("{}/{}: {}", problem_index + 1, benchmark_set.size(),
                 problem_name);
    auto problem = get_archived(problem_name);
    replace_inequalities(problem);

    const auto [n, d] = problem.A.shape();
    std::println("  size: {} x {} (nz = {})", n, d, problem.A.nonzero_count());

    for (size_t groups_count = 3; groups_count <= 6; ++groups_count) {
      seekers::AllRowsParameters tables_params{
          .max_diff = 2,
          .groups_count = groups_count,
          .threshold = 8,
          .entries_reduction = true,
      };

      auto [result, stats] = seekers::AllRows<double, DoubleHasher>::seek(
          problem.A, tables_params);

      std::println(os, "{},{},{},{}", problem_name, groups_squared(problem.A),
                   groups_count, stats.pairs_considered);
      os.flush();
    }
  }
}
