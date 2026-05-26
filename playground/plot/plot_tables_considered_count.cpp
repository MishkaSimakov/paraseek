#include <cassert>
#include <chrono>
#include <map>
#include <print>

#include "../../src/variables/ExpressionDisjointSet.h"
#include "../../src/variables/Reducer.h"
#include "problems/ProblemMatrix.h"
#include "problems/ProblemStatistics.h"
#include "problems/ProblemsNames.h"
#include "utils/Printing.h"
#include "variables/BruteForce.h"
#include "variables/Tables.h"

using namespace std::chrono_literals;

int main() {
  std::ofstream os(paths::log("tables_considered_count.csv"));
  std::println(os, "problem,groups_squared,groups_count,considered_count");

  for (size_t problem_index = 0; problem_index < benchmark_set.size();
       ++problem_index) {
    const auto& problem_name = benchmark_set[problem_index];

    std::println("{}/{}: {}", problem_index + 1, benchmark_set.size(),
                 problem_name);
    auto problem = get_problem(problem_name, true);

    const auto [n, d] = problem.A.shape();
    std::println("  size: {} x {} (nz = {})", n, d, problem.A.nonzero_count());

    for (size_t groups_count = 3; groups_count <= 6; ++groups_count) {
      seekers::TablesParameters tables_params{
          .groups_count = groups_count,
          .max_small_row_size = 8,
          .entries_reduction = true,
          .log_prefix = "",
          .log_entries_growth = false,
      };

      auto seeker =
          seekers::Tables<double, seekers::DoubleHasher>(2, tables_params);

      seeker.seek(problem.A);

      auto stats = seeker.get_stats();

      std::println(os, "{},{},{},{}", problem_name, groups_squared(problem.A),
                   groups_count, stats.pairs_considered);
      os.flush();
    }
  }
}
