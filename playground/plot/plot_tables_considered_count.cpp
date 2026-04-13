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
  std::ofstream os("tables_considered_count.csv");
  std::println(os, "problem,groups_count,considered_count");

  for (size_t problem_index = 0; problem_index < problems_names.size();
       ++problem_index) {
    const auto& problem_name = problems_names[problem_index];

    std::println("{}/{}: {}", problem_index + 1, problems_names.size(),
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

      std::println(os, "{},{},{}", problem_name, groups_count,
                   stats.pairs_considered);
    }
  }
}
