#include <cassert>
#include <chrono>
#include <map>
#include <print>

#include "ExpressionDisjointSet.h"
#include "Reducer.h"
#include "problems/ProblemMatrix.h"
#include "seekers/Tables.h"

using namespace std::chrono_literals;

int main() {
  const std::vector<std::string> problems = {
      "app1-2",
      "square41",
  };

  for (size_t problem_index = 0; problem_index < problems.size();
       ++problem_index) {
    const auto& problem_name = problems[problem_index];

    std::println("{}/{}: {}", problem_index + 1, problems.size(), problem_name);
    auto problem = get_problem(problem_name, true);

    const auto [n, d] = problem.A.shape();
    std::println("  size: {} x {} (nz = {})", n, d, problem.A.nonzero_count());

    // solve using entries reduction
    {
      seekers::TablesParameters tables_params{
          .groups_count = 4,
          .max_small_row_size = 8,
          .entries_reduction = true,
          .log_prefix = problem_name + "_with_reduction_",
          .log_entries_growth = true,
      };

      seekers::Tables<double, seekers::DoubleHasher>(2, tables_params)
          .seek(problem.A);
    }

    // solve without entries reduction
    {
      seekers::TablesParameters tables_params{
          .groups_count = 4,
          .max_small_row_size = 8,
          .entries_reduction = false,
          .log_prefix = problem_name + "_no_reduction_",
          .log_entries_growth = true,
      };

      seekers::Tables<double, seekers::DoubleHasher>(2, tables_params)
          .seek(problem.A);
    }
  }
}
