#include <chrono>
#include <print>

#include "../src/variables/Reducer.h"
#include "problems/ProblemMatrix.h"
#include "problems/ProblemsNames.h"
#include "splitters/Evaluate.h"
#include "splitters/RandomSplitter.h"
#include "utils/Printing.h"
#include "variables/BruteForce.h"
#include "variables/Tables.h"

int main() {
  const auto problem_name = "presolved_square47";
  std::println("{}", problem_name);

  const auto problem = get_problem(problem_name, true);

  const auto [n, d] = problem.A.shape();
  std::println("  size: {} x {} (nz = {})", n, d, problem.A.nonzero_count());

  for (size_t groups_count = 3; groups_count <= 10; ++groups_count) {
    std::println("groups_count = {}", groups_count);

    // with small
    {
      seekers::TablesParameters params{
          .groups_count = groups_count,
          .max_small_row_size = 0,
          .log_prefix = "with_small_",
      };

      auto seeker =
          seekers::Tables<double, seekers::DoubleHasher,
                          splitters::RandomSplitter<double>>(2, params);
      seeker.seek(problem.A);
    }

    // without small
    {
      seekers::TablesParameters params{
          .groups_count = groups_count,
          .max_small_row_size = 8,
          .log_prefix = "without_small_",
      };

      auto seeker =
          seekers::Tables<double, seekers::DoubleHasher,
                          splitters::RandomSplitter<double>>(2, params);
      seeker.seek(problem.A);
    }

    // without small + greedy
    {
      seekers::TablesParameters params{
        .groups_count = groups_count,
        .max_small_row_size = 8,
        .log_prefix = "without_small_greedy_",
    };

      auto seeker =
          seekers::Tables<double, seekers::DoubleHasher,
                          splitters::GreedySplitter<double>>(2, params);
      seeker.seek(problem.A);
    }
  }
}
