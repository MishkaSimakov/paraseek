#include <chrono>
#include <print>

#include "Reducer.h"
#include "problems/ProblemMatrix.h"
#include "problems/ProblemsNames.h"
#include "seekers/BruteForce.h"
#include "seekers/Tables.h"
#include "seekers/TablesSimpleHashing.h"
#include "utils/Printing.h"

int main() {
  const auto problem_name = "app1-2";
  std::println("{}", problem_name);

  auto problem = get_problem_matrices(problem_name, true);

  seekers::TablesParameters params{
      .groups_count = 4,
      .max_small_row_size = 8,
  };

  for (size_t i = 0; i < 5; ++i) {
    const auto [n, d] = problem.A.shape();

    std::println("  size: {} x {} (nz = {})", n, d,
                 problem.A.nonzero_count());

    auto seeker = seekers::Tables(2, params);
    auto result = seeker.seek(problem.A);

    std::println("  done!");

    auto for_reducer = result_for_reducer(result);
    std::println("  size = {}", for_reducer.size());

    auto [new_problem, mapping] = Reducer().apply(problem, for_reducer);

    std::println("  {} x {} -> {} x {}", problem.A.shape().first,
                 problem.A.shape().second, new_problem.A.shape().first,
                 new_problem.A.shape().second);

    std::println("  ----------------------");
    problem = std::move(new_problem);
  }
}
