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
#include "seekers/TablesSimpleHashing.h"
#include "utils/Printing.h"

using namespace std::chrono_literals;

int main() {
  std::ofstream os(paths::log("runtimes.csv"));
  std::println(os, "problem_name,tables_time,bf_time");

  for (size_t problem_index = 0; problem_index < problems_names.size();
       ++problem_index) {
    std::println("{}/{}: {}", problem_index + 1, problems_names.size(),
                 problems_names[problem_index]);

    auto problem = get_problem_matrices(problems_names[problem_index], true);

    seekers::TablesParameters params{
        .groups_count = 4,
        .max_small_row_size = 8,
    };

    // solve using tables
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

      const size_t dn = problem.A.shape().first - new_problem.A.shape().first;
      const size_t dd = problem.A.shape().second - new_problem.A.shape().second;
      std::println("  dn = {}, dd = {}", dn, dd);

      if (dn == 0 && dd == 0) {
        break;
      }

      std::println("  ----------------------");
      problem = std::move(new_problem);
    }
  }
}
