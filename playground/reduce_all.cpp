#include <cassert>
#include <chrono>
#include <map>
#include <print>

#include "../src/variables/ExpressionDisjointSet.h"
#include "../src/variables/Reducer.h"
#include "problems/ProblemMatrix.h"
#include "problems/ProblemsNames.h"
#include "utils/Printing.h"
#include "variables/BruteForce.h"
#include "variables/Tables.h"

using namespace std::chrono_literals;

int main() {
  constexpr size_t max_diff = 2;

  const auto filename = std::format("problems_reduction_{}.csv", max_diff);
  std::ofstream os(paths::log(filename));
  std::println(os,
               "problem_name,rows_count,cols_count,new_rows_count,new_cols_"
               "count,reduction_iterations");

  for (size_t problem_index = 0; problem_index < benchmark_set.size();
       ++problem_index) {
    std::println("{}/{}: {}", problem_index + 1, benchmark_set.size(),
                 benchmark_set[problem_index]);

    auto problem =
        get_problem("presolved_" + benchmark_set[problem_index], true);
    const auto [init_n, init_d] = problem.A.shape();

    seekers::TablesParameters params{
        .groups_count = 4,
        .max_small_row_size = 8,
    };

    size_t iterations = 0;

    // solve using tables
    for (size_t i = 0; i < 5; ++i) {
      ++iterations;
      const auto [n, d] = problem.A.shape();

      std::println("  size: {} x {} (nz = {})", n, d,
                   problem.A.nonzero_count());

      auto seeker =
          seekers::Tables<double, seekers::DoubleHasher>(max_diff, params);
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

    std::println(os, "{},{},{},{},{},{}", benchmark_set[problem_index], init_n,
                 init_d, problem.A.shape().first, problem.A.shape().second,
                 iterations);
    os.flush();
  }
}
