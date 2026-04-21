#include <cassert>
#include <chrono>
#include <map>
#include <print>

#include "ExpressionDisjointSet.h"
#include "ReduceNonzeros.h"
#include "Reducer.h"
#include "problems/ProblemMatrix.h"
#include "problems/ProblemsNames.h"
#include "seekers/Tables.h"
#include "utils/Printing.h"

using namespace std::chrono_literals;

int main() {
  constexpr size_t max_diff = 2;

  const auto filename = std::format("problems_nz_reduction_{}.csv", max_diff);
  std::ofstream os(paths::log(filename));
  std::println(
      os, "problem_name,rows_count,cols_count,new_rows_count,new_cols_count");

  for (size_t problem_index = 0; problem_index < benchmark_set.size();
       ++problem_index) {
    std::println("{}/{}: {}", problem_index + 1, benchmark_set.size(),
                 benchmark_set[problem_index]);

    if (benchmark_set[problem_index] != "neos-3402454-bohle") {
      continue;
    }

    auto problem =
        get_problem("presolved_" + benchmark_set[problem_index], true);
    const auto [n, d] = problem.A.shape();

    std::println("  size: {} x {} (nz = {})", n, d, problem.A.nonzero_count());

    std::optional<Problem<double>> new_problem = std::nullopt;

    auto duration = timing::timeit([&] {
      new_problem =
          ReduceNonzeros<double, seekers::DoubleHasher>(6, 3).apply(problem);
    });

    std::println("  nz count: {} -> {}", problem.A.nonzero_count(),
                 new_problem->A.nonzero_count());
    std::println("  duration: {}", duration);

    std::println(os, "{},{},{},{},{}", benchmark_set[problem_index], n, d,
                 problem.A.shape().first, problem.A.shape().second);
    os.flush();
  }
}

