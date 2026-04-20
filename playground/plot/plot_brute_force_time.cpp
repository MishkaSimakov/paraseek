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

  const auto filename = std::format("brute_force_time_{}.csv", max_diff);
  std::ofstream os(paths::log(filename));
  std::println(
      os,
      "problem,rows_count,cols_count,nonzeros_count,groups_squared,bf_time");

  for (size_t problem_index = 0; problem_index < benchmark_set.size();
       ++problem_index) {
    const auto& problem_name = benchmark_set[problem_index];

    std::println("{}/{}: {}", problem_index + 1, benchmark_set.size(),
                 problem_name);
    auto problem = get_problem(problem_name, true);

    const auto [n, d] = problem.A.shape();
    std::println("  size: {} x {} (nz = {})", n, d, problem.A.nonzero_count());

    seekers::BruteForceParameters params{
        .size_limit = 1'000'000,
        .deadline = timing::Deadline::after(120s),
    };

    auto bf_time = timing::timeit(
        [&] { seekers::BruteForce<double>(max_diff, params).seek(problem.A); });

    std::println(os, "{},{},{},{},{},{}", problem_name, n, d,
                 problem.A.nonzero_count(), groups_squared(problem.A, max_diff),
                 bf_time.count());
    os.flush();
  }
}
