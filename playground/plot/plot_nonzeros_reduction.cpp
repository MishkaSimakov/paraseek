#include <cassert>
#include <chrono>
#include <map>
#include <print>

#include "../../src/nonzeros/ReduceNonzeros.h"
#include "../../src/variables/ExpressionDisjointSet.h"
#include "../../src/variables/Reducer.h"
#include "problems/ProblemMatrix.h"
#include "problems/ProblemsNames.h"
#include "utils/Printing.h"
#include "variables/Tables.h"

using namespace std::chrono_literals;

std::vector<std::string> get_presolved_benchmark_set() {
  std::vector<std::string> result;

  for (const auto& problem : benchmark_set) {
    result.push_back("presolved_" + problem);
  }

  return result;
}

int main() {
  constexpr size_t groups_count = 4;
  constexpr size_t selected_groups_count = 2;

  const auto filename = std::format("problems_nz_reduction_{}_{}.csv",
                                    groups_count, selected_groups_count);
  std::ofstream os(paths::log(filename));
  std::println(os,
               "problem_name,rows_count,old_nz_count,new_nz_count,duration");

  // const auto& problems = get_presolved_benchmark_set();
  std::vector<std::string> problems = {"model_presolved256961v2"};

  for (size_t problem_index = 0; problem_index < problems.size();
       ++problem_index) {
    std::println("{}/{}: {}", problem_index + 1, problems.size(),
                 problems[problem_index]);

    auto problem = get_problem(problems[problem_index], true);
    const auto [n, d] = problem.A.shape();

    std::println("  size: {} x {} (nz = {})", n, d, problem.A.nonzero_count());

    std::optional<Problem<double>> new_problem = std::nullopt;

    auto duration = timing::timeit([&] {
      new_problem = ReduceNonzeros<double, seekers::DoubleHasher>(
                        groups_count, selected_groups_count)
                        .apply(problem);
    });

    std::println("  nz count: {} -> {}", problem.A.nonzero_count(),
                 new_problem->A.nonzero_count());
    std::println("  duration: {}", duration);

    std::println(os, "{},{},{},{},{}", problems[problem_index],
                 problem.A.shape().first, problem.A.nonzero_count(),
                 new_problem->A.nonzero_count(), duration.count());
    os.flush();
  }
}
