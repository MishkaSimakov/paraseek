#include <cassert>
#include <chrono>
#include <map>
#include <print>

#include "nonzeros/ReduceNonzeros.h"
#include "problems/ArchivedProblem.h"
#include "problems/ProblemStatistics.h"
#include "problems/ProblemsNames.h"
#include "problems/ReplaceInequalities.h"
#include "utils/CsvWriter.h"
#include "utils/Printing.h"
#include "utils/Time.h"

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

  const auto filename =
      std::format("{}_problems_nz_reduction_{}_{}.csv", paths::time_prefix(),
                  groups_count, selected_groups_count);
  std::ofstream os(paths::log(filename));
  CsvWriter csv{os};

  csv << std::tuple{"problem_name",
                    "rows_count",
                    "cols_count",
                    "old_nz_count",
                    "new_nz_count",
                    "old_equalities_count",
                    "new_equalities_count",
                    "old_trivial_equalities_count",
                    "new_trivial_equalities_count",
                    "duration_ns"};

  const auto& problems = get_presolved_benchmark_set();

  for (size_t problem_index = 0; problem_index < problems.size();
       ++problem_index) {
    std::println("{}/{}: {}", problem_index + 1, problems.size(),
                 problems[problem_index]);

    std::optional<Problem<double>> problem;

    try {
      problem = get_archived(problems[problem_index]);
    } catch (const std::exception& err) {
      std::cerr << err.what() << std::endl;
      continue;
    }

    const auto [n, d] = problem->A.shape();

    std::println("  size: {} x {} (nz = {})", n, d, problem->A.nonzero_count());

    std::optional<Problem<double>> new_problem = std::nullopt;

    auto duration = timing::timeit([&] {
      new_problem = ReduceNonzeros<double, DoubleHasher>(groups_count,
                                                         selected_groups_count)
                        .apply(*problem);
    });

    std::println("  nz count: {} -> {}", problem->A.nonzero_count(),
                 new_problem->A.nonzero_count());
    std::println("  duration: {}", duration);

    csv << std::tuple{problems[problem_index],
                      problem->A.shape().first,
                      problem->A.shape().second,
                      problem->A.nonzero_count(),
                      new_problem->A.nonzero_count(),
                      equalities_count(*problem),
                      equalities_count(*new_problem),
                      trivial_equalities_count(*problem),
                      trivial_equalities_count(*new_problem),
                      duration.count()};
    csv.flush();
  }
}
