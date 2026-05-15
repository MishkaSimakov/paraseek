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
#include "utils/Accumulators.h"
#include "utils/Printing.h"

using namespace std::chrono_literals;

int main() {
  constexpr size_t max_diff = 2;
  constexpr size_t max_small_row_size = 8;

  const auto filename =
      std::format("tables_time_{}_{}.csv", max_diff, max_small_row_size);
  std::ofstream os(paths::log(filename));
  std::println(os,
               "problem,rows_count,cols_count,nonzeros_count,groups_squared,"
               "entries_considered,small_rows_count,big_rows_time,small_rows_"
               "time,total_time");

  const auto& names = benchmark_set;

  for (size_t problem_index = 0; problem_index < names.size();
       ++problem_index) {
    const auto& problem_name = names[problem_index];
    std::println("{}/{}: {}", problem_index + 1, names.size(), problem_name);

    try {
      auto problem = get_problem(problem_name, true);

      const auto [n, d] = problem.A.shape();
      std::println("  size: {} x {} (nz = {})", n, d,
                   problem.A.nonzero_count());

      seekers::TablesParameters params{
          .groups_count = 4,
          .max_small_row_size = max_small_row_size,
          .small_column_limit = 2,
          .entries_reduction = true,
          .log_prefix = "",
          .log_entries_growth = false,
      };

      ArithmeticMean<double> average_duration;
      ArithmeticMean<double> average_small_rows_duration;
      ArithmeticMean<double> average_big_rows_duration;

      seekers::Statistics stats;

      for (size_t i = 0; i < 10; ++i) {
        auto seeker =
            seekers::Tables<double, seekers::DoubleHasher>(max_diff, params);

        auto duration =
            timing::timeit([&] { auto result = seeker.seek(problem.A); });

        stats = seeker.get_stats();

        average_duration.record(duration.count());
        average_small_rows_duration.record(stats.small_rows_duration.count());
        average_big_rows_duration.record(stats.big_rows_duration.count());
      }

      size_t small_rows_count = 0;
      const auto transposed = problem.A.get_transposed();
      for (size_t i = 0; i < n; ++i) {
        if (transposed[i].size() <= max_small_row_size + max_diff) {
          ++small_rows_count;
        }
      }

      std::println(os, "{},{},{},{},{},{},{},{},{},{}", problem_name, n, d,
                   problem.A.nonzero_count(),
                   groups_squared(problem.A, max_diff),
                   stats.entries_considered, small_rows_count,
                   average_big_rows_duration.mean(),
                   average_small_rows_duration.mean(), average_duration.mean());
      os.flush();
    } catch (const std::exception& exception) {
      std::println("{}", exception.what());
    }
  }
}
