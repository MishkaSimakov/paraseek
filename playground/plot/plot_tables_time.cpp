#include <chrono>
#include <fstream>
#include <print>

#include "problems/ArchivedProblem.h"
#include "problems/ProblemStatistics.h"
#include "problems/ProblemsNames.h"
#include "problems/ReplaceInequalities.h"
#include "utils/Accumulators.h"
#include "vars/AllRows.h"
#include "vars/Reducer.h"

using namespace std::chrono_literals;

int main() {
  constexpr size_t max_diff = 2;
  constexpr size_t max_small_row_size = 8;

  const auto filename =
      std::format("tables_time_{}_{}.csv", max_diff, max_small_row_size);
  std::ofstream os(paths::log(filename));

  std::cout << paths::log(filename) << std::endl;
  if (!os) {
    throw std::runtime_error("Failed to open output file.");
  }

  std::println(os,
               "problem,rows_count,cols_count,nonzeros_count,groups_squared,"
               "entries_considered,small_rows_count,big_rows_time,small_rows_"
               "time");

  const auto& names = benchmark_set;

  for (size_t problem_index = 0; problem_index < names.size();
       ++problem_index) {
    const auto& problem_name = names[problem_index];
    std::println("{}/{}: {}", problem_index + 1, names.size(), problem_name);

    try {
      auto problem = get_archived(problem_name);
      replace_inequalities(problem);

      const auto [n, d] = problem.A.shape();
      std::println("  size: {} x {} (nz = {})", n, d,
                   problem.A.nonzero_count());

      seekers::AllRowsParameters params{
          .max_diff = max_diff,
          .groups_count = 4,
          .threshold = max_small_row_size,
          .small_column_limit = 2,
          .entries_reduction = true,
      };

      ArithmeticMean<double> average_small_rows_duration;
      ArithmeticMean<double> average_big_rows_duration;

      size_t entries_considered;

      for (size_t i = 0; i < 10; ++i) {
        auto [_, stats] =
            seekers::AllRows<double, DoubleHasher>::seek(problem.A, params);

        entries_considered = stats.entries_considered;

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

      std::println(os, "{},{},{},{},{},{},{},{},{}", problem_name, n, d,
                   problem.A.nonzero_count(),
                   groups_squared(problem.A, max_diff), entries_considered,
                   small_rows_count, average_big_rows_duration.mean(),
                   average_small_rows_duration.mean());
      os.flush();
    } catch (const std::exception& exception) {
      std::println("{}", exception.what());
    }
  }
}
