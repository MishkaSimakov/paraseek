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

const std::vector<std::string> bad_problems_names = {
    "nursesched-medium04",
    "neos-4647032-veleka",
    "n3seq24",
    "bab3",
    "hypothyroid-k1",
    "neos-4647027-thurso",
    "neos-4647030-tutaki",
    "8div-n59k12",
    "neos-2991472-kalu",
    "neos-5157194-moruya",
};

int main() {
  constexpr size_t max_diff = 2;

  const auto filename = std::format("tables_time_{}_new.csv", max_diff);
  std::ofstream os(paths::log(filename));
  std::println(os,
               "problem,rows_count,cols_count,nonzeros_count,groups_squared,"
               "big_rows_time,small_rows_time,total_time");

  const auto& names = collection_set;

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
          .groups_count = 6,
          .max_small_row_size = 8,
          .small_column_limit = 2,
          .entries_reduction = true,
          .log_prefix = "",
          .log_entries_growth = false,
      };

      auto seeker =
          seekers::Tables<double, seekers::DoubleHasher>(max_diff, params);

      auto duration = timing::timeit([&] {
        auto result = seeker.seek(problem.A);

        std::println("  singular: {}", result.singular.size());
        std::println("  bipartite: {}", result.bipartite.size());
      });

      auto stats = seeker.get_stats();

      std::println(os, "{},{},{},{},{},{},{},{}", problem_name, n, d,
                   problem.A.nonzero_count(),
                   groups_squared(problem.A, max_diff),
                   stats.big_rows_duration.count(),
                   stats.small_rows_duration.count(), duration.count());
      os.flush();
    } catch (const std::exception& exception) {
      std::println("{}", exception.what());
    }
  }
}

// 21'676'706'333 (4)
// 68'023'620'666 (3)
// 5'326'576'875

// 89'621'815'792
// 190'866'383'125