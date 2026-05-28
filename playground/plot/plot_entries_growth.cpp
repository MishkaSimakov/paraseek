#include <chrono>
#include <fstream>
#include <print>

#include "problems/ArchivedProblem.h"
#include "problems/ReplaceInequalities.h"
#include "utils/Hashers.h"
#include "vars/AllRows.h"

using namespace std::chrono_literals;

struct EntriesCountLogger {
  std::string filename;

  void operator()(std::vector<size_t> entries_count) const {
    std::ofstream os(paths::log(filename + ".csv"));

    std::println(os, "count");

    for (size_t value : entries_count) {
      std::println(os, "{}", value);
    }
  }
};

int main() {
  const std::vector<std::string> problems = {
      "app1-2",
      "square41",
  };

  for (size_t problem_index = 0; problem_index < problems.size();
       ++problem_index) {
    const auto& problem_name = problems[problem_index];

    std::println("{}/{}: {}", problem_index + 1, problems.size(), problem_name);
    auto problem = get_archived(problem_name);
    replace_inequalities(problem);

    const auto [n, d] = problem.A.shape();
    std::println("  size: {} x {} (nz = {})", n, d, problem.A.nonzero_count());

    // solve using entries reduction
    {
      seekers::AllRowsParameters params{
          .max_diff = 2,
          .groups_count = 4,
          .threshold = 8,
          .entries_reduction = true,
          .entries_count_logger =
              EntriesCountLogger{problem_name +
                                 "_with_reduction_entries_growth"},
      };

      seekers::AllRows<double, DoubleHasher>::seek(problem.A, params);
    }

    // solve without entries reduction
    {
      seekers::AllRowsParameters params{
          .max_diff = 2,
          .groups_count = 4,
          .threshold = 8,
          .entries_reduction = false,
          .entries_count_logger =
              EntriesCountLogger{problem_name + "_no_reduction_entries_growth"},
      };

      seekers::AllRows<double, DoubleHasher>::seek(problem.A, params);
    }
  }
}
