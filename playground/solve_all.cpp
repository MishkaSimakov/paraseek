#include <chrono>
#include <print>

#include "problems/ProblemMatrix.h"
#include "problems/ProblemsNames.h"
#include "seekers/BruteForce.h"
#include "seekers/Tables.h"
#include "utils/Printing.h"

int main() {
  std::ofstream os("output.csv");
  std::println(os, "problem_name,time,small_rows_time,big_rows_time");

  for (size_t i = 0; i < problems_names.size(); ++i) {
    const auto problem_name = problems_names[i];
    std::println("{}/{}: {}", i + 1, problems_names.size(), problem_name);

    auto matrix = get_problem_matrix(problem_name);
    std::println("  size: {} x {}", matrix.shape().first,
                 matrix.shape().second);

    seekers::TablesParameters params{
        .groups_count = 4,
        .max_small_row_size = 8,
    };

    auto seeker = seekers::Tables(2, params);
    seekers::Result result;

    auto duration = timing::timeit([&] { result = seeker.seek(matrix); });

    auto stats = seeker.get_stats();

    std::println(os, "{},{},{},{}", problem_name, duration.count(),
                 stats.small_rows_duration.count(),
                 stats.big_rows_duration.count());
  }
}
