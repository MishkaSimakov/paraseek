#include <chrono>
#include <print>

#include "../src/DisjointSet.h"
#include "problems/ProblemMatrix.h"
#include "problems/ProblemsNames.h"
#include "seekers/BruteForce.h"
#include "seekers/Tables.h"
#include "utils/Printing.h"

int main() {
  std::ofstream os(paths::log("runtimes.csv"));
  std::println(os, "problem_name,time,small_rows_time,big_rows_time");

  for (size_t i = 0; i < problems_names.size(); ++i) {
    const auto problem_name = problems_names[i];
    std::println("{}/{}: {}", i + 1, problems_names.size(), problem_name);

    auto matrix = get_problem_matrix(problem_name, true);
    const auto [n, d] = matrix.shape();

    std::println("  size: {} x {} (nz = {})", n, d, matrix.nonzero_count());

    seekers::TablesParameters params{
        .groups_count = 4,
        .max_small_row_size = 8,
        .log_entries_growth = false,
        .log_entries_per_row = false,
    };

    auto seeker = seekers::Tables(2, params);
    seekers::Result result;

    auto duration = timing::timeit([&] { result = seeker.seek(matrix); });
    std::println("  done!");

    auto stats = seeker.get_stats();

    std::println(os, "{},{},{},{}", problem_name, duration.count(),
                 stats.small_rows_duration.count(),
                 stats.big_rows_duration.count());

    // const auto transposed = matrix.get_transposed();
    // auto ds = DisjointSet(d);
    //
    // for (auto [i, j] : result.singular) {
    //   auto [distance, ratio] = similarity::hamming(transposed[i], transposed[j]);
    // }

    // auto normalized = seekers::normalize_result(result);
    // std::println("  size: {}", normalized.singular.size());
    // auto chosen = seekers::only_one(matrix.shape().first, result);
    //
    // size_t removed_count = 0;
    // for (const size_t i : chosen) {
    //   if (i != -1) {
    //     ++removed_count;
    //   }
    // }
    //
    // std::println("  removed count = {}", removed_count);
  }
}
