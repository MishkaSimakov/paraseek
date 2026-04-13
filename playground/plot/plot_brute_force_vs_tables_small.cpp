#include <cassert>
#include <chrono>
#include <map>
#include <print>
#include <random>

#include "../../tests/helpers/RandomProblem.h"
#include "seekers/BruteForce.h"
#include "seekers/Tables.h"

using namespace std::chrono_literals;

int main() {
  constexpr size_t problems_count = 10'000;

  std::ofstream os(paths::log("brute_force_vs_tables_small.csv"));
  std::println(os, "rows_count,cols_count,nonzeros_count,tables_time,bf_time");

  std::uniform_int_distribution<size_t> size_distribution(100, 2000);
  std::default_random_engine random;

  for (size_t problem_index = 0; problem_index < problems_count;
       ++problem_index) {
    const size_t n = size_distribution(random);
    const size_t d = size_distribution(random);

    const auto problem =
        generate_random_problem<double>(n, d, random, true, 0.004);

    std::println("{}/{}:", problem_index + 1, problems_count);
    std::println("  size: {} x {} (nz = {})", n, d, problem.A.nonzero_count());

    // solve using tables
    seekers::TablesParameters tables_params{
        .groups_count = 4,
        .max_small_row_size = 8,
        .entries_reduction = true,
        .log_prefix = "",
        .log_entries_growth = false,
    };

    auto tables_time = timing::timeit([&] {
      seekers::Tables<double, seekers::DoubleHasher>(2, tables_params)
          .seek(problem.A);
    });

    // solve using brute force
    auto bf_time =
        timing::timeit([&] { seekers::BruteForce<double>(2).seek(problem.A); });

    std::println(os, "{},{},{},{},{}", n, d, problem.A.nonzero_count(),
                 tables_time.count(), bf_time.count());
    os.flush();
  }
}
