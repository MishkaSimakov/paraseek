#include <cassert>
#include <chrono>
#include <map>
#include <print>

#include "ExpressionDisjointSet.h"
#include "Reducer.h"
#include "problems/ProblemMatrix.h"
#include "problems/ProblemsNames.h"
#include "seekers/BruteForce.h"
#include "seekers/Tables.h"
#include "seekers/TablesSimpleHashing.h"
#include "utils/Printing.h"

using namespace std::chrono_literals;

int main() {
  std::ofstream os(paths::log("runtimes.csv"));
  std::println(os, "problem_name,tables_time,bf_time");

  for (size_t problem_index = 0; problem_index < problems_names.size();
       ++problem_index) {
    std::println("{}/{}: {}", problem_index + 1, problems_names.size(),
                 problems_names[problem_index]);

    auto problem = get_problem_matrices(problems_names[problem_index], true);

    const auto [n, d] = problem.A.shape();

    std::println("  size: {} x {} (nz = {})", n, d, problem.A.nonzero_count());

    // solve using tables
    seekers::TablesParameters tables_params{
        .groups_count = 4,
        .max_small_row_size = 8,
    };

    const auto tables_duration = timing::timeit([&] {
      seekers::Tables(2, tables_params).seek(problem.A);
    });

    // solve using brute force
    // seekers::BruteForceParameters bf_params{
    // .size_limit = 1'000'000, .deadline = timing::Deadline::after(60s)};

    // const auto bf_duration = timing::timeit(
    // [&] { seekers::BruteForce(2, bf_params).seek(problem.A); });

    // std::println(os, "{},{},{}", problems_names[problem_index],
    // tables_duration, bf_duration);

    std::println(os, "{},{}", problems_names[problem_index], tables_duration);
    os.flush();
  }
}
