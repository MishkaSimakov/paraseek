#include <chrono>
#include <print>

#include "problems/ProblemMatrix.h"
#include "problems/ProblemsNames.h"
#include "seekers/BruteForce.h"
#include "seekers/Tables.h"
#include "utils/Printing.h"

int main() {
  std::ofstream os("output.csv");
  std::println(os, "problem_name,time");

  for (size_t i = 0; i < problems_names.size(); ++i) {
    const auto problem_name = problems_names[i];
    std::println("{}/{}: {}", i + 1, problems_names.size(), problem_name);

    auto matrix = get_problem_matrix(problem_name);
    std::println("  size: {} x {}", matrix.shape().first,
                 matrix.shape().second);

    seekers::TablesParameters params{
        .groups_count = 4,
        .max_small_row_size = 10,
    };

    auto start = std::chrono::steady_clock::now();
    auto [singular, bipartite] = seekers::Tables(2, params).seek(matrix);
    auto end = std::chrono::steady_clock::now();

    std::println(os, "{},{}", problem_name, (end - start).count());

    std::println("  singular part: {}", singular.size());
    std::println("  bipartite part: {}", bipartite.size());
    std::println("  duration: {}", end - start);
  }
}
