#include <chrono>
#include <print>

#include "problems/ProblemMatrix.h"
#include "problems/ProblemsNames.h"
#include "seekers/BruteForce.h"
#include "seekers/Tables.h"
#include "utils/Printing.h"

int main() {
  std::ofstream os("brute_force_vs_tables.csv");
  std::println(os, "problem_name,brute_force_time,tables_time");

  for (size_t i = 0; i < std::min(150uz, problems_names.size()); ++i) {
    const auto problem_name = problems_names[i];
    std::println("{}/{}: {}", i + 1, problems_names.size(), problem_name);

    auto matrix = get_problem_matrix(problem_name);
    std::println("  size: {} x {}", matrix.shape().first,
                 matrix.shape().second);

    size_t brute_force_time;
    {
      auto start = std::chrono::steady_clock::now();
      auto result = seekers::BruteForce(2).seek(matrix);
      auto end = std::chrono::steady_clock::now();

      brute_force_time = (end - start).count();
    }

    size_t tables_time;
    {
      auto start = std::chrono::steady_clock::now();
      auto [singular, bipartite] = seekers::Tables(2).seek(matrix);
      auto end = std::chrono::steady_clock::now();

      tables_time = (end - start).count();
    }

    std::println(os, "{},{},{}", problem_name, brute_force_time, tables_time);
  }
}
