#include <chrono>
#include <print>

#include "matrix/NPY.h"
#include "problems/ProblemMatrix.h"
#include "problems/ProblemsNames.h"
#include "seekers/BruteForce.h"
#include "seekers/Tables.h"
#include "utils/Printing.h"

int main() {
  const auto problem_name = "lectsched-5-obj";
  std::println("{}", problem_name);

  auto matrix = get_problem_matrix(problem_name);
  std::println("  size: {} x {}", matrix.shape().first, matrix.shape().second);

  seekers::TablesParameters params{
      .groups_count = 4,
      .max_small_row_size = 8,
  };

  auto start = std::chrono::steady_clock::now();
  auto result = seekers::Tables(2, params).seek(matrix);
  auto end = std::chrono::steady_clock::now();

  std::println("  singular part: {}", result.singular.size());
  std::println("  bipartite part: {}", result.bipartite.size());
  std::println("  duration: {}", end - start);

  auto normalized = seekers::normalize_result(result);
  std::println("  normalized: {}", normalized.singular.size());
}
