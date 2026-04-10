#include <chrono>
#include <print>

#include "problems/ProblemMatrix.h"
#include "problems/ProblemsNames.h"
#include "seekers/BruteForce.h"
#include "seekers/Tables.h"
#include "utils/Printing.h"

int main() {
  // const auto problem_name = "gmu-35-40";
  // std::println("{}", problem_name);
  //
  // auto matrix = get_problem_matrices(problem_name, false).A;

  const CSCMatrix<double> matrix = {
      {0, 0, 0, 0, -0.853554, 0, 0, 0, -0.677405, 0},
      {0, 2.80367, 0.791415, -0.880154, 0, -1.65152, -3.25476, -0.325852, 0,
       -1.10218},
      {1.72653, -26.8861, -7.58936, 4.94937, 0, 15.8374, 31.2119, 3.12479, 0,
       10.5694},
      {0, 0, 0, 0, 5.0967, 0, -1.03632, 0, 5.37045, 0},
      {2.35089, 2.74656, 0.775293, -0.862225, 0, -1.61788, -3.18846, -0.319214,
       -4.92745, -1.07972},
  };
  std::println("  size: {} x {}", matrix.shape().first, matrix.shape().second);

  seekers::TablesParameters params{
      .groups_count = 4,
      .max_small_row_size = 4,
  };

  auto start = std::chrono::steady_clock::now();
  auto result = seekers::Tables(2, params).seek(matrix);
  auto end = std::chrono::steady_clock::now();

  auto bf_result = seekers::BruteForce(2).seek(matrix);

  printing::print_result(matrix, seekers::normalize_result(result).singular);

  std::println("  bf_size = {}", bf_result.size());
  std::println("  singular part: {}", seekers::normalize_result(result).singular.size());
  std::println("  bipartite part: {}", result.bipartite.size());
  std::println("  duration: {}", end - start);
}
