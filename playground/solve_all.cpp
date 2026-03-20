#include <print>

#include "problems/ProblemMatrix.h"
#include "problems/ProblemsNames.h"
#include "seekers/BruteForce.h"
#include "seekers/Tables.h"

int main() {
  for (size_t i = 0; i < problems_names.size(); ++i) {
    const auto problem_name = problems_names[i];
    std::println("{}/{}: {}", i + 1, problems_names.size(), problem_name);

    auto matrix = get_problem_matrix(problem_name);
    std::println("  size: {} x {}", matrix.shape().first,
                 matrix.shape().second);

    auto [singular, bipartite] = seekers::Tables(2).seek(matrix);

    std::println("  singular part: {}", singular.size());
    std::println("  bipartite part: {}", bipartite.size());
  }
}
