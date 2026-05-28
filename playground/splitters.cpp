#include <chrono>
#include <fstream>
#include <print>

#include "problems/ArchivedProblem.h"
#include "problems/ProblemsNames.h"
#include "problems/ReplaceInequalities.h"
#include "splitters/Evaluate.h"
#include "splitters/Greedy.h"
#include "utils/Paths.h"

int main() {
  std::ofstream os(paths::log("splitter_score_greedy_improved.csv"));
  std::println(os, "problem_name,score");

  for (size_t i = 0; i < benchmark_set.size(); ++i) {
    const auto problem_name = benchmark_set[i];

    std::println("{}/{}: {}", i + 1, benchmark_set.size(), problem_name);

    auto problem = get_archived(problem_name);
    replace_inequalities(problem);

    const auto& matrix = problem.A;

    std::println("  size: {} x {} (nz = {})", matrix.shape().first,
                 matrix.shape().second, matrix.nonzero_count());

    auto result = splitters::Greedy<double>().split(matrix, 4);
    auto evaluation = splitters::evaluate(matrix, result, 2);

    if (!evaluation.is_valid) {
      throw std::runtime_error("Invalid split.");
    }

    std::println("  score = {}", evaluation.score);
    std::println(os, "{},{}", problem_name, evaluation.score);
  }
}
