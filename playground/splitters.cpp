#include <chrono>
#include <print>

#include "problems/ProblemMatrix.h"
#include "problems/ProblemsNames.h"
#include "seekers/BruteForce.h"
#include "seekers/Tables.h"
#include "splitters/Evaluate.h"
#include "splitters/GreedySplitter.h"
#include "splitters/LocalSearchSplitter.h"
#include "splitters/RandomSplitter.h"
#include "utils/Printing.h"

int main() {
  std::ofstream os(paths::log("splitter_score_greedy_improved.csv"));
  std::println(os, "problem_name,score");

  for (size_t i = 0; i < problems_names.size(); ++i) {
    const auto problem_name = problems_names[i];
    if (problem_name != "neos-3402454-bohle") {
      continue;
    }

    std::println("{}/{}: {}", i + 1, problems_names.size(), problem_name);

    auto matrix = get_problem_matrix(problem_name, true);
    std::println("  size: {} x {} (nz = {})", matrix.shape().first,
                 matrix.shape().second, matrix.nonzero_count());

    auto result = splitters::GreedySplitter().split(matrix, 4);
    auto evaluation = splitters::evaluate(matrix, result, 2);

    if (!evaluation.is_valid) {
      throw std::runtime_error("Invalid split.");
    }

    std::println("  score = {}", evaluation.score);
    std::println(os, "{},{}", problem_name, evaluation.score);

    // form matrix out of zero rows
    auto [n, d] = matrix.shape();
    std::vector<bool> is_zero(n, true);

    for (size_t col : result[0]) {
      for (const auto [row, value] : matrix.get_column(col)) {
        is_zero[row] = false;
      }
    }

    std::unordered_map<size_t, size_t> zero_rows_enumeration;
    for (size_t row = 0; row < n; ++row) {
      if (is_zero[row]) {
        zero_rows_enumeration.emplace(row, zero_rows_enumeration.size());
      }
    }

    CSCMatrix<double> matrix2(zero_rows_enumeration.size());

    for (size_t group = 1; group < 4; ++group) {
      for (size_t col : result[group]) {
        matrix2.add_column();

        for (const auto [row, value] : matrix.get_column(col)) {
          if (is_zero[row]) {
            matrix2.push_to_last_column(zero_rows_enumeration.at(row), value);
          }
        }
      }
    }

    {
      auto result = splitters::GreedySplitter().split(matrix2, 4);
      auto evaluation = splitters::evaluate(matrix2, result, 2);

      if (!evaluation.is_valid) {
        throw std::runtime_error("Invalid split.");
      }

      std::println("  score = {}", evaluation.score);
      std::println(os, "{},{}", problem_name, evaluation.score);
    }
  }
}

// 3321330228380 (random)
// 2008473144618 (columns sorted asc)
// 1282282151750 (columns sorted desc)
// 645261033342
// 645261033342
// 480418685786
// 439912266088
// 439892174401
// 414038540365