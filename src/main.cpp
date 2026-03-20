#include <numeric>
#include <print>

#include "ProblemMatrix.h"
#include "ProblemsNames.h"
#include "seekers/BruteForce.h"
#include "seekers/Tables.h"
#include "similarity/ZipRows.h"

void print_result(const CSCMatrix<double>& matrix,
                  const std::vector<std::pair<size_t, size_t>>& result,
                  std::optional<size_t> max_cnt = std::nullopt) {
  size_t cnt = 0;

  for (auto [i, j] : result) {
    if (max_cnt && cnt > *max_cnt) {
      return;
    }

    if (i > j) {
      std::swap(i, j);
    }

    std::println("({}, {}):", i, j);

    auto xs = matrix.get_row(i);
    auto ys = matrix.get_row(j);

    for (auto [i, x, y] : SparseZipRange{xs, ys}) {
      std::print("{:20.5f}", x);
    }
    std::print("\n");

    for (auto [i, x, y] : SparseZipRange{xs, ys}) {
      std::print("{:20.5f}", y);
    }
    std::print("\n\n");

    ++cnt;
  }
}

void print_diff(const CSCMatrix<double>& matrix,
                const std::vector<std::pair<size_t, size_t>>& left,
                const std::vector<std::pair<size_t, size_t>>& right) {
  std::unordered_set left_set(left.begin(), left.end());
  for (auto [a, b] : right) {
    left_set.erase({a, b});
    left_set.erase({b, a});
  }

  std::println("+ left:");
  std::vector left_diff(left_set.begin(), left_set.end());
  print_result(matrix, left_diff);
}

void print_small_rows_cnt(const CSCMatrix<double>& matrix, size_t count) {
  auto [n, d] = matrix.shape();

  std::vector<size_t> counts(n, 0);
  for (size_t col = 0; col < d; ++col) {
    for (auto [row, value] : matrix.get_column(col)) {
      ++counts[row];
    }
  }

  std::vector<size_t> counts_counts(d + 1, 0);
  for (size_t row = 0; row < n; ++row) {
    ++counts_counts[counts[row]];
  }

  for (size_t i = 1; i <= count && i < counts_counts.size(); ++i) {
    std::println("  #{} - {}", i, counts_counts[i]);
  }
}

int main() {
  // std::ofstream os("output.csv");
  // std::println(os, "problem,found,considered");

  for (size_t i = 0; i < problems_names.size(); ++i) {
    const auto problem_name = problems_names[i];
    // const auto problem_name = "gmu-35-40";
    std::println("problem: {}", problem_name);

    auto matrix = get_problem_matrix(problem_name);

    std::println("  starting seeking");
    auto tbls_result = seekers::Tables(2).seek(matrix);
    std::println("  finished seeking");

    // std::println("----------   tables   ----------");
    // print_result(matrix, {tbls_normalized.begin(), tbls_normalized.end()});
    //
    // std::println("---------- brute force ----------");
    // print_result(matrix, bf_result);
  }

  // std::println("----------   tables   ----------");
  // print_result(matrix, result);
  //
  // std::println("---------- brute force ----------");
  // print_result(matrix, bf_result);

  // auto sh_result = seekers::SimHash().seek(matrix, bf_result);
  // std::println("found: {} pairs", sh_result.size());

  // for (const auto& problem_name : problems) {
  //   auto matrix = get_problem_matrix(problem_name);
  //
  //   auto [n, d] = matrix.shape();
  //
  //   std::vector<size_t> permutation(d);
  //   std::iota(permutation.begin(), permutation.end(), 0);
  //   std::ranges::shuffle(permutation, std::default_random_engine{});
  //
  //   std::vector<size_t> counts(n);
  //   for (size_t i = 0; 2 * i < d; ++i) {
  //     for (auto [row, _] : matrix.get_column(permutation[i])) {
  //       ++counts[row];
  //     }
  //   }
  //
  //   size_t zero_count = 0;
  //   for (size_t cnt : counts) {
  //     if (cnt == 0) {
  //       ++zero_count;
  //     }
  //   }
  //
  //   std::println("{} - {}/{} = {}", problem_name, zero_count, n,
  //   static_cast<double>(zero_count) / static_cast<double>(n));
  // }
}
