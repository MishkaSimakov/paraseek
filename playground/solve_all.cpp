#include <cassert>
#include <chrono>
#include <map>
#include <print>

#include "ExpressionDisjointSet.h"
#include "problems/ProblemMatrix.h"
#include "problems/ProblemsNames.h"
#include "seekers/BruteForce.h"
#include "seekers/Tables.h"
#include "utils/Printing.h"

Problem simple_problem() {
  CSCMatrix<double> A = {
      {1, 2, 0, 0},
      {2, 4, 1, 2},
      {0, 0, 3, 3},
  };

  Matrix<double> b = {{1}, {2}, {3}};

  std::vector<Bound<double>> bounds(3);

  return {A, b, bounds};
}

int main() {
  std::ofstream os(paths::log("runtimes.csv"));
  std::println(os, "problem_name,time,small_rows_time,big_rows_time");

  // for (size_t problem_index = 0; problem_index < problems_names.size();
  // ++problem_index) {
  auto [A, b, bounds] = simple_problem();
  const auto [n, d] = A.shape();

  std::println("  size: {} x {} (nz = {})", n, d, A.nonzero_count());

  seekers::TablesParameters params{
      .groups_count = 4,
      .max_small_row_size = 8,
  };

  auto seeker = seekers::Tables(2, params);
  auto result = seeker.seek(A);

  std::println("  done!");

  const auto transposed = A.get_transposed();
  auto ds = ExpressionDisjointSet(d);

  auto normalized_result = seekers::normalize_result(result);
  std::println("  size = {}", normalized_result.singular.size());

  for (auto [i, j] : normalized_result.singular) {
    auto ratio = similarity::hamming_leq(transposed[i], transposed[j], 2);

    assert(ratio.has_value());

    // subtract vectors
    double rhs_diff = b[i, 0] - *ratio * b[j, 0];

    SparseVector<double> diff;
    for (auto [i, x, y] : SparseZipRange{transposed[i], transposed[j]}) {
      const double value = x - *ratio * y;

      if (FieldTraits<double>::is_nonzero(value)) {
        diff.emplace_back(i, value);
      }
    }

    if (diff.size() == 2) {
      ds.unite(diff[0].first, diff[1].first,
               {-diff[1].second / diff[0].second, rhs_diff / diff[0].second});
      std::println("diff.size() = {}", diff.size());
    }
  }

  std::unordered_map<size_t, size_t> classes_enumeration;
  std::vector<std::vector<size_t>> classes(d);

  for (size_t col = 0; col < d; ++col) {
    auto [leader, expr] = ds.leader(col);

    auto [itr, _] =
        classes_enumeration.emplace(leader, classes_enumeration.size());
    classes[itr->second].push_back(col);
  }

  CSCMatrix<double> updated_A(n);
  auto updated_b = b;

  for (const auto [leader, index] : classes_enumeration) {
    std::map<size_t, double> column;

    for (const size_t col : classes[index]) {
      const auto [_, expr] = ds.leader(col);

      for (const auto [row, value] : A.get_column(col)) {
        column[row] += value * expr.first;
        updated_b[row, 0] -= value * expr.second;
      }
    }

    updated_A.add_column();
    for (const auto [row, value] : column) {
      updated_A.push_to_last_column(row, value);
    }
  }

  std::cout << updated_A << std::endl;
  std::cout << updated_b << std::endl;
}
