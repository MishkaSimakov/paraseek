#pragma once

#include <algorithm>
#include <vector>

#include "matrix/CSCMatrix.h"
#include "seekers/Statistics.h"
#include "utils/Hamming.h"

namespace seekers {

class BruteForce {
  const size_t max_diff;

  Statistics statistics_;

  static void add_to_result(std::vector<std::pair<size_t, size_t>>& result,
                            size_t i, size_t j) {
    if (i > j) {
      std::swap(i, j);
    }

    result.emplace_back(i, j);
  }

 public:
  explicit BruteForce(size_t max_diff) : max_diff(max_diff) {}

  std::vector<std::pair<size_t, size_t>> seek(const CSCMatrix<double>& matrix) {
    auto [n, d] = matrix.shape();

    std::vector<std::pair<size_t, size_t>> result;
    std::vector<std::pair<size_t, size_t>> counts(n);

    for (size_t i = 0; i < n; ++i) {
      counts[i] = {i, 0};
    }
    for (size_t i = 0; i < d; ++i) {
      for (auto [row, _] : matrix.get_column(i)) {
        ++counts[row].second;
      }
    }

    std::ranges::sort(counts, {}, [](auto p) { return p.second; });

    // precalculate transposed matrix
    std::vector<SparseVector<double>> rows(n);
    for (size_t col = 0; col < d; ++col) {
      for (auto [row, value] : matrix.get_column(col)) {
        rows[row].emplace_back(col, value);
      }
    }

    for (size_t i = 0; i < n; ++i) {
      for (size_t j = i + 1;
           j < n && counts[j].second <= counts[i].second + max_diff; ++j) {
        ++statistics_.pairs_considered;

        size_t diff =
            similarity::hamming(rows[counts[i].first], rows[counts[j].first]);

        if (diff <= max_diff) {
          add_to_result(result, counts[i].first, counts[j].first);
        }
      }
    }

    return result;
  }

  Statistics get_stats() const { return statistics_; }
};

}  // namespace seekers
