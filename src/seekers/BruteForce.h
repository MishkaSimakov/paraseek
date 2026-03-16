#pragma once

#include <algorithm>
#include <vector>

#include "Checker.h"
#include "matrix/CSCMatrix.h"
#include "similarity/Jaccard.h"

namespace seekers {

struct BruteForceJaccard {
  double min_similarity;
};

struct BruteForceHamming {
  size_t max_diff;
};

using BruteForceMode = std::variant<BruteForceJaccard, BruteForceHamming>;

class BruteForce {
  const BruteForceMode mode;

  static std::vector<std::pair<size_t, size_t>> seek_hamming(
      const CSCMatrix<double>& matrix, size_t max_diff) {
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
        size_t diff =
            similarity::hamming(rows[counts[i].first], rows[counts[j].first]);

        if (diff <= max_diff) {
          result.emplace_back(counts[i].first, counts[j].first);
        }
      }
    }

    return result;
  }

  static std::vector<std::pair<size_t, size_t>> seek_jaccard(
      const CSCMatrix<double>& matrix, double min_similarity) {
    auto [n, d] = matrix.shape();

    std::vector<std::pair<size_t, size_t>> result;

    // precalculate transposed matrix
    std::vector<SparseVector<double>> rows(n);
    for (size_t col = 0; col < d; ++col) {
      for (auto [row, value] : matrix.get_column(col)) {
        rows[row].emplace_back(col, value);
      }
    }

    for (size_t i = 0; i < n; ++i) {
      for (size_t j = i + 1; j < n; ++j) {
        double sim = similarity::jaccard(rows[i], rows[j]);

        if (sim >= min_similarity) {
          result.emplace_back(i, j);
        }
      }
    }

    return result;
  }

 public:
  explicit BruteForce(BruteForceMode mode) : mode(mode) {}

  std::vector<std::pair<size_t, size_t>> seek(const CSCMatrix<double>& matrix) {
    if (std::holds_alternative<BruteForceHamming>(mode)) {
      return seek_hamming(matrix, std::get<BruteForceHamming>(mode).max_diff);
    }

    return seek_jaccard(matrix,
                        std::get<BruteForceJaccard>(mode).min_similarity);
  }
};

}  // namespace seekers
