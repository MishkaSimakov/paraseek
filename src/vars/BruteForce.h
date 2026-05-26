#pragma once

#include <algorithm>
#include <vector>

#include "Hamming.h"
#include "matrix/CSCMatrix.h"
#include "utils/Time.h"

namespace seekers {

struct BruteForceParameters {
  // If answer size exceeds size_limit pairs, then new pairs will not be added
  // to it. But algorithm won't stop working.
  std::optional<size_t> size_limit = std::nullopt;

  std::optional<timing::Deadline> deadline = std::nullopt;
};

struct BruteForceStatistics {
  size_t pairs_considered;
  timing::Duration duration;
};

template <typename Field>
class BruteForce {
  const size_t max_diff;
  const BruteForceParameters params;

  BruteForceStatistics statistics_;

  void add_to_result(std::vector<std::pair<size_t, size_t>>& result, size_t i,
                     size_t j) {
    if (params.size_limit && result.size() >= *params.size_limit) {
      return;
    }

    if (i > j) {
      std::swap(i, j);
    }

    result.emplace_back(i, j);
  }

 public:
  explicit BruteForce(size_t max_diff, BruteForceParameters params = {})
      : max_diff(max_diff), params(params) {}

  std::vector<std::pair<size_t, size_t>> seek(const CSCMatrix<Field>& matrix) {
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
    const auto transposed = matrix.get_transposed();

    for (size_t i = 0; i < n; ++i) {
      for (size_t j = i + 1;
           j < n && counts[j].second <= counts[i].second + max_diff; ++j) {
        ++statistics_.pairs_considered;

        if (similarity::hamming_leq(transposed[counts[i].first],
                                    transposed[counts[j].first], max_diff)) {
          add_to_result(result, counts[i].first, counts[j].first);
        }
      }

      if (params.deadline && params.deadline->is_over()) {
        break;
      }
    }

    return result;
  }

  BruteForceStatistics get_stats() const { return statistics_; }
};

}  // namespace seekers
