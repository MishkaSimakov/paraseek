#pragma once

#include "matrix/CSCMatrix.h"

namespace splitters {

template <typename Field>
class GreedySplitter {
 public:
  GreedySplitter() = default;

  std::vector<std::vector<size_t>> split(const CSCMatrix<Field>& matrix,
                                         size_t groups_count) {
    const auto [n, d] = matrix.shape();

    std::vector<std::vector<size_t>> groups(groups_count);
    std::vector<size_t> groups_total_zeros(groups_count, n);
    std::vector<std::vector<bool>> groups_zero_rows_mask(n);

    for (size_t row = 0; row < n; ++row) {
      groups_zero_rows_mask[row].resize(groups_count, false);
    }

    for (size_t col = 0; col < d; ++col) {
      std::vector<size_t> increment_per_group(groups_count, 0);

      for (const auto [row, value] : matrix.get_column(col)) {
        for (size_t i = 0; i < groups_count; ++i) {
          if (!groups_zero_rows_mask[row][i]) {
            ++increment_per_group[i];
          }
        }
      }

      const size_t shift = std::hash<size_t>()(col) % groups_count;
      size_t max_increment_group = shift;

      for (size_t i = 0; i < groups_count; ++i) {
        size_t group_id = (i + shift) % groups_count;

        if (increment_per_group[group_id] >
                increment_per_group[max_increment_group] ||
            increment_per_group[group_id] ==
                    increment_per_group[max_increment_group] &&
                groups_total_zeros[group_id] >
                    groups_total_zeros[max_increment_group]) {
          max_increment_group = group_id;
        }
      }

      groups[max_increment_group].push_back(col);

      for (const auto [row, value] : matrix.get_column(col)) {
        if (!groups_zero_rows_mask[row][max_increment_group]) {
          groups_zero_rows_mask[row][max_increment_group] = true;
          --groups_total_zeros[max_increment_group];
        }
      }
    }

    return groups;
  }
};

}  // namespace splitters
