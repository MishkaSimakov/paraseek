#pragma once

#include "Greedy.h"
#include "matrix/CSCMatrix.h"
#include "utils/Random.h"

namespace splitters {

// Doesn't work. Experimental attempt to implement local search for splitting.
class LocalSearch {
  std::pair<size_t, size_t> find_group(
      const std::vector<std::vector<size_t>>& groups, size_t col) {
    for (size_t group_id = 0; group_id < groups.size(); ++group_id) {
      for (size_t i = 0; i < groups[group_id].size(); ++i) {
        if (col == groups[group_id][i]) {
          return {group_id, i};
        }
      }
    }

    std::unreachable();
  }

  std::optional<size_t> find_replace_column(
      const std::vector<size_t>& groups,
      const std::vector<std::vector<size_t>>& row_coverage,
      const std::vector<std::vector<std::pair<size_t, double>>>& transposed,
      const size_t row, const size_t unconnected_group,
      const CSCMatrix<double>& matrix) {
    for (const size_t col : transposed[row] | std::views::keys) {
      // try to change groups[col] to unconnected_group
      bool can_replace = true;

      for (const size_t col_row : matrix.get_column(col) | std::views::keys) {
        if (row_coverage[col_row][groups[col]] == 1) {
          can_replace = false;
          break;
        }
      }

      if (can_replace) {
        return col;
      }
    }

    return std::nullopt;
  }

 public:
  LocalSearch() = default;

  std::vector<std::vector<size_t>> split(const CSCMatrix<double>& matrix,
                                         size_t groups_count) {
    const auto [n, d] = matrix.shape();

    auto transposed = matrix.get_transposed();

    const auto greedy_solution = Greedy<double>().split(matrix, groups_count);

    std::vector<size_t> groups(d, 0);
    for (size_t i = 0; i < groups_count; ++i) {
      for (const size_t col : greedy_solution[i]) {
        groups[col] = i;
      }
    }

    // row_coverage[row][group] = # of nonzero columns in the row from the group
    std::vector<std::vector<size_t>> row_coverage(n);
    for (size_t row = 0; row < n; ++row) {
      row_coverage[row].resize(groups_count, 0);

      for (const size_t col : transposed[row] | std::views::keys) {
        ++row_coverage[row][groups[col]];
      }
    }

    for (size_t row = 0; row < n; ++row) {
      size_t unconnected_group = -1;

      for (size_t group = 0; group < groups_count; ++group) {
        if (row_coverage[row][group] == 0) {
          unconnected_group = group;
        }
      }

      if (unconnected_group == -1) {
        continue;
      }

      // one-level group change algorithm
      auto replace_col = find_replace_column(groups, row_coverage, transposed,
                                             row, unconnected_group, matrix);

      if (!replace_col) {
        continue;
      }

      for (const size_t col_row :
           matrix.get_column(*replace_col) | std::views::keys) {
        --row_coverage[col_row][groups[*replace_col]];
        ++row_coverage[col_row][unconnected_group];
      }

      groups[*replace_col] = unconnected_group;
    }

    std::vector<std::vector<size_t>> result(groups_count);

    for (size_t col = 0; col < d; ++col) {
      result[groups[col]].push_back(col);
    }

    return result;
  }
};

}  // namespace splitters
