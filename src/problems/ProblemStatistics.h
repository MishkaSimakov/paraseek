#pragma once

#include <algorithm>
#include <unordered_map>

#include "matrix/CSCMatrix.h"

template <typename Field>
size_t groups_squared(const CSCMatrix<Field>& matrix) {
  const auto [n, d] = matrix.shape();
  const auto rows_sizes = matrix.get_rows_sizes();

  std::unordered_map<size_t, size_t> sizes_counts;
  for (size_t size : rows_sizes) {
    ++sizes_counts[size];
  }

  size_t result = 0;
  for (const size_t count : sizes_counts | std::views::values) {
    result += count * count;
  }

  return result;
}

template <typename Field>
size_t groups_squared(const CSCMatrix<Field>& matrix, size_t max_diff) {
  const auto [n, d] = matrix.shape();
  auto rows_sizes = matrix.get_rows_sizes();

  std::ranges::sort(rows_sizes);

  size_t current_size = -1;
  size_t current_group_end = 0;

  size_t result = 0;

  for (size_t row = 0; row < n; ++row) {
    if (rows_sizes[row] != current_size) {
      // update current_group_end
      while (current_group_end < n &&
             rows_sizes[current_group_end] <= rows_sizes[row] + max_diff) {
        ++current_group_end;
      }

      current_size = rows_sizes[row];
    }

    result += current_group_end - row - 1;
  }

  return result;
}
