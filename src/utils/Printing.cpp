#include "Printing.h"

#include <algorithm>
#include <unordered_map>

#include "utils/Hashers.h"

void printing::print_result(
    const CSCMatrix<double>& matrix,
    const std::vector<std::pair<size_t, size_t>>& result,
    std::optional<size_t> max_count) {
  size_t cnt = 0;

  for (auto [i, j] : result) {
    if (max_count && cnt > *max_count) {
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

void printing::print_diff(const CSCMatrix<double>& matrix,
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

void printing::print_small_rows_cnt(const CSCMatrix<double>& matrix,
                                    size_t count) {
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

void printing::print_most_frequent(const std::vector<size_t>& classes,
                                   std::optional<size_t> max_count) {
  std::unordered_map<size_t, size_t> counts;
  for (size_t c : classes) {
    ++counts[c];
  }

  std::vector<std::pair<size_t, size_t>> sizes(counts.begin(), counts.end());
  std::ranges::sort(sizes, {},
                    [&](auto p) { return -static_cast<int>(p.second); });

  max_count = max_count ? std::min(*max_count, sizes.size()) : sizes.size();

  for (size_t i = 0; i < *max_count; ++i) {
    std::println("  {} {}", sizes[i].second, sizes[i].first);
  }
}
