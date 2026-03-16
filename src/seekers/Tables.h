#pragma once

#include <algorithm>
#include <numeric>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "matrix/CSCMatrix.h"
#include "utils/Hashers.h"

namespace seekers {

class Tables {
  struct Node {
    size_t prev;
    size_t next;
  };

  static void print_top_classes(const std::vector<size_t>& classes,
                                size_t count) {
    std::unordered_map<size_t, size_t> counts;
    for (size_t c : classes) {
      ++counts[c];
    }

    std::vector<std::pair<size_t, size_t>> sizes(counts.begin(), counts.end());
    std::ranges::sort(sizes, {},
                      [&](auto p) { return -static_cast<int>(p.second); });

    for (size_t i = 0; i < count; ++i) {
      std::println("  {} {}", sizes[i].second, sizes[i].first);
    }
  }

  void seek_table(const CSCMatrix<double>& matrix,
                  const std::vector<SparseVector<double>>& transposed,
                  const std::vector<std::vector<size_t>>& groups,
                  const std::vector<bool>& groups_mask,
                  std::unordered_set<std::pair<size_t, size_t>>& result) {
    auto [n, d] = matrix.shape();

    size_t classes_cnt = 1;
    std::vector<size_t> classes(n, 0);
    std::vector<double> coef(n, 0);

    std::unordered_map<std::pair<size_t, int>, size_t> map;

    for (size_t group_id = 0; group_id < groups.size(); ++group_id) {
      if (!groups_mask[group_id]) {
        continue;
      }

      for (size_t col : groups[group_id]) {
        map.clear();

        for (auto [row, value] : matrix.get_column(col)) {
          if (coef[row] == 0) {
            coef[row] = value;
          }

          int normalized = static_cast<int>(value / coef[row]);
          auto [itr, new_class] =
              map.emplace(std::pair{classes[row], normalized}, classes_cnt);

          if (new_class) {
            ++classes_cnt;
          }

          classes[row] = itr->second;
        }
      }
    }

    std::vector<size_t> indices(n);
    std::iota(indices.begin(), indices.end(), 0);

    std::ranges::sort(indices, {}, [&](size_t i) { return classes[i]; });

    // for each class remaining positions must be checked
    print_top_classes(classes, 5);

    for (size_t i = 0; i < n; ++i) {
      if (classes[indices[i]] == 0) {
        continue;
      }

      for (size_t j = i + 1;
           j < n && classes[indices[j]] == classes[indices[i]]; ++j) {
        if (similarity::hamming(transposed[indices[i]],
                                transposed[indices[j]]) <= 2) {
          result.emplace(indices[i], indices[j]);
        }
      }
    }
  }

 public:
  // max_diff = 2
  Tables() = default;

  std::vector<std::pair<size_t, size_t>> seek(const CSCMatrix<double>& matrix) {
    constexpr size_t groups_count = 4;
    constexpr size_t selected_groups_count = groups_count - 2;

    auto [n, d] = matrix.shape();
    std::vector<size_t> permutation(d);

    std::vector<SparseVector<double>> transposed(n);
    for (size_t col = 0; col < d; ++col) {
      for (auto [row, value] : matrix.get_column(col)) {
        transposed[row].emplace_back(col, value);
      }
    }

    // find groups
    // for each rows stores a bitmask:
    // bit i is zero for row j if j-th row is empty in i-th group
    std::vector<uint32_t> zero_rows(n, 0);
    std::vector<std::vector<size_t>> groups(groups_count);

    static_assert(groups_count <= std::numeric_limits<uint32_t>::digits);

    for (size_t col = 0; col < d; ++col) {
      bool found_any = false;

      for (auto [row, _] : matrix.get_column(col)) {
        // if there is a group without elements in this row, add this column
        // into this group
        if (zero_rows[row] != (1 << groups_count) - 1) {
          for (size_t i = 0; i < groups_count; ++i) {
            size_t group = (i + col) % groups_count;

            if ((zero_rows[row] >> group & 1) == 0) {
              // add this column to group
              found_any = true;
              groups[group].push_back(col);

              for (auto [j, _] : matrix.get_column(col)) {
                zero_rows[j] |= 1 << group;
              }

              break;
            }
          }

          break;
        }
      }

      if (!found_any) {
        groups[col % groups_count].push_back(col);
      }
    }

    std::vector<bool> mask(groups_count, false);
    std::fill_n(mask.begin(), selected_groups_count, true);

    std::unordered_set<std::pair<size_t, size_t>> result;

    size_t i = 0;

    do {
      seek_table(matrix, transposed, groups, mask, result);
    } while (std::ranges::prev_permutation(mask).found);

    return {result.begin(), result.end()};
  }
};

}  // namespace seekers
