#pragma once

#include <algorithm>
#include <numeric>
#include <unordered_map>
#include <vector>

#include "matrix/CSCMatrix.h"

namespace seekers {

class Tables {
  struct Node {
    size_t prev;
    size_t next;
  };

  void seek_table(const CSCMatrix<double>& matrix,
                  const std::vector<size_t>& permutation) {
    auto [n, d] = matrix.shape();

    size_t classes_cnt = 1;
    std::vector<size_t> classes(n, 0);
    std::unordered_map<size_t, size_t> map;

    for (size_t i = 0; 2 * i < d; ++i) {
      map.clear();

      for (auto [row, _] : matrix.get_column(permutation[i])) {
        auto [itr, new_class] = map.emplace(classes[row], classes_cnt);

        if (new_class) {
          ++classes_cnt;
        }

        classes[row] = itr->second;
      }
    }

    std::vector<size_t> indices(n);
    std::iota(indices.begin(), indices.end(), 0);

    std::ranges::sort(indices, {}, [&](size_t i) { return classes[i]; });

    std::unordered_map<size_t, size_t> counts;
    for (size_t i = 0; i < n; ++i) {
      ++counts[classes[i]];
    }

    std::vector<std::pair<size_t, size_t>> sizes;
    for (auto [i, j] : counts) {
      sizes.emplace_back(i, j);
    }
    std::ranges::sort(sizes, {},
                      [&](auto p) { return -static_cast<int>(p.second); });

    for (size_t i = 0; i < 1; ++i) {
      std::println("  {} {}", sizes[i].second, sizes[i].first);
    }
  }

 public:
  // max_diff = 2
  Tables() = default;

  std::vector<std::pair<size_t, size_t>> seek(const CSCMatrix<double>& matrix) {
    // 4 permutations
    auto [n, d] = matrix.shape();
    std::vector<size_t> permutation(d);
    std::iota(permutation.begin(), permutation.end(), 0);

    seek_table(matrix, permutation);

    //
    std::vector<uint8_t> parts(n, 0);

    for (size_t col = 0; col < d; ++col) {
      for (auto [row, _] : matrix.get_column(col)) {
        // if there is a group without elements in this row, add this column
        // into this group
        if (parts[row] != 0b1111) {
          for (size_t i = 0; i < 4; ++i) {
            size_t group = (i + col) % 4;

            if ((parts[row] >> group & 1) == 0) {
              // add this column to group
              for (auto [j, _] : matrix.get_column(col)) {
                parts[j] |= 1 << group;
              }

              break;
            }
          }

          break;
        }
      }
    }

    std::vector<size_t> zero_cnt(16, 0);

    for (uint8_t part : parts) {
      for (size_t group1 = 0; group1 < 4; ++group1) {
        for (size_t group2 = 0; group2 < 4; ++group2) {
          if ((part >> group1 & 1) == 0 && (part >> group2 & 1) == 0) {
            ++zero_cnt[group1 * 4 + group2];
          }
        }
      }
    }

    for (size_t group1 = 0; group1 < 4; ++group1) {
      for (size_t group2 = 0; group2 < 4; ++group2) {
        std::cout << "  " << zero_cnt[group1 * 4 + group2];
      }
      std::cout << std::endl;
    }

    std::cout << std::endl;

    return {};
  }
};

}  // namespace seekers
