#pragma once

#include <algorithm>
#include <cmath>
#include <numeric>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "matrix/CSCMatrix.h"
#include "seekers/Statistics.h"
#include "similarity/Hamming.h"
#include "similarity/ZipRows.h"
#include "utils/Hashers.h"

namespace seekers {

// https://stackoverflow.com/questions/20511347/a-good-hash-function-for-a-vector/72073933
class RowHasher {
  uint32_t hash_;

 public:
  explicit RowHasher(uint32_t seed) : hash_(seed) {}

  template <typename T>
  RowHasher& operator<<(T value) {
    uint32_t x = std::hash<T>()(value);

    x = ((x >> 16) ^ x) * 0x45d9f3b;
    x = ((x >> 16) ^ x) * 0x45d9f3b;
    x = (x >> 16) ^ x;
    hash_ ^= x + 0x9e3779b9 + (hash_ << 6) + (hash_ >> 2);

    return *this;
  }

  uint32_t get_hash() const { return hash_; }
};

class Tables {
  struct Block {
    size_t class_id;
    double front;
  };

  constexpr static size_t max_diff = 2;
  const size_t groups_count;

  Statistics statistics_;

  std::vector<SparseVector<double>> transposed_;
  std::unordered_set<std::pair<size_t, size_t>> result_;

  // Prepares double for hashing
  static int normalize_double(double value) { return std::round(value * 1000); }

  // Returns vector of (class hash, class size) sorted by size (desc)
  static std::vector<std::pair<size_t, size_t>> get_classes_sizes(
      const std::vector<size_t>& classes) {
    std::unordered_map<size_t, size_t> counts;
    for (size_t c : classes) {
      ++counts[c];
    }

    std::vector<std::pair<size_t, size_t>> sizes(counts.begin(), counts.end());
    std::ranges::sort(sizes, {},
                      [&](auto p) { return -static_cast<int>(p.second); });

    return sizes;
  }

  static void print_top_classes(const std::vector<size_t>& classes,
                                size_t count) {
    auto sizes = get_classes_sizes(classes);

    for (size_t i = 0; i < count && i < sizes.size(); ++i) {
      std::println("  {} {}", sizes[i].second, sizes[i].first);
    }
  }

  // For each row a list of blocks is returned
  static std::vector<std::vector<Block>> get_blocks(
      const CSCMatrix<double>& matrix,
      const std::vector<std::vector<size_t>>& groups) {
    auto [n, d] = matrix.shape();

    std::vector<std::vector<Block>> blocks(n);
    for (size_t i = 0; i < n; ++i) {
      blocks[i].resize(groups.size(), Block{0, 0});
    }

    std::unordered_map<std::pair<size_t, int>, size_t> map;

    for (size_t group_id = 0; group_id < groups.size(); ++group_id) {
      size_t classes_cnt = 1;

      for (const size_t col : groups[group_id]) {
        map.clear();

        for (auto [row, value] : matrix.get_column(col)) {
          if (blocks[row][group_id].front == 0) {
            blocks[row][group_id].front = value;
          }

          int normalized =
              normalize_double(value / blocks[row][group_id].front);

          auto [itr, new_class] =
              map.emplace(std::pair{blocks[row][group_id].class_id, normalized},
                          classes_cnt);

          if (new_class) {
            ++classes_cnt;
          }

          blocks[row][group_id].class_id = itr->second;
        }
      }
    }

    return blocks;
  }

  void consider_pair(size_t i, size_t j) {
    ++statistics_.pairs_considered;

    if (similarity::hamming(transposed_[i], transposed_[j]) <= max_diff) {
      result_.emplace(i, j);
    } else {
      // std::println("({}, {}):", i, j);
      //
      // auto xs = transposed_[i];
      // auto ys = transposed_[j];
      //
      // for (auto [i, x, y] : SparseZipRange{xs, ys}) {
      //   std::print("{:8} ", x);
      // }
      // std::print("\n");
      //
      // for (auto [i, x, y] : SparseZipRange{xs, ys}) {
      //   std::print("{:8} ", y);
      // }
      // std::print("\n\n");
    }
  }

  void consider_pair_fast(size_t i, size_t j,
                          const std::vector<double>& front) {
    ++statistics_.pairs_considered;

    if (front[i] != 0) {
      double ratio = front[i] / front[j];
      if (similarity::hamming_fixed_ratio(transposed_[i], transposed_[j],
                                          ratio) <= max_diff) {
        result_.emplace(i, j);
      }
    } else {
      if (similarity::fast_hamming(transposed_[i], transposed_[j]) <=
          max_diff) {
        result_.emplace(i, j);
      }
    }
  }

  void seek_table(const CSCMatrix<double>& matrix,
                  const std::vector<bool>& groups_mask,
                  const std::vector<std::vector<Block>>& blocks,
                  const std::vector<std::vector<size_t>>& groups) {
    auto [n, d] = matrix.shape();

    std::vector<size_t> merged_classes(n);
    std::vector<double> front(n, 0);
    std::unordered_map<size_t, size_t> classes_sizes;

    for (size_t i = 0; i < n; ++i) {
      if (transposed_[i].size() <= max_diff * 2) {
        // small rows are processed separately
        merged_classes[i] = 0;
      } else {
        // TODO: correct seed? RowHasher?
        StreamHasher hasher;

        for (size_t group_id = 0; group_id < groups_mask.size(); ++group_id) {
          if (groups_mask[group_id]) {
            if (front[i] == 0) {
              front[i] = blocks[i][group_id].front;
            }

            hasher << blocks[i][group_id].class_id;
          }

          merged_classes[i] = hasher.get_hash();
        }
      }

      ++classes_sizes[merged_classes[i]];
    }

    size_t prev = 0;
    for (auto& [key, value] : classes_sizes) {
      value += prev;
      prev = value;
    }

    std::vector<size_t> indices(n);
    for (size_t i = 0; i < n; ++i) {
      indices[--classes_sizes[merged_classes[i]]] = i;
    }

    // for each class remaining positions must be checked
    // print_top_classes(merged_classes, 5);

    for (size_t i = 0; i < n; ++i) {
      if (merged_classes[indices[i]] == 0) {
        continue;
      }

      // process class using full pairwise search
      for (size_t j = i + 1;
           j < n && merged_classes[indices[j]] == merged_classes[indices[i]];
           ++j) {
        consider_pair_fast(indices[i], indices[j], front);
      }
    }
  }

  void process_row(
      size_t row_index, const SparseVector<double>& row, size_t removed,
      std::vector<std::unordered_multimap<size_t, size_t>>& hashes) {
    RowHasher hasher(row.size());

    for (const size_t index : row | std::views::keys) {
      hasher << index;
    }

    size_t hash = hasher.get_hash();
    auto [begin, end] = hashes[max_diff - removed].equal_range(hash);

    for (auto itr = begin; itr != end; ++itr) {
      consider_pair(row_index, itr->second);
    }

    for (size_t i = removed; i <= max_diff; ++i) {
      hashes[i].emplace(hash, row_index);
    }
  }

  void traverse_row_combinations(
      size_t row, size_t i, SparseVector<double>& current, size_t removed,
      std::vector<std::unordered_multimap<size_t, size_t>>& hashes) {
    if (i == transposed_[row].size()) {
      process_row(row, current, removed, hashes);
      return;
    }

    current.push_back(transposed_[row][i]);
    traverse_row_combinations(row, i + 1, current, removed, hashes);
    current.pop_back();

    if (removed < max_diff) {
      traverse_row_combinations(row, i + 1, current, removed + 1, hashes);
    }
  }

  void seek_small() {
    // process small rows using hashing
    size_t n = transposed_.size();

    size_t small_rows_cnt = 0;
    for (size_t i = 0; i < n; ++i) {
      if (transposed_[i].size() > 3 * max_diff) {
        ++small_rows_cnt;
      }
    }

    // stores pairs of (hash, row index)
    std::vector<std::unordered_multimap<size_t, size_t>> hashes(max_diff + 1);
    for (size_t i = 0; i <= max_diff; ++i) {
      hashes[i].reserve(small_rows_cnt * 3);
    }

    SparseVector<double> buffer;

    for (size_t i = 0; i < n; ++i) {
      if (transposed_[i].size() > 3 * max_diff) {
        continue;
      }

      traverse_row_combinations(i, 0, buffer, 0, hashes);
    }
  }

 public:
  // max_diff = 2
  explicit Tables(size_t groups_count) : groups_count(groups_count) {}

  std::vector<std::pair<size_t, size_t>> seek(const CSCMatrix<double>& matrix) {
    const size_t selected_groups_count = groups_count - max_diff;

    auto [n, d] = matrix.shape();
    std::vector<size_t> permutation(d);

    transposed_.resize(n);
    for (size_t col = 0; col < d; ++col) {
      for (auto [row, value] : matrix.get_column(col)) {
        transposed_[row].emplace_back(col, value);
      }
    }

    // find groups
    // for each rows stores a bitmask:
    // bit i is zero for row j if j-th row is empty in i-th group
    std::vector<uint32_t> zero_rows(n, 0);
    std::vector<std::vector<size_t>> groups(groups_count);

    for (size_t col = 0; col < d; ++col) {
      bool found_any = false;

      // for (auto [row, _] : matrix.get_column(col)) {
      //   // if there is a group without elements in this row, add this column
      //   // into this group
      //   if (zero_rows[row] != (1 << groups_count) - 1) {
      //     for (size_t i = 0; i < groups_count; ++i) {
      //       size_t group = (i + col) % groups_count;
      //
      //       if ((zero_rows[row] >> group & 1) == 0) {
      //         // add this column to group
      //         found_any = true;
      //         groups[group].push_back(col);
      //
      //         for (auto [j, _] : matrix.get_column(col)) {
      //           zero_rows[j] |= 1 << group;
      //         }
      //
      //         break;
      //       }
      //     }
      //
      //     break;
      //   }
      // }

      if (!found_any) {
        groups[col % groups_count].push_back(col);
      }
    }

    auto blocks = get_blocks(matrix, groups);

    std::vector<bool> mask(groups_count, false);
    std::fill_n(mask.begin(), selected_groups_count, true);

    // do {
    //   seek_table(matrix, mask, blocks, groups);
    // } while (std::ranges::prev_permutation(mask).found);

    // now process small rows
    seek_small();

    return {result_.begin(), result_.end()};
  }

  Statistics get_stats() const { return statistics_; }
};

}  // namespace seekers

// 1992169627
