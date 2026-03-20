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

std::unordered_set<std::pair<size_t, size_t>> normalize_tables_result(
    const std::vector<std::pair<size_t, size_t>>& singular,
    const std::vector<std::pair<std::vector<size_t>, std::vector<size_t>>>&
        bipartite) {
  std::unordered_set<std::pair<size_t, size_t>> result;

  for (auto [i, j] : singular) {
    if (i > j) {
      result.emplace(j, i);
    } else {
      result.emplace(i, j);
    }
  }

  for (const auto& [left, right] : bipartite) {
    for (size_t i : left) {
      for (size_t j : right) {
        if (i == j) {
          continue;
        }

        if (i > j) {
          result.emplace(j, i);
        } else {
          result.emplace(i, j);
        }
      }
    }
  }

  return result;
}

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

struct TablesParameters {
  size_t groups_count;
  size_t max_small_row_size;
};

class Tables {
  struct Block {
    size_t class_id;
    double front;
  };

  const size_t max_diff;
  const size_t groups_count;
  const size_t max_small_row_size;

  Statistics statistics_;

  std::vector<SparseVector<double>> transposed_;

  std::vector<std::pair<size_t, size_t>> result_singular_;
  std::vector<std::pair<std::vector<size_t>, std::vector<size_t>>>
      result_bipartite_;

  // Prepares double for hashing
  static int64_t normalize_double(double value) {
    return std::round(value * 1e10);
  }

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

    std::unordered_map<std::pair<size_t, int64_t>, size_t> map;

    for (size_t group_id = 0; group_id < groups.size(); ++group_id) {
      size_t classes_cnt = 1;

      for (const size_t col : groups[group_id]) {
        map.clear();

        for (auto [row, value] : matrix.get_column(col)) {
          if (blocks[row][group_id].front == 0) {
            blocks[row][group_id].front = value;
          }

          auto normalized =
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

  void add_to_answer(size_t i, size_t j) {
    if (i > j) {
      std::swap(i, j);
    }

    result_singular_.emplace_back(i, j);
  }

  void consider_pair(size_t i, size_t j) {
    ++statistics_.pairs_considered;

    if (similarity::hamming(transposed_[i], transposed_[j]) <= max_diff) {
      add_to_answer(i, j);
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
        add_to_answer(i, j);
      }
    } else {
      if (similarity::fast_hamming(transposed_[i], transposed_[j]) <=
          max_diff) {
        add_to_answer(i, j);
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
      if (transposed_[i].size() <= max_small_row_size) {
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
      std::vector<std::vector<std::pair<size_t, size_t>>>& hashes) {
    RowHasher hasher(0);

    for (const auto [index, value] : row) {
      hasher << index << normalize_double(value / row[0].second);
    }

    // if (row_index == 25265 || row_index == 80466) {
    //   std::println("row index: {}, hash: {}", row_index, hasher.get_hash());
    //
    //   for (const auto [index, value] : row) {
    //     std::print("  ({}, {:.6f})", index, value);
    //   }
    //   std::print("\n");
    //   for (const auto [index, value] : row) {
    //     std::print("  ({}, {:.6f})", index, value / row[0].second);
    //   }
    //   std::print("\n\n");
    // }

    hashes[row.size()].emplace_back(hasher.get_hash(), row_index);
  }

  void traverse_row_combinations(
      size_t row, size_t i, SparseVector<double>& current, size_t removed,
      std::vector<std::vector<std::pair<size_t, size_t>>>& hashes) {
    if (i == transposed_[row].size()) {
      process_row(row, current, removed, hashes);
      return;
    }

    current.push_back(transposed_[row][i]);
    traverse_row_combinations(row, i + 1, current, removed, hashes);
    current.pop_back();

    if (removed < max_diff && removed + 2 < transposed_[row].size()) {
      traverse_row_combinations(row, i + 1, current, removed + 1, hashes);
    }
  }

  void compare_hashes(const std::vector<std::pair<size_t, size_t>>& left,
                      const std::vector<std::pair<size_t, size_t>>& right) {
    size_t right_hash = 0;
    size_t right_begin = 0;
    size_t right_end = 0;

    size_t left_begin = 0;
    size_t left_end = 0;

    while (left_begin < left.size()) {
      size_t left_hash = left[left_begin].first;

      while (left_end < left.size() && left[left_end].first == left_hash) {
        ++left_end;
      }

      if (left_hash != right_hash || left_begin == 0) {
        // update right_begin and right_end
        right_begin = right_end;

        while (right_begin < right.size() &&
               right[right_begin].first < left_hash) {
          ++right_begin;
        }

        right_end = right_begin;
        while (right_end < right.size() &&
               right[right_end].first <= left_hash) {
          ++right_end;
        }

        right_hash = left_hash;
      }

      if (right_begin == right_end) {
        left_begin = left_end;
        continue;
      }

      std::vector<size_t> left_elements(left_end - left_begin);
      std::vector<size_t> right_elements(right_end - right_begin);

      for (size_t i = left_begin; i < left_end; ++i) {
        left_elements[i - left_begin] = left[i].second;
      }

      for (size_t i = right_begin; i < right_end; ++i) {
        right_elements[i - right_begin] = right[i].second;
      }

      result_bipartite_.emplace_back(std::move(left_elements),
                                     std::move(right_elements));

      left_begin = left_end;
    }
  }

  // TODO: cnt_0 and cnt_1 are small, maybe transform them into uint8_t
  // or even merge into one uint8_t
  // for every entry cnt_0 + cnt_1 <= max_diff
  struct SmallRowEntry {
    size_t cnt_0 = 0;
    size_t cnt_1 = 0;

    size_t class_id = 0;

    double front = 0;
  };

  void seek_small(const CSCMatrix<double>& matrix) {
    const auto [n, d] = matrix.shape();
    const size_t max_size = max_small_row_size + max_diff;

    std::vector<std::vector<SmallRowEntry>> rows(n);
    for (size_t i = 0; i < n; ++i) {
      rows[i].push_back(SmallRowEntry{
          .cnt_0 = 0,
          .cnt_1 = 0,
          .class_id = 0,
          .front = 0,
      });
    }

    size_t class_cnt = 1;
    std::unordered_map<std::pair<size_t, int64_t>, size_t> map;

    auto get_class = [&](size_t prev_class, int64_t norm_value) {
      const auto [itr, inserted] =
          map.emplace(std::pair{prev_class, norm_value}, class_cnt);

      if (inserted) {
        ++class_cnt;
      }

      return itr->second;
    };

    for (size_t col = 0; col < d; ++col) {
      map.clear();

      for (auto [row, value] : matrix.get_column(col)) {
        if (transposed_[row].size() > max_size) {
          continue;
        }

        auto& entries = rows[row];
        std::vector<SmallRowEntry> new_entries;

        for (auto& entry : entries) {
          if (entry.cnt_0 + entry.cnt_1 < max_diff) {
            // class 0 (create new entry)
            {
              SmallRowEntry new_entry = entry;
              ++new_entry.cnt_0;

              new_entries.push_back(new_entry);
            }

            // class 1 (create new entry)
            {
              SmallRowEntry new_entry = entry;
              ++new_entry.cnt_1;
              new_entry.class_id = get_class(entry.class_id, 1);

              new_entries.push_back(new_entry);
            }
          }

          // class 2 (update current entry)
          if (entry.front == 0) {
            entry.front = value;
          }

          auto normalized = normalize_double(value / entry.front);
          entry.class_id = get_class(entry.class_id, normalized);
        }

        entries.append_range(new_entries);
      }
    }

    std::vector<std::vector<std::vector<std::pair<size_t, size_t>>>> m(
        max_diff + 1);
    for (size_t i = 0; i <= max_diff; ++i) {
      m[i].resize(max_diff + 1);
    }

    for (size_t row = 0; row < n; ++row) {
      if (transposed_[row].size() > max_size) {
        continue;
      }

      for (const auto entry : rows[row]) {
        m[entry.cnt_0][entry.cnt_1].emplace_back(entry.class_id, row);
      }
    }

    for (size_t i = 0; i <= max_diff; ++i) {
      for (size_t j = 0; j <= max_diff; ++j) {
        std::ranges::sort(m[i][j], {}, [](auto p) { return p.first; });
      }
    }

    for (size_t i = 0; i <= max_diff; ++i) {
      for (size_t j = i; i + j <= max_diff; ++j) {
        for (size_t k = 0; i + j + k <= max_diff; ++k) {
          compare_hashes(m[i][k], m[j][k]);
        }
      }
    }
  }

 public:
  explicit Tables(size_t max_diff)
      : max_diff(max_diff),
        groups_count(max_diff * 2),
        max_small_row_size(max_diff * 2) {
    if (max_diff != 2) {
      throw std::runtime_error("Not implemented");
    }
  }

  // max_diff = 2
  explicit Tables(size_t max_diff, TablesParameters params)
      : max_diff(max_diff),
        groups_count(params.groups_count),
        max_small_row_size(params.max_small_row_size) {
    if (max_diff != 2) {
      throw std::runtime_error("Not implemented");
    }
  }

  auto seek(const CSCMatrix<double>& matrix) {
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

    do {
      seek_table(matrix, mask, blocks, groups);
    } while (std::ranges::prev_permutation(mask).found);

    // std::println("  finished big rows");
    // now process small rows
    // seek_small();

    // std::println("  finished small rows");

    // remove duplicates from the result
    std::ranges::sort(result_singular_);

    std::vector<std::pair<size_t, size_t>> unique;
    std::ranges::unique_copy(result_singular_, std::back_inserter(unique));

    // TODO: add all rows such that:
    // transposed_[i].size() + transposed_[j].size() <= max_diff

    // TODO: add all rows that intersect in exactly one element
    // this can be done by traversing columns
    seek_small(matrix);

    return std::pair{std::move(unique), std::move(result_bipartite_)};
  }

  Statistics get_stats() const { return statistics_; }
};

}  // namespace seekers
