#pragma once

#include <algorithm>
#include <cmath>
#include <numeric>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "Result.h"
#include "matrix/CSCMatrix.h"
#include "seekers/Statistics.h"
#include "splitters/GreedySplitter.h"
#include "utils/Hamming.h"
#include "utils/Hashers.h"
#include "utils/Logging.h"

namespace seekers {

/*
 * If you can look into the seeds of time
 * And say which grain will grow and which will not,
 * Speak, then, to me, who neither beg nor fear
 * Your favors nor your hate.
 */

struct TablesSimpleHashingParameters {
  size_t groups_count;
  size_t max_small_row_size;
};

class TablesSimpleHashing {
  struct Block {
    size_t hash;
    double front;
  };

  const size_t max_diff;
  const size_t groups_count;
  const size_t max_small_row_size;

  // stores indices of rows sorted by size (e.g. sorted_[0] is the smallest row)
  std::vector<size_t> sorted_;
  size_t big_rows_start_;

  Statistics statistics_;

  std::vector<SparseVector<double>> transposed_;

  std::vector<std::pair<size_t, size_t>> result_singular_;
  std::vector<std::pair<std::vector<size_t>, std::vector<size_t>>>
      result_bipartite_;

  // Prepares double for hashing
  static int64_t normalize_double(double value) {
    return std::round(value * 1e10);
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

    for (size_t group_id = 0; group_id < groups.size(); ++group_id) {
      for (const size_t col : groups[group_id]) {
        for (auto [row, value] : matrix.get_column(col)) {
          if (blocks[row][group_id].front == 0) {
            blocks[row][group_id].front = value;
          }

          RowHasher hasher(blocks[row][group_id].hash);

          // hash position and value
          hasher << col;
          hasher << normalize_double(value / blocks[row][group_id].front);

          blocks[row][group_id].hash = hasher.get_hash();
        }
      }
    }

    return blocks;
  }

  void add_to_answer(size_t i, size_t j) {
    if (i < j) {
      result_singular_.emplace_back(i, j);
    } else if (i > j) {
      result_singular_.emplace_back(j, i);
    }
  }

  void consider_pair(size_t i, size_t j, const std::vector<double>& front) {
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
                  const std::vector<std::vector<Block>>& blocks) {
    auto [n, d] = matrix.shape();

    std::vector<size_t> merged_classes(n);
    std::vector<double> front(n, 0);
    std::unordered_map<size_t, size_t> classes_sizes;

    for (size_t i = big_rows_start_; i < n; ++i) {
      // TODO: correct seed? RowHasher?
      const size_t row = sorted_[i];

      StreamHasher hasher;

      for (size_t group_id = 0; group_id < groups_mask.size(); ++group_id) {
        if (groups_mask[group_id]) {
          if (front[row] == 0) {
            front[row] = blocks[row][group_id].front;
          }

          hasher << blocks[row][group_id].hash;
        }
      }

      merged_classes[row] = hasher.get_hash();
      ++classes_sizes[merged_classes[row]];
    }

    size_t prev = 0;
    for (size_t& value : classes_sizes | std::views::values) {
      const size_t old = value;

      value = prev;
      prev += old;
    }

    // total count of big rows is stored in prev after previous loop
    // use stable sort, so that rows are remain sorted by size inside classes
    std::vector<size_t> indices(n - big_rows_start_);

    for (size_t i = big_rows_start_; i < n; ++i) {
      const size_t row = sorted_[i];
      indices[classes_sizes[merged_classes[row]]++] = row;
    }

    for (size_t i = 0; i < indices.size(); ++i) {
      for (size_t j = i + 1;
           j < indices.size() &&
           merged_classes[indices[j]] == merged_classes[indices[i]] &&
           transposed_[indices[j]].size() <=
               transposed_[indices[i]].size() + max_diff;
           ++j) {
        consider_pair(indices[i], indices[j], front);
      }
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

  struct SmallRowEntry {
    size_t cnt_0 = 0;
    size_t cnt_1 = 0;

    size_t hash = 0;

    double front = 0;
  };

  void seek_small(const CSCMatrix<double>& matrix) {
    const auto [n, d] = matrix.shape();
    const size_t max_size = max_small_row_size + max_diff;

    size_t small_rows_cnt = 0;
    for (size_t i = 0; i < n; ++i) {
      if (transposed_[i].size() <= max_size) {
        ++small_rows_cnt;
      }
    }

    if (small_rows_cnt == 0) {
      return;
    }

    constexpr size_t hashtable_multiplier = 17;
    const size_t hashtable_size = hashtable_multiplier * small_rows_cnt;

    std::vector<size_t> counts(hashtable_size * (max_diff + 1) * (max_diff + 1),
                               0);

    // a + b <= 2
    // 0 0
    // 0 1
    // 1 0
    // 1 1
    // 2 0
    // 0 2

    // TODO: counts size can be two times smaller
    // max_diff + (max_diff - 1) + ... + (max_diff - k - 1) = max_diff * k - (k
    // - 1) * k / 2

    const auto get_index = [&](SmallRowEntry entry) {
      return (entry.hash % hashtable_size) * (max_diff + 1) * (max_diff + 1) +
             entry.cnt_0 * (max_diff + 1) + entry.cnt_1;
    };

    // std::vector<std::unordered_map<int64_t, size_t>> real_counts(9);

    std::vector<std::vector<SmallRowEntry>> rows(n);
    for (size_t i = 0; i < n; ++i) {
      if (transposed_[i].size() <= max_size) {
        rows[i].push_back(SmallRowEntry{});
      }
    }

    for (size_t col = 0; col < d; ++col) {
      std::println("column: {}", col);
      const size_t column_size = matrix.get_column(col).size();

      if (column_size <= 2) {
        // compare and add rows manually
        for (size_t i = 0; i < column_size; ++i) {
          for (size_t j = i + 1; j < column_size; ++j) {
            if (similarity::hamming_leq(transposed_[i], transposed_[j],
                                        max_diff)) {
              add_to_answer(i, j);
            }
          }
        }

        for (const auto [row, value] : matrix.get_column(col)) {
          auto& entries = rows[row];

          for (auto& entry : entries) {
            --counts[get_index(entry)];

            ++entry.cnt_0;

            if (entry.cnt_0 + entry.cnt_1 <= max_diff) {
              ++counts[get_index(entry)];
            }
          }

          std::erase_if(entries, [&](const SmallRowEntry& entry) {
            return entry.cnt_0 + entry.cnt_1 > max_diff;
          });
        }

        continue;
      }

      std::println("  adding...");
      for (auto [row, value] : matrix.get_column(col)) {
        auto& entries = rows[row];
        // DO NOT REMOVE THIS: entries.size() does change during loop execution
        const size_t entries_size = entries.size();

        for (size_t i = 0; i < entries_size; ++i) {
          --counts[get_index(entries[i])];

          if (entries[i].cnt_0 + entries[i].cnt_1 < max_diff) {
            // class 0 (create new entry)
            {
              SmallRowEntry new_entry = entries[i];
              ++new_entry.cnt_0;

              ++counts[get_index(new_entry)];
              entries.push_back(new_entry);
            }
          }

          if (entries[i].cnt_0 + entries[i].cnt_1 < max_diff) {
            // class 1 (create new entry)
            {
              SmallRowEntry new_entry = entries[i];
              ++new_entry.cnt_1;

              RowHasher hasher(entries[i].hash);

              hasher << col;
              hasher << 0;

              new_entry.hash = hasher.get_hash();

              ++counts[get_index(new_entry)];
              entries.push_back(new_entry);
            }
          }

          // class 2 (update current entry)
          {
            if (entries[i].front == 0) {
              entries[i].front = value;
            }

            RowHasher hasher(entries[i].hash);

            hasher << col;
            hasher << normalize_double(value / entries[i].front);

            entries[i].hash = hasher.get_hash();

            ++counts[get_index(entries[i])];
          }
        }
      }

      // erase unique entries
      std::println("  removing unique...");
      for (auto [row, value] : matrix.get_column(col)) {
        std::erase_if(rows[row], [&](const SmallRowEntry& entry) {
          if (counts[get_index(entry)] == 1) {
            counts[get_index(entry)] = 0;
            return true;
          }

          return false;
        });
      }

      // reset counts
      // for (auto [row, value] : matrix.get_column(col)) {
      //   for (const auto& entry : rows[row]) {
      //     if (counts[get_index(entry.hash, col, entry.cnt_0, entry.cnt_1)] -
      //             real_counts[entry.cnt_0 * (max_diff + 1) + entry.cnt_1]
      //                        [entry.hash] >
      //         10) {
      //       std::println(
      //           "collision, {}, {}",
      //           counts[get_index(entry.hash, col, entry.cnt_0, entry.cnt_1)],
      //           real_counts[entry.cnt_0 * (max_diff + 1) + entry.cnt_1]
      //                      [entry.hash]);
      //     }
      //   }
      // }

      // for (size_t i = 0; i < 9; ++i) {
      // real_counts[i].clear();
      // }

      // calculate and log total entries count
      size_t total_count = 0;
      for (size_t row = 0; row < n; ++row) {
        total_count += rows[row].size();
      }

      logging::log_csv<size_t>({{"count", total_count}},
                               "entries_count_simple.csv");
    }

    // calculate prefix sums
    for (size_t i = 1; i < hashtable_size; ++i) {
      for (size_t j = 0; j < (max_diff + 1) * (max_diff + 1); ++j) {
        counts[i * (max_diff + 1) * (max_diff + 1) + j] +=
            counts[(i - 1) * (max_diff + 1) * (max_diff + 1) + j];
      }
    }

    std::vector<std::vector<std::vector<std::pair<size_t, size_t>>>> m(
        max_diff + 1);
    for (size_t i = 0; i <= max_diff; ++i) {
      m[i].resize(max_diff + 1);

      for (size_t j = 0; j <= max_diff; ++j) {
        m[i][j].resize(
            counts[get_index(SmallRowEntry{hashtable_size - 1, i, j})]);
      }
    }

    for (size_t row = 0; row < n; ++row) {
      for (const auto entry : rows[row]) {
        m[entry.cnt_0][entry.cnt_1][--counts[get_index(entry)]] = {entry.hash,
                                                                   row};
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
  TablesSimpleHashing(size_t max_diff, TablesSimpleHashingParameters params)
      : max_diff(max_diff),
        groups_count(params.groups_count),
        max_small_row_size(params.max_small_row_size) {
    if (max_diff != 2) {
      throw std::invalid_argument("Not implemented");
    }

    if (params.max_small_row_size < max_diff * 2) {
      throw std::invalid_argument(
          "max_small_row_size must be at least max_diff * 2");
    }
  }

  explicit TablesSimpleHashing(size_t max_diff)
      : TablesSimpleHashing(max_diff, {.groups_count = max_diff * 2,
                                       .max_small_row_size = max_diff * 2}) {}

  Result seek(const CSCMatrix<double>& matrix) {
    const size_t selected_groups_count = groups_count - max_diff;

    auto [n, d] = matrix.shape();
    transposed_ = matrix.get_transposed();

    sorted_ = std::vector<size_t>(n);
    std::iota(sorted_.begin(), sorted_.end(), 0);

    std::ranges::sort(sorted_, {},
                      [&](size_t row) { return transposed_[row].size(); });

    big_rows_start_ = n;

    for (size_t i = 0; i < n; ++i) {
      if (transposed_[sorted_[i]].size() > max_small_row_size) {
        big_rows_start_ = i;
        break;
      }
    }

    //
    auto groups = splitters::GreedySplitter().split(matrix, groups_count);

    auto blocks = get_blocks(matrix, groups);

    std::vector<bool> mask(groups_count, false);
    std::fill_n(mask.begin(), selected_groups_count, true);

    do {
      seek_table(matrix, mask, blocks);
    } while (std::ranges::prev_permutation(mask).found);

    // remove duplicates from the singular part of the result
    std::ranges::sort(result_singular_);

    std::vector<std::pair<size_t, size_t>> unique;
    std::ranges::unique_copy(result_singular_, std::back_inserter(unique));

    // process all small rows
    seek_small(matrix);

    return Result{
        .singular = std::move(unique),
        .bipartite = std::move(result_bipartite_),
    };
  }

  Statistics get_stats() const { return statistics_; }
};

}  // namespace seekers
