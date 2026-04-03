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
#include "utils/Hamming.h"
#include "utils/Hashers.h"
#include "utils/Logging.h"
#include "utils/Printing.h"

namespace seekers {

/*
 * If you can look into the seeds of time
 * And say which grain will grow and which will not,
 * Speak, then, to me, who neither beg nor fear
 * Your favors nor your hate.
 */

// TODO: cnt_0 and cnt_1 are small, maybe transform them into uint8_t
// or even merge into one uint8_t
// for every entry cnt_0 + cnt_1 <= max_diff
struct SmallRowEntry {
  size_t cnt_0 = 0;
  size_t cnt_1 = 0;

  size_t class_id = 0;

  double front = 0;
};

class ClassesStorage {
  struct ClassInfo {
    size_t map_column{0};
    std::unordered_map<int64_t, size_t> map;
    std::vector<size_t> counts;
  };

  std::vector<size_t> free_classes_;
  std::vector<ClassInfo> classes_;

  const size_t max_diff;

  void extend_classes_storage(size_t new_size) {
    size_t old_size = classes_.size();
    classes_.resize(new_size);

    for (size_t i = old_size; i < new_size; ++i) {
      classes_[i].counts.resize((max_diff + 1) * (max_diff + 1), 0);
      free_classes_.push_back(i);
    }
  }

 public:
  ClassesStorage(size_t max_diff, size_t rows_count) : max_diff(max_diff) {
    extend_classes_storage(std::max(rows_count, 64uz));
  }

  size_t& get_rows_count(size_t class_id, size_t cnt_0, size_t cnt_1) {
    return classes_[class_id].counts[cnt_0 * (max_diff + 1) + cnt_1];
  }

  size_t& get_rows_count(const SmallRowEntry entry) {
    return get_rows_count(entry.class_id, entry.cnt_0, entry.cnt_1);
  }

  void try_free_class(size_t class_id) {
    bool is_zero = true;

    for (const size_t count : classes_[class_id].counts) {
      if (count != 0) {
        is_zero = false;
        break;
      }
    }

    if (is_zero) {
      free_classes_.push_back(class_id);
    }
  }

  size_t get_class(size_t prev_class, int64_t value, size_t column) {
    if (free_classes_.empty()) {
      extend_classes_storage(classes_.size() * 2);
    }

    if (classes_[prev_class].map_column != column) {
      classes_[prev_class].map.clear();
      classes_[prev_class].map_column = column;
    }

    auto [itr, inserted] =
        classes_[prev_class].map.emplace(value, free_classes_.back());

    if (inserted) {
      free_classes_.pop_back();
    }

    return itr->second;
  }

  size_t pop_free_class() {
    if (free_classes_.empty()) {
      extend_classes_storage(classes_.size() * 2);
    }

    size_t free_class = free_classes_.back();
    free_classes_.pop_back();
    return free_class;
  }

  std::vector<std::vector<size_t>> accumulate_counts() {
    const size_t counts_size = (max_diff + 1) * (max_diff + 1);

    for (size_t i = 1; i < classes_.size(); ++i) {
      for (size_t j = 0; j < counts_size; ++j) {
        classes_[i].counts[j] += classes_[i - 1].counts[j];
      }
    }

    const size_t last_class = classes_.size() - 1;

    std::vector<std::vector<size_t>> total(max_diff + 1);

    for (size_t i = 0; i <= max_diff; ++i) {
      total[i].resize(max_diff + 1, 0);

      for (size_t j = 0; j <= max_diff; ++j) {
        total[i][j] = get_rows_count(last_class, i, j);
      }
    }

    return total;
  }
};

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

  size_t total_considered_big_{0};

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

    for (size_t i = 0; i < n; ++i) {
      if (transposed_[i].size() <= max_small_row_size) {
        continue;
      }

      // TODO: correct seed? RowHasher?
      StreamHasher hasher;

      for (size_t group_id = 0; group_id < groups_mask.size(); ++group_id) {
        if (groups_mask[group_id]) {
          if (front[i] == 0) {
            front[i] = blocks[i][group_id].front;
          }

          hasher << blocks[i][group_id].class_id;
        }
      }

      merged_classes[i] = hasher.get_hash();
      ++classes_sizes[merged_classes[i]];
    }

    size_t prev = 0;
    for (auto& [key, value] : classes_sizes) {
      value += prev;
      prev = value;
    }

    // total count of big rows is stored in prev after previous loop
    std::vector<size_t> indices(prev);
    for (size_t i = 0; i < n; ++i) {
      if (transposed_[i].size() > max_small_row_size) {
        indices[--classes_sizes[merged_classes[i]]] = i;
      }
    }

    for (size_t i = 0; i < prev; ++i) {
      if (transposed_[indices[i]].size() <= max_small_row_size) {
        continue;
      }

      for (size_t j = i + 1;
           j < prev && merged_classes[indices[j]] == merged_classes[indices[i]];
           ++j) {
        ++total_considered_big_;
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

  void seek_small(const CSCMatrix<double>& matrix) {
    const auto [n, d] = matrix.shape();
    const size_t max_size = max_small_row_size + max_diff;

    size_t small_rows_cnt = 0;
    for (size_t i = 0; i < n; ++i) {
      if (transposed_[i].size() <= max_size) {
        ++small_rows_cnt;
      }
    }

    ClassesStorage classes(max_diff, small_rows_cnt);

    size_t init_class = classes.pop_free_class();
    classes.get_rows_count(init_class, 0, 0) = small_rows_cnt;

    std::vector<std::vector<SmallRowEntry>> rows(n);
    for (size_t i = 0; i < n; ++i) {
      if (transposed_[i].size() <= max_size) {
        rows[i].push_back(SmallRowEntry{
            .cnt_0 = 0,
            .cnt_1 = 0,
            .class_id = init_class,
            .front = 0,
        });
      }
    }

    for (size_t col = 0; col < d; ++col) {
      for (auto [row, value] : matrix.get_column(col)) {
        if (transposed_[row].size() > max_size) {
          continue;
        }

        auto& entries = rows[row];
        size_t entries_size = entries.size();

        for (size_t i = 0; i < entries_size; ++i) {
          --classes.get_rows_count(entries[i]);

          if (entries[i].cnt_0 + entries[i].cnt_1 < max_diff) {
            // class 0 (create new entry)
            {
              SmallRowEntry new_entry = entries[i];
              ++new_entry.cnt_0;

              entries.push_back(new_entry);
              ++classes.get_rows_count(new_entry);
            }

            // class 1 (create new entry)
            {
              SmallRowEntry new_entry = entries[i];
              ++new_entry.cnt_1;
              new_entry.class_id =
                  classes.get_class(entries[i].class_id, 0, col);

              entries.push_back(new_entry);
              ++classes.get_rows_count(new_entry);
            }
          }

          // class 2 (update current entry)
          {
            if (entries[i].front == 0) {
              entries[i].front = value;
            }

            auto normalized = normalize_double(value / entries[i].front);
            entries[i].class_id =
                classes.get_class(entries[i].class_id, normalized, col);

            ++classes.get_rows_count(entries[i]);
          }
        }
      }

      for (auto [row, value] : matrix.get_column(col)) {
        if (transposed_[row].size() > max_size) {
          continue;
        }

        auto& entries = rows[row];
        std::vector<SmallRowEntry> new_entries;

        for (auto& entry : entries) {
          bool is_unique = true;

          for (size_t i = 0; entry.cnt_0 + entry.cnt_1 + i <= max_diff; ++i) {
            const size_t count =
                classes.get_rows_count(entry.class_id, i, entry.cnt_1);

            if ((i != entry.cnt_0 && count > 0) ||
                (i == entry.cnt_0 && count > 1)) {
              is_unique = false;
              break;
            }
          }

          if (is_unique) {
            // the row is removed
            --classes.get_rows_count(entry);
            classes.try_free_class(entry.class_id);
          } else {
            new_entries.push_back(entry);
          }
        }

        entries = std::move(new_entries);
      }

      // std::unordered_set<size_t> current_classes;
      // for (const auto& row : rows) {
      //   for (const auto& entry : row) {
      //     current_classes.insert(entry.class_id);
      //   }
      // }
      //
      // logging::log_value(current_classes.size(), "classes_cnt.csv");
    }

    auto total = classes.accumulate_counts();

    std::vector<std::vector<std::vector<std::pair<size_t, size_t>>>> m(
        max_diff + 1);
    for (size_t i = 0; i <= max_diff; ++i) {
      m[i].resize(max_diff + 1);

      for (size_t j = 0; j <= max_diff; ++j) {
        m[i][j].resize(total[i][j]);
      }
    }

    for (size_t row = 0; row < n; ++row) {
      if (transposed_[row].size() > max_size) {
        continue;
      }

      for (const auto entry : rows[row]) {
        m[entry.cnt_0][entry.cnt_1][--classes.get_rows_count(entry)] = {
            entry.class_id, row};
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
  explicit Tables(size_t max_diff, TablesParameters params)
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

  explicit Tables(size_t max_diff)
      : Tables(max_diff, {.groups_count = max_diff * 2,
                          .max_small_row_size = max_diff * 2}) {}

  Result seek(const CSCMatrix<double>& matrix) {
    const size_t selected_groups_count = groups_count - max_diff;

    auto [n, d] = matrix.shape();
    transposed_ = matrix.get_transposed();

    //
    std::vector<std::vector<size_t>> groups(groups_count);

    // TODO: maybe implement some algorithm for choosing groups
    for (size_t col = 0; col < d; ++col) {
      groups[std::hash<size_t>()(col) % groups_count].push_back(col);
    }

    // greedy algorithm
    // std::vector<std::vector<bool>> groups_zero_rows_mask(n);
    // for (size_t row = 0; row < n; ++row) {
    //   groups_zero_rows_mask[row].resize(groups_count, false);
    // }
    //
    // for (size_t col = 0; col < d; ++col) {
    //   std::vector<size_t> increment_per_group(groups_count, 0);
    //
    //   for (const auto [row, value] : matrix.get_column(col)) {
    //     if (transposed_[row].size() <= max_small_row_size) {
    //       continue;
    //     }
    //
    //     for (size_t i = 0; i < groups_count; ++i) {
    //       if (!groups_zero_rows_mask[row][i]) {
    //         ++increment_per_group[i];
    //       }
    //     }
    //   }
    //
    //   size_t max_increment_group = 0;
    //   for (size_t i = 0; i < groups_count; ++i) {
    //     if (increment_per_group[i] >
    //     increment_per_group[max_increment_group]) {
    //       max_increment_group = i;
    //     }
    //   }
    //
    //   groups[max_increment_group].push_back(col);
    //
    //   for (const auto [row, value] : matrix.get_column(col)) {
    //     if (transposed_[row].size() <= max_small_row_size) {
    //       continue;
    //     }
    //
    //     groups_zero_rows_mask[row][max_increment_group] = true;
    //   }
    // }

    auto blocks = get_blocks(matrix, groups);

    std::vector<bool> mask(groups_count, false);
    std::fill_n(mask.begin(), selected_groups_count, true);

    statistics_.big_rows_duration = timing::timeit([&] {
      do {
        seek_table(matrix, mask, blocks);
      } while (std::ranges::prev_permutation(mask).found);
    });

    // remove duplicates from the singular part of the result
    std::ranges::sort(result_singular_);

    std::vector<std::pair<size_t, size_t>> unique;
    std::ranges::unique_copy(result_singular_, std::back_inserter(unique));

    // process all small rows
    statistics_.small_rows_duration =
        timing::timeit([&] { seek_small(matrix); });

    return Result{
        .singular = std::move(unique),
        .bipartite = std::move(result_bipartite_),
    };
  }

  Statistics get_stats() const { return statistics_; }
};

}  // namespace seekers
