#pragma once

#include <algorithm>
#include <cmath>
#include <numeric>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "Hamming.h"
#include "Result.h"
#include "matrix/CSCMatrix.h"
#include "seekers/Statistics.h"
#include "splitters/GreedySplitter.h"
#include "splitters/RandomSplitter.h"
#include "utils/Hashers.h"
#include "utils/Logging.h"

namespace seekers {

// Stores information about patterns for small rows
class ClassesStorage {
  struct ClassInfo {
    size_t map_column{0};
    std::unordered_map<size_t, size_t> map;
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

  size_t get_class(size_t prev_class, size_t hash, size_t column) {
    if (free_classes_.empty()) {
      extend_classes_storage(classes_.size() * 2);
    }

    if (classes_[prev_class].map_column != column) {
      classes_[prev_class].map.clear();
      classes_[prev_class].map_column = column;
    }

    auto [itr, inserted] =
        classes_[prev_class].map.emplace(hash, free_classes_.back());

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

struct TablesParameters {
  size_t groups_count;
  size_t max_small_row_size;

  bool entries_reduction = true;
  size_t small_column_limit = 2;

  std::string log_prefix;
  bool log_entries_growth;
};

// Hashes doubles, but not too precisely
struct DoubleHasher {
  size_t operator()(double value) const {
    return std::hash<double>()(std::round(value * 1e10));
  }
};

template <typename Field, typename FieldHasher,
          typename Splitter = splitters::RandomSplitter<Field>>
class Tables {
  // TODO: cnt_0 and cnt_1 are small, maybe transform them into uint8_t
  // or even merge into one uint8_t
  // for every entry cnt_0 + cnt_1 <= max_diff
  struct SmallRowEntry {
    size_t cnt_0 = 0;
    size_t cnt_1 = 0;

    size_t class_id = 0;
    Field front = 0;
  };

  struct Block {
    size_t class_id;
    Field front;
  };

  [[no_unique_address]]
  FieldHasher hasher_;

  const size_t max_diff;
  const TablesParameters params;

  // stores indices of rows sorted by size (e.g. sorted_[0] is the smallest row)
  std::vector<size_t> sorted_;
  size_t big_rows_start_{};

  Statistics statistics_;

  std::vector<SparseVector<Field>> transposed_;

  std::vector<std::pair<size_t, size_t>> result_singular_;
  std::vector<std::pair<std::vector<size_t>, std::vector<size_t>>>
      result_bipartite_;

  // For each row a list of blocks is returned
  // Uses smart hashing that avoids collisions
  std::vector<std::vector<Block>> get_blocks_smart(
      const CSCMatrix<Field>& matrix,
      const std::vector<std::vector<size_t>>& groups) const {
    auto [n, d] = matrix.shape();

    std::vector<std::vector<Block>> blocks(n);
    for (size_t i = 0; i < n; ++i) {
      blocks[i].resize(groups.size(), Block{0, 0});
    }

    std::unordered_map<std::pair<size_t, size_t>, size_t> map;

    for (size_t group_id = 0; group_id < groups.size(); ++group_id) {
      size_t classes_cnt = 1;

      for (const size_t col : groups[group_id]) {
        map.clear();

        for (const auto [row, value] : matrix.get_column(col)) {
          if (blocks[row][group_id].front == 0) {
            blocks[row][group_id].front = value;
          }

          const size_t hash = hasher_(value / blocks[row][group_id].front);

          auto [itr, new_class] = map.emplace(
              std::pair{blocks[row][group_id].class_id, hash}, classes_cnt);

          if (new_class) {
            ++classes_cnt;
          }

          blocks[row][group_id].class_id = itr->second;
        }
      }
    }

    return blocks;
  }

  // Same as get_blocks_smart, but uses simple hashing that allows collisions
  std::vector<std::vector<Block>> get_blocks_simple(
      const CSCMatrix<Field>& matrix,
      const std::vector<std::vector<size_t>>& groups) const {
    auto [n, d] = matrix.shape();

    std::vector<std::vector<Block>> blocks(n);
    for (size_t i = 0; i < n; ++i) {
      blocks[i].resize(groups.size(), Block{0, 0});
    }

    for (size_t group_id = 0; group_id < groups.size(); ++group_id) {
      for (const size_t col : groups[group_id]) {
        for (const auto [row, value] : matrix.get_column(col)) {
          if (blocks[row][group_id].front == 0) {
            blocks[row][group_id].front = value;
          }

          const size_t value_hash =
              hasher_(value / blocks[row][group_id].front);

          RowHasher hasher(blocks[row][group_id].class_id);

          hasher << col << value_hash;

          blocks[row][group_id].class_id = hasher.get_hash();
        }
      }
    }

    return blocks;
  }

  void add_to_answer(size_t i, size_t j) {
    if (i > j) {
      result_singular_.emplace_back(j, i);
    } else if (i < j) {
      result_singular_.emplace_back(i, j);
    }
  }

  void consider_pair(size_t i, size_t j, const std::vector<Field>& front) {
    if (front[i] != 0) {
      const Field ratio = front[i] / front[j];

      if (similarity::hamming_fixed_ratio_leq(transposed_[i], transposed_[j],
                                              ratio, max_diff)) {
        add_to_answer(i, j);
      }
    } else {
      if (similarity::hamming_leq(transposed_[i], transposed_[j], max_diff)) {
        add_to_answer(i, j);
      }
    }
  }

  void seek_table(const CSCMatrix<Field>& matrix,
                  const std::vector<bool>& groups_mask,
                  const std::vector<std::vector<Block>>& blocks) {
    auto [n, d] = matrix.shape();

    std::vector<size_t> merged_classes(n);
    std::vector<Field> front(n, 0);
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

          hasher << blocks[row][group_id].class_id;
        }
      }

      merged_classes[row] = hasher.get_hash();
      ++classes_sizes[merged_classes[row]];
    }

    // TODO: try simple std::sort here
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
        ++statistics_.pairs_considered;
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

      if (&left == &right && left_end - left_begin == 1 &&
          right_end - right_begin == 1) {
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

  void seek_small(const CSCMatrix<Field>& matrix) {
    const auto [n, d] = matrix.shape();
    const size_t max_size = params.max_small_row_size + max_diff;

    // sort columns by size
    std::vector<size_t> columns_order(d);
    std::iota(columns_order.begin(), columns_order.end(), 0);

    std::ranges::sort(columns_order, {}, [&](size_t col) {
      return matrix.get_column(col).size();
    });

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

    // for logging
    std::vector<size_t> total_entries_count;
    size_t current_entries_count = small_rows_cnt;

    for (const size_t col : columns_order) {
      if (params.log_entries_growth) {
        total_entries_count.push_back(current_entries_count);
      }

      const size_t column_size = matrix.get_column(col).size();

      size_t col_small_rows = 0;
      for (const auto [row, value] : matrix.get_column(col)) {
        if (transposed_[row].size() <= max_size) {
          ++col_small_rows;
        }
      }

      // TODO: for small problems, 1 is optimal (and big loop can be removed)
      // for big problems, 2-3 looks optimal
      // maybe should add a compile-time parameter
      if (col_small_rows <= params.small_column_limit &&
          params.entries_reduction) {
        // compare and add rows manually
        for (size_t i = 0; i < column_size; ++i) {
          const size_t row1 = matrix.get_column(col)[i].first;

          if (transposed_[row1].size() > max_size) {
            continue;
          }

          for (size_t j = i + 1; j < column_size; ++j) {
            const size_t row2 = matrix.get_column(col)[j].first;

            if (transposed_[row2].size() > max_size) {
              continue;
            }

            if (similarity::hamming_leq(transposed_[row1], transposed_[row2],
                                        max_diff)) {
              add_to_answer(row1, row2);
            }
          }
        }

        for (const auto [row, value] : matrix.get_column(col)) {
          auto& entries = rows[row];
          // DO NOT REMOVE THIS: entries.size() does change during loop
          // execution
          const size_t entries_size = entries.size();

          for (size_t i = 0; i < entries_size; ++i) {
            --classes.get_rows_count(entries[i].class_id, entries[i].cnt_0,
                                     entries[i].cnt_1);
            ++entries[i].cnt_0;

            // otherwise this entry will be removed later
            if (entries[i].cnt_0 + entries[i].cnt_1 <= max_diff) {
              ++classes.get_rows_count(entries[i].class_id, entries[i].cnt_0,
                                       entries[i].cnt_1);
            }
          }

          size_t erased_cnt =
              std::erase_if(entries, [&](const SmallRowEntry& entry) {
                return entry.cnt_0 + entry.cnt_1 > max_diff;
              });

          current_entries_count -= erased_cnt;
        }

        continue;
      }

      for (auto [row, value] : matrix.get_column(col)) {
        auto& entries = rows[row];
        // DO NOT REMOVE THIS: entries.size() does change during loop execution
        const size_t entries_size = entries.size();

        for (size_t i = 0; i < entries_size; ++i) {
          ++statistics_.entries_considered;

          --classes.get_rows_count(entries[i].class_id, entries[i].cnt_0,
                                   entries[i].cnt_1);
          --current_entries_count;

          if (entries[i].cnt_0 + entries[i].cnt_1 < max_diff) {
            // class 0 (create new entry)
            {
              SmallRowEntry new_entry = entries[i];
              ++new_entry.cnt_0;

              entries.push_back(new_entry);
              ++classes.get_rows_count(new_entry.class_id, new_entry.cnt_0,
                                       new_entry.cnt_1);
              ++current_entries_count;
            }
          }

          if (entries[i].cnt_0 + entries[i].cnt_1 < max_diff) {
            // class 1 (create new entry)
            {
              SmallRowEntry new_entry = entries[i];
              ++new_entry.cnt_1;
              new_entry.class_id =
                  classes.get_class(entries[i].class_id, 0, col);

              entries.push_back(new_entry);
              ++classes.get_rows_count(new_entry.class_id, new_entry.cnt_0,
                                       new_entry.cnt_1);
              ++current_entries_count;
            }
          }

          // class 2 (update current entry)
          {
            if (entries[i].front == 0) {
              entries[i].front = value;
            }

            const size_t hash = hasher_(value / entries[i].front);
            entries[i].class_id =
                classes.get_class(entries[i].class_id, hash, col);
            ++classes.get_rows_count(entries[i].class_id, entries[i].cnt_0,
                                     entries[i].cnt_1);
            ++current_entries_count;
          }
        }
      }

      if (!params.entries_reduction) {
        continue;
      }

      for (auto [row, value] : matrix.get_column(col)) {
        size_t erased_cnt =
            std::erase_if(rows[row], [&](const SmallRowEntry& entry) {
              for (size_t i = 0; entry.cnt_0 + entry.cnt_1 + i <= max_diff;
                   ++i) {
                const size_t count =
                    classes.get_rows_count(entry.class_id, i, entry.cnt_1);

                if ((i != entry.cnt_0 && count > 0) ||
                    (i == entry.cnt_0 && count > 1)) {
                  return false;
                }
              }

              // the entry has unique class => it is removed
              --classes.get_rows_count(entry.class_id, entry.cnt_0,
                                       entry.cnt_1);
              classes.try_free_class(entry.class_id);

              return true;
            });

        current_entries_count -= erased_cnt;
      }
    }

    if (params.log_entries_growth) {
      std::ofstream os(paths::log(params.log_prefix + "entries_growth.csv"));

      std::println(os, "count");

      for (size_t value : total_entries_count) {
        std::println(os, "{}", value);
      }
    }

    const auto total = classes.accumulate_counts();

    std::vector<std::vector<std::vector<std::pair<size_t, size_t>>>> m(
        max_diff + 1);
    for (size_t i = 0; i <= max_diff; ++i) {
      m[i].resize(max_diff + 1);

      for (size_t j = 0; j <= max_diff; ++j) {
        m[i][j].resize(total[i][j]);
      }
    }

    for (size_t row = 0; row < n; ++row) {
      for (const auto entry : rows[row]) {
        m[entry.cnt_0][entry.cnt_1][--classes.get_rows_count(
            entry.class_id, entry.cnt_0, entry.cnt_1)] = {entry.class_id, row};
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
  Tables(size_t max_diff, TablesParameters params)
      : max_diff(max_diff), params(params) {
    if (params.max_small_row_size < max_diff * 2) {
      std::cerr << "Warning: max_small_row_size must be at least max_diff * 2. "
                   "Otherwise correctness is not guaranteed!"
                << std::endl;
    }
  }

  explicit Tables(size_t max_diff)
      : Tables(max_diff, {.groups_count = max_diff * 2,
                          .max_small_row_size = max_diff * 2}) {}

  Result seek(const CSCMatrix<Field>& matrix) {
    const size_t selected_groups_count = params.groups_count - max_diff;

    auto [n, d] = matrix.shape();
    transposed_ = matrix.get_transposed();

    sorted_ = std::vector<size_t>(n);
    std::iota(sorted_.begin(), sorted_.end(), 0);

    std::ranges::sort(sorted_, {},
                      [&](size_t row) { return transposed_[row].size(); });

    big_rows_start_ = n;

    for (size_t i = 0; i < n; ++i) {
      if (transposed_[sorted_[i]].size() > params.max_small_row_size) {
        big_rows_start_ = i;
        break;
      }
    }

    //
    statistics_.big_rows_duration = timing::timeit([&] {
      auto groups = Splitter().split(matrix, params.groups_count);

      auto blocks = get_blocks_simple(matrix, groups);

      std::vector<bool> mask(params.groups_count, false);
      std::fill_n(mask.begin(), selected_groups_count, true);

      do {
        seek_table(matrix, mask, blocks);
      } while (std::ranges::prev_permutation(mask).found);
    });

    // process all small rows
    statistics_.small_rows_duration =
        timing::timeit([&] { seek_small(matrix); });

    // remove duplicates from the singular part of the result
    std::ranges::sort(result_singular_);

    std::vector<std::pair<size_t, size_t>> unique;
    std::ranges::unique_copy(result_singular_, std::back_inserter(unique));

    return Result{
        .singular = std::move(unique),
        .bipartite = std::move(result_bipartite_),
    };
  }

  Statistics get_stats() const { return statistics_; }
};

}  // namespace seekers
