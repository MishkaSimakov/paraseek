#pragma once

#include <algorithm>
#include <cmath>
#include <map>
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

  size_t get_used_classes_count() const {
    return classes_.size() - free_classes_.size();
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

  bool log_entries_growth = false;
  bool log_entries_per_row = false;
};

class Tables {
  struct Block {
    size_t class_id;
    double front;
  };

  const size_t max_diff;
  const TablesParameters params;

  Statistics statistics_{};

  // rows that should be compared using methods for small rows
  std::vector<bool> small_rows_mask_;

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
  std::vector<std::vector<Block>> get_blocks(
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
      const double ratio = front[i] / front[j];

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

  void seek_table(const CSCMatrix<double>& matrix,
                  const std::vector<bool>& groups_mask,
                  const std::vector<std::vector<Block>>& blocks) {
    auto [n, d] = matrix.shape();

    std::vector<size_t> merged_classes(n);
    std::vector<double> front(n, 0);
    std::map<size_t, size_t> classes_sizes;
    std::set<size_t> big_classes;

    for (size_t row = 0; row < n; ++row) {
      // TODO: correct seed? RowHasher?
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

    // remove small rows from big classes
    for (const auto [class_id, size] : classes_sizes) {
      if (size > 1000) {
        big_classes.insert(class_id);
      }
    }

    std::println("{}", big_classes.size());

    for (size_t row = 0; row < n; ++row) {
      if (big_classes.contains(merged_classes[row]) &&
          transposed_[row].size() <= params.max_small_row_size + max_diff) {
        small_rows_mask_[row] = true;
        --classes_sizes[merged_classes[row]];
      }
    }

    {
      std::vector<std::pair<size_t, size_t>> by_size;
      for (auto [cls, size] : classes_sizes) {
        by_size.emplace_back(cls, size);
      }

      std::ranges::sort(by_size, {}, [](auto p) { return p.second; });

      std::println("  largest classes:");
      for (size_t i = 0; i < 5; ++i) {
        const auto p = by_size[by_size.size() - i - 1];
        std::println("    {} (size = {})", p.first, p.second);
      }
    }

    size_t prev = 0;
    for (auto& class_size : classes_sizes | std::views::values) {
      if (class_size == 1) {
        continue;
      }

      class_size += prev;
      prev = class_size;
    }

    std::println("prev = {}", prev);
    // total count of big rows is stored in prev after previous loop
    std::vector<size_t> indices(prev);
    for (size_t row = 0; row < n; ++row) {
      if (big_classes.contains(merged_classes[row]) &&
          transposed_[row].size() <= params.max_small_row_size + max_diff) {
        continue;
      }

      if (classes_sizes[merged_classes[row]] == 1) {
        continue;
      }

      indices[--classes_sizes[merged_classes[row]]] = row;
    }

    for (size_t i = 0; i < prev; ++i) {
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

      if (&left == &right && left_end == left_begin + 1 &&
          right_end == right_begin + 1) {
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

    // tried sorting columns by size with no success
    const size_t max_size = params.max_small_row_size + max_diff;

    size_t small_rows_cnt = 0;
    for (size_t row = 0; row < n; ++row) {
      if (small_rows_mask_[row]) {
        ++small_rows_cnt;
      }
    }

    ClassesStorage classes(max_diff, small_rows_cnt);

    size_t init_class = classes.pop_free_class();
    classes.get_rows_count(init_class, 0, 0) = small_rows_cnt;

    std::vector<std::vector<SmallRowEntry>> rows(n);
    for (size_t row = 0; row < n; ++row) {
      if (small_rows_mask_[row]) {
        rows[row].push_back(SmallRowEntry{
            .cnt_0 = 0,
            .cnt_1 = 0,
            .class_id = init_class,
            .front = 0,
        });
      }
    }

    for (size_t col = 0; col < d; ++col) {
      for (auto [row, value] : matrix.get_column(col)) {
        if (!small_rows_mask_[row]) {
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

      // classes pruning
      for (const size_t row : matrix.get_column(col) | std::views::keys) {
        if (!small_rows_mask_[row]) {
          continue;
        }

        std::vector<SmallRowEntry> saved_entries;

        for (auto& entry : rows[row]) {
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
            saved_entries.push_back(entry);
          }
        }

        rows[row] = std::move(saved_entries);
      }

      // log entries growth
      if (params.log_entries_growth) {
        // classes count
        size_t total_entries = 0;
        for (const auto& entries : rows) {
          total_entries += entries.size();
        }

        // considered entries
        size_t considered_entries = 0;

        for (auto [row, value] : matrix.get_column(col)) {
          if (small_rows_mask_[row]) {
            considered_entries += rows[row].size();
          }
        }

        logging::log_csv<size_t>(
            {
                {"column", col},
                {"classes_cnt", classes.get_used_classes_count()},
                {"considered_entries", considered_entries},
            },
            "classes_cnt.csv");
      }
    }

    if (params.log_entries_per_row) {
      std::unordered_set<size_t> logged_classes;

      for (size_t row = 0; row < n; ++row) {
        for (const auto& entry : rows[row]) {
          if (logged_classes.contains(entry.class_id)) {
            continue;
          }

          logged_classes.insert(entry.class_id);

          for (size_t c0 = 0; c0 <= 2; ++c0) {
            for (size_t c1 = 0; c0 + c1 <= 2; ++c1) {
              logging::log_csv<size_t>(
                  {{"class_id", entry.class_id},
                   {"c0", c0},
                   {"c1", c1},
                   {"size", classes.get_rows_count(entry.class_id, c0, c1)}},
                  "classes_sizes.csv");
            }
          }
        }
      }
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
      if (!small_rows_mask_[row]) {
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
      : max_diff(max_diff), params(params) {
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
    const size_t selected_groups_count = params.groups_count - max_diff;

    auto [n, d] = matrix.shape();

    transposed_ = matrix.get_transposed();
    small_rows_mask_.resize(n, false);

    const auto groups =
        splitters::GreedySplitter().split(matrix, params.groups_count);

    auto blocks = get_blocks(matrix, groups);

    std::vector<bool> mask(params.groups_count, false);
    std::fill_n(mask.begin(), selected_groups_count, true);

    statistics_.big_rows_duration = timing::timeit([&] {
      do {
        seek_table(matrix, mask, blocks);
      } while (std::ranges::prev_permutation(mask).found);
    });

    {
      std::vector<size_t> col_group(d);

      for (size_t i = 0; i < groups.size(); ++i) {
        for (const size_t col : groups[i]) {
          col_group[col] = i;
        }
      }

      size_t small_rows_cnt = 0;
      size_t nonzero_small_rows_cnt = 0;

      for (size_t row = 0; row < n; ++row) {
        if (transposed_[row].size() > params.max_small_row_size + max_diff) {
          continue;
        }

        ++small_rows_cnt;

        std::vector<bool> mask(params.groups_count, false);
        std::fill_n(mask.begin(), selected_groups_count, true);

        bool all_nonzero = true;

        do {
          bool has_nonzero = false;
          for (const auto col : transposed_[row] | std::views::keys) {
            if (mask[col_group[col]]) {
              has_nonzero = true;
              break;
            }
          }

          if (!has_nonzero) {
            all_nonzero = false;
            break;
          }
        } while (std::ranges::prev_permutation(mask).found);

        if (all_nonzero) {
          ++nonzero_small_rows_cnt;
        }
      }

      logging::log_csv<size_t>({{"small_rows", small_rows_cnt},
                                {"nonzero_small_rows", nonzero_small_rows_cnt}},
                               "nonzero_rows.csv");
    }

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

// 73
//   largest classes:
//     128885107609 (size = 496)
//     128885087214 (size = 493)
//     128885108496 (size = 484)
//     128885108273 (size = 483)
//     128885047942 (size = 478)
// prev = 1453618
// 44
//   largest classes:
//     128885105188 (size = 490)
//     128885108367 (size = 485)
//     128915575830 (size = 483)
//     128885083439 (size = 479)
//     128885108281 (size = 475)
// prev = 1525718
// 50
//   largest classes:
//     128884897303 (size = 497)
//     128914326917 (size = 496)
//     128885108247 (size = 496)
//     128885108451 (size = 493)
//     128885023931 (size = 492)
// prev = 1581870
// 108
//   largest classes:
//     128885108421 (size = 1178)
//     128885057069 (size = 542)
//     128885081214 (size = 500)
//     128885108332 (size = 500)
//     128885045515 (size = 490)
// prev = 1357329
// 54
//   largest classes:
//     128885108421 (size = 669)
//     128885078693 (size = 493)
//     128884688458 (size = 486)
//     128885053350 (size = 482)
//     128885108451 (size = 467)
// prev = 1425828
// 47
//   largest classes:
//     128885105934 (size = 486)
//     128885108322 (size = 476)
//     128915804287 (size = 467)
//     128885085323 (size = 456)
//     128885108553 (size = 446)
// prev = 1455672
