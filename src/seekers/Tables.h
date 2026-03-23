#pragma once

#include <algorithm>
#include <cmath>
#include <numeric>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "../utils/Hamming.h"
#include "../utils/ZipRows.h"
#include "matrix/CSCMatrix.h"
#include "seekers/Statistics.h"
#include "utils/Hashers.h"
#include "utils/Logging.h"

namespace seekers {

std::vector<std::pair<size_t, size_t>> normalize_tables_result(
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

  return {result.begin(), result.end()};
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

  // TODO: cnt_0 and cnt_1 are small, maybe transform them into uint8_t
  // or even merge into one uint8_t
  // for every entry cnt_0 + cnt_1 <= max_diff
  struct SmallRowEntry {
    size_t cnt_0 = 0;
    size_t cnt_1 = 0;

    size_t class_id = 0;
    bool* is_unique = nullptr;

    double front = 0;
  };

  void seek_small(const CSCMatrix<double>& matrix) {
    const auto [n, d] = matrix.shape();
    const size_t max_size = max_small_row_size + max_diff;

    size_t small_rows_cnt = 0;
    std::vector<std::vector<SmallRowEntry>> rows(n);
    for (size_t i = 0; i < n; ++i) {
      if (transposed_[i].size() > max_size) {
        continue;
      }

      rows[i].push_back(SmallRowEntry{
          .cnt_0 = 0,
          .cnt_1 = 0,
          .class_id = 0,
          .front = 0,
      });
      ++small_rows_cnt;
    }

    size_t class_cnt = 1;
    std::unordered_map<std::pair<size_t, int64_t>, std::pair<size_t, bool>> map;

    // class_id -> (i, j) -> rows_count
    std::unordered_map<size_t, std::vector<size_t>> rows_cnt;

    {
      std::vector<size_t> counts((max_diff + 1) * (max_diff + 1), 0);
      counts[0] = small_rows_cnt;

      rows_cnt.emplace(0, std::move(counts));
    }

    auto get_class = [&](size_t prev_class,
                         int64_t norm_value) -> std::pair<size_t, bool>& {
      const auto [itr, inserted] = map.emplace(
          std::pair{prev_class, norm_value}, std::pair{class_cnt, true});

      if (inserted) {
        // new class is created, insert entry in rows_cnt for it
        std::vector<size_t> counts((max_diff + 1) * (max_diff + 1), 0);
        rows_cnt.emplace(class_cnt, std::move(counts));

        ++class_cnt;
      } else {
        itr->second.second = false;
      }

      return itr->second;
    };

    auto get_index = [&](size_t cnt_0, size_t cnt_1) -> size_t {
      return cnt_0 * (max_diff + 1) + cnt_1;
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
          --rows_cnt.at(entry.class_id)[get_index(entry.cnt_0, entry.cnt_1)];

          if (entry.cnt_0 + entry.cnt_1 < max_diff) {
            // class 0 (create new entry)
            {
              SmallRowEntry new_entry = entry;
              ++new_entry.cnt_0;
              new_entry.is_unique = nullptr;

              new_entries.push_back(new_entry);

              ++rows_cnt.at(
                  entry.class_id)[get_index(new_entry.cnt_0, new_entry.cnt_1)];
            }

            // class 1 (create new entry)
            {
              SmallRowEntry new_entry = entry;
              ++new_entry.cnt_1;

              auto& new_class = get_class(entry.class_id, 0);
              new_entry.class_id = new_class.first;
              new_entry.is_unique = &new_class.second;

              new_entries.push_back(new_entry);

              ++rows_cnt.at(new_entry.class_id)[get_index(new_entry.cnt_0,
                                                          new_entry.cnt_1)];
            }
          }

          // class 2 (update current entry)
          {
            auto new_entry = entry;

            if (new_entry.front == 0) {
              new_entry.front = value;
            }

            auto normalized = normalize_double(value / new_entry.front);
            auto& new_class = get_class(new_entry.class_id, normalized);
            new_entry.class_id = new_class.first;
            new_entry.is_unique = &new_class.second;

            new_entries.push_back(new_entry);

            ++rows_cnt.at(new_entry.class_id)[get_index(new_entry.cnt_0,
                                                        new_entry.cnt_1)];
          }
        }

        entries = std::move(new_entries);
      }

      for (auto [row, value] : matrix.get_column(col)) {
        if (transposed_[row].size() > max_size) {
          continue;
        }

        auto& entries = rows[row];
        std::vector<SmallRowEntry> new_entries;

        for (auto& entry : entries) {
          bool is_unique = true;

          auto& counts = rows_cnt.at(entry.class_id);
          for (size_t i = 0; entry.cnt_0 + entry.cnt_1 + i <= max_diff; ++i) {
            if (counts[get_index(i, entry.cnt_1)] != 0) {
              is_unique = false;
              break;
            }
          }

          if (is_unique) {
            // the row is removed
            --counts[get_index(entry.cnt_0, entry.cnt_1)];
          } else {
            new_entries.push_back(entry);
          }
        }

        entries = std::move(new_entries);
      }
    }

    // std::vector<size_t> entries_sizes;
    // for (auto& entries : rows) {
    // entries_sizes.push_back(entries.size());
    // }

    size_t total_entries_count = 0;
    for (auto& entries : rows) {
      total_entries_count += entries.size();
    }
    logging::log_value(total_entries_count, "entries_count.csv");
    // std::ranges::sort(entries_sizes);
    // for (size_t i = 0; i < 10; ++i) {
    //   std::println("    {}", entries_sizes[entries_sizes.size() - i - 1]);
    // }

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

  auto seek(const CSCMatrix<double>& matrix) {
    const size_t selected_groups_count = groups_count - max_diff;

    auto [n, d] = matrix.shape();
    transposed_ = matrix.get_transposed();

    //
    std::vector<std::vector<size_t>> groups(groups_count);

    // TODO: maybe implement some algorithm for choosing groups
    for (size_t col = 0; col < d; ++col) {
      groups[std::hash<size_t>()(col) % groups_count].push_back(col);
    }

    auto blocks = get_blocks(matrix, groups);

    std::vector<bool> mask(groups_count, false);
    std::fill_n(mask.begin(), selected_groups_count, true);

    do {
      seek_table(matrix, mask, blocks);
    } while (std::ranges::prev_permutation(mask).found);

    std::println("  finished big rows");

    // remove duplicates from the singular part of the result
    std::ranges::sort(result_singular_);

    std::vector<std::pair<size_t, size_t>> unique;
    std::ranges::unique_copy(result_singular_, std::back_inserter(unique));

    // process all small rows
    seek_small(matrix);

    std::println("  finished small rows");

    return std::pair{std::move(unique), std::move(result_bipartite_)};
  }

  Statistics get_stats() const { return statistics_; }
};

}  // namespace seekers
