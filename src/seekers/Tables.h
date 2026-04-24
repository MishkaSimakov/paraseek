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
#include "utils/Accumulators.h"
#include "utils/Hashers.h"
#include "utils/Logging.h"
#include "utils/Primes.h"

namespace seekers {

/*
 * If you can look into the seeds of time
 * And say which grain will grow and which will not,
 * Speak, then, to me, who neither beg nor fear
 * Your favors nor your hate.
 */

static ArithmeticMean<double> mean_attempts;

// TODO: cnt_0 and cnt_1 are small, maybe transform them into uint8_t
// or even merge into one uint8_t
// for every entry cnt_0 + cnt_1 <= max_diff
template <typename Field>
struct SmallRowEntry {
  size_t cnt_0 = 0;
  size_t cnt_1 = 0;

  size_t class_id = 0;
  Field front = 0;
};

template <typename Field>
class ClassesStorage {
  struct ClassRecord {
    bool occupied{false};

    size_t parent_class = -1;
    size_t column{0};
    Field value{0};

    std::vector<size_t> counts;
  };

  const size_t max_diff;
  std::vector<ClassRecord> storage_;

  size_t used_classes_count_{0};

  static size_t find_prime_at_least(size_t value) {
    for (const size_t prime : primes::prime_list) {
      if (prime >= value) {
        return prime;
      }
    }

    throw std::runtime_error("value is too big!");
  }

  size_t get_index(size_t cnt_0, size_t cnt_1) const {
    return (max_diff + 1) * cnt_0 + cnt_1 - (cnt_0 - 1) * cnt_0 / 2;
  }

  size_t counts_size() const { return (max_diff + 1) * (max_diff + 2) / 2; }

 public:
  explicit ClassesStorage(size_t max_diff) : max_diff(max_diff) {}

  void occupy_class(size_t class_id) {
    assert(!storage_[class_id].occupied);

    storage_[class_id].occupied = true;
    ++used_classes_count_;
  }

  void try_free_class(size_t class_id) {
    if (!storage_[class_id].occupied) {
      return;
    }

    for (size_t i = 0; i < counts_size(); ++i) {
      if (storage_[class_id].counts[i] != 0) {
        return;
      }
    }

    --used_classes_count_;
    storage_[class_id].occupied = false;
  }

  void try_extend() {
    if (used_classes_count_ * 4 >= storage_.size()) {
      extend(used_classes_count_ * 8);
    }
  }

  void extend(size_t size) {
    if (storage_.size() >= size) {
      return;
    }

    const size_t new_size = find_prime_at_least(size);
    const size_t old_size = storage_.size();

    storage_.resize(new_size);

    for (size_t i = old_size; i < new_size; ++i) {
      storage_[i].counts.resize(counts_size(), 0);
    }
  }

  size_t get_size() const { return storage_.size(); }

  void add_count(SmallRowEntry<Field> entry, size_t delta) {
    storage_[entry.class_id].counts[get_index(entry.cnt_0, entry.cnt_1)] +=
        delta;
  }

  void sub_count(SmallRowEntry<Field> entry, size_t delta) {
    storage_[entry.class_id].counts[get_index(entry.cnt_0, entry.cnt_1)] -=
        delta;
  }

  const size_t& get_count(size_t class_id, size_t cnt_0, size_t cnt_1) const {
    return storage_[class_id].counts[get_index(cnt_0, cnt_1)];
  }

  const size_t& get_count(SmallRowEntry<Field> entry) const {
    return get_count(entry.class_id, entry.cnt_0, entry.cnt_1);
  }

  size_t get_class(size_t prev_class, size_t column, Field value,
                   size_t value_hash) {
    size_t current_class = prev_class;

    const size_t stride =
        tuple_hasher_fn(value_hash, column) % (storage_.size() - 1) + 1;

    size_t total_attempts = 0;

    while (true) {
      current_class = (current_class + stride) % storage_.size();

      ClassRecord& record = storage_[current_class];

      if (!record.occupied) {
        ++used_classes_count_;

        record.occupied = true;

        record.parent_class = prev_class;
        record.column = column;
        record.value = value;

        mean_attempts.record(total_attempts);
        return current_class;
      }

      if (record.parent_class == prev_class && record.column == column &&
          !FieldTraits<Field>::is_nonzero(record.value - value)) {
        mean_attempts.record(total_attempts);
        return current_class;
      }

      ++total_attempts;
    }

    std::unreachable();
  }

  std::vector<std::vector<size_t>> accumulate_counts() {
    for (size_t i = 1; i < storage_.size(); ++i) {
      for (size_t j = 0; j < counts_size(); ++j) {
        storage_[i].counts[j] += storage_[i - 1].counts[j];
      }
    }

    std::vector<std::vector<size_t>> total(max_diff + 1);

    for (size_t i = 0; i <= max_diff; ++i) {
      total[i].resize(max_diff + 1);

      for (size_t j = 0; j <= max_diff; ++j) {
        total[i][j] = storage_.back().counts[get_index(i, j)];
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
          typename Splitter = splitters::GreedySplitter<Field>>
class Tables {
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

    size_t prev = 0;
    for (size_t& value : classes_sizes | std::views::values) {
      const size_t old = value;

      value = prev;
      prev += old;
    }

    // std::vector<std::pair<size_t, size_t>> sizes;
    // for (auto [class_id, size] : classes_sizes) {
    //   sizes.emplace_back(class_id, size);
    // }
    //
    // std::ranges::sort(sizes, {}, [](auto p) {
    //   return -static_cast<double>(p.second);
    // });
    //

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

    ClassesStorage<Field> classes(max_diff);
    classes.extend(n * 3);

    classes.occupy_class(0);

    std::vector<std::vector<SmallRowEntry<Field>>> rows(n);
    for (size_t i = 0; i < n; ++i) {
      if (transposed_[i].size() <= max_size) {
        rows[i].push_back(SmallRowEntry<Field>{
            .cnt_0 = 0,
            .cnt_1 = 0,
            .class_id = 0,
            .front = 0,
        });

        classes.add_count(rows[i].front(), 1);
      }
    }

    // for logging
    std::vector<size_t> total_entries_count;
    size_t current_entries_count = small_rows_cnt;

    for (const size_t col : columns_order) {
      classes.extend(current_entries_count * 3);

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

          for (SmallRowEntry<Field>& entry : entries) {
            classes.sub_count(entry, 1);
            ++entry.cnt_0;

            // otherwise this entry will be removed later
            if (entry.cnt_0 + entry.cnt_1 <= max_diff) {
              classes.add_count(entry, 1);
            }
          }

          size_t erased_cnt =
              std::erase_if(entries, [&](const SmallRowEntry<Field>& entry) {
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
          const SmallRowEntry entry = entries[i];

          classes.sub_count(entry, 1);
          --current_entries_count;

          if (entry.cnt_0 + entry.cnt_1 < max_diff) {
            // class 0 (create new entry)
            {
              SmallRowEntry new_entry = entry;
              ++new_entry.cnt_0;

              entries.push_back(new_entry);
              classes.add_count(new_entry, 1);

              ++current_entries_count;
            }

            // class 1 (create new entry)
            {
              SmallRowEntry new_entry = entries[i];
              ++new_entry.cnt_1;
              new_entry.class_id =
                  classes.get_class(entry.class_id, col, 0, 41);

              classes.add_count(new_entry, 1);
              entries.push_back(new_entry);

              ++current_entries_count;
            }
          }

          // class 2 (update current entry)
          {
            if (entries[i].front == 0) {
              entries[i].front = value;
            }

            const Field normalized = value / entries[i].front;

            const size_t hash = hasher_(normalized);
            entries[i].class_id =
                classes.get_class(entry.class_id, col, normalized, hash);

            classes.add_count(entries[i], 1);
            ++current_entries_count;
          }
        }
      }

      if (!params.entries_reduction) {
        continue;
      }

      for (auto [row, value] : matrix.get_column(col)) {
        size_t erased_cnt =
            std::erase_if(rows[row], [&](const SmallRowEntry<Field>& entry) {
              for (size_t i = 0; entry.cnt_0 + entry.cnt_1 + i <= max_diff;
                   ++i) {
                const size_t count =
                    classes.get_count(entry.class_id, i, entry.cnt_1);

                if ((i != entry.cnt_0 && count > 0) ||
                    (i == entry.cnt_0 && count > 1)) {
                  return false;
                }
              }

              // the entry has unique class => it is removed
              classes.sub_count(entry, 1);
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

      for (size_t j = 0; i + j <= max_diff; ++j) {
        m[i][j].resize(total[i][j]);
      }
    }

    for (size_t row = 0; row < n; ++row) {
      for (const auto entry : rows[row]) {
        classes.sub_count(entry, 1);

        m[entry.cnt_0][entry.cnt_1][classes.get_count(entry)] = {entry.class_id,
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

// /Users/mihailsimakov/Documents/Programs/raspisator/paraseek/cmake-build-release/playground/paraseek_playground_solve_all
// native-profiler-starter: waiting for profiler...
// native-profiler-starter: starting target executable itself...
// 221/240: savsched1
//   size: 295989 x 588093 (nz = 2030025)
//  duration = 439328917ns
// 222/240: s100
//   size: 14733 x 364632 (nz = 1778132)
//  duration = 51048958ns
// 223/240: ns1760995
//   size: 615388 x 633079 (nz = 2469135)
//  duration = 1006750042ns
// 224/240: co-100
//   size: 2187 x 50431 (nz = 1997831)
//  duration = 36972333ns
// 225/240: bab2
//   size: 17245 x 160600 (nz = 2040414)
//  duration = 57086166ns
// 226/240: neos-3402294-bobin
//   size: 591076 x 593772 (nz = 2625756)
//  duration = 1013255083ns
// 227/240: neos-5104907-jarama
//   size: 489818 x 688870 (nz = 2396562)
//  duration = 1221269709ns
// 228/240: ns1644855
//   size: 40698 x 70598 (nz = 2151094)
//  duration = 135638291ns
// 229/240: supportcase22
//   size: 260602 x 267091 (nz = 2488790)
//  duration = 574233417ns
// 230/240: supportcase12
//   size: 166781 x 819983 (nz = 2354804)
//  duration = 253621667ns
// 231/240: roi5alpha10n8
//   size: 4665 x 110813 (nz = 2374887)
//  duration = 37637709ns
// 232/240: supportcase7
//   size: 6532 x 145236 (nz = 2851937)
//  duration = 66970250ns
// 233/240: neos-4647030-tutaki
//   size: 8382 x 18182 (nz = 3958970)
//  duration = 92286750ns
// 234/240: neos-5114902-kasavu
//   size: 961170 x 1416338 (nz = 4946550)
//  duration = 2251279833ns
// 235/240: supportcase19
//   size: 10713 x 1429460 (nz = 4287456)
//  duration = 123960375ns
// 236/240: neos-5052403-cygnet
//   size: 38268 x 71136 (nz = 4936572)
//  duration = 196265916ns
// 237/240: neos-2075418-temuka
//   size: 349602 x 471906 (nz = 7959863)
//  duration = 638160625ns
// 238/240: neos-3402454-bohle
//   size: 2897380 x 2900076 (nz = 11850972)
//  duration = 5079115208ns
// 239/240: square41
//   size: 40160 x 62270 (nz = 13566462)
//  duration = 260913208ns
// 240/240: square47
//   size: 61591 x 95072 (nz = 27329898)
//  duration = 464252583ns
//
// Process finished with exit code 0
