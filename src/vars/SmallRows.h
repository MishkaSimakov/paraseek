#pragma once

#include <algorithm>
#include <fstream>
#include <functional>
#include <numeric>
#include <string>
#include <vector>

#include "ClassesStorage.h"
#include "Hamming.h"
#include "Result.h"
#include "matrix/CSCMatrix.h"
#include "utils/Paths.h"
#include "utils/Time.h"

namespace seekers {

struct SmallRowsParameters {
  size_t max_row_size;
  size_t max_diff;

  bool entries_reduction = true;
  size_t small_column_limit = 2;

  std::function<void(std::vector<size_t>)> entries_count_logger;
};

struct SmallRowsStatistics {
  // SmallRowEntries considered in seek_small
  size_t entries_considered{0};

  timing::Duration duration;
};

template <typename Field, typename Hasher>
class SmallRows {
  // TODO: cnt_0 and cnt_1 are small, maybe transform them into uint8_t
  // or even merge into one uint8_t
  // for every entry cnt_0 + cnt_1 <= max_diff
  struct SmallRowEntry {
    size_t cnt_0 = 0;
    size_t cnt_1 = 0;

    size_t class_id = 0;
    Field front = 0;
  };

  [[no_unique_address]]
  Hasher hasher_;

  SmallRowsStatistics stats_{};

  Result result_;

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

      result_.bipartite.emplace_back(std::move(left_elements),
                                     std::move(right_elements));

      left_begin = left_end;
    }
  }

  // Calculates final result based on patterns
  void calculate_result(ClassesStorage classes,
                        std::vector<std::vector<SmallRowEntry>> rows,
                        const SmallRowsParameters& params) {
    const auto total = classes.accumulate_counts();

    std::vector<std::vector<std::vector<std::pair<size_t, size_t>>>> m(
        params.max_diff + 1);
    for (size_t i = 0; i <= params.max_diff; ++i) {
      m[i].resize(params.max_diff + 1);

      for (size_t j = 0; j <= params.max_diff; ++j) {
        m[i][j].resize(total[i][j]);
      }
    }

    for (size_t row = 0; row < rows.size(); ++row) {
      for (const auto entry : rows[row]) {
        m[entry.cnt_0][entry.cnt_1][--classes.get_rows_count(
            entry.class_id, entry.cnt_0, entry.cnt_1)] = {entry.class_id, row};
      }
    }

    for (size_t i = 0; i <= params.max_diff; ++i) {
      for (size_t j = i; i + j <= params.max_diff; ++j) {
        for (size_t k = 0; i + j + k <= params.max_diff; ++k) {
          compare_hashes(m[i][k], m[j][k]);
        }
      }
    }
  }

  // If column has few nonzero values, then we can manually check all pairs
  // without increasing entries count.
  ssize_t update_entries_small_column(
      std::span<const std::pair<size_t, Field>> column,
      const std::vector<SparseVector<Field>>& transposed,
      ClassesStorage& classes, std::vector<std::vector<SmallRowEntry>>& rows,
      const SmallRowsParameters& params) {
    ssize_t entries_count_change = 0;

    // compare and add rows manually
    for (size_t i = 0; i < column.size(); ++i) {
      const size_t row1 = column[i].first;

      if (transposed[row1].size() > params.max_row_size) {
        continue;
      }

      for (size_t j = i + 1; j < column.size(); ++j) {
        const size_t row2 = column[j].first;

        if (transposed[row2].size() > params.max_row_size) {
          continue;
        }

        if (similarity::hamming_leq(transposed[row1], transposed[row2],
                                    params.max_diff)) {
          result_.add(row1, row2);
        }
      }
    }

    for (const auto [row, value] : column) {
      auto& entries = rows[row];
      // DO NOT REMOVE THIS: entries.size() does change during loop
      // execution
      const size_t entries_size = entries.size();

      for (size_t i = 0; i < entries_size; ++i) {
        --classes.get_rows_count(entries[i].class_id, entries[i].cnt_0,
                                 entries[i].cnt_1);
        ++entries[i].cnt_0;

        // otherwise this entry will be removed later
        if (entries[i].cnt_0 + entries[i].cnt_1 <= params.max_diff) {
          ++classes.get_rows_count(entries[i].class_id, entries[i].cnt_0,
                                   entries[i].cnt_1);
        }
      }

      entries_count_change -=
          std::erase_if(entries, [&](const SmallRowEntry& entry) {
            return entry.cnt_0 + entry.cnt_1 > params.max_diff;
          });
    }

    return entries_count_change;
  }

  // Updates entries using values from the @column. Returns change in entries
  // count: new_entries_count - old_entries_count
  ssize_t update_entries(std::span<const std::pair<size_t, Field>> column,
                         size_t column_index, ClassesStorage& classes,
                         std::vector<std::vector<SmallRowEntry>>& rows,
                         const SmallRowsParameters& params) {
    ssize_t entries_count_change = 0;

    for (auto [row, value] : column) {
      auto& entries = rows[row];
      // DO NOT REMOVE THIS: entries.size() does change during loop execution
      const size_t entries_size = entries.size();

      for (size_t i = 0; i < entries_size; ++i) {
        ++stats_.entries_considered;

        --classes.get_rows_count(entries[i].class_id, entries[i].cnt_0,
                                 entries[i].cnt_1);
        --entries_count_change;

        if (entries[i].cnt_0 + entries[i].cnt_1 < params.max_diff) {
          // class 0 (create new entry)
          {
            SmallRowEntry new_entry = entries[i];
            ++new_entry.cnt_0;

            entries.push_back(new_entry);
            ++classes.get_rows_count(new_entry.class_id, new_entry.cnt_0,
                                     new_entry.cnt_1);
            ++entries_count_change;
          }
        }

        if (entries[i].cnt_0 + entries[i].cnt_1 < params.max_diff) {
          // class 1 (create new entry)
          {
            SmallRowEntry new_entry = entries[i];
            ++new_entry.cnt_1;
            new_entry.class_id =
                classes.get_class(entries[i].class_id, 0, column_index);

            entries.push_back(new_entry);
            ++classes.get_rows_count(new_entry.class_id, new_entry.cnt_0,
                                     new_entry.cnt_1);
            ++entries_count_change;
          }
        }

        // class 2 (update current entry)
        {
          if (entries[i].front == 0) {
            entries[i].front = value;
          }

          const size_t hash = hasher_(value / entries[i].front);
          entries[i].class_id =
              classes.get_class(entries[i].class_id, hash, column_index);
          ++classes.get_rows_count(entries[i].class_id, entries[i].cnt_0,
                                   entries[i].cnt_1);
          ++entries_count_change;
        }
      }
    }

    return entries_count_change;
  }

  // Filters entries after column is processed. Returns erased entries count.
  size_t filter_entries(std::span<const std::pair<size_t, Field>> column,
                        ClassesStorage& classes,
                        std::vector<std::vector<SmallRowEntry>>& rows,
                        const SmallRowsParameters& params) {
    size_t erased_count = 0;

    for (auto [row, value] : column) {
      erased_count += std::erase_if(rows[row], [&](const SmallRowEntry& entry) {
        for (size_t i = 0; entry.cnt_0 + entry.cnt_1 + i <= params.max_diff;
             ++i) {
          const size_t count =
              classes.get_rows_count(entry.class_id, i, entry.cnt_1);

          if ((i != entry.cnt_0 && count > 0) ||
              (i == entry.cnt_0 && count > 1)) {
            return false;
          }
        }

        // the entry has unique class => it is removed
        --classes.get_rows_count(entry.class_id, entry.cnt_0, entry.cnt_1);
        classes.try_free_class(entry.class_id);

        return true;
      });
    }

    return erased_count;
  }

  void seek_impl(const CSCMatrix<Field>& matrix,
                 const SmallRowsParameters& params) {
    const auto [n, d] = matrix.shape();
    const auto transposed = matrix.get_transposed();

    // sort columns by size
    std::vector<size_t> columns_order(d);
    std::iota(columns_order.begin(), columns_order.end(), 0);

    std::ranges::sort(columns_order, {}, [&](size_t col) {
      return matrix.get_column(col).size();
    });

    size_t small_rows_cnt = 0;
    for (size_t i = 0; i < n; ++i) {
      if (transposed[i].size() <= params.max_row_size) {
        ++small_rows_cnt;
      }
    }

    ClassesStorage classes(params.max_diff, small_rows_cnt);

    size_t init_class = classes.pop_free_class();
    classes.get_rows_count(init_class, 0, 0) = small_rows_cnt;

    std::vector<std::vector<SmallRowEntry>> rows(n);
    for (size_t i = 0; i < n; ++i) {
      if (transposed[i].size() <= params.max_row_size) {
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
      if (params.entries_count_logger) {
        total_entries_count.push_back(current_entries_count);
      }

      size_t small_rows_in_column_cnt = 0;
      for (const auto [row, value] : matrix.get_column(col)) {
        if (transposed[row].size() <= params.max_row_size) {
          ++small_rows_in_column_cnt;
        }
      }

      if (small_rows_in_column_cnt <= params.small_column_limit &&
          params.entries_reduction) {
        current_entries_count += update_entries_small_column(
            matrix.get_column(col), transposed, classes, rows, params);
      } else {
        current_entries_count +=
            update_entries(matrix.get_column(col), col, classes, rows, params);
      }

      if (params.entries_reduction) {
        current_entries_count -=
            filter_entries(matrix.get_column(col), classes, rows, params);
      }
    }

    if (params.entries_count_logger) {
      params.entries_count_logger(std::move(total_entries_count));
    }

    calculate_result(std::move(classes), std::move(rows), params);
  }

 public:
  static std::pair<Result, SmallRowsStatistics> seek(
      const CSCMatrix<Field>& matrix, const SmallRowsParameters& params) {
    auto seeker = SmallRows();

    seeker.stats_.duration =
        timing::timeit([&] { seeker.seek_impl(matrix, params); });

    return {std::move(seeker.result_), seeker.stats_};
  }
};

}  // namespace seekers
