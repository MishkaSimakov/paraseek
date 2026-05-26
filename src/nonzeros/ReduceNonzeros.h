#pragma once

#include <algorithm>
#include <cmath>
#include <numeric>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "Hamming.h"
#include "matrix/CSCMatrix.h"
#include "problems/Problem.h"
#include "splitters/RandomSplitter.h"
#include "utils/Hashers.h"
#include "utils/Logging.h"
#include "vars/Statistics.h"

template <typename Field, typename FieldHasher,
          typename Splitter = splitters::RandomSplitter<Field>>
class ReduceNonzeros {
  struct Block {
    size_t class_id{0};
    Field front{0};
    size_t nz_count{0};
  };

  [[no_unique_address]]
  FieldHasher hasher_;

  [[no_unique_address]]
  Splitter splitter_;

  const size_t groups_count;
  const size_t selected_groups_count;

  // subtrahends_[i] = (j, \alpha) means that row i should be changed to:
  // r_i + \alpha * r_j
  // subtrahends_[i] = (-1, \alpha) means that row i should not be changed
  std::vector<std::pair<size_t, Field>> subtrahends_;

  std::vector<size_t> rows_sizes_;

  // For each row a list of blocks is returned
  std::vector<std::vector<Block>> get_blocks(
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
          ++blocks[row][group_id].nz_count;

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

  struct MergedBlock {
    size_t row;
    size_t class_id;
    Field front;
    size_t nz_count;
  };

  // Returns subspan of merged that contains merged classes info for active rows
  // Uses std::ranges::sort for sorting.
  std::span<MergedBlock> get_merged_classes_sort(
      const std::vector<bool>& groups_mask,
      const std::vector<std::vector<Block>>& blocks,
      const Problem<Field>& problem, std::vector<MergedBlock>& merged) {
    auto [n, d] = problem.A.shape();

    size_t active_count = 0;
    for (size_t row = 0; row < n; ++row) {
      if (subtrahends_[row].first != -1) {
        continue;
      }

      MergedBlock result{
          .row = row,
          .class_id = 0,
          .front = 0,
          .nz_count = 0,
      };

      StreamHasher hasher;

      for (size_t group_id = 0; group_id < groups_mask.size(); ++group_id) {
        if (groups_mask[group_id]) {
          result.nz_count += blocks[row][group_id].nz_count;

          if (result.front == 0) {
            result.front = blocks[row][group_id].front;
          }

          hasher << blocks[row][group_id].class_id;
        }
      }

      if (result.nz_count == 0) {
        continue;
      }

      result.class_id = hasher.get_hash();

      merged[active_count] = result;
      ++active_count;
    }

    auto active_span = std::span{merged.begin(), merged.begin() + active_count};

    std::ranges::sort(active_span, {},
                      [](const MergedBlock& block) { return block.class_id; });

    return active_span;
  }

  void seek_table(const std::vector<bool>& groups_mask,
                  const std::vector<std::vector<Block>>& blocks,
                  const Problem<Field>& problem) {
    auto [n, d] = problem.A.shape();

    std::vector<MergedBlock> merged(n);
    auto active = get_merged_classes_sort(groups_mask, blocks, problem, merged);

    auto class_begin = active.begin();
    auto class_end = active.begin();

    while (class_begin != active.end()) {
      while (class_end != active.end() &&
             class_end->class_id == class_begin->class_id) {
        ++class_end;
      }

      // find optimal row for the current class
      // it should satisfy the following criteria:
      // 1. it must be equality
      // 2. number of non-zeros in unselected groups should be as small as
      // possible

      std::optional<MergedBlock> candidate = std::nullopt;
      for (const MergedBlock& block : std::span{class_begin, class_end}) {
        if (!problem.rhs_bounds[block.row].is_fixed()) {
          continue;
        }

        if (!candidate ||
            rows_sizes_[block.row] - block.nz_count <
                rows_sizes_[candidate->row] - candidate->nz_count) {
          candidate = block;
        }
      }

      if (!candidate) {
        class_begin = class_end;
        continue;
      }

      const MergedBlock subtracted = *candidate;

      for (const MergedBlock& block : std::span{class_begin, class_end}) {
        if (block.row == subtracted.row) {
          continue;
        }

        subtrahends_[block.row] = {subtracted.row,
                                   -block.front / subtracted.front};
      }

      class_begin = class_end;
    }
  }

 public:
  explicit ReduceNonzeros(size_t groups_count, size_t selected_group_count,
                          Splitter splitter = {})
      : splitter_(std::move(splitter)),
        groups_count(groups_count),
        selected_groups_count(selected_group_count) {}

  Problem<Field> apply(Problem<Field> problem) {
    auto [n, d] = problem.A.shape();

    subtrahends_.resize(n, std::pair{-1, 0});
    rows_sizes_ = problem.A.get_rows_sizes();

    auto groups = splitter_.split(problem.A, groups_count);

    auto blocks = get_blocks(problem.A, groups);

    std::vector<bool> mask(groups_count, false);
    std::fill_n(mask.begin(), selected_groups_count, true);

    do {
      seek_table(mask, blocks, problem);
    } while (std::ranges::prev_permutation(mask).found);

    size_t found_count = 0;
    for (size_t i = 0; i < n; ++i) {
      if (subtrahends_[i].first != -1) {
        ++found_count;
      }
    }

    std::println("  found count = {}", found_count);

    std::vector<Field> dense(n, 0);

    // calculate hamming distance for subtracted rows, and discard subtraction
    // if non-zeros count would increase

    // if subtrahend_[i] = (j, \alpha), then
    // intersection[i] holds the number of coefficients that are equal in
    // transposed_[i] and -\alpha * transposed_[j]
    std::vector<size_t> intersection(n, 0);

    for (size_t col = 0; col < d; ++col) {
      for (const auto [row, value] : problem.A.get_column(col)) {
        dense[row] = value;
      }

      for (const auto [row, value] : problem.A.get_column(col)) {
        if (subtrahends_[row].first == -1) {
          continue;
        }

        if (FieldTraits<Field>::is_nonzero(dense[subtrahends_[row].first])) {
          ++intersection[row];
        }
        if (!FieldTraits<Field>::is_nonzero(
                dense[row] +
                subtrahends_[row].second * dense[subtrahends_[row].first])) {
          ++intersection[row];
        }
      }

      for (const auto [row, value] : problem.A.get_column(col)) {
        dense[row] = 0;
      }
    }

    for (size_t row = 0; row < n; ++row) {
      if (subtrahends_[row].first == -1) {
        continue;
      }

      const size_t new_nz_count = rows_sizes_[row] +
                                  rows_sizes_[subtrahends_[row].first] -
                                  intersection[row];

      if (rows_sizes_[row] <= new_nz_count) {
        subtrahends_[row].first = -1;
      }
    }

    // subtract rows
    std::vector<std::vector<size_t>> minuends(n);
    for (size_t row = 0; row < n; ++row) {
      if (subtrahends_[row].first != -1) {
        minuends[subtrahends_[row].first].push_back(row);
      }
    }

    CSCMatrix<Field> new_matrix(n);

    for (size_t col = 0; col < d; ++col) {
      for (const auto [row, value] : problem.A.get_column(col)) {
        dense[row] = value;
      }

      for (const auto [row, value] : problem.A.get_column(col)) {
        for (const size_t minuend : minuends[row]) {
          dense[minuend] += subtrahends_[minuend].second * value;
        }
      }

      new_matrix.add_column();

      for (const auto [row, value] : problem.A.get_column(col)) {
        new_matrix.push_to_last_column(row, dense[row]);
        dense[row] = 0;
      }

      for (const auto [row, value] : problem.A.get_column(col)) {
        for (const size_t minuend : minuends[row]) {
          if (dense[minuend] != 0) {
            new_matrix.push_to_last_column(minuend, dense[minuend]);
            dense[minuend] = 0;
          }
        }
      }
    }

    assert(new_matrix.nonzero_count() <= problem.A.nonzero_count());

    problem.A = std::move(new_matrix);

    auto new_bounds = problem.rhs_bounds;

    for (size_t row = 0; row < n; ++row) {
      if (subtrahends_[row].first != -1) {
        new_bounds[row] += subtrahends_[row].second *
                           problem.rhs_bounds[subtrahends_[row].first];
      }
    }

    problem.rhs_bounds = std::move(new_bounds);

    return std::move(problem);
  }
};
