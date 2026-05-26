#pragma once

#include <algorithm>
#include <numeric>
#include <unordered_map>

#include "Hamming.h"
#include "Result.h"
#include "matrix/CSCMatrix.h"
#include "splitters/Random.h"
#include "utils/Hashers.h"
#include "utils/Time.h"

namespace seekers {

struct BigRowsParameters {
  size_t min_row_size;
  size_t groups_count;
  size_t max_diff;
};

struct BigRowsStatistics {
  size_t pairs_considered{0};

  timing::Duration duration;
};

template <typename Field, typename Hasher,
          typename Splitter = splitters::Random<Field>>
class BigRows {
  struct Block {
    size_t class_id;
    Field front;
  };

  [[no_unique_address]]
  Hasher hasher_;

  // stores indices of rows sorted by size (e.g. sorted_[0] is the smallest row)
  std::vector<size_t> sorted_;

  // index of the first row with size >= min_row_size
  size_t big_rows_start_{};

  const CSCMatrix<Field>& matrix_;
  std::vector<SparseVector<Field>> transposed_;

  // buffers for seek_impl
  std::vector<size_t> merged_classes_;
  std::vector<Field> front_;
  std::unordered_map<size_t, size_t> classes_sizes_;
  std::vector<size_t> indices_;

  BigRowsParameters params_;

  Result result_;

  BigRowsStatistics stats_{};

  // For each row returns a list of blocks that aggregate information in column
  // groups. Uses simple hashing to calculate class_id and allows collisions.
  std::vector<std::vector<Block>> get_blocks(
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

  // seeks row pairs that are parallel in groups where @groups_mask is true, and
  // adds them to the result
  void seek_in_groups(const std::vector<bool>& groups_mask,
                      const std::vector<std::vector<Block>>& blocks) {
    const size_t n = transposed_.size();

    std::fill(front_.begin(), front_.end(), 0);
    classes_sizes_.clear();

    // calculate merged classes for each row
    for (size_t i = big_rows_start_; i < n; ++i) {
      const size_t row = sorted_[i];

      StreamHasher hasher;

      for (size_t group_id = 0; group_id < groups_mask.size(); ++group_id) {
        if (groups_mask[group_id]) {
          if (front_[row] == 0) {
            front_[row] = blocks[row][group_id].front;
          }

          hasher << blocks[row][group_id].class_id;
        }
      }

      merged_classes_[row] = hasher.get_hash();
      ++classes_sizes_[merged_classes_[row]];
    }

    // counting sort
    // total count of big rows is stored in prev after the next loop
    // use stable sort, so that rows are remain sorted by size inside classes
    size_t prev = 0;
    for (size_t& value : classes_sizes_ | std::views::values) {
      const size_t old = value;

      value = prev;
      prev += old;
    }

    for (size_t i = big_rows_start_; i < n; ++i) {
      const size_t row = sorted_[i];
      indices_[classes_sizes_[merged_classes_[row]]++] = row;
    }

    // check all pairs inside group with same hash
    for (size_t i = 0; i < n - big_rows_start_; ++i) {
      for (size_t j = i + 1;
           j < n - big_rows_start_ &&
           merged_classes_[indices_[j]] == merged_classes_[indices_[i]] &&
           transposed_[indices_[j]].size() <=
               transposed_[indices_[i]].size() + params_.max_diff;
           ++j) {
        consider_pair(indices_[i], indices_[j], front_);
      }
    }
  }

  void consider_pair(size_t i, size_t j, const std::vector<Field>& front) {
    ++stats_.pairs_considered;

    if (front[i] != 0) {
      const Field ratio = front[i] / front[j];

      if (similarity::hamming_fixed_ratio_leq(transposed_[i], transposed_[j],
                                              ratio, params_.max_diff)) {
        result_.add(i, j);
      }
    } else {
      if (similarity::hamming_leq(transposed_[i], transposed_[j],
                                  params_.max_diff)) {
        result_.add(i, j);
      }
    }
  }

  void seek_impl() {
    const size_t selected_groups_count =
        params_.groups_count - params_.max_diff;

    const size_t n = transposed_.size();

    std::iota(sorted_.begin(), sorted_.end(), 0);
    std::ranges::sort(sorted_, {},
                      [&](size_t row) { return transposed_[row].size(); });

    big_rows_start_ = n;

    for (size_t i = 0; i < n; ++i) {
      if (transposed_[sorted_[i]].size() >= params_.min_row_size) {
        big_rows_start_ = i;
        break;
      }
    }

    //
    const auto groups = Splitter().split(matrix_, params_.groups_count);
    const auto blocks = get_blocks(matrix_, groups);

    std::vector<bool> mask(params_.groups_count, false);
    std::fill_n(mask.begin(), selected_groups_count, true);

    do {
      seek_in_groups(mask, blocks);
    } while (std::ranges::prev_permutation(mask).found);
  }

  BigRows(const CSCMatrix<Field>& matrix, const BigRowsParameters& params)
      : sorted_(matrix.shape().first),
        matrix_(matrix),
        transposed_(matrix.get_transposed()),
        merged_classes_(matrix.shape().first),
        front_(matrix.shape().first),
        indices_(matrix.shape().first),
        params_(params) {}

 public:
  static std::pair<Result, BigRowsStatistics> seek(
      const CSCMatrix<Field>& matrix, const BigRowsParameters& params) {
    BigRows seeker(matrix, params);

    seeker.stats_.duration = timing::timeit([&] { seeker.seek_impl(); });

    return {std::move(seeker.result_), seeker.stats_};
  }
};

}  // namespace seekers
