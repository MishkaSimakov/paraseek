#pragma once

#include <algorithm>
#include <map>
#include <ranges>
#include <unordered_map>

#include "Splitter.h"
#include "matrix/CSCMatrix.h"

namespace splitters {

inline int64_t normalize_double(double value) {
  return static_cast<int64_t>(std::round(value * 1e10));
}

struct Block {
  size_t class_id;
  double front;
};

struct EvaluationResult {
  size_t score;
  bool is_valid;
};

EvaluationResult evaluate(const CSCMatrix<double>& matrix,
                          const std::vector<std::vector<size_t>>& groups,
                          size_t max_diff) {
  const auto [n, d] = matrix.shape();

  size_t score = 0;
  std::unordered_map<std::pair<size_t, int64_t>, size_t> map;

  std::vector<size_t> rows_order(n);
  std::iota(rows_order.begin(), rows_order.end(), 0);

  std::vector<bool> mask(groups.size(), false);
  std::fill_n(mask.begin(), groups.size() - max_diff, true);

  do {
    std::vector<Block> blocks(n);
    size_t classes_cnt = 1;

    for (size_t group_id = 0; group_id < groups.size(); ++group_id) {
      if (!mask[group_id]) {
        continue;
      }

      for (const size_t col : groups[group_id]) {
        map.clear();

        for (auto [row, value] : matrix.get_column(col)) {
          if (blocks[row].front == 0) {
            blocks[row].front = value;
          }

          auto normalized = normalize_double(value / blocks[row].front);

          auto [itr, new_class] = map.emplace(
              std::pair{blocks[row].class_id, normalized}, classes_cnt);

          if (new_class) {
            ++classes_cnt;
          }

          blocks[row].class_id = itr->second;
        }
      }
    }

    std::map<size_t, size_t> classes_sizes;
    for (const Block& block : blocks) {
      ++classes_sizes[block.class_id];
    }

    // std::vector<std::pair<size_t, size_t>> by_size;
    // for (auto [cls, size] : classes_sizes) {
    //   by_size.emplace_back(cls, size);
    // }
    //
    // std::ranges::sort(by_size, {}, [](auto p) { return p.second; });
    //
    // std::println("  largest classes:");
    // for (size_t i = 0; i < 5; ++i) {
    //   const auto p = by_size[by_size.size() - i - 1];
    //   std::println("    {} (size = {})", p.first, p.second);
    // }

    for (size_t value : classes_sizes | std::views::values) {
      if (value > 1) {
        score += value * value;
      }
    }
  } while (std::ranges::prev_permutation(mask).found);

  // check correctness
  std::vector<bool> columns(d, false);

  for (size_t group_id = 0; group_id < groups.size(); ++group_id) {
    for (const size_t col : groups[group_id]) {
      if (columns[col]) {
        return {.score = 0, .is_valid = false};
      }

      columns[col] = true;
    }
  }

  for (bool col : columns) {
    if (!col) {
      return {.score = 0, .is_valid = false};
    }
  }

  return {
      .score = score,
      .is_valid = true,
  };
}

}  // namespace splitters
