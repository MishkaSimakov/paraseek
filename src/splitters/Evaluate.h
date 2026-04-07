#pragma once

#include <algorithm>
#include <map>
#include <ranges>
#include <unordered_map>

#include "Splitter.h"
#include "matrix/CSCMatrix.h"

namespace splitters {

inline int64_t normalize_double(double value) {
  return std::round(value * 1e10);
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
                          const std::vector<std::vector<size_t>>& groups) {
  const auto [n, d] = matrix.shape();
  const auto transposed = matrix.get_transposed();

  size_t score = 0;
  std::unordered_map<std::pair<size_t, int64_t>, size_t> map;

  std::vector<size_t> rows_order(n);
  std::iota(rows_order.begin(), rows_order.end(), 0);

  for (size_t group_id = 0; group_id < groups.size(); ++group_id) {
    std::vector<Block> blocks(n);
    size_t classes_cnt = 1;

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

    // std::map<size_t, size_t> classes_sizes;
    // for (const Block& block : blocks) {
    //   ++classes_sizes[block.class_id];
    // }
    //
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

    // for (size_t value : classes_sizes | std::views::values) {
    //   if (value > 1) {
    //     score += value * value;
    //   }
    // }

    std::ranges::sort(rows_order, {}, [&](size_t row) {
      return std::pair{blocks[row].class_id, transposed[row].size()};
    });

    size_t i_start = 0;
    while (i_start < n) {
      const size_t j_start = i_start;
      size_t j_end = i_start + 1;

      for (; j_end < n &&
             blocks[rows_order[j_end]].class_id ==
                 blocks[rows_order[i_start]].class_id &&
             transposed[rows_order[j_end]].size() <=
                 transposed[rows_order[j_start]].size() + 2;
           ++j_end) {
      }

      size_t i_end = i_start + 1;
      while (j_end < n &&
             blocks[rows_order[i_end]].class_id ==
                 blocks[rows_order[i_start]].class_id &&
             transposed[rows_order[i_start]].size() ==
                 transposed[rows_order[i_end]].size()) {
        ++i_end;
      }

      score += (i_end - i_start - 1) * (j_end - j_start - 1);

      i_start = i_end;
    }
  }

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
