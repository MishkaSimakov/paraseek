#pragma once

#include <vector>

#include "matrix/CSCMatrix.h"
#include "utils/String.h"

namespace seekers {

class SimHash {
  static uint32_t hash32(size_t value) {
    return static_cast<uint64_t>(value + 1) * 1299827;
  }

  static void apply(std::array<int32_t, 32>& counts, uint32_t hash) {
    for (size_t i = 0; i < 32; ++i) {
      counts[i] += (hash >> i & 1) == 0 ? -1 : 1;
    }
  }

  static uint32_t counts_to_hash(const std::array<int32_t, 32>& counts) {
    uint32_t result = 0;

    for (size_t i = 0; i < 32; ++i) {
      if (counts[i] >= 0) {
        result |= 1 << i;
      }
    }

    return result;
  }

  static size_t hamming_distance(uint32_t x, uint32_t y) {
    size_t result = 0;

    for (size_t i = 0; i < 32; ++i) {
      if ((x >> i & 1) != (y >> i & 1)) {
        ++result;
      }
    }

    return result;
  }

 public:
  SimHash() = default;

  std::vector<std::pair<size_t, size_t>> seek(
      const CSCMatrix<double>& matrix,
      const std::vector<std::pair<size_t, size_t>>& bf_result) {
    auto [n, d] = matrix.shape();

    std::vector<std::array<int32_t, 32>> counts(n);
    for (size_t col = 0; col < d; ++col) {
      uint32_t col_hash = hash32(col);

      for (auto [row, _] : matrix.get_column(col)) {
        apply(counts[row], col_hash);
      }
    }

    std::vector<uint32_t> hashes(n);
    for (size_t i = 0; i < n; ++i) {
      hashes[i] = counts_to_hash(counts[i]);
    }

    // stupid solution for now
    std::vector<std::pair<size_t, size_t>> result;
    size_t total_checks = 0;

    // for (size_t i = 0; i < n; ++i) {
    //   for (size_t j = i + 1; j < n; ++j) {
    //     if (hamming_distance(hashes[i], hashes[j]) <= 3) {
    //       ++total_checks;
    //       auto x = matrix.get_row(i);
    //       auto y = matrix.get_row(j);
    //
    //       if (check_rows(x, y, 3)) {
    //         result.emplace_back(i, j);
    //       }
    //     }
    //   }
    // }

    // std::println("total checks = {}", total_checks);

    for (auto [i, j] : bf_result) {
      std::println("{}", hamming_distance(hashes[i], hashes[j]));
    }

    return result;
  }
};

}  // namespace seekers
