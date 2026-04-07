#pragma once

#include "matrix/CSCMatrix.h"

namespace splitters {

class RandomSplitter {
 public:
  RandomSplitter() = default;

  std::vector<std::vector<size_t>> split(const CSCMatrix<double>& matrix,
                                         size_t groups_count) {
    const auto [n, d] = matrix.shape();
    std::vector<std::vector<size_t>> groups(groups_count);

    for (size_t col = 0; col < d; ++col) {
      groups[std::hash<size_t>()(col) % groups_count].push_back(col);
    }

    return groups;
  }
};

}  // namespace splitters
