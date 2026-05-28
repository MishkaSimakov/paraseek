#pragma once

#include <concepts>
#include <vector>

#include "matrix/CSCMatrix.h"

namespace splitters {

template <typename T>
concept Splitter =
    requires(T splitter, CSCMatrix<double> matrix, size_t groups_count) {
      {
        splitter.split(matrix, groups_count)
      } -> std::same_as<std::vector<std::vector<size_t>>>;
    };

}  // namespace splitters
