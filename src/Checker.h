#pragma once

#include "matrix/CSCMatrix.h"
#include "similarity/Hamming.h"

// Checks that there exists alpha != 0 s.t.:
// matrix[i] - alpha * matrix[j] has no more that max_diff nonzero elements
inline bool check_rows(const CSCMatrix<double>& matrix, size_t i, size_t j,
                       size_t max_diff) {
  auto first_row = matrix.get_row(i);
  auto second_row = matrix.get_row(j);

  return similarity::hamming(first_row, second_row) <= max_diff;
}
