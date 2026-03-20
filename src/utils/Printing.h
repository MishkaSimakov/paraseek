#pragma once

#include <optional>
#include <print>
#include <unordered_set>
#include <vector>

#include "ZipRows.h"
#include "matrix/CSCMatrix.h"

namespace printing {

void print_result(const CSCMatrix<double>& matrix,
                  const std::vector<std::pair<size_t, size_t>>& result,
                  std::optional<size_t> max_count = std::nullopt);

void print_diff(const CSCMatrix<double>& matrix,
                const std::vector<std::pair<size_t, size_t>>& left,
                const std::vector<std::pair<size_t, size_t>>& right);

void print_small_rows_cnt(const CSCMatrix<double>& matrix, size_t count);

void print_most_frequent(const std::vector<size_t>& classes,
                              std::optional<size_t> max_count = std::nullopt);

}  // namespace printing
