#pragma once

#include <iostream>
#include <ranges>
#include <span>
#include <sstream>
#include <vector>

#include "FieldTraits.h"
#include "Matrix.h"

template <typename Field>
using SparseVector = std::vector<std::pair<size_t, Field>>;

template <typename Field>
class CSCMatrix {
  std::vector<std::pair<size_t, Field>> entries_;
  std::vector<size_t> index_pointers_;

  size_t rows_cnt_;

 public:
  // creates (height, 0) sparse matrix
  explicit CSCMatrix(size_t height)
      : index_pointers_(1, 0), rows_cnt_(height) {}

  explicit CSCMatrix(const Matrix<Field>& matrix)
      : index_pointers_(matrix.get_width() + 1),
        rows_cnt_(matrix.get_height()) {
    auto [n, d] = matrix.shape();

    index_pointers_[0] = 0;
    size_t nonzero_cnt = 0;

    for (size_t col = 0; col < d; ++col) {
      for (size_t row = 0; row < n; ++row) {
        if (FieldTraits<Field>::is_nonzero(matrix[row, col])) {
          entries_.emplace_back(row, matrix[row, col]);
          ++nonzero_cnt;
        }
      }

      index_pointers_[col + 1] = nonzero_cnt;
    }
  }

  CSCMatrix(const std::initializer_list<std::initializer_list<Field>>& values)
      : CSCMatrix(Matrix<Field>(values)) {}

  std::pair<size_t, size_t> shape() const {
    return {rows_cnt_, index_pointers_.size() - 1};
  }

  std::span<const std::pair<size_t, Field>> get_column(size_t col) const {
    return std::span(entries_.begin() + index_pointers_[col],
                     entries_.begin() + index_pointers_[col + 1]);
  }

  std::span<std::pair<size_t, Field>> get_column(size_t col) {
    return std::span(entries_.begin() + index_pointers_[col],
                     entries_.begin() + index_pointers_[col + 1]);
  }

  void add_column(std::span<const Field> dense, size_t shift = 0) {
    size_t height = rows_cnt_;

    if (dense.size() + shift > height) {
      throw std::invalid_argument("Column is too large for the matrix.");
    }

    size_t nonzero_cnt = index_pointers_.back();
    for (size_t i = 0; i < dense.size(); ++i) {
      if (FieldTraits<Field>::is_nonzero(dense[i])) {
        entries_.emplace_back(i + shift, dense[i]);

        ++nonzero_cnt;
      }
    }

    index_pointers_.push_back(nonzero_cnt);
  }

  void add_column(std::span<const std::pair<size_t, Field>> sparse) {
    entries_.append_range(sparse);
    index_pointers_.push_back(entries_.size());
  }

  // adds zero column
  void add_column() { index_pointers_.push_back(index_pointers_.back()); }

  void push_to_last_column(size_t row, const Field& value) {
    if (!FieldTraits<Field>::is_nonzero(value)) {
      return;
    }

    ++index_pointers_.back();

    entries_.emplace_back(row, value);
  }

  size_t nonzero_count() const { return entries_.size(); }

  double density() const {
    auto [n, d] = shape();
    return static_cast<double>(nonzero_count()) / static_cast<double>(n * d);
  }

  void clear() {
    entries_.clear();
    index_pointers_.clear();

    index_pointers_.push_back(0);
  }

  void pop_back_column() {
    index_pointers_.pop_back();
    entries_.resize(index_pointers_.back());
  }

  // Row elements indices are guaranteed to be in ascending order
  SparseVector<Field> get_row(size_t index) const {
    SparseVector<Field> result;
    get_row(index, result);
    return result;
  }

  // Row elements indices are guaranteed to be in ascending order
  void get_row(size_t index, SparseVector<Field>& result) const {
    auto [n, d] = shape();

    if (index > n) {
      throw std::invalid_argument("Row index is outside of bounds.");
    }

    for (size_t i = 0; i < d; ++i) {
      for (auto [row, value] : get_column(i)) {
        if (row == index) {
          result.emplace_back(i, value);
        }
      }
    }
  }
};

template <typename Field>
std::ostream& operator<<(std::ostream& os, const CSCMatrix<Field>& matrix) {
  auto [n, m] = matrix.shape();

  Matrix<std::string> result(n, m, "-");

  for (size_t col = 0; col < m; ++col) {
    for (const auto& [row, value] : matrix.get_column(col)) {
      std::stringstream ss;
      ss << value;
      result[row, col] = ss.str();
    }
  }

  os << result;

  return os;
}
