#pragma once

#include <Highs.h>

#include <cassert>
#include <stdexcept>
#include <string>
#include <vector>

#include "Problem.h"

namespace mps {

namespace detail {

template <typename Field>
std::optional<Field> highs_lower(double v) {
  return v <= -kHighsInf ? std::nullopt
                         : std::optional<Field>(static_cast<Field>(v));
}

template <typename Field>
std::optional<Field> highs_upper(double v) {
  return v >= kHighsInf ? std::nullopt
                        : std::optional<Field>(static_cast<Field>(v));
}

}  // namespace detail

template <typename Field>
Problem<Field> read(const std::string& path) {
  Highs highs;
  highs.setOptionValue("output_flag", false);

  if (highs.readModel(path) != HighsStatus::kOk) {
    throw std::runtime_error("HiGHS failed to read: " + path);
  }

  const HighsLp& lp = highs.getLp();

  assert(lp.a_matrix_.isColwise());

  const size_t n = lp.num_row_;
  const size_t d = lp.num_col_;

  CSCMatrix<Field> A(static_cast<size_t>(n));
  for (size_t col = 0; col < d; ++col) {
    const size_t start = lp.a_matrix_.start_[col];
    const size_t end = lp.a_matrix_.start_[col + 1];

    std::vector<std::pair<size_t, Field>> column;
    column.reserve(end - start);
    for (size_t i = start; i < end; ++i) {
      column.emplace_back(static_cast<size_t>(lp.a_matrix_.index_[i]),
                          static_cast<Field>(lp.a_matrix_.value_[i]));
    }
    A.add_column(column);
  }

  std::vector<Bound<Field>> rhs_bounds(n);
  for (size_t i = 0; i < n; ++i) {
    rhs_bounds[i] = {detail::highs_lower<Field>(lp.row_lower_[i]),
                     detail::highs_upper<Field>(lp.row_upper_[i])};
  }

  const bool maximize = lp.sense_ == ObjSense::kMaximize;
  const Field sign = maximize ? Field(-1) : Field(1);

  std::vector<Field> c(d);
  for (int j = 0; j < d; ++j) {
    c[j] = sign * static_cast<Field>(lp.col_cost_[j]);
  }
  const Field shift = sign * static_cast<Field>(lp.offset_);

  std::vector<Bound<Field>> bounds(d);
  for (int j = 0; j < d; ++j) {
    bounds[j] = {detail::highs_lower<Field>(lp.col_lower_[j]),
                 detail::highs_upper<Field>(lp.col_upper_[j])};
  }

  std::vector<bool> is_integer(d, false);

  for (size_t j = 0; j < d; ++j) {
    if (!lp.integrality_.empty()) {
      if (lp.integrality_[j] != HighsVarType::kContinuous &&
          lp.integrality_[j] != HighsVarType::kInteger) {
        throw std::runtime_error("Unsupported type of integrality.");
      }

      is_integer[j] = lp.integrality_[j] == HighsVarType::kInteger;
    }
  }

  return Problem<Field>(std::move(A), std::move(rhs_bounds), std::move(c),
                        std::move(bounds), std::move(is_integer), shift);
}

}  // namespace mps
