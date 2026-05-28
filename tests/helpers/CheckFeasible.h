#pragma once

#include <gtest/gtest.h>

#include "problems/Problem.h"

template <typename Field>
void check_feasible(const Problem<Field>& problem,
                    const std::vector<double>& x) {
  const auto [n, d] = problem.A.shape();

  ASSERT_EQ(x.size(), d);

  std::vector<double> real(n, 0);

  for (size_t col = 0; col < d; ++col) {
    for (const auto [row, value] : problem.A.get_column(col)) {
      real[row] += value * x[col];
    }
  }

  for (size_t row = 0; row < n; ++row) {
    ASSERT_TRUE(problem.rhs_bounds[row].is_inside(real[row], 1e-5));
  }
}
