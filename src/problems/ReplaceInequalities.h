#pragma once

#include "Problem.h"

// Transforms all inequalities to equalities by adding slack variables.
template <typename Field>
void replace_inequalities(Problem<Field>& problem) {
  const auto [n, d] = problem.A.shape();

  for (size_t i = 0; i < n; ++i) {
    if (problem.rhs_bounds[i].is_fixed()) {
      continue;
    }

    problem.c.push_back(0);
    problem.bounds.push_back(-problem.rhs_bounds[i]);

    problem.rhs_bounds[i] = Bound<Field>{0, 0};

    problem.A.add_column();
    problem.A.push_to_last_column(i, 1);
  }
}
