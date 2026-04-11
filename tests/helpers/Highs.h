#pragma once

#include <Highs.h>

#include <vector>

#include "problems/Problem.h"

struct Solution {
  HighsModelStatus status;
  std::vector<double> x;
  double objective;
};

template <std::convertible_to<double> Field>
HighsLp to_highs(const Problem<Field>& problem) {
  const auto [n, d] = problem.A.shape();

  HighsLp lp;
  lp.num_col_ = static_cast<int>(d);
  lp.num_row_ = static_cast<int>(n);

  // Objective
  lp.col_cost_.resize(d);
  for (size_t i = 0; i < d; ++i) {
    lp.col_cost_[i] = static_cast<double>(problem.c[i]);
  }

  // Bounds
  std::vector<double> lower(d);
  std::vector<double> upper(d);

  for (size_t i = 0; i < d; ++i) {
    lower[i] =
        problem.bounds[i]
            .lower
            .transform([](Field value) { return static_cast<double>(value); })
            .value_or(-kHighsInf);

    upper[i] =
        problem.bounds[i]
            .upper
            .transform([](Field value) { return static_cast<double>(value); })
            .value_or(kHighsInf);
  }

  lp.col_lower_ = std::move(lower);
  lp.col_upper_ = std::move(upper);

  // Matrix
  lp.row_lower_.resize(n);
  lp.row_upper_.resize(n);

  for (size_t row = 0; row < n; ++row) {
    lp.row_lower_[row] = lp.row_upper_[row] =
        static_cast<double>(problem.b[row]);
  }

  lp.a_matrix_.start_.resize(d + 1);
  lp.a_matrix_.index_.clear();
  lp.a_matrix_.value_.clear();

  size_t nnz = 0;
  lp.a_matrix_.start_[0] = 0;

  for (size_t col = 0; col < d; ++col) {
    for (auto [row, value] : problem.A.get_column(col)) {
      lp.a_matrix_.index_.push_back(static_cast<int>(row));
      lp.a_matrix_.value_.push_back(static_cast<double>(value));
      ++nnz;
    }

    lp.a_matrix_.start_[col + 1] = static_cast<int>(nnz);
  }

  return lp;
}

template <std::convertible_to<double> Field>
Solution solve(const Problem<Field>& problem) {
  if (problem.proven_unfeasible) {
    return {
        .status = HighsModelStatus::kInfeasible,
        .objective = 0,
        .x = std::vector<double>(problem.A.shape().second, 0),
    };
  }

  Highs highs;
  highs.setOptionValue("output_flag", false);

  auto lp = to_highs(problem);

  highs.passModel(lp);
  highs.writeModel("model.lp");
  highs.run();

  Solution sol;
  sol.status = highs.getModelStatus();

  const auto& highs_sol = highs.getSolution();
  sol.x = highs_sol.col_value;
  sol.objective = highs.getInfo().objective_function_value + problem.shift;

  return sol;
}
