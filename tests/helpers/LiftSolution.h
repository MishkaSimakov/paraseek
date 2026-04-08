#pragma once

#include <vector>

#include "LinearExpression.h"
#include "Reducer.h"

inline std::vector<double> lift_solution(
    const std::vector<double>& x_new,
    const std::vector<VariableExpression>& mapping) {
  std::vector<double> x_old(mapping.size());

  for (size_t i = 0; i < mapping.size(); ++i) {
    const auto& [variable, expr] = mapping[i];

    x_old[i] = expr.value_at(x_new.at(variable));
  }

  return x_old;
}
