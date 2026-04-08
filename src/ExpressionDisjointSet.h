#pragma once

#include <vector>

#include "LinearExpression.h"

class ExpressionDisjointSet {
  // Для каждой переменной храню ту переменную, через которую она выражается
  // в parent. В expression храню само выражение. Например, если:
  // x_i = a x_j + b, то parent[i] = j, expression[i] = {a, b}
  std::vector<size_t> parent_;
  std::vector<LinearExpression<double>> expression_;
  std::vector<size_t> size_;

 public:
  explicit ExpressionDisjointSet(size_t count)
      : parent_(count), expression_(count), size_(count) {
    for (size_t i = 0; i < count; ++i) {
      parent_[i] = i;
      expression_[i] = LinearExpression<double>::id();
      size_[i] = 1;
    }
  }

  // Для переменной v возвращает u - представителя её класса эквивалентности, а
  // также выражение v через u, то есть (a, b) такие, что x_v = a x_u + b
  std::pair<size_t, LinearExpression<double>> leader(size_t v) {
    if (parent_[v] == v) {
      return {v, LinearExpression<double>::id()};
    }

    auto [u, expr] = leader(parent_[v]);

    parent_[v] = u;
    expression_[v] = expression_[v].composed_with(expr);

    return {parent_[v], expression_[v]};
  }

  // Если expression = (a, b), то v = a u + b
  void unite(size_t v, size_t u, LinearExpression<double> expr) {
    auto [leader_v, expr_v] = leader(v);
    auto [leader_u, expr_u] = leader(u);

    // TODO: in this case we may get rid of variable if expressions are
    // different
    if (leader_u == leader_v) {
      return;
    }

    if (size_[leader_v] > size_[leader_u]) {
      std::swap(leader_v, leader_u);
      std::swap(expr_v, expr_u);
      std::swap(u, v);

      expr = expr.inversed();
    }

    size_[leader_u] += size_[leader_v];

    parent_[leader_v] = leader_u;
    expression_[leader_v] =
        expr_v.inversed().composed_with(expr).composed_with(expr_u);
  }
};
