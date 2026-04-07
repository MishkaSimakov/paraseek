#pragma once

#include <vector>

class DisjointSet {
  // Для каждой переменной храню ту переменную, через которую она выражается
  // в parent. В expression храню само выражение. Например, если:
  // x_i = a x_j + b, то parent[i] = j, expression[i] = {a, b}
  std::vector<size_t> parent;
  std::vector<std::pair<double, double>> expression;
  std::vector<size_t> size;

  static std::pair<double, double> compose_expr(
      std::pair<double, double> left, std::pair<double, double> right) {
    return {left.first * right.first, left.first * right.second + left.second};
  }

  static std::pair<double, double> inverse_expr(
      std::pair<double, double> expr) {
    return {1. / expr.first, -expr.second / expr.first};
  }

 public:
  explicit DisjointSet(size_t count)
      : parent(count), expression(count), size(count) {
    for (size_t i = 0; i < count; ++i) {
      parent[i] = i;
      expression[i] = {1, 0};
      size[i] = 1;
    }
  }

  // Для переменной v возвращает u - представителя её класса эквивалентности, а
  // также выражение v через u, то есть (a, b) такие, что x_v = a x_u + b
  std::pair<size_t, std::pair<double, double>> leader(size_t v) {
    if (parent[v] == v) {
      return {v, {1, 0}};
    }

    auto [u, expr] = leader(parent[v]);

    parent[v] = u;
    expression[v] = compose_expr(expression[v], expr);

    return {parent[v], expression[v]};
  }

  // Если expression = (a, b), то v = a u + b
  void unite(size_t v, size_t u, std::pair<double, double> expr) {
    auto [leader_v, expr_v] = leader(v);
    auto [leader_u, expr_u] = leader(u);

    // TODO: in this case we may get rid of variable if expressions are
    // different
    if (leader_u == leader_v) {
      return;
    }

    if (size[leader_v] > size[leader_u]) {
      std::swap(leader_v, leader_u);
      std::swap(expr_v, expr_u);
      std::swap(u, v);

      expr = inverse_expr(expr);
    }

    size[leader_u] += size[leader_v];

    parent[leader_v] = leader_u;
    expression[leader_v] =
        compose_expr(compose_expr(inverse_expr(expr_v), expr), expr_u);
  }
};
