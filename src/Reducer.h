#pragma once

#include <cassert>
#include <map>
#include <unordered_map>

#include "DisjointSet.h"
#include "ExpressionDisjointSet.h"
#include "matrix/CSCMatrix.h"
#include "problems/Bound.h"
#include "problems/Problem.h"
#include "utils/Hamming.h"
#include "utils/ZipRows.h"

struct VariableExpression {
  size_t variable;
  LinearExpression<double> expression;
};

class Reducer {
 public:
  Reducer() = default;

  // Returns new problem and expression of old variables in terms of new ones.
  std::pair<Problem, std::vector<VariableExpression>> apply(
      const Problem& problem,
      const std::vector<std::pair<size_t, size_t>>& rows) {
    const auto [n, d] = problem.A.shape();

    const auto transposed = problem.A.get_transposed();

    auto cols_ds = ExpressionDisjointSet(d);
    auto rows_ds = DisjointSet(n);

    for (const auto [i, j] : rows) {
      auto ratio = similarity::hamming_leq(transposed[i], transposed[j], 2);

      if (!ratio) {
        continue;
      }

      // subtract vectors
      double rhs_diff = problem.b[i] - *ratio * problem.b[j];

      SparseVector<double> diff;
      for (auto [i, x, y] : SparseZipRange{transposed[i], transposed[j]}) {
        const double value = x - *ratio * y;

        if (FieldTraits<double>::is_nonzero(value)) {
          diff.emplace_back(i, value);
        }
      }

      if (diff.size() == 2) {
        // rows_ds.unite(i, j);

        cols_ds.unite(diff[0].first, diff[1].first,
                      LinearExpression<double>(-diff[1].second / diff[0].second,
                                               rhs_diff / diff[0].second));
      }
    }

    // 0 - saved, -1 - removed
    std::vector<size_t> saved_rows(n, -1);
    size_t saved_rows_count = 0;

    for (size_t i = 0; i < n; ++i) {
      size_t leader = rows_ds.leader(i);

      if (saved_rows[leader] == -1) {
        saved_rows[leader] = saved_rows_count;
        ++saved_rows_count;
      }
    }

    std::unordered_map<size_t, size_t> classes_enumeration;
    std::vector<std::vector<size_t>> classes(d);

    for (size_t col = 0; col < d; ++col) {
      auto [leader, expr] = cols_ds.leader(col);

      auto [itr, _] =
          classes_enumeration.emplace(leader, classes_enumeration.size());
      classes[itr->second].push_back(col);
    }

    // TODO: here problem may be proven to be unfeasible?
    std::vector<double> updated_b(saved_rows_count);
    for (size_t i = 0; i < n; ++i) {
      if (saved_rows[i] != -1) {
        updated_b[saved_rows[i]] = problem.b[i];
      }
    }

    CSCMatrix<double> updated_A(saved_rows_count);
    std::vector<double> updated_c(classes_enumeration.size(), 0);
    std::vector<Bound<double>> updated_bounds(classes_enumeration.size());
    double updated_shift = problem.shift;

    std::vector<VariableExpression> mapping(d);

    size_t i = 0;

    for (const size_t index : classes_enumeration | std::views::values) {
      std::map<size_t, double> column;
      double objective_coef = 0;

      for (const size_t col : classes[index]) {
        const auto [_, expr] = cols_ds.leader(col);

        mapping[col] = {
            .variable = i,
            .expression = expr,
        };

        for (const auto [row, value] : problem.A.get_column(col)) {
          if (saved_rows[row] != -1) {
            column[saved_rows[row]] += value * expr.alpha();
            updated_b[saved_rows[row]] -= value * expr.beta();
          }
        }

        objective_coef += problem.c[col] * expr.alpha();
        updated_shift += problem.c[col] * expr.beta();

        updated_bounds[i] ^= expr.inversed().value_at(problem.bounds[col]);
      }

      updated_A.add_column();
      for (const auto [row, value] : column) {
        updated_A.push_to_last_column(row, value);
      }

      updated_c[i] = objective_coef;

      ++i;
    }

    return std::pair{
        Problem(updated_A, updated_b, updated_c, updated_bounds, updated_shift),
        mapping};
  }
};
