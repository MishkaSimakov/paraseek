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

    // cols_ds[d] is a special constant variable
    // if it can be proven that x_i = \beta, then x_i would be in the same
    // set with cols_ds[d]
    auto cols_ds = ExpressionDisjointSet(d + 1);
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

      rows_ds.unite(i, j);
      if (diff.size() == 2) {
        cols_ds.unite(diff[0].first, diff[1].first,
                      LinearExpression<double>(-diff[1].second / diff[0].second,
                                               rhs_diff / diff[0].second));
      } else if (diff.size() == 1) {
        cols_ds.unite(diff[0].first, d,
                      LinearExpression<double>(0, rhs_diff / diff[0].second));
      } else if (diff.size() == 0) {
        // TODO: problem may be proven to be unfeasible here
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
    std::vector<std::vector<size_t>> classes(d + 1);

    for (size_t col = 0; col <= d; ++col) {
      auto [leader, expr] = cols_ds.leader(col);

      auto [itr, _] =
          classes_enumeration.emplace(leader, classes_enumeration.size());
      classes[itr->second].push_back(col);
    }

    // classes_enumeration.size() - 1 because of the constant variables class
    auto result =
        Problem::with_size(saved_rows_count, classes_enumeration.size() - 1);
    result.shift = problem.shift;

    // TODO: here problem may be proven to be unfeasible?
    for (size_t i = 0; i < n; ++i) {
      if (saved_rows[i] != -1) {
        result.b[saved_rows[i]] = problem.b[i];
      }
    }

    std::vector<VariableExpression> mapping(d);

    const size_t constant_class = cols_ds.leader(d).first;
    size_t i = 0;

    for (const auto [class_id, index] : classes_enumeration) {
      std::map<size_t, double> column;
      double objective_coef = 0;
      Bound<double> new_bound;

      for (const size_t col : classes[index]) {
        if (col == d) {
          continue;
        }

        const auto [_, expr] = cols_ds.leader(col);

        mapping[col] = {
            .variable = class_id != constant_class ? i : 0,
            .expression = expr,
        };

        for (const auto [row, value] : problem.A.get_column(col)) {
          if (saved_rows[row] != -1) {
            column[saved_rows[row]] += value * expr.alpha();
            result.b[saved_rows[row]] -= value * expr.beta();
          }
        }

        objective_coef += problem.c[col] * expr.alpha();
        result.shift += problem.c[col] * expr.beta();

        new_bound ^= expr.inversed().value_at(problem.bounds[col]);
      }

      if (constant_class == class_id) {
        // already updated rhs and cost, should check bounds
        for (const size_t col : classes[index]) {
          if (col == d) {
            continue;
          }

          const auto expr = cols_ds.leader(col).second;
          assert(expr.alpha() == 0);

          if (!problem.bounds[col].is_inside(expr.beta())) {
            result.proven_unfeasible = true;
          }
        }

        continue;
      }

      result.bounds[i] = new_bound;
      result.c[i] = objective_coef;

      result.A.add_column();
      for (const auto [row, value] : column) {
        result.A.push_to_last_column(row, value);
      }

      ++i;
    }

    return std::pair{result, mapping};
  }
};
