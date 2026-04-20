#pragma once

#include <cassert>
#include <map>
#include <unordered_map>
#include <unordered_set>

#include "DisjointSet.h"
#include "ExpressionDisjointSet.h"
#include "Hamming.h"
#include "matrix/CSCMatrix.h"
#include "problems/Bound.h"
#include "problems/Problem.h"
#include "seekers/Result.h"
#include "utils/Hashers.h"
#include "utils/ZipRows.h"

inline std::vector<std::pair<size_t, size_t>> result_for_reducer(
    const seekers::Result& result) {
  std::unordered_set<std::pair<size_t, size_t>> pairs;

  auto emplace_pair = [&](size_t i, size_t j) {
    if (i > j) {
      pairs.emplace(j, i);
    } else if (i < j) {
      pairs.emplace(i, j);
    }
  };

  for (auto [i, j] : result.singular) {
    emplace_pair(i, j);
  }

  for (const auto& [left, right] : result.bipartite) {
    if (left.empty() || right.empty()) {
      continue;
    }

    if (left.size() < right.size()) {
      const size_t j = left[0];

      for (const size_t i : right) {
        emplace_pair(i, j);
      }
    } else {
      const size_t j = right[0];

      for (const size_t i : left) {
        emplace_pair(i, j);
      }
    }
  }

  return {pairs.begin(), pairs.end()};
}

struct VariableExpression {
  size_t variable;
  LinearExpression<double> expression;
};

class Reducer {
 public:
  Reducer() = default;

  // Returns new problem and expression of old variables in terms of new ones.
  std::pair<Problem<double>, std::vector<VariableExpression>> apply(
      const Problem<double>& problem,
      const std::vector<std::pair<size_t, size_t>>& rows) {
    const auto [n, d] = problem.A.shape();

    const auto transposed = problem.A.get_transposed();

    // cols_ds[d] is a special constant variable
    // if it can be proven that x_i = \beta, then x_i would be in the same
    // set with cols_ds[d]
    auto cols_ds = ExpressionDisjointSet(d + 1);
    auto rows_ds = DisjointSet(n);

    for (const auto [i, j] : rows) {
      if (i == j) {
        continue;
      }

      auto ratio = similarity::hamming_leq(transposed[i], transposed[j], 2);

      if (!ratio) {
        std::println("error!!!!");
        continue;
      }

      // apply transformation only to equalities
      if (!problem.rhs_bounds[i].is_fixed() ||
          !problem.rhs_bounds[j].is_fixed()) {
        continue;
      }

      // subtract vectors
      double rhs_diff =
          *problem.rhs_bounds[i].lower - *ratio * *problem.rhs_bounds[j].lower;

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
    auto result = Problem<double>::with_size(saved_rows_count,
                                             classes_enumeration.size() - 1);
    result.shift = problem.shift;

    // TODO: here problem may be proven to be unfeasible?
    for (size_t i = 0; i < n; ++i) {
      if (saved_rows[i] != -1) {
        result.rhs_bounds[saved_rows[i]] = problem.rhs_bounds[i];
      }
    }

    std::vector<VariableExpression> mapping(d);

    const size_t constant_class = cols_ds.leader(d).first;
    size_t i = 0;

    for (const auto [class_id, index] : classes_enumeration) {
      if (constant_class == class_id) {
        for (const size_t col : classes[index]) {
          if (col == d) {
            continue;
          }

          // x_col = a x_u + b = f1(x_u)
          // x_d = a x_u + b = f2(x_u)
          // then: x_col = f1(f2^{-1}(x_d))
          const auto f1 = cols_ds.leader(col).second;
          const auto f2 = cols_ds.leader(d).second;

          // x_col = value
          const double value = f1.composed_with(f2.inversed()).beta();
          assert(f1.composed_with(f2.inversed()).alpha() == 0);

          mapping[col] = {
              .variable = 0,
              .expression = LinearExpression<double>(0, value),
          };

          if (!problem.bounds[col].is_inside(value)) {
            result.proven_unfeasible = true;
          }

          for (const auto [row, coef] : problem.A.get_column(col)) {
            if (saved_rows[row] != -1) {
              result.rhs_bounds[saved_rows[row]] -= value * coef;
            }
          }

          result.shift += problem.c[col] * value;
        }

        continue;
      }

      std::map<size_t, double> column;

      for (const size_t col : classes[index]) {
        if (col == d) {
          continue;
        }

        const auto [_, expr] = cols_ds.leader(col);

        mapping[col] = {
            .variable = i,
            .expression = expr,
        };

        for (const auto [row, value] : problem.A.get_column(col)) {
          if (saved_rows[row] != -1) {
            column[saved_rows[row]] += value * expr.alpha();
            result.rhs_bounds[saved_rows[row]] -= value * expr.beta();
          }
        }

        result.c[i] += problem.c[col] * expr.alpha();
        result.shift += problem.c[col] * expr.beta();

        result.bounds[i] ^= expr.inversed().value_at(problem.bounds[col]);
      }

      result.A.add_column();
      for (const auto [row, value] : column) {
        result.A.push_to_last_column(row, value);
      }

      ++i;
    }

    return std::pair{result, mapping};
  }
};
