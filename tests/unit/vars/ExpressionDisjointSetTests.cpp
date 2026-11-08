#include <gtest/gtest.h>

#include "vars/ExpressionDisjointSet.h"

// x = y
// x = y + 5
TEST(ExpressionDisjointTests, infeasible) {
  ExpressionDisjointSet set(2);

  ASSERT_EQ(set.unite(0, 1, LinearExpression<double>(1, 0)),
            ExpressionDisjointSet::UniteResult::OK);

  ASSERT_EQ(set.unite(0, 1, LinearExpression<double>(1, 5)),
            ExpressionDisjointSet::UniteResult::INFEASIBLE);
}

// x = y
// x = 2y
TEST(ExpressionDisjointTests, fixed) {
  ExpressionDisjointSet set(2);

  ASSERT_EQ(set.unite(0, 1, LinearExpression<double>(1, 0)),
            ExpressionDisjointSet::UniteResult::OK);

  ASSERT_EQ(set.unite(0, 1, LinearExpression<double>(2, 0)),
            ExpressionDisjointSet::UniteResult::OK);

  // now x = y = 0
  ASSERT_EQ(set.leader(0).first, set.leader(1).first);
  ASSERT_EQ(set.leader(1).first, set.leader(2).first);
}