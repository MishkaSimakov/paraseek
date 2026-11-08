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

// x = y
// x = 1
// x = 2y
TEST(ExpressionDisjointTests, small_1) {
  ExpressionDisjointSet set(2);

  set.unite(0, 1, LinearExpression<double>(1, 0));
  set.unite(0, 2, LinearExpression<double>(1, 1));
  const auto result = set.unite(0, 1, LinearExpression<double>(2, 0));

  ASSERT_EQ(result, ExpressionDisjointSet::UniteResult::INFEASIBLE);
}

// x = 2y
// x = 0
// x = 3y
TEST(ExpressionDisjointTests, small_2) {
  ExpressionDisjointSet set(2);

  set.unite(0, 1, LinearExpression<double>(2, 0));
  set.unite(0, 2, LinearExpression<double>(1, 0));
  set.unite(0, 1, LinearExpression<double>(3, 0));
}

// x = y
// x = 1
// x = 2y - 1
TEST(ExpressionDisjointTests, small_3) {
  ExpressionDisjointSet set(2);

  ASSERT_EQ(set.unite(0, 1, LinearExpression<double>(1, 0)),
            ExpressionDisjointSet::UniteResult::OK);

  ASSERT_EQ(set.unite(0, 2, LinearExpression<double>(1, 1)),
            ExpressionDisjointSet::UniteResult::OK);

  ASSERT_EQ(set.unite(0, 1, LinearExpression<double>(2, -1)),
            ExpressionDisjointSet::UniteResult::OK);
}

// x = y
// x = 0
// x = 2y + 1
TEST(ExpressionDisjointTests, small_4) {
  ExpressionDisjointSet set(2);

  ASSERT_EQ(set.unite(0, 1, LinearExpression<double>(1, 0)),
            ExpressionDisjointSet::UniteResult::OK);

  ASSERT_EQ(set.unite(0, 2, LinearExpression<double>(1, 0)),
            ExpressionDisjointSet::UniteResult::OK);

  ASSERT_EQ(set.unite(0, 1, LinearExpression<double>(2, 1)),
            ExpressionDisjointSet::UniteResult::INFEASIBLE);
}

// x = y + 1
// x = 2y
// => x = 2, y = 1
TEST(ExpressionDisjointTests, small_5) {
  ExpressionDisjointSet set(2);

  ASSERT_EQ(set.unite(0, 1, LinearExpression<double>(1, 1)),
            ExpressionDisjointSet::UniteResult::OK);

  ASSERT_EQ(set.unite(0, 1, LinearExpression<double>(2, 0)),
            ExpressionDisjointSet::UniteResult::OK);

  const auto [leader_x, expr_x] = set.leader(0);
  const auto [leader_y, expr_y] = set.leader(1);
  const auto [leader_zero, expr_zero] = set.leader(2);

  ASSERT_EQ(leader_x, leader_zero);
  ASSERT_EQ(leader_y, leader_zero);

  ASSERT_DOUBLE_EQ(expr_x.composed_with(expr_zero.inversed()).value_at(0), 2);
  ASSERT_DOUBLE_EQ(expr_y.composed_with(expr_zero.inversed()).value_at(0), 1);
}

// x = y + 1
// x = 1
// y = 2x - 2
TEST(ExpressionDisjointTests, small_6) {
  ExpressionDisjointSet set(2);

  ASSERT_EQ(set.unite(0, 1, LinearExpression<double>(1, 1)),
            ExpressionDisjointSet::UniteResult::OK);

  ASSERT_EQ(set.unite(0, 2, LinearExpression<double>(1, 1)),
            ExpressionDisjointSet::UniteResult::OK);

  ASSERT_EQ(set.unite(1, 0, LinearExpression<double>(2, -2)),
            ExpressionDisjointSet::UniteResult::OK);
}
