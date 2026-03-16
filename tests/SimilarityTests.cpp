#include <gtest/gtest.h>

#include "Checker.h"
#include "similarity/Jaccard.h"
#include "similarity/Union.h"

TEST(SimilarityTests, Hamming1) {
  const SparseVector<double> x = {{0, 1}, {1, 1}, {2, 1}};
  const SparseVector<double> y = {{0, 3}, {1, 3}, {2, 3}};

  ASSERT_TRUE(similarity::hamming(x, y) == 0);
  ASSERT_TRUE(similarity::hamming(y, x) == 0);
}

TEST(SimilarityTests, Hamming2) {
  const SparseVector<double> x = {{0, 1}, {1, 1}, {2, 1}, {3, 3}};
  const SparseVector<double> y = {{0, 3}, {1, 3}, {2, 3}};

  ASSERT_TRUE(similarity::hamming(x, y) == 1);
  ASSERT_TRUE(similarity::hamming(y, x) == 1);
}

TEST(SimilarityTests, Hamming3) {
  const SparseVector<double> x = {{0, 1}, {1, 1}, {2, 1}, {5, 1}};
  const SparseVector<double> y = {{0, 3}, {1, 3}, {3, 2}, {5, 3}};

  ASSERT_TRUE(similarity::hamming(x, y) == 2);
  ASSERT_TRUE(similarity::hamming(y, x) == 2);
}

TEST(SimilarityTests, Hamming4) {
  const SparseVector<double> x = {{0, 1}, {1, 1}};
  const SparseVector<double> y = {{3, 2}, {5, 3}};

  ASSERT_TRUE(similarity::hamming(x, y) == 4);
  ASSERT_TRUE(similarity::hamming(y, x) == 4);
}

TEST(SimilarityTests, Hamming5) {
  const SparseVector<double> x = {};
  const SparseVector<double> y = {{3, 2}, {5, 3}};

  ASSERT_TRUE(similarity::hamming(x, y) == 2);
  ASSERT_TRUE(similarity::hamming(y, x) == 2);
}

TEST(SimilarityTests, Hamming6) {
  const SparseVector<double> x = {{0, 1}, {2, 2}, {3, 3}};
  const SparseVector<double> y = {{0, 1}};

  ASSERT_TRUE(similarity::hamming(x, y) == 2);
  ASSERT_TRUE(similarity::hamming(y, x) == 2);
}

TEST(SimilarityTests, Union1) {
  const SparseVector<double> x = {};
  const SparseVector<double> y = {{0, 1}, {1, 2}, {2, 3}};

  ASSERT_EQ(similarity::union_size(x, y), 3);
  ASSERT_EQ(similarity::union_size(y, x), 3);
}

TEST(SimilarityTests, Union2) {
  const SparseVector<double> x = {{0, 1}};
  const SparseVector<double> y = {{0, 1}, {1, 2}, {2, 3}};

  ASSERT_EQ(similarity::union_size(x, y), 3);
  ASSERT_EQ(similarity::union_size(y, x), 3);
}

TEST(SimilarityTests, Union3) {
  const SparseVector<double> x = {{1, 2}, {3, 4}};
  const SparseVector<double> y = {{0, 1}, {1, 2}, {2, 3}};

  ASSERT_EQ(similarity::union_size(x, y), 4);
  ASSERT_EQ(similarity::union_size(y, x), 4);
}

TEST(SimilarityTests, Union4) {
  const SparseVector<double> x = {{3, 4}, {4, 5}, {5, 6}};
  const SparseVector<double> y = {{0, 1}, {1, 2}, {2, 3}};

  ASSERT_EQ(similarity::union_size(x, y), 6);
  ASSERT_EQ(similarity::union_size(y, x), 6);
}

TEST(SimilarityTests, Jaccard1) {
  const SparseVector<double> x = {{0, 1}, {1, 1}, {2, 1}};
  const SparseVector<double> y = {{1, 1}, {2, 1}, {3, 1}};

  // |A \/ B| = 4, |A /\ B| = 2 => J = 0.5
  ASSERT_DOUBLE_EQ(similarity::jaccard(x, y), 0.5);
  ASSERT_DOUBLE_EQ(similarity::jaccard(y, x), 0.5);
}
