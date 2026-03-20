#include <gtest/gtest.h>

#include "matrix/CSCMatrix.h"
#include "utils/Hamming.h"

TEST(HammingTests, Hamming1) {
  const SparseVector<double> x = {{0, 1}, {1, 1}, {2, 1}};
  const SparseVector<double> y = {{0, 3}, {1, 3}, {2, 3}};

  ASSERT_EQ(similarity::hamming(x, y), 0);
  ASSERT_EQ(similarity::hamming(y, x), 0);
}

TEST(HammingTests, Hamming2) {
  const SparseVector<double> x = {{0, 1}, {1, 1}, {2, 1}, {3, 3}};
  const SparseVector<double> y = {{0, 3}, {1, 3}, {2, 3}};

  ASSERT_TRUE(similarity::hamming(x, y) == 1);
  ASSERT_TRUE(similarity::hamming(y, x) == 1);
}

TEST(HammingTests, Hamming3) {
  const SparseVector<double> x = {{0, 1}, {1, 1}, {2, 1}, {5, 1}};
  const SparseVector<double> y = {{0, 3}, {1, 3}, {3, 2}, {5, 3}};

  ASSERT_TRUE(similarity::hamming(x, y) == 2);
  ASSERT_TRUE(similarity::hamming(y, x) == 2);
}

TEST(HammingTests, Hamming4) {
  const SparseVector<double> x = {{0, 1}, {1, 1}};
  const SparseVector<double> y = {{3, 2}, {5, 3}};

  ASSERT_TRUE(similarity::hamming(x, y) == 4);
  ASSERT_TRUE(similarity::hamming(y, x) == 4);
}

TEST(HammingTests, Hamming5) {
  const SparseVector<double> x = {};
  const SparseVector<double> y = {{3, 2}, {5, 3}};

  ASSERT_TRUE(similarity::hamming(x, y) == 2);
  ASSERT_TRUE(similarity::hamming(y, x) == 2);
}

TEST(HammingTests, Hamming6) {
  const SparseVector<double> x = {{0, 1}, {2, 2}, {3, 3}};
  const SparseVector<double> y = {{0, 1}};

  ASSERT_TRUE(similarity::hamming(x, y) == 2);
  ASSERT_TRUE(similarity::hamming(y, x) == 2);
}

TEST(HammingTests, FastHammingTest1) {
  const SparseVector<double> x = {{0, 1}, {1, 1}, {2, 1}};
  const SparseVector<double> y = {{0, 1}, {1, 1}, {2, 1}, {4, 1}};

  ASSERT_EQ(similarity::fast_hamming(x, y), 1);
}

TEST(HammingTests, FastHammingTest2) {
  const SparseVector<double> x = {{0, 1}, {1, 1}, {2, 1}, {5, 6}};
  const SparseVector<double> y = {{0, 1}, {1, 1}, {2, 1}, {4, 1}};

  ASSERT_EQ(similarity::fast_hamming(x, y), 2);
}