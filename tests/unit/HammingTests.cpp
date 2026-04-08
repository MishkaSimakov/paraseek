#include <gtest/gtest.h>

#include "matrix/CSCMatrix.h"
#include "utils/Hamming.h"

void check_hamming_leq(const SparseVector<double>& left,
                       const SparseVector<double>& right, size_t true_distance,
                       double true_ratio) {
  for (size_t i = 0; i < true_distance; ++i) {
    ASSERT_FALSE(similarity::hamming_leq(left, right, i));
    ASSERT_FALSE(similarity::hamming_leq(right, left, i));
  }

  for (size_t i = true_distance; i < true_distance + 3; ++i) {
    ASSERT_TRUE(similarity::hamming_leq(left, right, i));
    ASSERT_DOUBLE_EQ(*similarity::hamming_leq(left, right, i), true_ratio);

    ASSERT_TRUE(similarity::hamming_leq(right, left, i));
    ASSERT_DOUBLE_EQ(*similarity::hamming_leq(right, left, i), 1. / true_ratio);
  }
}

TEST(HammingTests, Hamming1) {
  const SparseVector<double> x = {{0, 1}, {1, 1}, {2, 1}};
  const SparseVector<double> y = {{0, 3}, {1, 3}, {2, 3}};

  ASSERT_EQ(similarity::hamming(x, y).first, 0);
  ASSERT_EQ(similarity::hamming(y, x).first, 0);

  check_hamming_leq(x, y, 0, 1. / 3.);
}

TEST(HammingTests, Hamming2) {
  const SparseVector<double> x = {{0, 1}, {1, 1}, {2, 1}, {3, 3}};
  const SparseVector<double> y = {{0, 3}, {1, 3}, {2, 3}};

  ASSERT_EQ(similarity::hamming(x, y).first, 1);
  ASSERT_EQ(similarity::hamming(y, x).first, 1);

  check_hamming_leq(x, y, 1, 1. / 3.);
}

TEST(HammingTests, Hamming3) {
  const SparseVector<double> x = {{0, 1}, {1, 1}, {2, 1}, {5, 1}};
  const SparseVector<double> y = {{0, 3}, {1, 3}, {3, 2}, {5, 3}};

  ASSERT_EQ(similarity::hamming(x, y).first, 2);
  ASSERT_EQ(similarity::hamming(y, x).first, 2);

  check_hamming_leq(x, y, 2, 1. / 3.);
}

TEST(HammingTests, Hamming4) {
  const SparseVector<double> x = {{0, 1}, {1, 1}};
  const SparseVector<double> y = {{3, 2}, {5, 3}};

  ASSERT_EQ(similarity::hamming(x, y).first, 4);
  ASSERT_EQ(similarity::hamming(y, x).first, 4);

  check_hamming_leq(x, y, 4, 1.);
}

TEST(HammingTests, Hamming5) {
  const SparseVector<double> x = {};
  const SparseVector<double> y = {{3, 2}, {5, 3}};

  ASSERT_EQ(similarity::hamming(x, y).first, 2);
  ASSERT_EQ(similarity::hamming(y, x).first, 2);

  check_hamming_leq(x, y, 2, 1.);
}

TEST(HammingTests, Hamming6) {
  const SparseVector<double> x = {{0, 1}, {2, 2}, {3, 3}};
  const SparseVector<double> y = {{0, 1}};

  ASSERT_EQ(similarity::hamming(x, y).first, 2);
  ASSERT_EQ(similarity::hamming(y, x).first, 2);

  check_hamming_leq(x, y, 2, 1.);
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
