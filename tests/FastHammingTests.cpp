#include <gtest/gtest.h>

#include "matrix/CSCMatrix.h"
#include "similarity/Hamming.h"

TEST(FastHammingTests, SmallTest1) {
  const SparseVector<double> x = {{0, 1}, {1, 1}, {2, 1}};
  const SparseVector<double> y = {{0, 1}, {1, 1}, {2, 1}, {4, 1}};

  ASSERT_EQ(similarity::fast_hamming(x, y), 1);
}

TEST(FastHammingTests, SmallTest2) {
  const SparseVector<double> x = {{0, 1}, {1, 1}, {2, 1}, {5, 6}};
  const SparseVector<double> y = {{0, 1}, {1, 1}, {2, 1}, {4, 1}};

  ASSERT_EQ(similarity::fast_hamming(x, y), 2);
}