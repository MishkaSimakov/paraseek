#include <gtest/gtest.h>

#include <random>
#include <utility>

#include "helpers/RandomProblem.h"
#include "helpers/Rational.h"
#include "problems/ArchivedProblem.h"
#include "problems/ProblemsNames.h"
#include "utils/Printing.h"
#include "vars/AllRows.h"
#include "vars/BruteForce.h"
#include "vars/Result.h"

void compare_all_with_brute_force(const CSCMatrix<Rational>& matrix,
                                  size_t max_diff) {
  seekers::AllRowsParameters params{
      .max_diff = max_diff,
      .groups_count = 2 * max_diff,
      .threshold = 10,
      .entries_reduction = true,
  };

  const auto expected = seekers::BruteForce<Rational>(max_diff).seek(matrix);
  const auto result =
      seekers::AllRows<Rational, RationalHasher>::seek(matrix, params).first;

  ASSERT_EQ(result.as_set(), std::set(expected.begin(), expected.end()));
}

TEST(AllRowsTests, CompareWithBruteForce1) {
  for (size_t i = 0; i < 50 && i < benchmark_set.size(); ++i) {
    SCOPED_TRACE(std::format("problem name: {}", benchmark_set[i]));

    auto matrix = get_archived<Rational>(benchmark_set[i]).A;
    ASSERT_NO_FATAL_FAILURE(compare_all_with_brute_force(matrix, 1));
  }
}

TEST(AllRowsTests, CompareWithBruteForce2) {
  for (size_t i = 0; i < 50 && i < benchmark_set.size(); ++i) {
    SCOPED_TRACE(std::format("problem name: {}", benchmark_set[i]));

    auto matrix = get_archived<Rational>(benchmark_set[i]).A;
    ASSERT_NO_FATAL_FAILURE(compare_all_with_brute_force(matrix, 2));
  }
}

TEST(AllRowsTests, CompareWithBruteForce3) {
  for (size_t i = 0; i < 50 && i < benchmark_set.size(); ++i) {
    SCOPED_TRACE(std::format("problem name: {}", benchmark_set[i]));

    auto matrix = get_archived<Rational>(benchmark_set[i]).A;
    ASSERT_NO_FATAL_FAILURE(compare_all_with_brute_force(matrix, 3));
  }
}

TEST(AllRowsTests, RandomizedTest1) {
  constexpr size_t rows_size = 4;

  std::default_random_engine engine;
  std::uniform_int_distribution<int> elements_distribution(-5, 5);

  for (size_t i = 0; i < 10000; ++i) {
    CSCMatrix<Rational> matrix(2);

    for (size_t j = 0; j < rows_size; ++j) {
      matrix.add_column();

      matrix.push_to_last_column(0, elements_distribution(engine));
      matrix.push_to_last_column(1, elements_distribution(engine));
    }

    seekers::AllRowsParameters params{
        .max_diff = 2,
        .groups_count = 4,
        .threshold = 2,
        .entries_reduction = true,
    };

    auto result =
        seekers::AllRows<Rational, RationalHasher>::seek(matrix, params)
            .first.to_singular();

    if (similarity::hamming(matrix.get_row(0), matrix.get_row(1)).first <= 2) {
      ASSERT_EQ(result.singular.size(), 1);
      ASSERT_EQ(result.singular[0], (std::pair{0, 1}));
    } else {
      ASSERT_EQ(result.singular.size(), 0);
    }
  }
}

TEST(AllRowsTests, RandomizedTest2) {
  std::default_random_engine random;

  for (size_t i = 0; i < 100'000; ++i) {
    SCOPED_TRACE(std::format("i = {}", i));

    auto problem = generate_random_problem<Rational>(5, 10, random);

    ASSERT_NO_FATAL_FAILURE(compare_all_with_brute_force(problem.A, 2));
  }
}

TEST(AllRowsTests, DifferentRowsSizes1) {
  for (size_t row_size = 0; row_size < 100; ++row_size) {
    CSCMatrix<Rational> matrix(2);

    for (size_t j = 0; j < row_size; ++j) {
      matrix.add_column();

      matrix.push_to_last_column(0, j);
      matrix.push_to_last_column(1, 2 * j);
    }

    seekers::AllRowsParameters params{
        .max_diff = 2,
        .groups_count = 4,
        .threshold = 10,
        .entries_reduction = true,
    };

    auto result =
        seekers::AllRows<Rational, RationalHasher>::seek(matrix, params)
            .first.to_singular();

    ASSERT_EQ(result.singular.size(), 1);
    ASSERT_EQ(result.singular[0], (std::pair{0, 1}));
  }
}

TEST(AllRowsTests, DifferentRowsSizes2) {
  for (size_t row_size = 0; row_size < 100; ++row_size) {
    CSCMatrix<Rational> matrix(2);

    // generate two rows with difference = 2
    for (size_t j = 0; j < row_size; ++j) {
      matrix.add_column();

      if (j >= 2) {
        matrix.push_to_last_column(0, j);
      }
      matrix.push_to_last_column(1, 2 * j);
    }

    seekers::AllRowsParameters params{
        .max_diff = 2,
        .groups_count = 4,
        .threshold = 10,
        .entries_reduction = true,
    };

    auto result =
        seekers::AllRows<Rational, RationalHasher>::seek(matrix, params)
            .first.to_singular();

    ASSERT_EQ(result.singular.size(), 1);
    ASSERT_EQ(result.singular[0], (std::pair{0, 1}));
  }
}

TEST(AllRowsTests, DifferentRowsSizes3) {
  for (size_t row_size = 0; row_size < 100; ++row_size) {
    CSCMatrix<Rational> matrix(2);

    // generate two rows with difference = 0
    for (size_t j = 0; j < row_size; ++j) {
      matrix.add_column();

      matrix.push_to_last_column(0, j);
      matrix.push_to_last_column(1, 2 * j);
    }

    seekers::AllRowsParameters params{
        .max_diff = 0,
        .groups_count = 4,
        .threshold = 10,
        .entries_reduction = true,
    };

    auto result =
        seekers::AllRows<Rational, RationalHasher>::seek(matrix, params)
            .first.to_singular();

    ASSERT_EQ(result.singular.size(), 1);
    ASSERT_EQ(result.singular[0], (std::pair{0, 1}));
  }
}
