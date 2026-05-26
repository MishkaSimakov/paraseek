#include <cassert>
#include <chrono>
#include <map>
#include <print>

#include "../src/variables/ExpressionDisjointSet.h"
#include "../src/variables/Reducer.h"
#include "problems/ProblemMatrix.h"
#include "problems/ProblemsNames.h"
#include "utils/Printing.h"
#include "variables/BruteForce.h"
#include "variables/Tables.h"

using namespace std::chrono_literals;

std::vector<size_t> get_groups(size_t d, size_t groups_count) {
  // d % groups_count групп размера d / groups_count + 1
  // d - d % groups_count групп размера d / groups_count

  std::vector<size_t> groups(d);
  for (size_t col = 0; col < d; ++col) {
    groups[col] = col / (d / groups_count + 1);
  }

  return groups;
}

int main() {
  const auto& problems = benchmark_set;

  // for (size_t problem_index = 0; problem_index < problems.size();
  //      ++problem_index) {
  //   const auto& problem_name = problems[problem_index];
  //
  //   std::println("{}/{}: {}", problem_index + 1, problems.size(),
  //   problem_name);
  //
  //   auto problem = get_problem(problem_name, true);
  // auto problem = get_problem("neos-3402454-bohle", true);
  auto problem = get_problem("problem1", true);

  const auto [n, d] = problem.A.shape();
  std::println("  size: {} x {} (nz = {})", n, d, problem.A.nonzero_count());

  auto sizes = problem.A.get_rows_sizes();
  auto transposed = problem.A.get_transposed();

  // calculate lsh
  std::vector<std::array<int32_t, 64>> walks(n);

  for (size_t col = 0; col < d; ++col) {
    uint64_t column_hash = static_cast<uint64_t>(col) * 16599994560425689207u;

    for (const auto [row, value] : problem.A.get_column(col)) {
      for (size_t bit = 0; bit < 64; ++bit) {
        walks[row][bit] += ((column_hash >> bit) & 1) == 0 ? 1 : -1;
      }
    }
  }

  std::vector<uint64_t> hashes(n, 0);

  size_t big_rows = 0;

  for (size_t row = 0; row < n; ++row) {
    if (sizes[row] < 10) {
      continue;
    }

    ++big_rows;
    for (size_t bit = 0; bit < 64; ++bit) {
      if (walks[row][bit] > 0) {
        hashes[row] |= static_cast<uint64_t>(1) << bit;
      }
    }
  }

  std::println("big rows count = {}", big_rows);

  // find hashes with diff <= 3
  std::vector<bool> mask(6, false);
  std::fill_n(mask.begin(), 3, true);

  auto groups = get_groups(64, 6);

  std::vector<size_t> order(n);
  std::iota(order.begin(), order.end(), 0);

  std::vector<std::pair<size_t, size_t>> result;

  do {
    uint64_t bitmask = 0;
    for (size_t bit = 0; bit < 64; ++bit) {
      if (mask[groups[bit]]) {
        bitmask |= static_cast<uint64_t>(1) << bit;
      }
    }

    std::cout << std::bitset<64>(bitmask) << std::endl;

    std::ranges::sort(
        order, {}, [&](size_t i) -> uint64_t { return hashes[i] & bitmask; });

    size_t current = 0;
    size_t total = 0;

    while (current < n) {
      size_t next = current + 1;

      while (next < n && (hashes[order[next]] & bitmask) ==
                             (hashes[order[current]] & bitmask)) {
        ++next;
      }

      if (hashes[order[current]] == 0) {
        current = next;
        continue;
      }

      for (size_t i = current; i < next; ++i) {
        for (size_t j = i + 1; j < next; ++j) {
          ++total;

          const size_t row1 = order[i];
          const size_t row2 = order[j];

          if (similarity::hamming_leq(transposed[row1], transposed[row2], 3)) {
            result.emplace_back(row1, row2);
          }
        }
      }

      current = next;
    }

    std::cout << total << std::endl;
  } while (std::ranges::prev_permutation(mask).found);

  auto normalized =
      seekers::normalize_result(seekers::Result{.singular = result});
  std::println("result size = {}", normalized.singular.size());
}
