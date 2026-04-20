#include <cassert>
#include <chrono>
#include <map>
#include <print>

#include "ExpressionDisjointSet.h"
#include "Reducer.h"
#include "problems/ProblemMatrix.h"
#include "problems/ProblemsNames.h"
#include "seekers/BruteForce.h"
#include "seekers/Tables.h"
#include "utils/Printing.h"

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
  auto problem = get_problem("app1-2", true);

  const auto [n, d] = problem.A.shape();
  std::println("  size: {} x {} (nz = {})", n, d, problem.A.nonzero_count());

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

  for (size_t row = 0; row < n; ++row) {
    for (size_t bit = 0; bit < 64; ++bit) {
      if (walks[row][bit] > 0) {
        hashes[row] |= static_cast<uint64_t>(1) << bit;
      }
    }
  }

  // find hashes with diff <= 3
  std::vector<bool> mask(6, false);
  std::fill_n(mask.begin(), 3, true);

  auto groups = get_groups(64, 6);

  do {
    uint64_t bitmask = 0;
    for (size_t bit = 0; bit < 64; ++bit) {
      if (mask[groups[bit]]) {
        bitmask |= static_cast<uint64_t>(1) << bit;
      }
    }

    std::cout << std::bitset<64>(bitmask) << std::endl;

    std::ranges::sort(hashes, {}, [bitmask](uint64_t hash) -> uint64_t {
      return hash & bitmask;
    });

    size_t current = 0;
    size_t total = 0;

    while (current < n) {
      size_t next = current + 1;

      while (next < n &&
             (hashes[next] & bitmask) == (hashes[current] & bitmask)) {
        ++next;
      }

      if (next > current + 10) {
        total += (next - current) * (next - current);
      }
      current = next;
    }

    std::cout << total << std::endl;
  } while (std::ranges::prev_permutation(mask).found);

  // }
}
