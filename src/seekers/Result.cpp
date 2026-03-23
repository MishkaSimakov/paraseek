#include "Result.h"

#include <unordered_set>

seekers::SingularResult seekers::normalize_result(const Result& result) {
  std::unordered_set<std::pair<size_t, size_t>> pairs;

  for (auto [i, j] : result.singular) {
    if (i > j) {
      pairs.emplace(j, i);
    } else {
      pairs.emplace(i, j);
    }
  }

  for (const auto& [left, right] : result.bipartite) {
    for (size_t i : left) {
      for (size_t j : right) {
        if (i == j) {
          continue;
        }

        if (i > j) {
          pairs.emplace(j, i);
        } else {
          pairs.emplace(i, j);
        }
      }
    }
  }

  return SingularResult{
      .singular = {pairs.begin(), pairs.end()},
  };
}
