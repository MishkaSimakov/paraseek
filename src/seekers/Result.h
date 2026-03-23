#pragma once

#include <unordered_set>
#include <vector>

#include "utils/Hashers.h"

namespace seekers {

struct SingularResult {
  std::vector<std::pair<size_t, size_t>> singular;

  std::unordered_set<std::pair<size_t, size_t>> as_set() const {
    return {singular.begin(), singular.end()};
  }
};

struct Result {
  std::vector<std::pair<size_t, size_t>> singular;
  std::vector<std::pair<std::vector<size_t>, std::vector<size_t>>> bipartite;
};

SingularResult normalize_result(const Result& result);

}  // namespace seekers
