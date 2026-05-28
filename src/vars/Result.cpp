#include "Result.h"

#include <algorithm>
#include <set>

void seekers::Result::add(size_t i, size_t j) {
  if (i < j) {
    singular.emplace_back(i, j);
  } else {
    singular.emplace_back(j, i);
  }
}

void seekers::Result::purge_singular() {
  std::ranges::sort(singular);

  std::vector<std::pair<size_t, size_t>> unique;
  std::ranges::unique_copy(singular, std::back_inserter(unique));

  singular = std::move(unique);
}

std::set<std::pair<size_t, size_t>> seekers::Result::as_set() const {
  std::set<std::pair<size_t, size_t>> result;

  const auto emplace_pair = [&](size_t i, size_t j) {
    if (i > j) {
      result.emplace(j, i);
    } else if (i < j) {
      result.emplace(i, j);
    }
  };

  for (auto [i, j] : singular) {
    emplace_pair(i, j);
  }

  for (const auto& [left, right] : bipartite) {
    for (size_t i : left) {
      for (size_t j : right) {
        emplace_pair(i, j);
      }
    }
  }

  return result;
}

seekers::SingularResult seekers::Result::to_singular() const {
  const auto set = as_set();

  return SingularResult{
      .singular = {set.begin(), set.end()},
  };
}

std::vector<size_t> seekers::only_one(size_t rows_count, const Result& result) {
  std::vector<size_t> chosen(rows_count, -1);

  for (const auto [i, j] : result.singular) {
    if (chosen[j] == -1) {
      chosen[j] = i;
    } else if (chosen[i] == -1) {
      chosen[i] = j;
    }
  }

  for (const auto& [left, right] : result.bipartite) {
    for (size_t i : left) {
      for (size_t j : right) {
        if (i == j) {
          continue;
        }

        if (chosen[j] == -1) {
          chosen[j] = i;
        } else if (chosen[i] == -1) {
          chosen[i] = j;
        }
      }
    }
  }

  return chosen;
}
