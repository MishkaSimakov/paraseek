#pragma once

#include <vector>

#include "LinearExpression.h"

class DisjointSet {
  std::vector<size_t> parent_;
  std::vector<size_t> size_;

 public:
  explicit DisjointSet(size_t count) : parent_(count), size_(count) {
    for (size_t i = 0; i < count; ++i) {
      parent_[i] = i;
      size_[i] = 1;
    }
  }

  size_t leader(size_t v) {
    if (parent_[v] == v) {
      return v;
    }

    return parent_[v] = leader(parent_[v]);
  }

  void unite(size_t v, size_t u) {
    v = leader(v);
    u = leader(u);

    if (size_[v] > size_[u]) {
      std::swap(v, u);
    }

    size_[u] += size_[v];

    parent_[v] = u;
  }
};
