#pragma once

#include <tuple>

#include "matrix/CSCMatrix.h"

template <class T>
class SparseZipIterator {
 public:
  using value_type = std::tuple<size_t, T, T>;
  using difference_type = std::ptrdiff_t;
  using iterator_category = std::forward_iterator_tag;

  SparseZipIterator(const SparseVector<T>& x, const SparseVector<T>& y,
                    size_t xi, size_t yi)
      : x_(x), y_(y), xi_(xi), yi_(yi) {}

  value_type operator*() const {
    if (yi_ >= y_.size() ||
        (xi_ < x_.size() && x_[xi_].first < y_[yi_].first)) {
      return {x_[xi_].first, x_[xi_].second, 0};
    }

    if (xi_ >= x_.size() || y_[yi_].first < x_[xi_].first) {
      return {y_[yi_].first, 0, y_[yi_].second};
    }

    return {x_[xi_].first, x_[xi_].second, y_[yi_].second};
  }

  SparseZipIterator& operator++() {
    if (yi_ >= y_.size() ||
        (xi_ < x_.size() && x_[xi_].first < y_[yi_].first)) {
      ++xi_;
    } else if (xi_ >= x_.size() || y_[yi_].first < x_[xi_].first) {
      ++yi_;
    } else {
      ++xi_;
      ++yi_;
    }

    return *this;
  }

  bool operator==(const SparseZipIterator& other) const {
    return xi_ == other.xi_ && yi_ == other.yi_;
  }

  bool operator!=(const SparseZipIterator& other) const {
    return !(*this == other);
  }

 private:
  const SparseVector<T>& x_;
  const SparseVector<T>& y_;

  size_t xi_;
  size_t yi_;
};

template <class T>
class SparseZipRange {
 public:
  SparseZipRange(const SparseVector<T>& x, const SparseVector<T>& y)
      : x_(x), y_(y) {}

  auto begin() const { return SparseZipIterator<T>(x_, y_, 0, 0); }

  auto end() const {
    return SparseZipIterator<T>(x_, y_, x_.size(), y_.size());
  }

 private:
  const SparseVector<T>& x_;
  const SparseVector<T>& y_;
};
