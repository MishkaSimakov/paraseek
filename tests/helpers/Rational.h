#pragma once

#include <boost/multiprecision/cpp_int.hpp>

#include "FieldTraits.h"
#include "utils/Hashers.h"

using Rational = boost::multiprecision::cpp_rational;

struct RationalHasher {
  size_t operator()(const Rational& value) const {
    StreamHasher hasher;

    hasher << hash_value(numerator(value));
    hasher << hash_value(denominator(value));

    return hasher.get_hash();
  }
};

template <>
struct FieldTraits<Rational> {
  static bool is_strictly_positive(const Rational& value) { return value > 0; }
  static bool is_strictly_negative(const Rational& value) { return value < 0; }

  static bool is_nonzero(const Rational& value) { return value != 0; }
};
