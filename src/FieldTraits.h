#pragma once

template <typename Field>
struct FieldTraits;

template <>
struct FieldTraits<double> {
  constexpr static double kEpsilon = 1e-10;

  static bool is_strictly_negative(double value) { return value < -kEpsilon; }
  static bool is_strictly_positive(double value) { return value > kEpsilon; }

  static bool is_nonzero(double value) {
    return is_strictly_negative(value) || is_strictly_positive(value);
  }
};
