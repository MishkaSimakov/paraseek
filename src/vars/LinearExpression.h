#pragma once

// Stores y = alpha * x + beta
template <typename Field>
class LinearExpression {
  Field alpha_;
  Field beta_;

 public:
  explicit LinearExpression(Field alpha = 1, Field beta = 0)
      : alpha_(alpha), beta_(beta) {}

  static LinearExpression id() { return LinearExpression(1, 0); }

  // If `this` is f1(x) and `other` is f2(x), then returns f1(f2(x))
  LinearExpression composed_with(LinearExpression other) const {
    return LinearExpression(alpha_ * other.alpha_,
                            alpha_ * other.beta_ + beta_);
  }

  LinearExpression inversed() const {
    return LinearExpression(1. / alpha_, -beta_ / alpha_);
  }

  Field alpha() const { return alpha_; }

  Field beta() const { return beta_; }

  template<typename T>
  T value_at(T x) const { return alpha_ * x + beta_; }
};
