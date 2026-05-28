#pragma once

#include <chrono>

namespace timing {

using Duration = std::chrono::nanoseconds;

template <typename F>
Duration timeit(F&& function) {
  auto t1 = std::chrono::high_resolution_clock::now();
  function();
  auto t2 = std::chrono::high_resolution_clock::now();

  return duration_cast<Duration>(t2 - t1);
}

class Deadline {
  using Clock = std::chrono::steady_clock;

  Clock::time_point deadline_;

 public:
  explicit Deadline(Clock::time_point deadline) : deadline_(deadline) {}

  static Deadline after(Clock::duration duration) {
    return Deadline(Clock::now() + duration);
  }

  bool is_over() const { return Clock::now() >= deadline_; }
};

}  // namespace timing
