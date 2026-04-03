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

}  // namespace timing
