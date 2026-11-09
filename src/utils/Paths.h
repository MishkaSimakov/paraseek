#pragma once

#include <filesystem>
#include <string_view>

#define QQUOTE(X) #X
#define QUOTE(X) QQUOTE(X)

namespace paths {

inline std::filesystem::path log(std::string_view filename) {
  return std::filesystem::path(LOG_PATH) / filename;
}

inline std::filesystem::path problem(std::string_view filename) {
  return std::filesystem::path(PROBLEMS_PATH) / filename;
}

// Generates a prefix for experiments files based on the current time
inline std::string time_prefix() {
  return std::to_string(
      std::chrono::steady_clock::now().time_since_epoch().count());
}

}  // namespace paths
