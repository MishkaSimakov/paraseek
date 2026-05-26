#pragma once

#include <filesystem>
#include <string_view>

#define QQUOTE(X) #X
#define QUOTE(X) QQUOTE(X)

namespace paths {

inline std::filesystem::path log(std::string_view filename) {
  return std::filesystem::path(QUOTE(LOG_PATH)) / filename;
}

inline std::filesystem::path problem(std::string_view filename) {
  return std::filesystem::path(QUOTE(PROBLEMS_PATH)) / filename;
}

}  // namespace paths
