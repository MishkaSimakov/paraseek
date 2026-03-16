#pragma once

#include <filesystem>
#include <string_view>

namespace paths {

inline std::filesystem::path problem(std::string_view filename) {
  return std::filesystem::path(PROBLEMS_PATH) / filename;
}

}  // namespace paths
