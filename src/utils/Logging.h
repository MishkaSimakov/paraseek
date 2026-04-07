#pragma once

#include <filesystem>
#include <fstream>
#include <memory>
#include <mutex>
#include <string_view>
#include <unordered_map>

#include "Paths.h"

namespace logging {

inline std::ofstream& get_log_fstream(std::string_view filename) {
  static std::unordered_map<std::string, std::unique_ptr<std::ofstream>>
      streams;
  static std::mutex mutex;

  std::lock_guard lock(mutex);

  auto it = streams.find(std::string(filename));
  if (it != streams.end()) {
    return *(it->second);
  }

  auto path = paths::log(filename);
  std::filesystem::create_directories(path.parent_path());

  auto os = std::make_unique<std::ofstream>(path, std::ofstream::app);
  if (!*os) {
    throw std::runtime_error("Failed to write to log file.");
  }

  auto& ref = *os;
  streams.emplace(std::string(filename), std::move(os));
  return ref;
}

template <typename T>
void log_value(const T& value, std::string_view filename) {
  get_log_fstream(filename) << value << "\n";
}

inline void log(std::string_view text, std::string_view filename) {
  get_log_fstream(filename) << text;
}

// --- CSV helpers ---

inline bool is_file_empty(const std::filesystem::path& path) {
  return !std::filesystem::exists(path) ||
         std::filesystem::file_size(path) == 0;
}

inline std::string escape_csv(std::string value) {
  bool need_quotes = value.find_first_of(",\"\n") != std::string::npos;

  if (value.find('"') != std::string::npos) {
    std::string escaped;
    escaped.reserve(value.size());
    for (char c : value) {
      if (c == '"')
        escaped += "\"\"";
      else
        escaped += c;
    }
    value = std::move(escaped);
  }

  if (need_quotes) {
    return "\"" + value + "\"";
  }

  return value;
}

template <typename T>
std::string to_string_csv(const T& value) {
  std::ostringstream oss;
  oss << value;
  return escape_csv(oss.str());
}

// --- main CSV logger ---

template <typename Field>
void log_csv(const std::vector<std::pair<std::string, Field>>& fields,
             std::string_view filename) {
  static std::unordered_map<std::string, bool> header_written;
  static std::mutex mutex;

  std::lock_guard lock(mutex);

  auto path = paths::log(filename);
  bool need_header = false;

  auto key = std::string(filename);

  // Check cached state first
  auto it = header_written.find(key);
  if (it == header_written.end()) {
    need_header = is_file_empty(path);
    header_written[key] = true;  // mark as initialized
  }

  auto& os = get_log_fstream(filename);

  // Write header if needed
  if (need_header) {
    for (size_t i = 0; i < fields.size(); ++i) {
      os << escape_csv(fields[i].first);
      if (i + 1 < fields.size()) os << ",";
    }
    os << "\n";
  }

  // Write values
  for (size_t i = 0; i < fields.size(); ++i) {
    os << to_string_csv(fields[i].second);
    if (i + 1 < fields.size()) os << ",";
  }
  os << "\n";

  if (!os) {
    throw std::runtime_error("Error occurred while logging");
  }
}

}  // namespace logging
