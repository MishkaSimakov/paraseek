#pragma once

#include <format>

#include "MPS.h"
#include "matrix/CSCMatrix.h"
#include "utils/Paths.h"

// First tries to find a file with the name @name.mps. Then tries to unarchive
// it from @name.mps.gz. Throws if it fails to find these files.
// Returns path to the unarchived file.
inline std::filesystem::path unarchive(const std::string& name) {
  auto archive_path = paths::problem(name + ".mps.gz");
  auto path = paths::problem(name + ".mps");

  if (!std::filesystem::exists(path)) {
    std::string command = std::format("gzip -dk {}", archive_path.string());
    int rs = std::system(command.c_str());

    if (rs != 0) {
      throw std::runtime_error(std::format(
          "Failed to decompress problem file: \"{}\"", archive_path.string()));
    }
  }

  return path;
}

// Reads the problem with the given name. First tries to find file with the name
// @name.mps, if this attempt fails, then tries to find archived version with
// the name @name.mps.gz. If archived version is found, unarchives and reads it.
template <typename Field = double>
Problem<Field> get_archived(const std::string& name) {
  return mps::read<Field>(unarchive(name));
}
