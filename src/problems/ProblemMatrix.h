#pragma once

#include <format>

#include "matrix/CSCMatrix.h"
#include "mps/MPS.h"
#include "utils/Paths.h"

inline Problem get_problem_matrices(std::string name,
                                         bool replace_inequalities = false) {
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

  auto reader = mps::MPSReader<double>(mps::Format::FREE);
  reader.read(path);

  return reader.get_problem(replace_inequalities);
}
