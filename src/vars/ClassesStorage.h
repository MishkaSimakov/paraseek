#pragma once

#include <unordered_map>
#include <vector>

namespace seekers {

// Stores information about patterns for small rows
class ClassesStorage {
  struct ClassInfo {
    size_t map_column{0};
    std::unordered_map<size_t, size_t> map;
    std::vector<size_t> counts;
  };

  std::vector<size_t> free_classes_;
  std::vector<ClassInfo> classes_;

  const size_t max_diff;

  void extend_classes_storage(size_t new_size) {
    size_t old_size = classes_.size();
    classes_.resize(new_size);

    for (size_t i = old_size; i < new_size; ++i) {
      classes_[i].counts.resize((max_diff + 1) * (max_diff + 1), 0);
      free_classes_.push_back(i);
    }
  }

 public:
  ClassesStorage(size_t max_diff, size_t rows_count) : max_diff(max_diff) {
    extend_classes_storage(std::max(rows_count, 64uz));
  }

  size_t& get_rows_count(size_t class_id, size_t cnt_0, size_t cnt_1) {
    return classes_[class_id].counts[cnt_0 * (max_diff + 1) + cnt_1];
  }

  void try_free_class(size_t class_id) {
    bool is_zero = true;

    for (const size_t count : classes_[class_id].counts) {
      if (count != 0) {
        is_zero = false;
        break;
      }
    }

    if (is_zero) {
      free_classes_.push_back(class_id);
    }
  }

  size_t get_class(size_t prev_class, size_t hash, size_t column) {
    if (free_classes_.empty()) {
      extend_classes_storage(classes_.size() * 2);
    }

    if (classes_[prev_class].map_column != column) {
      classes_[prev_class].map.clear();
      classes_[prev_class].map_column = column;
    }

    auto [itr, inserted] =
        classes_[prev_class].map.emplace(hash, free_classes_.back());

    if (inserted) {
      free_classes_.pop_back();
    }

    return itr->second;
  }

  size_t pop_free_class() {
    if (free_classes_.empty()) {
      extend_classes_storage(classes_.size() * 2);
    }

    size_t free_class = free_classes_.back();
    free_classes_.pop_back();
    return free_class;
  }

  std::vector<std::vector<size_t>> accumulate_counts() {
    const size_t counts_size = (max_diff + 1) * (max_diff + 1);

    for (size_t i = 1; i < classes_.size(); ++i) {
      for (size_t j = 0; j < counts_size; ++j) {
        classes_[i].counts[j] += classes_[i - 1].counts[j];
      }
    }

    const size_t last_class = classes_.size() - 1;

    std::vector<std::vector<size_t>> total(max_diff + 1);

    for (size_t i = 0; i <= max_diff; ++i) {
      total[i].resize(max_diff + 1, 0);

      for (size_t j = 0; j <= max_diff; ++j) {
        total[i][j] = get_rows_count(last_class, i, j);
      }
    }

    return total;
  }
};

}  // namespace seekers
