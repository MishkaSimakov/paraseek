#pragma once

#include "utils/Time.h"

namespace seekers {

struct Statistics {
  // false positives + true positives
  size_t pairs_considered{0};

  timing::Duration big_rows_duration;
  timing::Duration small_rows_duration;
};

}  // namespace seekers
