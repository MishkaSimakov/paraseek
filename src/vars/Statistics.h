#pragma once

#include "utils/Time.h"

namespace seekers {

struct Statistics {
  // pairs of rows considered in seek_tables
  size_t pairs_considered{0};

  // SmallRowEntries considered in seek_small
  size_t entries_considered{0};

  timing::Duration big_rows_duration;
  timing::Duration small_rows_duration;
};

}  // namespace seekers
